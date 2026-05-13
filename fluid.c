#include "fluid.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SIZE(N) ((N+2)*(N+2))
#define IX(i, j, N) ((i) + (N+2)*(j))
#define SWAP(x0,x) { float *tmp = x0; x0 = x; x = tmp; }
#define ITERATIONS 10

FluidGrid *fluid_create(int N, float dt, float diff, float visc, float vorticity) {
    FluidGrid *f = (FluidGrid *)malloc(sizeof(FluidGrid));
    f->N = N;
    f->dt = dt;
    f->diff = diff;
    f->visc = visc;
    f->vorticity = vorticity;

    int size = SIZE(N);
    float *slab = (float *)calloc(12 * size, sizeof(float));
    f->u           = slab + 0 * size;
    f->v           = slab + 1 * size;
    f->u_prev      = slab + 2 * size;
    f->v_prev      = slab + 3 * size;
    f->dens_r      = slab + 4 * size;
    f->dens_r_prev = slab + 5 * size;
    f->dens_g      = slab + 6 * size;
    f->dens_g_prev = slab + 7 * size;
    f->dens_b      = slab + 8 * size;
    f->dens_b_prev = slab + 9 * size;
    f->curl        = slab + 10 * size;
    f->obstacles   = slab + 11 * size;

    return f;
}

void fluid_free(FluidGrid *f) {
    if (f) {
        free(f->u); 
        free(f);
    }
}

void fluid_add_velocity(FluidGrid *f, int x, int y, float amtX, float amtY) {
    int N = f->N;
    if (x < 0 || x > N + 1 || y < 0 || y > N + 1) return;
    f->u_prev[IX(x, y, N)] += amtX;
    f->v_prev[IX(x, y, N)] += amtY;
}

void fluid_add_density(FluidGrid *f, int x, int y, float r, float g, float b) {
    int N = f->N;
    if (x < 0 || x > N + 1 || y < 0 || y > N + 1) return;
    f->dens_r_prev[IX(x, y, N)] += r;
    f->dens_g_prev[IX(x, y, N)] += g;
    f->dens_b_prev[IX(x, y, N)] += b;
}

void fluid_add_obstacle(FluidGrid *f, int x, int y, float radius, float amount) {
    int N = f->N;
    int r = (int)ceil(radius);
    for (int j = -r; j <= r; j++) {
        for (int i = -r; i <= r; i++) {
            if (i * i + j * j <= radius * radius) {
                int cx = x + i;
                int cy = y + j;
                if (cx >= 1 && cx <= N && cy >= 1 && cy <= N) {
                    f->obstacles[IX(cx, cy, N)] += amount;
                    if (f->obstacles[IX(cx, cy, N)] > 1.0f) f->obstacles[IX(cx, cy, N)] = 1.0f;
                }
            }
        }
    }
}

static void add_source(int N, float *x, float *s, float dt) {
    int size = SIZE(N);
    #pragma omp parallel for
    for (int i = 0; i < size; i++) {
        x[i] += dt * s[i];
    }
}

static void set_bnd(int N, int b, float *x, float *obstacles) {
    for (int i = 1; i <= N; i++) {
        x[IX(0, i, N)]     = b == 1 ? -x[IX(1, i, N)] : x[IX(1, i, N)];
        x[IX(N + 1, i, N)] = b == 1 ? -x[IX(N, i, N)] : x[IX(N, i, N)];
        x[IX(i, 0, N)]     = b == 2 ? -x[IX(i, 1, N)] : x[IX(i, 1, N)];
        x[IX(i, N + 1, N)] = b == 2 ? -x[IX(i, N, N)] : x[IX(i, N, N)];
    }
    x[IX(0, 0, N)]         = 0.5f * (x[IX(1, 0, N)] + x[IX(0, 1, N)]);
    x[IX(0, N + 1, N)]     = 0.5f * (x[IX(1, N + 1, N)] + x[IX(0, N, N)]);
    x[IX(N + 1, 0, N)]     = 0.5f * (x[IX(N, 0, N)] + x[IX(N + 1, 1, N)]);
    x[IX(N + 1, N + 1, N)] = 0.5f * (x[IX(N, N + 1, N)] + x[IX(N + 1, N, N)]);

    if (obstacles) {
        #pragma omp parallel for collapse(2)
        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                if (obstacles[IX(i, j, N)] > 0.5f) {
                    x[IX(i, j, N)] = 0.0f;
                }
            }
        }
    }
}

static void lin_solve(int N, int b, float *x, float *x0, float a, float c, float *obstacles) {
    float inv_c = 1.0f / c;
    for (int k = 0; k < ITERATIONS; k++) {
        #pragma omp parallel for collapse(2)
        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                x[IX(i, j, N)] = (x0[IX(i, j, N)] + a * (x[IX(i - 1, j, N)] + x[IX(i + 1, j, N)] + x[IX(i, j - 1, N)] + x[IX(i, j + 1, N)])) * inv_c;
            }
        }
        set_bnd(N, b, x, obstacles);
    }
}

static void diffuse(int N, int b, float *x, float *x0, float diff, float dt, float *obstacles) {
    float a = dt * diff * N * N;
    lin_solve(N, b, x, x0, a, 1.0f + 4.0f * a, obstacles);
}

static void advect(int N, int b, float *d, float *d0, float *u, float *v, float dt, float *obstacles) {
    int i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;

    dt0 = dt * N;
    #pragma omp parallel for collapse(2) private(x, y, i0, j0, i1, j1, s0, t0, s1, t1)
    for (int j = 1; j <= N; j++) {
        for (int i = 1; i <= N; i++) {
            x = i - dt0 * u[IX(i, j, N)];
            y = j - dt0 * v[IX(i, j, N)];
            if (isnan(x)) x = 0.5f;
            if (isnan(y)) y = 0.5f;
            if (x < 0.5f) x = 0.5f;
            if (x > N + 0.5f) x = N + 0.5f;
            i0 = (int)x;
            i1 = i0 + 1;
            if (y < 0.5f) y = 0.5f;
            if (y > N + 0.5f) y = N + 0.5f;
            j0 = (int)y;
            j1 = j0 + 1;
            s1 = x - i0;
            s0 = 1.0f - s1;
            t1 = y - j0;
            t0 = 1.0f - t1;
            d[IX(i, j, N)] = s0 * (t0 * d0[IX(i0, j0, N)] + t1 * d0[IX(i0, j1, N)]) +
                          s1 * (t0 * d0[IX(i1, j0, N)] + t1 * d0[IX(i1, j1, N)]);
        }
    }
    set_bnd(N, b, d, obstacles);
}

static void project(int N, float *u, float *v, float *p, float *div, float *obstacles) {
    #pragma omp parallel for collapse(2)
    for (int j = 1; j <= N; j++) {
        for (int i = 1; i <= N; i++) {
            div[IX(i, j, N)] = -0.5f * (u[IX(i + 1, j, N)] - u[IX(i - 1, j, N)] + v[IX(i, j + 1, N)] - v[IX(i, j - 1, N)]) / N;
            p[IX(i, j, N)] = 0;
        }
    }
    set_bnd(N, 0, div, obstacles);
    set_bnd(N, 0, p, obstacles);
    lin_solve(N, 0, p, div, 1.0f, 4.0f, obstacles);

    #pragma omp parallel for collapse(2)
    for (int j = 1; j <= N; j++) {
        for (int i = 1; i <= N; i++) {
            u[IX(i, j, N)] -= 0.5f * N * (p[IX(i + 1, j, N)] - p[IX(i - 1, j, N)]);
            v[IX(i, j, N)] -= 0.5f * N * (p[IX(i, j + 1, N)] - p[IX(i, j - 1, N)]);
        }
    }
    set_bnd(N, 1, u, obstacles);
    set_bnd(N, 2, v, obstacles);
}

static void vorticity_confinement(int N, float *u, float *v, float *curl, float vorticity, float dt, float *obstacles) {
    #pragma omp parallel for collapse(2)
    for (int j = 1; j < N; j++) {
        for (int i = 1; i < N; i++) {
            float dw_dy = (u[IX(i, j + 1, N)] - u[IX(i, j - 1, N)]) * 0.5f * N;
            float dv_dx = (v[IX(i + 1, j, N)] - v[IX(i - 1, j, N)]) * 0.5f * N;
            curl[IX(i, j, N)] = fabs(dv_dx - dw_dy);
        }
    }

    #pragma omp parallel for collapse(2)
    for (int j = 2; j < N - 1; j++) {
        for (int i = 2; i < N - 1; i++) {
            float dw_dx = (curl[IX(i + 1, j, N)] - curl[IX(i - 1, j, N)]) * 0.5f * N;
            float dw_dy = (curl[IX(i, j + 1, N)] - curl[IX(i, j - 1, N)]) * 0.5f * N;
            
            float length = sqrtf(dw_dx * dw_dx + dw_dy * dw_dy) + 0.000001f;
            dw_dx /= length;
            dw_dy /= length;

            float v_curl = (v[IX(i + 1, j, N)] - v[IX(i - 1, j, N)]) * 0.5f * N - (u[IX(i, j + 1, N)] - u[IX(i, j - 1, N)]) * 0.5f * N;

            u[IX(i, j, N)] += dt * vorticity * ( dw_dy * v_curl);
            v[IX(i, j, N)] += dt * vorticity * (-dw_dx * v_curl);
            
            float MAX_VEL = 200.0f;
            u[IX(i, j, N)] = fmaxf(-MAX_VEL, fminf(MAX_VEL, u[IX(i, j, N)]));
            v[IX(i, j, N)] = fmaxf(-MAX_VEL, fminf(MAX_VEL, v[IX(i, j, N)]));
        }
    }
    set_bnd(N, 1, u, obstacles);
    set_bnd(N, 2, v, obstacles);
}

static void vel_step(int N, float *u, float *v, float *u0, float *v0, float visc, float dt, float vorticity, float *curl, float *obstacles) {
    add_source(N, u, u0, dt);
    add_source(N, v, v0, dt);
    
    if (vorticity > 0.0f) {
        vorticity_confinement(N, u, v, curl, vorticity, dt, obstacles);
    }

    SWAP(u0, u);
    diffuse(N, 1, u, u0, visc, dt, obstacles);
    SWAP(v0, v);
    diffuse(N, 2, v, v0, visc, dt, obstacles);
    project(N, u, v, u0, v0, obstacles);
    SWAP(u0, u);
    SWAP(v0, v);
    advect(N, 1, u, u0, u0, v0, dt, obstacles);
    advect(N, 2, v, v0, u0, v0, dt, obstacles);
    project(N, u, v, u0, v0, obstacles);
}

static void dens_step(int N, float *x, float *x0, float *u, float *v, float diff, float dt, float *obstacles) {
    add_source(N, x, x0, dt);
    SWAP(x0, x);
    diffuse(N, 0, x, x0, diff, dt, obstacles);
    SWAP(x0, x);
    advect(N, 0, x, x0, u, v, dt, obstacles);
}

void fluid_step(FluidGrid *f) {
    int N = f->N;
    vel_step(N, f->u, f->v, f->u_prev, f->v_prev, f->visc, f->dt, f->vorticity, f->curl, f->obstacles);
    
    dens_step(N, f->dens_r, f->dens_r_prev, f->u, f->v, f->diff, f->dt, f->obstacles);
    dens_step(N, f->dens_g, f->dens_g_prev, f->u, f->v, f->diff, f->dt, f->obstacles);
    dens_step(N, f->dens_b, f->dens_b_prev, f->u, f->v, f->diff, f->dt, f->obstacles);

    int size = SIZE(N);
    memset(f->u_prev, 0, size * sizeof(float));
    memset(f->v_prev, 0, size * sizeof(float));
    memset(f->dens_r_prev, 0, size * sizeof(float));
    memset(f->dens_g_prev, 0, size * sizeof(float));
    memset(f->dens_b_prev, 0, size * sizeof(float));
}
