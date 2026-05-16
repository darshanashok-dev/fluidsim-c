#ifndef FLUID_H
#define FLUID_H

typedef struct {
    int    N;
    float  dt, diff, visc, vorticity;
    float *u, *v, *u_prev, *v_prev;
    float *dens_r, *dens_r_prev;
    float *dens_g, *dens_g_prev;
    float *dens_b, *dens_b_prev;
    float *curl;
    float *obstacles;
} FluidGrid;

FluidGrid *fluid_create(int N, float dt, float diff, float visc, float vorticity);
void       fluid_free(FluidGrid *f);
void       fluid_step(FluidGrid *f);
void       fluid_add_velocity(FluidGrid *f, int x, int y, float amtX, float amtY);
void       fluid_add_density(FluidGrid *f, int x, int y, float r, float g, float b);
void       fluid_add_obstacle(FluidGrid *f, int x, int y, float radius, float amount);

#endif // FLUID_H
