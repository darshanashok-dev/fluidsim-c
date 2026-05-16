#include "raylib.h"
#include "fluid.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define IX(i, j, N) ((i) + (N+2)*(j))

int main(void) {
    InitWindow(800, 800, "Fluid Sim");
    SetTargetFPS(60);

    int N = 128;
    float dt = 0.1f;
    float diffusion = 0.0001f;
    float viscosity = 0.0000001f;
    float vorticity = 0.0f; // Disabled for beautiful smooth aesthetic

    FluidGrid *f = fluid_create(N, dt, diffusion, viscosity, vorticity);

    Color *pixels = (Color *)malloc(N * N * sizeof(Color));
    Image img = {
        .data = pixels,
        .width = N,
        .height = N,
        .mipmaps = 1,
        .format = PIXELFORMAT_UNCOMPRESSED_R8G8B8A8
    };
    Texture2D tex = LoadTextureFromImage(img);

    float prevMouseX = GetMouseX();
    float prevMouseY = GetMouseY();

    while (!WindowShouldClose()) {
        // Hot-Swapping Resolution
        int newN = N;
        if (IsKeyPressed(KEY_ONE)) newN = 64;
        if (IsKeyPressed(KEY_TWO)) newN = 128;
        if (IsKeyPressed(KEY_THREE)) newN = 256;
        
        if (newN != N) {
            fluid_free(f);
            free(pixels);
            UnloadTexture(tex);
            N = newN;
            f = fluid_create(N, dt, diffusion, viscosity, vorticity);
            pixels = (Color *)malloc(N * N * sizeof(Color));
            img.data = pixels;
            img.width = N;
            img.height = N;
            tex = LoadTextureFromImage(img);
        }

        float mouseX = GetMouseX();
        float mouseY = GetMouseY();

        int cx = (int)(mouseX * N / 800.0f);
        int cy = (int)(mouseY * N / 800.0f);
        float dx = mouseX - prevMouseX;
        float dy = mouseY - prevMouseY;

        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            float t = GetTime() * 2.0f;
            float r = (sinf(t) * 0.5f + 0.5f) * 2000.0f;
            float g = (sinf(t + 2.094f) * 0.5f + 0.5f) * 2000.0f;
            float b = (sinf(t + 4.188f) * 0.5f + 0.5f) * 2000.0f;
            
            int radius = (N == 256) ? 8 : ((N == 128) ? 4 : 2);
            for (int j = -radius; j <= radius; j++) {
                for (int i = -radius; i <= radius; i++) {
                    float d = sqrtf(i*i + j*j);
                    if (d <= radius) {
                        float amt = 1.0f - (d / radius);
                        int nx = cx + i;
                        int ny = cy + j;
                        if (nx >= 1 && nx <= N && ny >= 1 && ny <= N) {
                            fluid_add_velocity(f, nx, ny, dx * 5.0f * amt, dy * 5.0f * amt);
                            fluid_add_density(f, nx, ny, r * amt, g * amt, b * amt);
                        }
                    }
                }
            }
        }

        if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
            float rad = (N == 256) ? 6.0f : ((N == 128) ? 3.0f : 1.5f);
            fluid_add_obstacle(f, cx, cy, rad, 1.0f);
        }

        if (IsKeyPressed(KEY_R)) {
            int size = (N+2)*(N+2);
            memset(f->u, 0, size*sizeof(float));
            memset(f->v, 0, size*sizeof(float));
            memset(f->u_prev, 0, size*sizeof(float));
            memset(f->v_prev, 0, size*sizeof(float));
            memset(f->dens_r, 0, size*sizeof(float));
            memset(f->dens_g, 0, size*sizeof(float));
            memset(f->dens_b, 0, size*sizeof(float));
            memset(f->dens_r_prev, 0, size*sizeof(float));
            memset(f->dens_g_prev, 0, size*sizeof(float));
            memset(f->dens_b_prev, 0, size*sizeof(float));
            memset(f->obstacles, 0, size*sizeof(float));
        }

        fluid_step(f);

        #pragma omp parallel for collapse(2)
        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                int idx = IX(i, j, N);
                float r = f->dens_r[idx];
                float g = f->dens_g[idx];
                float b = f->dens_b[idx];
                
                if (r > 255.0f) r = 255.0f;
                if (g > 255.0f) g = 255.0f;
                if (b > 255.0f) b = 255.0f;
                
                Color c = { (unsigned char)r, (unsigned char)g, (unsigned char)b, 255 };

                // Draw obstacles
                if (f->obstacles[idx] > 0.1f) {
                    unsigned char obs = (unsigned char)(f->obstacles[idx] * 100.0f);
                    c = (Color){ obs, obs, obs, 255 };
                }

                pixels[(j - 1) * N + (i - 1)] = c;

                // Dissipation
                f->dens_r[idx] *= 0.995f;
                f->dens_g[idx] *= 0.995f;
                f->dens_b[idx] *= 0.995f;
                f->obstacles[idx] *= 0.99f;
            }
        }

        UpdateTexture(tex, pixels);

        prevMouseX = mouseX;
        prevMouseY = mouseY;

        BeginDrawing();
            ClearBackground(BLACK);
            DrawTexturePro(tex, 
                (Rectangle){0, 0, (float)N, (float)N}, 
                (Rectangle){0, 0, 800.0f, 800.0f}, 
                (Vector2){0, 0}, 0.0f, WHITE);
            
            DrawText(TextFormat("FPS: %i | Res: %ix%i (Keys 1/2/3)", GetFPS(), N, N), 10, 10, 20, RAYWHITE);
        EndDrawing();
    }

    free(pixels);
    fluid_free(f);
    CloseWindow();

    return 0;
}
