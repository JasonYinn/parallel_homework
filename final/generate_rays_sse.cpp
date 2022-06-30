#include <stdio.h>
#include <time.h>
#include <emmintrin.h>

int main(int argc, char* argv[]) {

    float frac = 1.0;
    if (argc > 1) {
        frac = atof(argv[1]);
    }

    int W = (int)(1296 * frac);
    int H = (int)(968 * frac);
    float fx = 749.2899503552605 / 2.0 * frac;
    float fy = 749.2899503552605 / 2.0 * frac;
    float cx = W * 0.5;
    float cy = H * 0.5;
    float c2w[4][4] = {{0.9923269748687744, 0.04875849187374115, -0.11362089961767197, -0.3664160668849945},
                    {-0.1236410066485405, 0.39132946729660034, -0.911906898021698, -2.359870433807373},
                    {0.0, 0.9189580678939819, 0.39435532689094543, 1.2530806064605713},
                    {0.0, 0.0, 0.0, 1.0}};
    
    
    float *rays = new float[H * W * 8];
    float* rays_ptr = rays;
    

    struct timespec sts, ets;
    timespec_get(&sts, TIME_UTC);
    
    // float cam_loc[4] = {c2w[0][0], c2w[0][1], c2w[0][2], c2w[0][3]};
    __m128 cam_pos = _mm_set_ps(c2w[0][3], c2w[0][2], c2w[0][1], c2w[0][0]);
    float c2w_T[4][4];
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            c2w_T[i][j] = c2w[j][i];
        }
    }

    __m128 dir_at_ij;
    __m128 temp_dir;
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            _mm_store_ps(rays_ptr + i * W * 8 + j * 8, cam_pos);

            float dirs[4] = {((float)j - cx) / fx, -((float)i - cy) / fy, -1.0, 0.0};
            dir_at_ij = _mm_set_ps1(0);
            for (int k = 0; k < 3; ++k) {
                __m128 c2w_row = _mm_loadu_ps(&c2w_T[0][0] + k * 4);
                temp_dir = _mm_set_ps1(dirs[k]);
                __m128 temp = _mm_mul_ps(temp_dir, c2w_row);
                dir_at_ij = _mm_add_ps(dir_at_ij, temp);

            }
            _mm_store_ps(rays_ptr + i * W * 8 + j * 8 + 4, dir_at_ij);
        }
    }

    timespec_get(&ets, TIME_UTC);
    time_t dsec = ets.tv_sec - sts.tv_sec;
    long dnsec = ets.tv_nsec - sts.tv_nsec;
    printf ("x86 sse, matrix size (%d, %d), time consume: %ld.%09lds\n", H, W, dsec, dnsec);
    
    // printf ("%f, %f, %f, %f\n", rays[0][0][4], rays[0][0][5], rays[0][0][6], rays[0][0][7]);

    return 0;
}