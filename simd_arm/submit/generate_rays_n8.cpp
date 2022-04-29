#include <stdio.h>
#include <time.h>
#include <arm_neon.h>

int main(void) {
    
    float frac = 0.125;

    int W = (int)(400 * frac);
    int H = (int)(400 * frac);
    float fx = 749.2899503552605 / 2.0 * frac;
    float fy = 749.2899503552605 / 2.0 * frac;
    float cx = W * 0.5;
    float cy = H * 0.5;
    float32_t c2w[4][4] = {{0.9923269748687744, 0.04875849187374115, -0.11362089961767197, -0.3664160668849945},
                    {-0.1236410066485405, 0.39132946729660034, -0.911906898021698, -2.359870433807373},
                    {0.0, 0.9189580678939819, 0.39435532689094543, 1.2530806064605713},
                    {0.0, 0.0, 0.0, 1.0}};
    
    
    float32_t rays[H][W][8];
    float32_t* rays_ptr = &rays[0][0][0];
    

    struct timespec sts, ets;
    timespec_get(&sts, TIME_UTC);
    
    float32_t cam_loc[4] = {c2w[0][0], c2w[0][1], c2w[0][2], c2w[0][3]};
    float32x4_t cam_pos = vld1q_f32(cam_loc);
    float32x4_t all_zeros = vmovq_n_f32(0);
    float32_t c2w_T[4][4];
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            c2w_T[i][j] = c2w[j][i];
        }
    }

    float32x4_t dir_at_ij;
    float32x4_t temp_dir;
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            vst1q_f32(rays_ptr + i * W * 8 + j * 8, cam_pos);

            float32_t dirs[4] = {((float32_t)j - cx) / fx, -((float32_t)i - cy) / fy, -1.0, 0.0};
            dir_at_ij = vmovq_n_f32(0);
            for (int k = 0; k < 3; ++k) {
                float32x4_t c2w_row = vld1q_f32((float32_t*)&c2w_T[0][0] + k * 4);
                temp_dir = vmovq_n_f32(dirs[k]);
                float32x4_t temp = vmulq_f32(temp_dir, c2w_row);
                dir_at_ij = vaddq_f32(dir_at_ij, temp);

            }
            vst1q_f32(rays_ptr + i * W * 8 + j * 8 + 4, dir_at_ij);
        }
    }

    timespec_get(&ets, TIME_UTC);
    time_t dsec = ets.tv_sec - sts.tv_sec;
    long dnsec = ets.tv_nsec - sts.tv_nsec;
    printf ("%ld.%09lds\n", dsec, dnsec);
    
    // printf ("%f, %f, %f, %f\n", rays[0][0][4], rays[0][0][5], rays[0][0][6], rays[0][0][7]);

    return 0;
}