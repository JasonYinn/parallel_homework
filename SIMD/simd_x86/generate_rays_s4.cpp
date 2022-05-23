#include <stdio.h>
#include <time.h>

int main(int argc, char* argv[]) {

    float frac = 0.25;

    int W = (int)(400 * frac);
    int H = (int)(400 * frac);
    float fx = 749.2899503552605 / 2.0 * frac;
    float fy = 749.2899503552605 / 2.0 * frac;
    float cx = W * 0.5;
    float cy = H * 0.5;
    float c2w[4][4] = {{0.9923269748687744, 0.04875849187374115, -0.11362089961767197, -0.3664160668849945},
                    {-0.1236410066485405, 0.39132946729660034, -0.911906898021698, -2.359870433807373},
                    {0.0, 0.9189580678939819, 0.39435532689094543, 1.2530806064605713},
                    {0.0, 0.0, 0.0, 1.0}};
    
    float rays[H][W][8];

    struct timespec sts, ets;
    timespec_get(&sts, TIME_UTC);
    
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            for (int k = 0; k < 4; ++k) {
                rays[i][j][k] = c2w[3][k];
                rays[i][j][k + 4] = 0.0;
            }
            float dirs[4] = {((float)j - cx) / fx, -((float)i - cy) / fy, -1.0, 0.0};
            for (int k = 0; k < 4; ++k) {
                for (int z = 0; z < 4; ++z) {
                    rays[i][j][k + 4] += dirs[z] * c2w[k][z];
                }
            }
        }
    }

    timespec_get(&ets, TIME_UTC);
    time_t dsec = ets.tv_sec - sts.tv_sec;
    long dnsec = ets.tv_nsec - sts.tv_nsec;
    printf ("%ld.%09lds\n", dsec, dnsec);
    // printf ("%f, %f, %f, %f\n", rays[0][0][4], rays[0][0][5], rays[0][0][6], rays[0][0][7]);
    return 0;
}