#include <stdio.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include <emmintrin.h>
#include <string.h>

struct render_params
{
    float cx;
    float cy;
    float fx;
    float fy;
    int H;
    int W;
};

struct pthread_params {
    int thread_id;
    int stride;
    float *rays;
    float *c2w;
    float *c2w_T;
    render_params renderParams;
};


void generate_ray_per_pixel(float *rays, float *c2w, render_params renderParams, int i, int j) {
    float dirs[4] = {((float)j - renderParams.cx) / renderParams.fx, -((float)i - renderParams.cy) / renderParams.fy, -1.0, 0.0};
    for (int k = 0; k < 4; ++k) {
        float temp = 0.;
        rays[i * renderParams.W * 8 + j * 8 + k] = c2w[3 * 4 + k];
        for (int z = 0; z < 4; ++z) {
            temp += dirs[z] * c2w[k * 4 + z];
        }
        rays[i * renderParams.W * 8 + j * 8 + 4 + k] = temp;
    }
}

void *ray_thread_row(void *pthreadParams_ptr) {
    struct pthread_params *pthreadParams = (struct pthread_params *)pthreadParams_ptr;

    int threadIdx = pthreadParams->thread_id;
    int stride = pthreadParams->stride;
    struct render_params renderParams = pthreadParams->renderParams;
    float *rays = pthreadParams->rays;
    float *c2w = pthreadParams->c2w;
    float *c2w_T = pthreadParams->c2w_T;
    
    __m128 dir_at_ij;
    __m128 temp_dir;
    __m128 cam_pos = _mm_set_ps(c2w[3], c2w[7], c2w[11], c2w[15]);
    for (int i = threadIdx; i < renderParams.H; i += stride) {
        for (int j = 0; j < renderParams.W; j += 1) {
            _mm_store_ps(rays + i * renderParams.W * 8 + j * 8, cam_pos);

            float dirs[4] = {((float)j - renderParams.cx) / renderParams.fx, -((float)i - renderParams.cy) / renderParams.fy, -1.0, 0.0};
            dir_at_ij = _mm_set_ps1(0);
            for (int k = 0; k < 3; ++k) {
                __m128 c2w_row = _mm_loadu_ps(c2w_T + k * 4);
                temp_dir = _mm_set_ps1(dirs[k]);
                __m128 temp = _mm_mul_ps(temp_dir, c2w_row);
                dir_at_ij = _mm_add_ps(dir_at_ij, temp);
                // printf ("Thread %d: %f, %f\n", pthreadParams->thread_id, rays[i * renderParams.W * 8 + j * 8 + 4 + k], temp);
            }
            _mm_store_ps(rays + i * renderParams.W * 8 + j * 8 + 4, dir_at_ij);
        }
    }
    // printf ("Thread %d finished\n", pthreadParams->thread_id);
    pthread_exit(NULL);

}

void *ray_thread_col(void *pthreadParams_ptr) {
    struct pthread_params *pthreadParams = (struct pthread_params *)pthreadParams_ptr;

    int threadIdx = pthreadParams->thread_id;
    int stride = pthreadParams->stride;
    struct render_params renderParams = pthreadParams->renderParams;
    float *rays = pthreadParams->rays;
    float *c2w = pthreadParams->c2w;
    float *c2w_T = pthreadParams->c2w_T;
    
    __m128 dir_at_ij;
    __m128 temp_dir;
    __m128 cam_pos = _mm_set_ps(c2w[3], c2w[7], c2w[11], c2w[15]);
    for (int i = 0; i < renderParams.H; i += 1) {
        for (int j = threadIdx; j < renderParams.W; j += stride) {
            _mm_store_ps(rays + i * renderParams.W * 8 + j * 8, cam_pos);

            float dirs[4] = {((float)j - renderParams.cx) / renderParams.fx, -((float)i - renderParams.cy) / renderParams.fy, -1.0, 0.0};
            dir_at_ij = _mm_set_ps1(0);
            for (int k = 0; k < 3; ++k) {
                __m128 c2w_row = _mm_loadu_ps(c2w_T + k * 4);
                temp_dir = _mm_set_ps1(dirs[k]);
                __m128 temp = _mm_mul_ps(temp_dir, c2w_row);
                dir_at_ij = _mm_add_ps(dir_at_ij, temp);
                // printf ("Thread %d: %f, %f\n", pthreadParams->thread_id, rays[i * renderParams.W * 8 + j * 8 + 4 + k], temp);
            }
            _mm_store_ps(rays + i * renderParams.W * 8 + j * 8 + 4, dir_at_ij);
        }
    }
    // printf ("Thread %d finished\n", pthreadParams->thread_id);
    pthread_exit(NULL);

}

void *ray_thread_onedim(void *pthreadParams_ptr) {
    struct pthread_params *pthreadParams = (struct pthread_params *)pthreadParams_ptr;

    int threadIdx = pthreadParams->thread_id;
    int stride = pthreadParams->stride;
    struct render_params renderParams = pthreadParams->renderParams;
    float *rays = pthreadParams->rays;
    float *c2w = pthreadParams->c2w;
    float *c2w_T = pthreadParams->c2w_T;
    int end = renderParams.H * renderParams.W;
    
    __m128 dir_at_ij;
    __m128 temp_dir;
    __m128 cam_pos = _mm_set_ps(c2w[3], c2w[7], c2w[11], c2w[15]);
    for (int i = threadIdx; i < end; i += stride) {
            int j = i % renderParams.H;
            _mm_store_ps(rays + i * 8, cam_pos);

            float dirs[4] = {((float)j - renderParams.cx) / renderParams.fx, -((float)i - renderParams.cy) / renderParams.fy, -1.0, 0.0};
            dir_at_ij = _mm_set_ps1(0);
            for (int k = 0; k < 3; ++k) {
                __m128 c2w_row = _mm_loadu_ps(c2w_T + k * 4);
                temp_dir = _mm_set_ps1(dirs[k]);
                __m128 temp = _mm_mul_ps(temp_dir, c2w_row);
                dir_at_ij = _mm_add_ps(dir_at_ij, temp);
                // printf ("Thread %d: %f, %f\n", pthreadParams->thread_id, rays[i * renderParams.W * 8 + j * 8 + 4 + k], temp);
            }
            _mm_store_ps(rays + i * 8 + 4, dir_at_ij);
    }
    // printf ("Thread %d finished\n", pthreadParams->thread_id);
    pthread_exit(NULL);

}

int main(int argc, char* argv[]) {

    float frac = 1.0;
    int NUM_THREADS = 8;
    void *(* ray_thread)(void *);
    ray_thread = ray_thread_row;
    const char *orders[3] = {"row", "col", "oneDim"};
    const char *order;

    if (argc > 1) {
        frac = atof(argv[1]);
    }
    if (argc > 2) {
        NUM_THREADS = atoi(argv[2]);
    }
    if (argc > 3) {
        if (strcmp(argv[3], "row") == 0) {
            ray_thread = ray_thread_row;
            order = orders[0];
        }
        else if (strcmp(argv[3], "col") == 0)
        {
            ray_thread = ray_thread_col;
            order = orders[1];
        }
        else if (strcmp(argv[3], "oneDim") == 0)
        {
            ray_thread = ray_thread_onedim;
            order = orders[2];
        }
        
    }
    else {
        ray_thread = ray_thread_row;
        order = orders[0];
    }

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
    
    float c2w_T[4][4];
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            c2w_T[i][j] = c2w[j][i];
        }
    }

    float rays[H][W][8];
    // float rays2[H][W][8];
    // float err = 0;

    render_params renderParams;
    renderParams.cx = cx;
    renderParams.cy = cy;
    renderParams.fx = fx;
    renderParams.fy = fy;
    renderParams.H  =  H;
    renderParams.W  =  W;

    pthread_params pthreadParams_array[NUM_THREADS];
    pthread_t threads[NUM_THREADS];
    int rc;

    for (int threadIdx = 0; threadIdx < NUM_THREADS; ++threadIdx) {
        pthreadParams_array[threadIdx].c2w = &c2w[0][0];
        pthreadParams_array[threadIdx].c2w_T = &c2w_T[0][0];
        pthreadParams_array[threadIdx].rays = &rays[0][0][0];
        pthreadParams_array[threadIdx].stride = NUM_THREADS;
        pthreadParams_array[threadIdx].thread_id = threadIdx;
        pthreadParams_array[threadIdx].renderParams = renderParams;

        rc = pthread_create(&threads[threadIdx], NULL, ray_thread, (void *)&pthreadParams_array[threadIdx]);
        if (rc) {
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }

    struct timespec sts, ets;
    timespec_get(&sts, TIME_UTC);
    

    for (int threadIdx = 0; threadIdx < NUM_THREADS; ++threadIdx) {
        pthread_join(threads[threadIdx], NULL);
    }

    timespec_get(&ets, TIME_UTC);
    time_t dsec = ets.tv_sec - sts.tv_sec;
    long dnsec = ets.tv_nsec - sts.tv_nsec;
    printf ("x86 sse pthread, num_threads %d, matrix size (%d, %d) in %s order, time consume: %ld.%09lds\n", NUM_THREADS, H, W, order, dsec, dnsec);
    pthread_exit(NULL);

    return 0;
}