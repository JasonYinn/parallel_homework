#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include "svjpeg.hpp"

struct render_params
{
    float cx;
    float cy;
    float fx;
    float fy;
    int H;
    int W;
    int dimx, dimy, dimz;
    float origins[3];
    float voxel_size;
};

template <typename scalar_t>
int readBin(const char* fname, scalar_t *arr, int len) {
    std::ifstream in(fname, std::ios::in | std::ios::binary);

    in.read((char *) arr, sizeof(scalar_t) * len);

    return in.gcount();
}

__device__ __forceinline__ float dot(float3 vec1, float3 vec2) {
    return vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z;
}

__device__ __forceinline__ float norm(float3 vec) {
    return sqrtf32(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}

__global__
void generate_ray_per_pixel(float *rays, float *c2w, render_params renderParams) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= 0 && i < renderParams.H && j >= 0 && j < renderParams.W){    
        float3 dirs = make_float3(((float)j - renderParams.cx) / renderParams.fx, -((float)i - renderParams.cy) / renderParams.fy, -1.0);
        float3 pose_line;
        for (int k = 0; k < 3; ++k) {
            rays[i * renderParams.W * 8 + j * 8 + k] = c2w[k * 4 + 3];
            pose_line = make_float3(c2w[k * 4 + 0], c2w[k * 4 + 1], c2w[k * 4 + 2]);
            // for (int z = 0; z < 3; ++z) {
            //     temp += dirs[z] * c2w[k * 4 + z];
            // }
            rays[i * renderParams.W * 8 + j * 8 + 4 + k] = dot(dirs, pose_line);
        }

        float norm = sqrt(rays[i * renderParams.W * 8 + j * 8 + 4 + 0] * rays[i * renderParams.W * 8 + j * 8 + 4 + 0] + 
                        rays[i * renderParams.W * 8 + j * 8 + 4 + 1] * rays[i * renderParams.W * 8 + j * 8 + 4 + 1] + 
                        rays[i * renderParams.W * 8 + j * 8 + 4 + 2] * rays[i * renderParams.W * 8 + j * 8 + 4 + 2]);
        rays[i * renderParams.W * 8 + j * 8 + 4 + 0] = rays[i * renderParams.W * 8 + j * 8 + 4 + 0] / norm;
        rays[i * renderParams.W * 8 + j * 8 + 4 + 1] = rays[i * renderParams.W * 8 + j * 8 + 4 + 1] / norm;
        rays[i * renderParams.W * 8 + j * 8 + 4 + 2] = rays[i * renderParams.W * 8 + j * 8 + 4 + 2] / norm;
        __syncthreads();
    }
}

__device__
float3 mul_scalar(float3 vec, float scalar) {
    return make_float3(vec.x * scalar, vec.y * scalar, vec.z * scalar);
}

__device__ __forceinline__
float3 add_vec(float3 vec1, float3 vec2) {
    return make_float3(vec1.x + vec2.x, vec1.y + vec2.y, vec1.z + vec2.z);
}

__device__ __forceinline__
int3 add_int3(int3 vec1, int3 vec2) {
    return make_int3(vec1.x + vec2.x, vec1.y + vec2.y, vec1.z + vec2.z);
}

__device__ __forceinline__ 
int3 pts2coords(float3 pts, float3 origins, float voxel_size) {
    return make_int3(int(floor(pts.x - origins.x) / voxel_size), int(floor(pts.y - origins.y) / voxel_size), int(floor(pts.z - origins.z) / voxel_size));
}

__device__ __forceinline__ 
float3 pts2weights(float3 pts, float3 origins, int3 coords, float voxel_size) {
    return make_float3((pts.x - origins.x) / voxel_size - float(coords.x), (pts.y - origins.y) / voxel_size - float(coords.y), (pts.z - origins.z) / voxel_size - float(coords.z));
}

__global__
void render_rays(const float *rays, const float *volume, uint8_t *imgs, render_params renderParams, float iso_val, float threshold, float max_depth, float stride, int samples) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (i >= 0 && i < renderParams.H && j >= 0 && j < renderParams.W) {
        float3 ray_o = make_float3(rays[i * renderParams.W * 8 + j * 8 + 0], rays[i * renderParams.W * 8 + j * 8 + 1], rays[i * renderParams.W * 8 + j * 8 + 2]);
        float3 ray_d = make_float3(rays[i * renderParams.W * 8 + j * 8 + 4], rays[i * renderParams.W * 8 + j * 8 + 5], rays[i * renderParams.W * 8 + j * 8 + 6]);
        float3 origin = make_float3(renderParams.origins[0], renderParams.origins[1], renderParams.origins[2]);

        const int3 offsets[8] = {make_int3(0, 0, 0), make_int3(0, 0, 1), make_int3(0, 1, 0), make_int3(0, 1, 1),
                                 make_int3(1, 0, 0), make_int3(1, 0, 1), make_int3(1, 1, 0), make_int3(1, 1, 1)};
        // const int offsets[8][3] = {{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1},
        //                           {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}};

        float3 pts, weights;
        int3 coords, tmp_coords;
        float ray_len = 0., tsdf = 0.;
        float3 color = make_float3(0., 0., 0.), tmp_color;
        float colors[3] = {0., 0., 0.};
        int oneDim_coord;
        float color_weight;
        
        for (int k = 1; k < samples; ++k) {
            tsdf = 0.;
            ray_len += stride;
            pts.x = ray_o.x + ray_len * ray_d.x;
            pts.y = ray_o.y + ray_len * ray_d.y;
            pts.z = ray_o.z + ray_len * ray_d.z;
            // ray_d = mul_scalar(ray_d, ray_len);
            // pts = add_vec(ray_o, ray_d);

            coords.x = int(floor((pts.x - origin.x) / renderParams.voxel_size));
            coords.y = int(floor((pts.y - origin.y) / renderParams.voxel_size));
            coords.z = int(floor((pts.z - origin.z) / renderParams.voxel_size));
            // coords = pts2coords(pts, origin, renderParams.voxel_size);
            // printf ("[%d, %d]: (%f, %f, %f), (%d, %d, %d)\n", i, j, pts.x, pts.y, pts.x, coords.x, coords.y, coords.z);
            
            weights.x = (pts.x - origin.x) / renderParams.voxel_size - float(coords.x);
            weights.y = (pts.y - origin.y) / renderParams.voxel_size - float(coords.y);
            weights.z = (pts.z - origin.z) / renderParams.voxel_size - float(coords.z);
            // weights = pts2weights(pts, origin, coords, renderParams.voxel_size);
            // printf ("[%d, %d]: (%f, %f, %f), (%f, %f, %f)\n", i, j, pts.x, pts.y, pts.x, weights.x, weights.y, weights.z);

            for (int offset = 0; offset < 8; ++offset) {
                tmp_coords.x = coords.x + offsets[offset].x;
                tmp_coords.y = coords.y + offsets[offset].y;
                tmp_coords.z = coords.z + offsets[offset].z;
                // tmp_coords = add_int3(coords, offsets[offset]);
                if (tmp_coords.x >= 0 && tmp_coords.x < renderParams.dimx && 
                    tmp_coords.y >= 0 && tmp_coords.y < renderParams.dimy &&
                    tmp_coords.z >= 0 && tmp_coords.z < renderParams.dimz) {
                        tsdf += abs(float(offsets[offset].x) - weights.x) * 
                                abs(float(offsets[offset].y) - weights.y) *
                                abs(float(offsets[offset].z) - weights.z) * 
                                volume[tmp_coords.x * renderParams.dimy * renderParams.dimz * 4 + tmp_coords.y * renderParams.dimz * 4 + tmp_coords.z * 4];
                    }
            }
            // printf ("[%d, %d]: %f, %f, %f, %f\n", i, j, pts.x, pts.y, pts.x, tsdf);
            if (abs(tsdf - iso_val) <= threshold) {
                color = make_float3(0., 0., 0.);
                for (int offset = 0; offset < 8; ++offset) {
                tmp_coords.x = coords.x + offsets[offset].x;
                tmp_coords.y = coords.y + offsets[offset].y;
                tmp_coords.z = coords.z + offsets[offset].z;
                // tmp_coords = add_int3(coords, offsets[offset]);
                if (tmp_coords.x >= 0 && tmp_coords.x < renderParams.dimx && 
                    tmp_coords.y >= 0 && tmp_coords.y < renderParams.dimy &&
                    tmp_coords.z >= 0 && tmp_coords.z < renderParams.dimz) {
                        oneDim_coord = tmp_coords.x * renderParams.dimy * renderParams.dimz * 4 + tmp_coords.y * renderParams.dimz * 4 + tmp_coords.z * 4;
                        
                        color_weight = abs(float(offsets[offset].x) - weights.x) * 
                                       abs(float(offsets[offset].y) - weights.y) *
                                       abs(float(offsets[offset].z) - weights.z);
                        // printf ("[%d, %d]: (%d, %d, %d), (%d, %d, %d), %f\n", i, j, tmp_coords.x, tmp_coords.y, tmp_coords.z, coords.x, coords.y, coords.z, color_weight);
                        tmp_color = make_float3(volume[oneDim_coord + 1], 
                                                volume[oneDim_coord + 2],
                                                volume[oneDim_coord + 3]);
                        // printf ("[%d, %d]: %f, %f, %f, %f\n", i, j, tmp_color.x, tmp_color.y, tmp_color.x, color_weight);
                        // color.x = color.x + volume[oneDim_coord + 1] * color_weight;
                        // color.y = color.y + volume[oneDim_coord + 2] * color_weight;
                        // color.z = color.z + volume[oneDim_coord + 3] * color_weight;
                        tmp_color = mul_scalar(tmp_color, color_weight);
                        color = add_vec(color, tmp_color);
                    }
                }
                // color.x = volume[coords.x * renderParams.dimy * renderParams.dimz * 4 + coords.y * renderParams.dimz * 4 + coords.z * 4 + 1];
                // color.y = volume[coords.x * renderParams.dimy * renderParams.dimz * 4 + coords.y * renderParams.dimz * 4 + coords.z * 4 + 2];
                // color.z = volume[coords.x * renderParams.dimy * renderParams.dimz * 4 + coords.y * renderParams.dimz * 4 + coords.z * 4 + 3];
                // printf ("[%d, %d]: (%f, %f, %f)\n", i, j, color.z, color.y, color.x);
                imgs[i * renderParams.W * 3 + j * 3 + 0] = uint8_t(color.z);
                imgs[i * renderParams.W * 3 + j * 3 + 1] = uint8_t(color.y);
                imgs[i * renderParams.W * 3 + j * 3 + 2] = uint8_t(color.x);
                break;
                // printf ("[%d, %d]: %f, %f, %f\n", i, j, color.z, color.y, color.x);
            }
        }
        // imgs[i * renderParams.W * 3 + j * 3 + 0] = int(color.z);
        // imgs[i * renderParams.W * 3 + j * 3 + 1] = int(color.y);
        // imgs[i * renderParams.W * 3 + j * 3 + 2] = int(color.x);
        __syncthreads();
    }
}

void scale_intr(float *intr, int origin_w, int origin_h, int target_w, int target_h) {
    float scale_w = float(origin_w) / float(target_w);
    float scale_h = float(origin_h) / float(target_h);

    for (int i = 0; i < 4; ++i) {
        intr[i] = intr[i] / scale_w;
        intr[i + 4] = intr[i + 4] / scale_h;
    }
}

int main(int argc, char* argv[]) {
    float frac = 1.0;
    int BLOCKSIZE = 8;

    if (argc > 1) {
        frac = atof(argv[1]);
    }
    if (argc > 2) {
        BLOCKSIZE = atoi(argv[2]);
    }

    const char pose_path[] = "./data_bin/pose.bin";
    const char intr_path[] = "./data_bin/intrinsic_color.bin";
    const char info_path[] = "./data_bin/render_info.bin";
    const char vol_path[]  = "./data_bin/volume.bin";
    const int origin_w = 1296;
    const int origin_h = 968;
    const int target_w = (int)(frac * origin_w);
    const int target_h = (int)(frac * origin_h);
    const int samples = 80;
    const float max_depth = 3.;
    const float stride = max_depth / float(samples);
    const float iso_val = 0.;
    const float threshold = 0.5;

    float *info;
    cudaMallocManaged((void **)&info, 9 * sizeof(float));
    int bytes = readBin(info_path, info, 9);
    // std::cout << bytes << " bytes of render info have been read" << std::endl;

    const int dims[3] = {(int)info[4], (int)info[5], (int)info[6]};
    const float origins[3] = {info[0], info[1], info[2]};
    const float voxel_size = info[3];
    const int img_size[2] = {(int)info[7], (int)info[8]};

    float *pose;
    cudaMallocManaged((void **)&pose, 16 * sizeof(float));
    bytes = readBin(pose_path, pose, 16);
    // std::cout << bytes << " bytes of pose have been read" << std::endl;

    float *intr;
    cudaMallocManaged((void **)&intr, 16 * sizeof(float));
    bytes = readBin(intr_path, intr, 16);
    // std::cout << bytes << " bytes of intrinsic have been read" << std::endl;
    scale_intr(intr, origin_w, origin_h, target_w, target_h);

    float *vol;
    cudaMallocManaged((void **)&vol, dims[0] * dims[1] * dims[2] * 4 * sizeof(float));
    bytes = readBin(vol_path, vol, dims[0] * dims[1] * dims[2] * 4);
    // std::cout << bytes << " bytes of volume have been read" << std::endl;
    
    float *rays;
    cudaMallocManaged((void **)&rays, target_h * target_w * 8 * sizeof(float));

    uint8_t *img, *img_cpu = new uint8_t[target_h * target_w * 3];
    cudaMallocManaged((void **)&img, target_h * target_w * 3 * sizeof(uint8_t));

    int device = 0;
    cudaMemPrefetchAsync((void *)info, 9 * sizeof(float), device, NULL);
    cudaMemPrefetchAsync((void *)pose, 16 * sizeof(float), device, NULL);
    cudaMemPrefetchAsync((void *)intr, 16 * sizeof(float), device, NULL);
    cudaMemPrefetchAsync((void *)vol, dims[0] * dims[1] * dims[2] * 4 * sizeof(float), device, NULL);
    cudaMemPrefetchAsync((void *)rays, target_h * target_w * 8 * sizeof(float), device, NULL);
    cudaMemPrefetchAsync((void *)img, target_h * target_w * 3 * sizeof(uint8_t), device, NULL);

    render_params renderParams;
    renderParams.cx = intr[2];
    renderParams.cy = intr[6];
    renderParams.fx = intr[0];
    renderParams.fy = intr[5];
    renderParams.H = target_h;
    renderParams.W = target_w;
    renderParams.dimx = dims[0];
    renderParams.dimy = dims[1];
    renderParams.dimz = dims[2];
    renderParams.origins[0] = origins[0];
    renderParams.origins[1] = origins[1];
    renderParams.origins[2] = origins[2];
    renderParams.voxel_size = voxel_size;

    const dim3 gridShape((target_h + BLOCKSIZE - 1) / BLOCKSIZE, (target_w + BLOCKSIZE - 1) / BLOCKSIZE);
    const dim3 blockShape(BLOCKSIZE, BLOCKSIZE);

    struct timespec sts, ets;
    timespec_get(&sts, TIME_UTC);

    generate_ray_per_pixel<<<gridShape, blockShape>>>(rays, pose, renderParams);
    cudaDeviceSynchronize();

    timespec_get(&ets, TIME_UTC);
    time_t dsec = ets.tv_sec - sts.tv_sec;
    long dnsec = ets.tv_nsec - sts.tv_nsec;
    printf ("image size [%d, %d], block size %d, generate rays: %ld.%09lds\n", target_h, target_w, BLOCKSIZE, dsec, dnsec);

    timespec_get(&sts, TIME_UTC);
    
    render_rays<<<gridShape, blockShape>>>(rays, vol, img, renderParams, iso_val, threshold, max_depth, stride, samples);
    cudaDeviceSynchronize();

    timespec_get(&ets, TIME_UTC);
    dsec = ets.tv_sec - sts.tv_sec;
    dnsec = ets.tv_nsec - sts.tv_nsec;
    printf ("image size: [%d, %d], block size %d, render rays: %ld.%09lds\n", target_h, target_w, BLOCKSIZE, dsec, dnsec);

    FILE *fp = fopen("./test_cuda.jpeg", "wb");
    svjpeg(fp, target_w, target_h, img);

    cudaFree(info);
    cudaFree(pose);
    cudaFree(intr);
    cudaFree(vol);
    cudaFree(rays);
    cudaFree(img);
    return 0;
}