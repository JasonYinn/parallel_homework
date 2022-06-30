#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include "svjpeg.hpp"

// using namespace std;

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
const int offsets[8][3] = {{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1},
                           {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}};

template <typename scalar_t>
int readBin(const char* fname, scalar_t *arr, int len) {
    std::ifstream in(fname, std::ios::in | std::ios::binary);

    in.read((char *) arr, sizeof(scalar_t) * len);

    return in.gcount();
}

void generate_ray_per_pixel(float *rays, float *c2w, render_params renderParams, int i, int j) {
    float dirs[4] = {((float)j - renderParams.cx) / renderParams.fx, -((float)i - renderParams.cy) / renderParams.fy, -1.0, 0.0};
    for (int k = 0; k < 4; ++k) {
        float temp = 0.;
        rays[i * renderParams.W * 8 + j * 8 + k] = c2w[k * 4 + 3];
        for (int z = 0; z < 4; ++z) {
            temp += dirs[z] * c2w[k * 4 + z];
        }
        rays[i * renderParams.W * 8 + j * 8 + 4 + k] = temp;
    }

    float norm = sqrt(rays[i * renderParams.W * 8 + j * 8 + 4 + 0] * rays[i * renderParams.W * 8 + j * 8 + 4 + 0] + 
                      rays[i * renderParams.W * 8 + j * 8 + 4 + 1] * rays[i * renderParams.W * 8 + j * 8 + 4 + 1] + 
                      rays[i * renderParams.W * 8 + j * 8 + 4 + 2] * rays[i * renderParams.W * 8 + j * 8 + 4 + 2]);
    rays[i * renderParams.W * 8 + j * 8 + 4 + 0] = rays[i * renderParams.W * 8 + j * 8 + 4 + 0] / norm;
    rays[i * renderParams.W * 8 + j * 8 + 4 + 1] = rays[i * renderParams.W * 8 + j * 8 + 4 + 1] / norm;
    rays[i * renderParams.W * 8 + j * 8 + 4 + 2] = rays[i * renderParams.W * 8 + j * 8 + 4 + 2] / norm;

    // printf ("fx: %f, fy: %f, cx: %f, cy: %f\n", renderParams.fx, renderParams.fy, renderParams.cx, renderParams.cy);
    // printf ("[%d, %d]: ray_o (%f, %f, %f), ray_d (%f, %f, %f)\n", i, j, rays[i * renderParams.W * 8 + j * 8 + 0], rays[i * renderParams.W * 8 + j * 8 + 1], rays[i * renderParams.W * 8 + j * 8 + 2], 
    //                                                                     rays[i * renderParams.W * 8 + j * 8 + 4], rays[i * renderParams.W * 8 + j * 8 + 5], rays[i * renderParams.W * 8 + j * 8 + 6]);
}

bool render_rays(const float *rays, const float *volume, uint8_t *imgs, render_params renderParams, float iso_val, float threshold, int i, int j, float max_depth, float stride, int samples) {
    float ray_o[3] = {rays[i * renderParams.W * 8 + j * 8 + 0], rays[i * renderParams.W * 8 + j * 8 + 1], rays[i * renderParams.W * 8 + j * 8 + 2]};
    float ray_d[3] = {rays[i * renderParams.W * 8 + j * 8 + 4], rays[i * renderParams.W * 8 + j * 8 + 5], rays[i * renderParams.W * 8 + j * 8 + 6]};
    float origin[3] = {renderParams.origins[0], renderParams.origins[1], renderParams.origins[2]};

    float pts[3], weights[3];
    int coords[3], tmp_coords[3];
    float ray_len = 0., tsdf;
    float color[3] = {0., 0., 0.};
    int oneDim_coord;
    float color_weight;
    
    for (int k = 1; k < samples; ++k) {
        tsdf = 0.;
        ray_len += stride;
        pts[0] = ray_o[0] + ray_len * ray_d[0];
        pts[1] = ray_o[1] + ray_len * ray_d[1];
        pts[2] = ray_o[2] + ray_len * ray_d[2];

        coords[0] = int(floor((pts[0] - origin[0]) / renderParams.voxel_size));
        coords[1] = int(floor((pts[1] - origin[1]) / renderParams.voxel_size));
        coords[2] = int(floor((pts[2] - origin[2]) / renderParams.voxel_size));
        
        weights[0] = (pts[0] - origin[0]) / renderParams.voxel_size - float(coords[0]);
        weights[1] = (pts[1] - origin[1]) / renderParams.voxel_size - float(coords[1]);
        weights[2] = (pts[2] - origin[2]) / renderParams.voxel_size - float(coords[2]);

        for (int offset = 0; offset < 8; ++offset) {
            tmp_coords[0] = coords[0] + offsets[offset][0];
            tmp_coords[1] = coords[1] + offsets[offset][1];
            tmp_coords[2] = coords[2] + offsets[offset][2];
            if (tmp_coords[0] >= 0 && tmp_coords[0] < renderParams.dimx && 
                tmp_coords[1] >= 0 && tmp_coords[1] < renderParams.dimy &&
                tmp_coords[2] >= 0 && tmp_coords[2] < renderParams.dimz) {
                    tsdf += abs(float(offsets[offset][0]) - weights[0]) * 
                            abs(float(offsets[offset][1]) - weights[1]) *
                            abs(float(offsets[offset][2]) - weights[2]) * 
                            volume[tmp_coords[0] * renderParams.dimy * renderParams.dimz * 4 + tmp_coords[1] * renderParams.dimz * 4 + tmp_coords[2] * 4];
                }
        }
        if (abs(tsdf - iso_val) <= threshold) {
            for (int offset = 0; offset < 8; ++offset) {
            tmp_coords[0] = coords[0] + offsets[offset][0];
            tmp_coords[1] = coords[1] + offsets[offset][1];
            tmp_coords[2] = coords[2] + offsets[offset][2];
            if (tmp_coords[0] >= 0 && tmp_coords[0] < renderParams.dimx && 
                tmp_coords[1] >= 0 && tmp_coords[1] < renderParams.dimy &&
                tmp_coords[2] >= 0 && tmp_coords[2] < renderParams.dimz) {
                    color_weight = abs(float(offsets[offset][0]) - weights[0]) * 
                                abs(float(offsets[offset][1]) - weights[1]) *
                                abs(float(offsets[offset][2]) - weights[2]);
                    color[0] += color_weight *
                                volume[tmp_coords[0] * renderParams.dimy * renderParams.dimz * 4 + tmp_coords[1] * renderParams.dimz * 4 + tmp_coords[2] * 4 + 1];
                    color[1] += color_weight * 
                                volume[tmp_coords[0] * renderParams.dimy * renderParams.dimz * 4 + tmp_coords[1] * renderParams.dimz * 4 + tmp_coords[2] * 4 + 2];
                    color[2] += color_weight * 
                                volume[tmp_coords[0] * renderParams.dimy * renderParams.dimz * 4 + tmp_coords[1] * renderParams.dimz * 4 + tmp_coords[2] * 4 + 3];
                }
            }
            imgs[i * renderParams.W * 3 + j * 3 + 0] = int(color[2]);
            imgs[i * renderParams.W * 3 + j * 3 + 1] = int(color[1]);
            imgs[i * renderParams.W * 3 + j * 3 + 2] = int(color[0]);
            return true;
        }
    }
    imgs[i * renderParams.W * 3 + j * 3 + 0] = int(color[0]);
    imgs[i * renderParams.W * 3 + j * 3 + 1] = int(color[1]);
    imgs[i * renderParams.W * 3 + j * 3 + 2] = int(color[2]);
    // printf ("[%d, %d]: (%d, %d, %d)\n", i, j, imgs[i * renderParams.W * 3 + j * 3 + 0], imgs[i * renderParams.W * 3 + j * 3 + 1], imgs[i * renderParams.W * 3 + j * 3 + 2]);
    // printf ("[%d, %d]: (%f, %f, %f)\n", i, j, color[0], color[1], color[2]);
    return false;
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

    float *info = new float[9];
    int bytes = readBin(info_path, info, 9);
    // std::cout << bytes << " bytes of render info have been read" << std::endl;

    const int dims[3] = {(int)info[4], (int)info[5], (int)info[6]};
    const float origins[3] = {info[0], info[1], info[2]};
    const float voxel_size = info[3];
    const int img_size[2] = {(int)info[7], (int)info[8]};

    float *pose = new float[16];
    bytes = readBin(pose_path, pose, 16);
    // std::cout << bytes << " bytes of pose have been read" << std::endl;

    float *intr = new float[16];
    bytes = readBin(intr_path, intr, 16);
    // std::cout << bytes << " bytes of intrinsic have been read" << std::endl;
    scale_intr(intr, origin_w, origin_h, target_w, target_h);

    float *vol = new float[dims[0] * dims[1] * dims[2] * 4];
    bytes = readBin(vol_path, vol, dims[0] * dims[1] * dims[2] * 4);
    // std::cout << bytes << " bytes of volume have been read" << std::endl;
    
    float *rays = new float[target_h * target_w * 8];

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

    struct timespec sts, ets;
    timespec_get(&sts, TIME_UTC);
    for (int i = 0; i < target_h; ++i) {
        for (int j = 0; j < target_w; ++j) {
            generate_ray_per_pixel(rays, pose, renderParams, i, j);
        }
    }
    timespec_get(&ets, TIME_UTC);
    time_t dsec = ets.tv_sec - sts.tv_sec;
    long dnsec = ets.tv_nsec - sts.tv_nsec;
    printf ("image size [%d, %d], generate rays: %ld.%09lds\n", target_h, target_w, dsec, dnsec);

    timespec_get(&sts, TIME_UTC);
    uint8_t *img = new uint8_t[target_h * target_w * 3];
    int cnt = 0;
    bool flag;
    for (int i = 0; i < target_h; ++i) {
        for (int j = 0; j < target_w; ++j) {
            flag = render_rays(rays, vol, img, renderParams, iso_val, threshold, i, j, max_depth, stride, samples);
            cnt += int(flag);
        }
    }

    timespec_get(&ets, TIME_UTC);
    dsec = ets.tv_sec - sts.tv_sec;
    dnsec = ets.tv_nsec - sts.tv_nsec;
    printf ("image size [%d, %d], render rays: %ld.%09lds\n", target_h, target_w, dsec, dnsec);

    FILE *fp = fopen("./test_serial.jpeg", "wb");
    svjpeg(fp, target_w, target_h, img);
    // WriteBMP((char *)img, "./test.bmp", target_w, target_h);

    free(info);
    free(pose);
    free(intr);
    free(vol);
    free(rays);
    free(img);
    return 0;
}