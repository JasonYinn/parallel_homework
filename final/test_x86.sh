export CUDA_VISIBLE_DEVICES=1

./build/generate_rays_pthread_sse_x86 1.0 8
./build/generate_rays_pthread_sse_x86 1.0 4 
./build/generate_rays_pthread_sse_x86 1.0 2 
./build/generate_rays_pthread_sse_x86 0.75 8 
./build/generate_rays_pthread_sse_x86 0.5 8 
./build/generate_rays_pthread_sse_x86 0.25 8 
./build/generate_rays_pthread_sse_x86 1.0 8 col 
./build/generate_rays_pthread_sse_x86 1.0 8 oneDim

./build/generate_rays_pthread_x86 1.0 8
./build/generate_rays_pthread_x86 1.0 4 
./build/generate_rays_pthread_x86 1.0 2 

./build/generate_rays_pthread_x86 0.75 8
./build/generate_rays_pthread_x86 0.50 8
./build/generate_rays_pthread_x86 0.25 8 

./build/generate_rays_serial_x86 1.0 
./build/generate_rays_serial_x86 0.75
./build/generate_rays_serial_x86 0.5 
./build/generate_rays_serial_x86 0.25 

./build/generate_rays_sse_x86 1.0 
./build/generate_rays_sse_x86 0.75
./build/generate_rays_sse_x86 0.5 
./build/generate_rays_sse_x86 0.25 

build/render_serial 0.25
build/render_serial 0.5
build/render_serial 0.75
build/render_serial 1.0

build/render_cuda 1.0 8
build/render_cuda 0.75 8
build/render_cuda 0.5 8
build/render_cuda 0.25 8

build/render_cuda 1.0 2
build/render_cuda 1.0 4
build/render_cuda 1.0 16
build/render_cuda 1.0 32