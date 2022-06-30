rm -r build
mkdir build

g++ -O2 -march=native -o ./build/generate_rays_serial_x86 generate_rays_serial.cpp

g++ -O2 -march=native -o ./build/generate_rays_sse_x86 generate_rays_sse.cpp

g++ -O2 -march=native -pthread -o ./build/generate_rays_pthread_x86 generate_rays_pthread.cpp

g++ -O2 -march=native -pthread -o ./build/generate_rays_pthread_sse_x86 generate_rays_pthread_sse.cpp

g++ render_serial.cpp -o build/render_serial

nvcc render_cuda.cu -o build/render_cuda