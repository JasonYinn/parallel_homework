g++ -O2 -march=native -o generate_rays_serial_x86 generate_rays_serial.cpp

g++ -O2 -march=native -o generate_rays_neon_x86 generate_rays_neon.cpp

g++ -O2 -march=native -o generate_rays_sse_x86 -march=native generate_rays_sse.cpp
# aarch64-linux-gnu-g++ -O2 -o generate_rays_sse -march=native -mcpu=cortex-a57 generate_rays_sse.cpp