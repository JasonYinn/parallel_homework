aarch64-linux-gnu-g++ -O2 -o generate_rays_serial -mcpu=cortex-a57 generate_rays_serial.cpp

aarch64-linux-gnu-g++ -O2 -o generate_rays_neon -mcpu=cortex-a57 generate_rays_neon.cpp

g++ -O2 -march=native -o generate_rays_sse -march=native generate_rays_sse.cpp