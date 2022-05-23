rm -r build
mkdir build

aarch64-linux-gnu-g++ -O2 -pthread -o ./build/generate_rays_pthread_neon_arm -mcpu=cortex-a73 generate_rays_pthread_neon.cpp

aarch64-linux-gnu-g++ -O2 -o ./build/generate_rays_serial_arm -mcpu=cortex-a73 generate_rays_serial.cpp

aarch64-linux-gnu-g++ -O2 -o ./build/generate_rays_neon_arm -mcpu=cortex-a73 generate_rays_neon.cpp

aarch64-linux-gnu-g++ -O2 -pthread -o ./build/generate_rays_pthread_arm -mcpu=cortex-a73 generate_rays_pthread.cpp
