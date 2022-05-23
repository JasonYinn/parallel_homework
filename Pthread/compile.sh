# g++ -O2 -march=native -o generate_rays_serial_x86 generate_rays_serial.cpp

# g++ -O2 -march=native -pthread -o generate_rays_pthread_x86 generate_rays_pthread.cpp

aarch64-linux-gnu-g++ -O2 -pthread -o generate_rays_pthread_neon -mcpu=cortex-a57 generate_rays_pthread_neon.cpp