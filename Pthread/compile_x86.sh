# g++ -O2 -march=native -o generate_rays_serial_x86 generate_rays_serial.cpp

# g++ -O2 -march=native -pthread -o generate_rays_pthread_x86 generate_rays_pthread.cpp

g++ -O2 -march=native -pthread -o generate_rays_pthread_sse generate_rays_pthread_sse.cpp