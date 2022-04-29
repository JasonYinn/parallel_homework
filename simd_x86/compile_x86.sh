g++ -O2 -march=native -o generate_rays_s1_x86 generate_rays_s1.cpp
g++ -O2 -march=native -o generate_rays_s2_x86 generate_rays_s2.cpp
g++ -O2 -march=native -o generate_rays_s4_x86 generate_rays_s4.cpp
g++ -O2 -march=native -o generate_rays_s8_x86 generate_rays_s8.cpp

g++ -O2 -march=native -o generate_rays_ss1_x86 -march=native generate_rays_ss1.cpp
g++ -O2 -march=native -o generate_rays_ss2_x86 -march=native generate_rays_ss2.cpp
g++ -O2 -march=native -o generate_rays_ss4_x86 -march=native generate_rays_ss4.cpp
g++ -O2 -march=native -o generate_rays_ss8_x86 -march=native generate_rays_ss8.cpp
# aarch64-linux-gnu-g++ -O2 -o generate_rays_sse -march=native -mcpu=cortex-a57 generate_rays_sse.cpp