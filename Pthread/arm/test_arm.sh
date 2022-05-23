qemu-aarch64 -L /usr/aarch64-linux-gnu/ ./build/generate_rays_pthread_neon_arm 1.0 8
qemu-aarch64 -L /usr/aarch64-linux-gnu/ ./build/generate_rays_pthread_neon_arm 1.0 4 
qemu-aarch64 -L /usr/aarch64-linux-gnu/ ./build/generate_rays_pthread_neon_arm 1.0 2 
qemu-aarch64 -L /usr/aarch64-linux-gnu/ ./build/generate_rays_pthread_neon_arm 0.5 8 
qemu-aarch64 -L /usr/aarch64-linux-gnu/ ./build/generate_rays_pthread_neon_arm 0.25 8 
qemu-aarch64 -L /usr/aarch64-linux-gnu/ ./build/generate_rays_pthread_neon_arm 0.125 8 
qemu-aarch64 -L /usr/aarch64-linux-gnu/ ./build/generate_rays_pthread_neon_arm 1.0 8 col
qemu-aarch64 -L /usr/aarch64-linux-gnu/ ./build/generate_rays_pthread_neon_arm 1.0 8 oneDim

qemu-aarch64 -L /usr/aarch64-linux-gnu/ ./build/generate_rays_pthread_arm 1.0 8
qemu-aarch64 -L /usr/aarch64-linux-gnu/ ./build/generate_rays_pthread_arm 1.0 4 
qemu-aarch64 -L /usr/aarch64-linux-gnu/ ./build/generate_rays_pthread_arm 1.0 2 

qemu-aarch64 -L /usr/aarch64-linux-gnu/ ./build/generate_rays_pthread_arm 0.5 8
qemu-aarch64 -L /usr/aarch64-linux-gnu/ ./build/generate_rays_pthread_arm 0.25 8 
qemu-aarch64 -L /usr/aarch64-linux-gnu/ ./build/generate_rays_pthread_arm 0.125 8 

qemu-aarch64 -L /usr/aarch64-linux-gnu/ ./build/generate_rays_serial_arm 1.0 
qemu-aarch64 -L /usr/aarch64-linux-gnu/ ./build/generate_rays_serial_arm 0.5 
qemu-aarch64 -L /usr/aarch64-linux-gnu/ ./build/generate_rays_serial_arm 0.25 
qemu-aarch64 -L /usr/aarch64-linux-gnu/ ./build/generate_rays_serial_arm 0.125 

qemu-aarch64 -L /usr/aarch64-linux-gnu/ ./build/generate_rays_neon_arm 1.0 
qemu-aarch64 -L /usr/aarch64-linux-gnu/ ./build/generate_rays_neon_arm 0.5 
qemu-aarch64 -L /usr/aarch64-linux-gnu/ ./build/generate_rays_neon_arm 0.25 
qemu-aarch64 -L /usr/aarch64-linux-gnu/ ./build/generate_rays_neon_arm 0.125 