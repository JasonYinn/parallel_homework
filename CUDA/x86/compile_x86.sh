rm -r build
mkdir build

g++ render_serial.cpp -o build/render_serial

nvcc render_cuda.cu -o build/render_cuda