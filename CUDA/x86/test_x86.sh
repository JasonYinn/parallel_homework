export CUDA_VISIBLE_DEVICES=1

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