MD on CPU
---------
This simple C code perfoms Lennard-Jones (12-16) liquid Argon molecular dynamics in 2D by a fast Verlet algorithm.

How to run:
g++ md_nve.cpp -O2 && octave movie.m  


MD on GPU
---------
This simple CUDA C++ code similarly simulates argon LJ (12-6) in 2D using CUDA thechnology on cuda enabled nvidia GPU.

How to run:
nvcc -O2 md_nve_gpu.cu -o md_gpu.out && ./md_gpu.out
