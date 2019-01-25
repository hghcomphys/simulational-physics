### MD on CPU
This simple C program perfoms Lennard-Jones (12-16) liquid Argon molecular dynamics in 2D by a fast Verlet algorithm.

How to run:
g++ -O2 md_nve.c && octave movie.m  


### MD on GPU
This simple CUDA C++ program similarly simulates argon LJ (12-6) in 2D using CUDA thechnology on cuda enabled nvidia GPU.

How to run:
nvcc -O2 md_nve_gpu.cu -o md_gpu.out && ./md_gpu.out
