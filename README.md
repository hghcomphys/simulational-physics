MD on CPU
---------
This simple C code perfoms Lennard-Jones (12-16) liquid Argon molecular dynamics in 2D by a fast Verlet algorithm.

How to run:
gcc md_nve.c -lm -O2 && octave movie.m  

MPI version:
mpicc mpi_md_nve.c -lm -O2 -o md_mpi.out && mpirun -np 2 md_mpi.out && rm -f md_mpi.out



MD on GPU
---------
This simple CUDA C++ code similarly simulates argon LJ (12-6) in 2D using CUDA thechnology on cuda enabled nvidia GPU.

How to run:
nvcc -O2 cuda_md_nve.cu -o md_gpu.out && ./md_gpu.out
