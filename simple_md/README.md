### Simple MD
It separately contains C and CUDA C++ programs in order to perform simple MD simulations for liquid Argon Lennard-Jones (12-6) in 2D using fast Verlet algorithm in NVE ensemble. Finally, a [MATLAB](https://nl.mathworks.com/products/matlab.html) script visualizes the output trajectory.

How to run (serial):
```  
$ gcc md_nve.c -lm -O2 
$ octave movie.m
```  
MPI:
```
$ mpicc mpi_md_nve.c -lm -O2 -o md_mpi.out
$ mpirun -np 2 md_mpi.out 
$ rm -f md_mpi.out
````
On GPU:
```
$ nvcc -O2 cuda_md_nve.cu -o md_gpu.out
$ ./md_gpu.out
```