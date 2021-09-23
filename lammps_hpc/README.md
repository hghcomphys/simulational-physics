### LAMMPS HPC
It contains [LAMMPS](https://lammps.sandia.gov/) scripts and carries out MD simulations of liquid water (SPC/E) on CPU and GPU in order to compare their computational performance. The CPU version includes parallel MPI and hybrid (MPI+OMP) simulations. Sample outputs presented in "out" folder.

How to run (MPI):
```
$ mpirun -np 8 lmp_mpi < mpi_spce.in 
```
Hybrid (MPI+OMP):
```
$ mpirun -np 4 lmp_omp -sf omp -pk omp 2 -in omp_spce.in  
```
On GPU:
```
$ lmp_gpu -sf gpu -pk gpu 1 -in gpu_spce.in 
```