# Molecular Simulations
This project contains a list of simulations (separate folders) each related to the topic of computational physics in particular molecular simulation.


---
### Water models  
It contains several [LAMMPS](https://lammps.sandia.gov/) scripts concern molecular dynamics (MD) simulations of liquid water at ambient conditions using different [water models](http://www1.lsbu.ac.uk/water/water_models.html) SPC/E, TIP3P, TIP4P, TIP4P, and ReaxFF.
Sample simulations outputs including mass density, radial distribution, and trajectory files (*xyz* format) located in the "out" folder.  

How to run (SPC/E):
```
lmp_serial < spce.in > spce.out 
```
MPI:
```
mpirun -np 4 lmp_mpi < spce.in > spce.out 
```
---
### LAMMPS HPC
It contains [LAMMPS](https://lammps.sandia.gov/) scripts and carries out MD simulations of liquid water (SPC/E) on CPU and GPU in order to compare their computational performance. The CPU version includes parallel MPI and hybrid (MPI+OMP) simulations. Sample outputs presented in "out" folder.

How to run (MPI):
```
mpirun -np 8 lmp_mpi < mpi_spce.in 
```
Hybrid (MPI+OMP):
```
mpirun -np 4 lmp_omp -sf omp -pk omp 2 -in omp_spce.in  
```
On GPU:
```
lmp_gpu -sf gpu -pk gpu 1 -in gpu_spce.in 
```
---
### Simple MD
It separately contains C and CUDA C++ programs in order to perform simple MD simulations for liquid Argon Lennard-Jones (12-6) in 2D using fast Verlet algorithm in NVE ensemble. Finally, a [MATLAB](https://nl.mathworks.com/products/matlab.html) script visualizes the output trajectory.

How to run (serial):
```  
gcc md_nve.c -lm -O2 
octave movie.m
```  
MPI:
```
mpicc mpi_md_nve.c -lm -O2 -o md_mpi.out
mpirun -np 2 md_mpi.out 
rm -f md_mpi.out
````
On GPU:
```
nvcc -O2 cuda_md_nve.cu -o md_gpu.out
./md_gpu.out
```
---
### Ising model
MATLAB scripts for Monte Carlo simulation of a 2D ising model in order to study the critical behavior of the system as a function of temperature.
  
How to run:
```
octave ising.m
```
