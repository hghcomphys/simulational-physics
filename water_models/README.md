### Water models  
It contains several [LAMMPS](https://lammps.sandia.gov/) scripts concern molecular dynamics (MD) simulations of liquid water at ambient conditions using different [water models](http://www1.lsbu.ac.uk/water/water_models.html) SPC/E, TIP3P, TIP4P, TIP4P, and ReaxFF.
Sample simulations outputs including mass density, radial distribution, and trajectory files (*xyz* format) located in the "out" folder.  

How to run (SPC/E):
```
$ lmp_serial < spce.in > spce.out 
```
MPI:
```
$ mpirun -np 4 lmp_mpi < spce.in > spce.out 
```