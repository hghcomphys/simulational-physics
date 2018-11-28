#!/bin/bash

echo 'It may take a while...'
for i in spce tip3p tip4p tip5p #reax
do
	echo $i
	mpirun -np 4 lmp_mpi < $i.in > $i.out
	
	cd out
	mv ../rdf.dat rdf-$i.dat
	mv ../dump.xyz dump-$i.xyz
	mv ../$i.out .
	cd ../	
done
rm -f log.lammps log.cite
echo 'done.'
echo 
