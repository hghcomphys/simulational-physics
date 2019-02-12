#!/bin/bash

echo 'It may take a while...'
for i in gpu omp mpi 
do
	echo $i

	if [ $i==mpi ]
	then
		mpirun -np 8 lmp_mpi < $i-spce.in > $i-spce.out
	fi

	if [ $i==gpu ]
	then
		lmp_serial -sf gpu -pk gpu 1  < $i-spce.in > $i-spce.out
	fi	
	
	if [ $i==omp ]
	then
		mpirun -np 4 lmp_mpi -sf omp -pk omp 2 -in $i-spce.in > $i-spce.out
    fi
	
	cd out
	mv ../rdf.dat rdf-$i.dat
	mv ../dump.xyz dump-$i.xyz
	mv ../$i-spce.out .
	cd ../	
done
rm -f log.lammps log.cite
echo 'done.'
echo 
