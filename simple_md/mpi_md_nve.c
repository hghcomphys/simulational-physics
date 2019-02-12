
/*  ---------------------------------------------------------------
	project name  : Molecular dynamics simulation (NVE ensemble)
	Written by    : Hosein Ghorbanfekr hgh.comphys@gmail.com
	Creation Date : July 2011
	Description   : This program simulates molecular dynamics
                    system in 2D using the fast-Verlet algorithem.
    Implementation: MPI                
    ---------------------------------------------------------------  */ 

// Includes

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mpi.h"

// Defs

#define N		100 // number of particles
#define L		(sqrt(N)*LatDis) // simulation length
#define h		.0001f // time step
#define LatDis  2.5f // lattice distance
#define vmax	10.f // maximum random velocity
#define MAX_N	500 // maximum number of particles
#define RAND	(2.f*((float)rand()/RAND_MAX-.5f))
#define root	0 // root processor
#define N_div	(N/procDim) // N for each processor 

// Functions

int  stop(int);
int  start(int);
void accel(void);
void setVCMzero(void);
void ShareCoordinates(void);
void InitializeParticles(void);
 
// Global variables

int   procIdx;
int   procDim;
float r[MAX_N][2];
float v[MAX_N][2];
float f[MAX_N][2];

// Implementation of functions

int main(int argc, char *argv[])
{
	int t;
	int i, d;
	FILE *fp;
	clock_t time0=clock();
	// Initialize MPI 
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &procIdx);
	MPI_Comm_size(MPI_COMM_WORLD, &procDim);
	// open output file on root processor
	if(procIdx==root) 
		fp = fopen("configuration.xyz","w");
	// initialize particles 
	InitializeParticles();
	accel(); 
	for(t=0;t<=1000;t++) // the main time loop
	{
		// update coordinates and velocites
		for(i=start(N);i<stop(N);i++) 
			for(d=0;d<2;d++)
			{
				r[i][d] = fmod( r[i][d] + v[i%N_div][d]*h + f[i%N_div][d]*h*h*.5f + L, L);
				v[i%N_div][d]+= f[i%N_div][d]*h*.5f;
			}
		// update velocities with new forces
		ShareCoordinates();
		accel();
		for(i=0;i<N_div;i++)
			for(d=0;d<2;d++)
				v[i][d]+= f[i][d]*h*.5f;			
		// save data	
		if(t%100==0 && procIdx==root) 
		{
			printf(" t --> %d\n", t);
			fprintf(fp, "%d\n\n", N);
			for(i=0;i<N;i++)
				fprintf(fp, "Ar\t%f\t%f\t%f\n", r[i][0], r[i][1], 0.f);
			fflush(fp);
		}
	} 
	if(procIdx==root) 
	{
		fclose(fp);
		printf("Finished after %.3f seconds.\n", (double)(clock() - time0) / CLOCKS_PER_SEC);
	}
	MPI_Finalize();
	return 0;
}

void accel(void)
{
	int	i, j;
	float tmpx, tmpy;
	float rij, Fijx, Fijy;
	for(i=0;i<N_div;i++) 
	{
		f[i][0] = 0.f;
		f[i][1] = 0.f;
	}
	for(i=start(N);i<stop(N);i++) 
		for(j=0;j<N;j++) 
			if(j!=i)
			{
				// PBC
				tmpx = (r[j][0]-r[i][0]);
				if(tmpx>=L/2.f)
					tmpx-=L;
				else if(tmpx<=-L/2.f)
					tmpx+=L;

				tmpy = (r[j][1]-r[i][1]);
				if(tmpy>=L/2.f)
					tmpy-=L;
				else if(tmpy<=-L/2.f)
					tmpy+=L;
				// LJ force field
				rij = sqrt(tmpx*tmpx + tmpy*tmpy);
				rij *=rij;
				Fijx = tmpx/rij;
				Fijy = tmpy/rij;
				rij *= rij*rij;
				Fijx*=(6.f/rij-12.f/rij/rij);
				Fijy*=(6.f/rij-12.f/rij/rij);
				// update forces
				f[i%N_div][0]+=Fijx;
				f[i%N_div][1]+=Fijy;
			}
}

void InitializeParticles(void)
{
	int i, j, k;
	// initialize coordinates
	for(j=start(sqrt(N)), k=start(N);j<stop(sqrt(N));j++) 
		for(i=0;i<sqrt(N);i++, k++) 
		{
			r[k][0] = i*LatDis;
			r[k][1] = j*LatDis;
		} 	
	ShareCoordinates();
	// initialize velocities
	for(i=0;i<N_div;i++) 
	{
		v[i][0] = RAND*vmax;
		v[i][1] = RAND*vmax;
	}
	setVCMzero();
}

void setVCMzero(void)
{   
	int i;
	float sumx=0.f,sumy=0.f;
	for(i=0;i<N_div;i++) {
		sumx+= v[i][0];
		sumy+= v[i][1];
	}
	for(i=0;i<N_div;i++) {
		v[i][0]-=sumx/N_div;
		v[i][1]-=sumy/N_div;
	}
}

void ShareCoordinates(void)
{
	int proc;
	MPI_Barrier(MPI_COMM_WORLD);
	// share coordinates between processors
	for(proc=0;proc<procDim;proc++)
		MPI_Bcast(r[proc*N_div], 2*N_div, MPI_FLOAT, proc, MPI_COMM_WORLD);
}

int start(int X) { return (procIdx*X/procDim); }

int stop(int X) { return ((procIdx+1)*X/procDim); }






