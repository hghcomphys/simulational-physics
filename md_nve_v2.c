
/*  ---------------------------------------------------------------
	project name  : Molecular dynamics simulation (NVE ensemble)
	Written by    : Hosein Ghorbanfekr hgh.comphys@gmail.com
	Creation Date : July 2011
	Description   : This program simulates molecular dynamics
                    system in 2D using the fast-Verlet algorithem.
    ---------------------------------------------------------------  */ 

#include <stdio.h>
#include <stdlib.h> 
#include <malloc.h>
#include <math.h>
#include <time.h>

// ------------------

#define N			100 // number of particles
#define d       	2.5f // lattice distance
#define L			(sqrtf(N)*d) // simulation box length
#define h			.0001f // time step
#define vmax		10.f // maximum random velocity
#define RAND		(2.f*((float)rand()/RAND_MAX-.5f))

// ------------------

typedef 
struct	{ 
	float x[2];
	float v[2];
	float a[2];
} Particle;

// ------------------

/*setting center of mass velocity to zero*/
void setVCMzero(Particle *P)
{   
	int		i;
	float	sumx=0, sumy=0;
	for(i=0;i<N;i++) {
		sumx+= (P+i)->v[0];
		sumy+= (P+i)->v[1];
	}
	for(i=0;i<N;i++) {
		(P+i)->v[0]-=sumx/N;
		(P+i)->v[1]-=sumy/N;
	}
}

/*compute force*/
void accel(Particle *P)
{
	int	i,j,k;
	float	rij, Fijx, Fijy, tmpx, tmpy;
	for(i=0;i<N;i++) 
	{
		(P+i)->a[0] = 0;
		(P+i)->a[1] = 0;
	}
	for(i=0;i<N;i++) 
	  for(j=0;j<N;j++) 
		if(i!=j)
		{
			tmpx = ((P+j)->x[0]-(P+i)->x[0]);
			if(tmpx>=L/2.)
				tmpx-=L;
			else if(tmpx<=-L/2)
				tmpx+=L;
			tmpy = ((P+j)->x[1]-(P+i)->x[1]);
			if(tmpy>=L/2.)
				tmpy-=L;
			else if(tmpy<=-L/2)
				tmpy+=L;

			rij = sqrt(tmpx*tmpx + tmpy*tmpy);
			rij *=rij;
			Fijx = tmpx/rij;
			Fijy = tmpy/rij;
			rij *= rij*rij;
			Fijx*=(6./rij-12./(rij*rij));
			Fijy*=(6./rij-12./(rij*rij));
			(P+i)->a[0]+=Fijx;
			(P+i)->a[1]+=Fijy;
		}	
}


/*calculate temperature*/
double calc_temp(Particle *P)
{   
	int		i;
	double  Ek = 0.0;
	for(i=0;i<N;i++) {
		Ek += ((P+i)->v[0]*(P+i)->v[0]+(P+i)->v[1]*(P+i)->v[1]);
	}
    return Ek/N; /*<Ek> = NKT*/
}



/*main function*/
int main(int argc, char *argv[])
{
	int	j, i, k;
	FILE *fp;
	clock_t time0=clock();
	Particle *P = (Particle *)calloc(sizeof(Particle), N);

	/*Initialize particles positions*/
	for(j=1, k=0;j<=sqrt((float)N);j++) 
		for(i=1;i<=sqrt((float)N);i++, k++) 
		{
			(P+k)->x[0] = (float)i*d;
			(P+k)->x[1] = (float)j*d;
		} 
	
    /*initialize particle velocities*/
	for(i=0;i<N;i++) 
	{
		(P+i)->v[0] = RAND * vmax;
		(P+i)->v[1] = RAND * vmax;
	}
	setVCMzero(P);
	
    /*open output file in .xyz format*/
    fp = fopen("configuration.xyz", "w");
    if(fp == NULL) {
      printf("ERROR: cannot open file!");
      exit(1);             
    }
       
	/*time integration*/
	accel(P);
	for(k=0;k<=1000;k++) 
	{
		for(i=0;i<N;i++) 
		{
			(P+i)->x[0] = (fmod( ((P+i)->x[0] + (P+i)->v[0]*h + (P+i)->a[0]*h*h*.5 )+ L , L));
			(P+i)->x[1] = (fmod( ((P+i)->x[1] + (P+i)->v[1]*h + (P+i)->a[1]*h*h*.5 )+ L , L));
			(P+i)->v[0] = (P+i)->v[0] + (P+i)->a[0]*h*.5;
			(P+i)->v[1] = (P+i)->v[1] + (P+i)->a[1]*h*.5;
		}
		accel(P);
		for(i=0;i<N;i++) 
		{
			(P+i)->v[0]+= (P+i)->a[0]*h*.5;
			(P+i)->v[1]+= (P+i)->a[1]*h*.5;
			
		}
		// save data	
		if(k%100==0) 
		{
            if (k==0) {
                printf(" Step Temp\n");
            }

			printf(" %10d %10.5f\n", k, calc_temp(P));
			fprintf(fp, "%d\n\n", N);
			for(i=0;i<N;i++)
				fprintf(fp, "Ar\t%f\t%f\t%f\n", (P+i)->x[0], (P+i)->x[1], 0.f);
			fflush(fp);
		}
	}
	free(P);
	fclose(fp);
	printf("Finished after %.3f seconds.\n", (double)(clock() - time0) / CLOCKS_PER_SEC);
	return 0;
}
