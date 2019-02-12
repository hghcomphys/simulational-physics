/*
---------------------------------------------------------------------
Project name  : Molecular dynamics simulation (NVE ensemble)
Written by    : Hosein Ghorbanfekr hgh.comphys@gmail.com
Creation Date : May 2009
Description   : This simple C program performs Lennard-Jones liquid 
                MD simulation in 2D using fast Verlet algorithem.
---------------------------------------------------------------------*/ 

//#include "stdafx.h"
#include "stdlib.h" 
#include "malloc.h"
#include "math.h"

// ------------------

#define N		50
#define Nh		100000
#define Nt		20
#define Nf		500
#define Nr		5000
#define h		.0001
#define Beta	.95
#define L		40.
#define Vmax	10.

// -------------------------------------------------------

#define Enable_restart	"yes"	// Enter "yes" or "no"
#define RAND			(2.*(rand()/(float)RAND_MAX - .5))

typedef 
struct	{ 
	double x[2];
	double v[2];
	double a[2];
}Vect;

// --------------------------

void Vcom(Vect *P)
{   
	int		i;
	double	sumx=0,sumy=0;

	for(i=0;i<N;i++) {
		sumx+= (P+i)->v[0];
		sumy+= (P+i)->v[1];
	}
	sumx/=(float)N;
	sumy/=(float)N;
	for(i=0;i<N;i++) {
		(P+i)->v[0]-=sumx;
		(P+i)->v[1]-=sumy;
	}
}

// --------------------------

double Kinetic(Vect *P)
{
	int i;
	double k=0;
	for(i=0;i<N;i++)
		k += (P+i)->v[0]*(P+i)->v[0] + (P+i)->v[1]*(P+i)->v[1];
	k *= .5;
	return k;
}

// --------------------------

double n_left(Vect *P) 
{
	int i,n=0;
	for(i=0;i<N;i++) {
		if((P+i)->x[0]<=L/2.)
			n++;
	}
	return (float)n/(float)N;
}

// --------------------------

void Dec_temprature(Vect *P)
{
	int i;
	for(i=0;i<N;i++) {
		(P+i)->v[0] *= Beta;
		(P+i)->v[1] *= Beta;
	}
}

// --------------------------

double accel(Vect *P,double *Pre)
{
	int				i,j,k;
	double			U=0;
	static double	rij,Fijx,Fijy,tmpx,tmpy;
	(*Pre) = 0;

	for(i=0;i<N;i++) {
		(P+i)->a[0] = 0;
		(P+i)->a[1] = 0;
	}
	for(i=0;i<N;i++) {
		for(j=0;j<i;j++) {
			tmpx = ((P+j)->x[0]-(P+i)->x[0]);
			tmpy = ((P+j)->x[1]-(P+i)->x[1]);
			if(tmpx>=L/2.)
				tmpx-=L;
			else if(tmpx<=-L/2)
				tmpx+=L;
			if(tmpy>=L/2.)
				tmpy-=L;
			else if(tmpy<=-L/2)
				tmpy+=L;
			rij = sqrt( tmpx*tmpx + tmpy*tmpy);
			rij *=rij;
			Fijx = tmpx/rij;
			Fijy = tmpy/rij;
			rij *= rij*rij;
			Fijx*=(6./rij-12./(rij*rij));
			Fijy*=(6./rij-12./(rij*rij));
			(P+i)->a[0]+=Fijx;
			(P+i)->a[1]+=Fijy;
			(P+j)->a[0]+=-Fijx;
			(P+j)->a[1]+=-Fijy;
			U+=-(1./rij-1./(rij*rij));
			*Pre +=-(Fijx*tmpx + Fijy*tmpy);
		}	
	}
	return U;
}

// ----------------------------------

int main(int argc, char *argv[])
{
	int		j,i=0,k=0,nf=1,nr=1;
	double  U,K,Pre;
	FILE	*fpRV,*fp_restart,*fp_data;
	Vect	*P;
	
	P = (Vect *)calloc(sizeof(Vect),N+1);
	if(P==NULL) {
		printf("System can not allocates memory.\n");
		return 1;
	}

	fp_restart = fopen("restart.dat","r");
	if((fp_restart!=NULL) & (Enable_restart=="yes")) { 
		for(i=0;i<N;i++)
			fscanf(fp_restart,"%lf%lf%lf%lf",&(P+i)->x[0],&(P+i)->x[1],&(P+i)->v[0],&(P+i)->v[1]);
		fclose(fp_restart);
		printf("restart.dat loaded ...\n");
		
	}
	else {
		for(j=1;j<=sqrt((float)N);j++) {
			for(i=1;i<=sqrt((float)N);i++) {
				(P+k)->x[0] = (float)i/sqrt((float)N)*L/2.;
				(P+k)->x[1] = (float)j/sqrt((float)N)*L;
				k++;
			} 
		}
		printf("N = %d\nAtoms are initiallized ...\n",k);
		for(i=0;i<N;i++) {
			(P+i)->v[0] = RAND * Vmax;
			(P+i)->v[1] = RAND * Vmax;
		}
		Vcom(P);	
	}

	fp_data = fopen("UKETPn.dat","w");
	fpRV = fopen("RV.dat","w");

	accel(P,&Pre);
	for(k=0;k<Nt;k++) {
		Dec_temprature(P);
		for(j=0;j<Nh;j++) {
			for(i=0;i<N;i++) {
				(P+i)->x[0] = (fmod( ((P+i)->x[0] + (P+i)->v[0]*h + (P+i)->a[0]*h*h*.5 )+100.*L , L));
				(P+i)->x[1] = (fmod( ((P+i)->x[1] + (P+i)->v[1]*h + (P+i)->a[1]*h*h*.5 )+100.*L , L));
				(P+i)->v[0] = (P+i)->v[0] + (P+i)->a[0]*h*.5;
				(P+i)->v[1] = (P+i)->v[1] + (P+i)->a[1]*h*.5;
			}
			U = accel(P,&Pre);
			for(i=0;i<N;i++) {
				(P+i)->v[0]+= (P+i)->a[0]*h*.5;
				(P+i)->v[1]+= (P+i)->a[1]*h*.5;
			}
			// --------------
			if(nf==Nf) {
				for(i=0;i<N;i++)
					fprintf(fpRV,"%1.10f\t%1.10f\t%1.10f\t%1.10f\n",(P+i)->x[0],(P+i)->x[1],(P+i)->v[0],(P+i)->v[1]);
				fflush(fpRV);
					K = Kinetic(P);
				fprintf(fp_data,"%1.10f\t%1.10f\t%1.10f\t%1.10f\t%1.10f\t%1.3f\n",U,K,U+K,K/N,(K + .5*Pre)/(L*L),n_left(P));
				fflush(fp_data);
				nf = 1;
			}
			else
				nf++;
			// ---------------
			if((nr==Nr) & (Enable_restart=="yes")) {
				fp_restart = fopen("restart.dat","w");
				for(i=0;i<N;i++)
					fprintf(fp_restart,"%1.10f\t%1.10f\t%1.10f\t%1.10f\n",(P+i)->x[0],(P+i)->x[1],(P+i)->v[0],(P+i)->v[1]);
				fclose(fp_restart);
				nr = 1;
			}
			else
				nr++;
			// ------------- 
			Vcom(P);
		}
	}
	free(P);
	fclose(fpRV);
	fclose(fp_data);
	return 0;
}
