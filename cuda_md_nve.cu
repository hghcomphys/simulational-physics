
/*  ----------------------------------------------------------
	Project name  : Molecular Dynamics Simulation on GPU
	Description   : This program simulates many particles 
					system using CUDA thechnology on
					cuda enabled nvidia GPU.
	
	What have done:
	- A simple (loop over all particles) & 2D code implemented
    ---------------------------------------------------------- */ 

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// Includes

#include<cstdio>
#include<iostream>
#include<fstream>
using namespace std;


// Macroes and Symbolic Constants

#define N		10000 // number of particles
#define d       	2.5f // lattice distance
#define L		(sqrtf(N)*d) // simulation box length
#define h		.0001f // time step
#define vmax		10.f // maximum random velocity
#define block_dim	32  // number of threads in each block
#define RAND		(2.f*((float)rand()/RAND_MAX-.5f))

// Typedefs

typedef struct
{
	float2 r;
	float2 v;
	float2 f;
} Particle;

// Functions

__host__   void setVCMzero(Particle *);
__global__ void kernel_ComputeForce(void *);
__host__   void InitializeParticle(Particle *);
__global__ void kernel_TimeIntegration_0(void *);
__global__ void kernel_TimeIntegration_1(void *);
__device__ float2 TowParticlesInteraction(float2, float2, float2);

// Implementation of functions

int main()
{
	// Multiple GPUs --------------------------
	int n_GPU;
	cudaGetDeviceCount(&n_GPU); 
	cudaDeviceProp prop;

	std::cout<<"Number of GPU is "<<n_GPU<<std::endl;
	for(int i=0;i<n_GPU;++i)
	{
		
		cudaGetDeviceProperties(&prop,i);
		std::cout<<"GPU No. "<<i<<" "<< prop.name <<std::endl;
	}
	cudaSetDevice(0);
	std::cout<<std::endl;

	Particle *dev_P; // device particles pointer
	Particle *hst_P = new Particle[N]; // host particles pointer allocation
	const int size = N*sizeof(Particle); // size of allocated memory
	ofstream out("configuration.xyz");
	InitializeParticle(hst_P);
	cudaMalloc((void**)&dev_P, size); // device particle pointer allocation
	cudaMemcpy(dev_P, hst_P, size, cudaMemcpyHostToDevice); // host to device memory copy
	kernel_ComputeForce<<<N/block_dim, block_dim>>>(dev_P); // compute force kernel function
	for(int t=0;t<=10000;t++)
	{
		if(t%100==0)
		{
			cudaMemcpy(hst_P, dev_P, size, cudaMemcpyDeviceToHost); // device to host memory copy
    			cout << " t --> " << t << endl;
			out << N << "\n\n";
			for(int i=0;i<N;i++)
				out <<"Ar\t"<<(hst_P+i)->r.x<<"\t"<<(hst_P+i)->r.y<<"\t"<<0.f<<endl;
		}
		// fast verlet time integration
		kernel_TimeIntegration_0<<<N/block_dim, block_dim>>>(dev_P);
		kernel_ComputeForce<<<N/block_dim, block_dim>>>(dev_P);
		kernel_TimeIntegration_1<<<N/block_dim, block_dim>>>(dev_P);
		
	}
	out.close();
	delete[] hst_P;
	cudaFree(dev_P);
	return 0;
}

__host__ void InitializeParticle(Particle *P)
{
	for (int j=0, k=0;j<(int)sqrt((double)N);j++) 
		for (int i=0;i<(int)sqrt((double)N);i++,k++) 
		{
			// initialize particle's coordinates
			(P+k)->r.x = (float)i*d; 
			(P+k)->r.y = (float)j*d; 	
			// initialize velocities
			(P+k)->v.x = vmax*RAND;
			(P+k)->v.y = vmax*RAND;
		}
	setVCMzero(P);
}

__host__ void setVCMzero(Particle *P)
{   
	int i;
	float2	Vcm={0.f, 0.f};
	for(i=0;i<N;i++) {
		Vcm.x+=(P+i)->v.x;
		Vcm.y+=(P+i)->v.y;
	}
	for(i=0;i<N;i++) {
		(P+i)->v.x-=Vcm.x/N;
		(P+i)->v.y-=Vcm.y/N;
	}
}
		
__device__ float2 TowParticlesInteraction(float2 ri, float2 rj, float2 fi)
{	
	float r;
	float2 rij, fij; 
    // periodic boundary condtion
	rij.x = rj.x - ri.x;
	if (rij.x>=L/2.f) 
		rij.x-=L;
	else if (rij.x<=-L/2.f) 
			rij.x+=L;

	rij.y = rj.y - ri.y;
	if (rij.y>=L/2.f) 
		rij.y-=L;
	else if (rij.y<=-L/2.f)
			rij.y+=L;	
	// LJ force field
	r = sqrtf(rij.x*rij.x + rij.y*rij.y);	
	if (r>0.000001f)
	{
		r*=r;
		fij.x = rij.x/r; 
		fij.y = rij.y/r;
		r*= r*r;
		fij.x*= (6.f/r-12.f/r/r); 
		fij.y*= (6.f/r-12.f/r/r);
		// update forces
		fi.x += fij.x;
		fi.y += fij.y;
	}
	return fi;
}

__global__ void kernel_ComputeForce(void *global_P)
{
	
	float2 force = {0.f, 0.f};
	Particle *P =(Particle *)global_P;
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	// compute force for each thread
	for (int i=0;i<N;i++)
		if(i != idx)
			force = TowParticlesInteraction((P+idx)->r, (P+i)->r, force);
	(P+idx)->f = force;
	__syncthreads();
}

__global__ void kernel_TimeIntegration_0(void *global_P)
{
	Particle *P =(Particle *)global_P;
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	// update coordinates and velocities
	(P+idx)->r.x=fmod((P+idx)->r.x + h * ( (P+idx)->v.x + 0.5f*h*(P+idx)->f.x ) + 10.f*L, L);
	(P+idx)->r.y=fmod((P+idx)->r.y + h * ( (P+idx)->v.y + 0.5f*h*(P+idx)->f.y ) + 10.f*L, L);
	(P+idx)->v.x+=0.5f*h*(P+idx)->f.x;
	(P+idx)->v.y+=0.5f*h*(P+idx)->f.y;	
	__syncthreads();
}

__global__ void kernel_TimeIntegration_1(void *global_P)
{
	Particle *P =(Particle *)global_P;
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	// ubdate velocities with new force
	(P+idx)->v.x+= 0.5f*h*(P+idx)->f.x;
	(P+idx)->v.y+= 0.5f*h*(P+idx)->f.y;			
	__syncthreads();
}

