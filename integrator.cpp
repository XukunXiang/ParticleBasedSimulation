#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "constants.h"
#include "integrator.h"
#include "forcepotential.h"
#include "mpi.h"

void VelocityVerlet(double R[N][dim],double P[N][dim],double F[N][dim],double dt){
	int i,j;
	double hdt = 0.5*dt,dtm = dt/m;
	//forwardVR
	for (i=0; i<N; i++){
		for (j=0; j<dim; j++){
			P[i][j] += F[i][j]*hdt;
			R[i][j] = R[i][j]+P[i][j]*dtm;
			R[i][j] = fmod(R[i][j],L);
		}
	}
	ForceCalculation(R,F);
	//forwardV
	for (i=0; i<N; i++){
		for (j=0; j<dim; j++){
			P[i][j] += F[i][j]*hdt;
		}
	}
}

void VelocityVerlet_MPI_3box(double r_d[],double p_d[],double f_d[],double dt,int numlocal){
	int i,j;
	double hdt = 0.5*dt,dtm = dt/m;

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double left = (rank-1)*L,right = rank*L;
	double lms = left-skin,lps = left+skin;
	double rms = right-skin,rps = right+skin;
	//forwardVR
	for (i=0; i<numlocal; i++){
		//x direction
		p_d[i*dim] += f_d[i*dim]*hdt;
		double rt = r_d[i*dim]+p_d[i*dim]*dtm;
//		if (rt < rms && rt >lps) 
//		r_d[i*dim] = rt;

		//r_d[i*dim] = fmod(r_d[i*dim],L);
		p_d[i*dim+1] += f_d[i*dim+1]*hdt;
		r_d[i*dim+1] = r_d[i*dim+1]+p_d[i*dim+1]*dtm;
		//r_d[i*dim+1] = fmod(r_d[i*dim+1],L);
	}
	ForceCalculation_MPI_3box(r_d,f_d,numlocal);
	//forwardV
	for (i=0; i<numlocal; i++){
		for (j=0; j<dim; j++){
			p_d[i*dim+j] += f_d[i*dim+j]*hdt;
		}
	}
}

