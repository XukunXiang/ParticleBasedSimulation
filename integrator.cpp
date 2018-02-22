#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "constants.h"
#include "integrator.h"
#include "forcepotential.h"

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
	//forwardVR
	for (i=0; i<N; i++){
		for (j=0; j<dim; j++){
			p_d[i*dim+j] += f_d[i*dim+j]*hdt;
			r_d[i*dim+j] = r_d[i*dim+j]+p_d[i*dim+j]*dtm;
			r_d[i*dim+j] = fmod(r_d[i*dim+j],L);
		}
	}
	ForceCalculation_MPI_3box(r_d,f_d,numlocal);
	//forwardV
	for (i=0; i<N; i++){
		for (j=0; j<dim; j++){
			p_d[i*dim+j] += f_d[i*dim+j]*hdt;
		}
	}
}

