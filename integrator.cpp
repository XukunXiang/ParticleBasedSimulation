#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "constants.h"
#include "integrator.h"
#include "forcepotential.h"

void VelocityVerlet(double R[N][3],double P[N][3],double F[N][3],double dt){
	int i,j;
	double hdt = 0.5*dt,dtm = dt/m;
	//forwardVR
	for (i=0; i<N; i++){
		for (j=0; j<3; j++){
			P[i][j] += F[i][j]*hdt;
			R[i][j] = R[i][j]+P[i][j]*dtm;
			R[i][j] = fmod(R[i][j],L);
		}
	}
	ForceCalculation(R,F);
	//forwardV
	for (i=0; i<N; i++){
		for (j=0; j<3; j++){
			P[i][j] += F[i][j]*hdt;
		}
	}
}

