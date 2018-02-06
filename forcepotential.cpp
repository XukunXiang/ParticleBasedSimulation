#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "constants.h"
#include "forcepotential.h"

void ForceCalculation(double R[N][3],double F[N][3]){
	int i,j,k;
	double rj[3],rel[3],rel_len2,fij,fijk;	
	
	//initialize F[N][3]
	for (i=0; i<N; i++){
		for (j=0; j<3; j++){
			F[i][j] = 0.0;
		}
	}
	for (j = 0; j < (N-1); j++) {
		for (k=0; k<3; k++){
			rj[k] = R[j][k];
		}
		for (i = j+1; i < N; i++){
			rel_len2 = 0.0;
			for (k = 0; k<3; k++) {
				rel[k] = R[i][k]-rj[k];
				rel[k] -=  L*nearbyint(rel[k]/L);
				rel_len2 += pow(rel[k],2.0);
			}
			//Poor man's neighbor list
			if (rel_len2>Fcutoff2) continue;
			fij = 4.0*(6.0*pow(rel_len2,-4.0)-12.0*pow(rel_len2,-7.0));
			for (k = 0; k<3; k++) {
				fijk = fij*rel[k];
				F[j][k] += fijk;
				F[i][k] -= fijk;
			}
		}
	}
}

double getPE(double R[N][3]) {
	double PE,rel_len2,rijk;
	int i,j,k;
	
	PE = 0.0;
	for (j = 0; j < (N-1); j++) {
		for (i = j+1; i < N; i++){
			rel_len2 = 0.0;
			for (k = 0; k<3; k++){
				rijk = R[i][k]-R[j][k];
				rijk -= L*nearbyint(rijk/L);
				rel_len2 += pow(rijk,2.0);		
			}
			PE += 4*(pow(rel_len2,-6.0) - pow(rel_len2,-3.0));
		}
	}
	return PE;
}
