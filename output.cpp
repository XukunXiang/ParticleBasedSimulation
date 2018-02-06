#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "forcepotential.h"

void output(double R[N][3],double P[N][3],double realt,FILE *frp,FILE *fenergy) {
	double KE = 0.0,PE,kei;
	int i,j;

	fprintf(frp,"%7d\n",N);
	fprintf(frp,"%11.3f \n",realt);
	for(i=0; i<N; i++){
		fprintf(frp,"%7d \t",1);
		for(j=0; j<3; j++){
			fprintf(frp,"%11.3f \t",R[i][j]);
		}
		kei = 0.0;
		for(j=0; j<3; j++){
			KE += P[i][j]*P[i][j];
			kei += P[i][j]*P[i][j];
			fprintf(frp,"%11.5f \t",P[i][j]);
		}
		fprintf(frp,"%11.5f \t",kei*0.5);
		fprintf(frp,"\n");
	}
	KE = 0.5*m*KE;
	PE = getPE(R);
	fprintf(fenergy,"%11.5f \t %11.5f \t %11.5f \t %11.5f \n", realt, PE, KE, PE+KE);
	printf("%11.5f \t %11.5f \t %11.5f \t %11.5f \n", realt, PE, KE, PE+KE);
}
