#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include "init.h"
#include "constants.h"

double getrand() {
	return (double)rand()/(double)RAND_MAX;
}

void Init_R(double R[N][3]){
	int numb, check, i, j;
	double r0[3],rel,s,thr;
	
	thr = cutoff*cutoff;	
	numb = 0;
	while (numb < N){
		for (j=0; j<3; j++)
			r0[j] = L*(2.0*getrand() - 1.0);		
		check = 0;
		//check for proximity
		for (i=0; i<numb; i++){
			s = 0.0;
			for (j=0; j<3; j++) {
				rel = R[i][j]-r0[j];
				rel -= L*nearbyint(rel/L);
				s += pow(rel,2.0);
			}
			if (s < thr){
				check = 1;
				break;
			}
		}
		if (check == 0){
			for (j=0; j<3; j++) {
				R[numb][j] = r0[j];
			}
			numb += 1;
		}
	}
}

void Init_R_FCC(double R[N][3]){
	int numb,i,j,k;
	
	numb = 0;
	for (i=0;i<nc;i++) {
		for(j=0;j<nc;j++) {
			for(k=0;k<nc;k++) {
				R[numb][0] = i*1.0*a0-L/2;
				R[numb][1] = j*1.0*a0-L/2;
				R[numb][2] = k*1.0*a0-L/2;
				R[numb+1][0] = (i*1.0+0.5)*a0-L/2;
				R[numb+1][1] = (j*1.0+0.5)*a0-L/2;
				R[numb+1][2] = k*1.0*a0-L/2;
				R[numb+2][0] = (i*1.0+0.5)*a0-L/2;
				R[numb+2][1] = j*1.0*a0-L/2;
				R[numb+2][2] = (k*1.0+0.5)*a0-L/2;
				R[numb+3][0] = i*1.0*a0-L/2;
				R[numb+3][1] = (j*1.0+0.5)*a0-L/2;
				R[numb+3][2] = (k*1.0+0.5)*a0-L/2;
				numb += 4;
			}
		}
	}
}

void Init_P(double P[N][3]) {
	int i,j;
	double temp,lambda = sqrt(T_init);
	double avgv[3]={0.0};

	for (i = 0; i < N; i++) {
		for (j = 0; j < 3; j++){
			temp = sqrt(-2.0*log(getrand()))*cos(2.0*PI*getrand());
			P[i][j] = temp;
			avgv[j] += temp; 
		}
	}
	for (j=0; j<3; j++){
		avgv[j] = avgv[j]/N;
	}
	for (i=0;i<N;i++){
		for (j = 0; j < 3; j++){
			P[i][j] = lambda*(P[i][j]-avgv[j]);
		}
	}
}

void Init_P_random(double P[N][3]) {
	int i,j;
	double temp;
	double avgv[3]={0.0};

	for (i = 0; i < N; i++) {
		for (j = 0; j < 3; j++){
			temp = getrand();
			P[i][j] = temp;
			avgv[j] += temp; 
		}
	}
	for (j=0; j<3; j++){
		avgv[j] = avgv[j]/N;
	}
	for (i=0;i<N;i++){
		for (j = 0; j < 3; j++){
			P[i][j] = P[i][j]-avgv[j];
		}
	}
}

double rescale(double P[N][3]){
	double T_current,s;
	int i,j;
	
	T_current = getT(P);
	s = sqrt(1+(T_init/T_current-1)/(double(Itime)/3.0));
	for (i=0;i<N;i++) {
		for (j=0;j<3;j++){
			P[i][j] = P[i][j]*s;
		}
	}
	return T_current;
}

double rescale_single(double P[N][3]){
	double T_current,lambda;
	int i,j;
	
	T_current = getT(P);
	lambda = sqrt(T_init/T_current);
	for (i=0;i<N;i++) {
		for (j=0;j<3;j++){
			P[i][j] = P[i][j]*lambda;
		}
	}
	return T_current;
}


double getT(double P[N][3]){
	double avgp[3]={0.0},p2=0.0,t_current;
	int i,j;
	//get <P>
	for (i=0;i<N;i++) {
		for (j=0;j<3;j++){
			avgp[j] += P[i][j];
		}
	}
	for (j=0;j<3;j++){
		avgp[j] = avgp[j]/N;
	}
	//substract <P> from every Pi
	for (i=0;i<N;i++) {
		for (j=0;j<3;j++){
			P[i][j] -= avgp[j];
			p2 += pow(P[i][j],2.0);
		}
	}
	t_current = sqrt(m*p2/(3.0*(N-1)));
	return t_current;
}
