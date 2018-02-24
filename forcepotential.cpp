#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "constants.h"
#include "forcepotential.h"
#include "mpi.h"
#include <stdio.h>

void ForceCalculation(double R[N][dim],double F[N][dim]){
	int i,j,k;
	double rj[dim],rel[dim],rel_len2,fij,fijk;	
	
	//initialize F[N][dim]
	for (i=0; i<N; i++){
		for (j=0; j<dim; j++){
			F[i][j] = 0.0;
		}
	}
	for (j = 0; j < (N-1); j++) {
		for (k=0; k<dim; k++){
			rj[k] = R[j][k];
		}
		for (i = j+1; i < N; i++){
			rel_len2 = 0.0;
			for (k = 0; k<dim; k++) {
				rel[k] = R[i][k]-rj[k];
				rel[k] -=  L*nearbyint(rel[k]/L);
				rel_len2 += pow(rel[k],2.0);
			}
			//Poor man's neighbor list
			if (rel_len2>Fcutoff2) continue;
			fij = 4.0*(6.0*pow(rel_len2,-4.0)-12.0*pow(rel_len2,-7.0));
			for (k = 0; k<dim; k++) {
				fijk = fij*rel[k];
				F[j][k] += fijk;
				F[i][k] -= fijk;
			}
		}
	}
}

void ForceCalculation_MPI_3box(double r_d[],double f_d[],int numlocal){
	int i,j,k;
	double rj[dim],rel[dim],rel_len2,fij,fijk;	 
	
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double lps = (rank-1)*L+skin,rms = rank*L-skin;
	int ngtl = 0, ngtr = 0;//number of ghost to neighbor
	int ngfl, ngfr; // number of ghost from neighbor
	int gtlist[numlocal],dummy=0;//outgoing ghost particle list
	double x;
	//initialize F[N][dim]
	for (i=0; i<numlocal; i++){
		for (j=0; j<dim; j++){
			f_d[i*dim+j] = 0.0;
		}
		//setup outgoing ghost particles
		x = r_d[i*dim];
		if (lps <= x and x <= rms) {
			gtlist[i] = 0;
			dummy +=1;
		}else if (x > rms) {
			gtlist[i] = 1;
			ngtr += 1;
		}else if (x < lps) {
			gtlist[i] = -1;
			ngtl += 1;
		}
	}
	printf("@%d:%d to left, %d to right; number check %d\n",rank,ngtl,ngtr,ngtl+dummy+ngtr);
	//packing outgoing ghost particles' position
	double *GhostToL = new double[ngtl*2];
	double *GhostToR = new double[ngtr*2];
	int idxL=0,idxR=0;
	for (i=0; i<numlocal; i++){
		if (gtlist[i] == 1){
			GhostToR[idxR*2] = r_d[i*2];
			GhostToR[idxR*2+1] = r_d[i*2+1];
			idxR += 1;
		} else if (gtlist[i] == -1){
			GhostToL[idxL*2] = r_d[i*2];
			GhostToL[idxL*2+1] = r_d[i*2+1];
			idxL += 1;
		}
	}
	//send/recv number first
	//send/recv ghost particle positions
	//check

	for (j = 0; j < (numlocal-1); j++) {
		for (k=0; k<dim; k++){
			rj[k] = r_d[j*dim+k];
		}
		for (i = j+1; i < numlocal; i++){
			rel_len2 = 0.0;
			for (k = 0; k<dim; k++) {
				rel[k] = r_d[i*dim+k]-rj[k];
				rel[k] -=  L*nearbyint(rel[k]/L);
				rel_len2 += pow(rel[k],2.0);
			}
			//Poor man's neighbor list
			if (rel_len2>Fcutoff2) continue;
			fij = 4.0*(6.0*pow(rel_len2,-4.0)-12.0*pow(rel_len2,-7.0));
			for (k = 0; k<dim; k++) {
				fijk = fij*rel[k];
				f_d[j*dim+k] += fijk;
				f_d[i*dim+k] -= fijk;
			}
		}
	}
}

double getPE(double R[N][dim]) {
	double PE,rel_len2,rijk;
	int i,j,k;
	
	PE = 0.0;
	for (j = 0; j < (N-1); j++) {
		for (i = j+1; i < N; i++){
			rel_len2 = 0.0;
			for (k = 0; k<dim; k++){
				rijk = R[i][k]-R[j][k];
				rijk -= L*nearbyint(rijk/L);
				rel_len2 += pow(rijk,2.0);		
			}
			PE += 4*(pow(rel_len2,-6.0) - pow(rel_len2,-3.0));
		}
	}
	return PE;
}
