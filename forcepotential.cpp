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

void ForceCalculation_MPI_3box(double *r_d,double *f_d,int numlocal){
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
	//packing outgoing ghost particles' position
	double *GhostToL = new double[ngtl*2];
	double *GhostToR = new double[ngtr*2];
	int *idGTL = new int[ngtl];
	int *idGTR = new int[ngtr];
	int idxL=0,idxR=0;
	for (i=0; i<numlocal; i++){
		if (gtlist[i] == 1){
			GhostToR[idxR*2] = r_d[i*2];
			GhostToR[idxR*2+1] = r_d[i*2+1];
			idGTR[idxR] = i;
			idxR += 1;
		} else if (gtlist[i] == -1){
			GhostToL[idxL*2] = r_d[i*2];
			GhostToL[idxL*2+1] = r_d[i*2+1];
			idGTL[idxL] = i;
			idxL += 1;
		}
	}
	
	MPI_Status sta[4];
	MPI_Request rq[4];
	//send/recv number first
	//**since we are on rank=1,2,3 and size=4
	//left neighbor(with periodic BC) for [0,1,2] in size=3 should be (rank-1+size)%size
	//so in our case, ((rank-1)-1+(size-1))%(size-1)+1
	int ln,rn;
	ln = (rank+size-3)%(size-1)+1;
	//(rank-1+1+size-1)%(size-1)+1
	rn = (rank+size-1)%(size-1)+1;
	//send/recv ghost particle number first, tag is 10*source+dest
	MPI_Irecv(&ngfl,1,MPI_INT,ln,10*ln+rank,MPI_COMM_WORLD,&rq[0]);
	MPI_Irecv(&ngfr,1,MPI_INT,rn,10*rn+rank,MPI_COMM_WORLD,&rq[1]);
	MPI_Isend(&ngtl,1,MPI_INT,ln,10*rank+ln,MPI_COMM_WORLD,&rq[2]);
	MPI_Isend(&ngtr,1,MPI_INT,rn,10*rank+rn,MPI_COMM_WORLD,&rq[3]);
	MPI_Waitall(4,rq,sta);
	//send/recv ghost particle positions
	double *GhostFromL = new double[ngfl*2];
	double *GhostFromR = new double[ngfr*2];
	MPI_Irecv(GhostFromL,2*ngfl,MPI_DOUBLE,ln,10*ln+rank,MPI_COMM_WORLD,&rq[0]);
	MPI_Irecv(GhostFromR,2*ngfr,MPI_DOUBLE,rn,10*rn+rank,MPI_COMM_WORLD,&rq[1]);
	MPI_Isend(GhostToL,2*ngtl,MPI_DOUBLE,ln,10*rank+ln,MPI_COMM_WORLD,&rq[2]);
	MPI_Isend(GhostToR,2*ngtr,MPI_DOUBLE,rn,10*rank+rn,MPI_COMM_WORLD,&rq[3]);
	MPI_Waitall(4,rq,sta);

	//check
	/*
	printf("rank%d: %d from ln%d, %d from rn%d\n",rank,ngfl,ln,ngfr,rn);
	for (i=0;i<3;i++){
		printf("#%d from %d to %d: %11.3f\t %11.3f \n",i,rank,ln,GhostToL[2*i],GhostToL[2*i+1]);
		printf("#%d from %d to %d: %11.3f\t %11.3f \n",i,rank,rn,GhostToR[2*i],GhostToR[2*i+1]);
		printf("#%d from %d to %d: %11.3f\t %11.3f \n",i,ln,rank,GhostFromL[2*i],GhostFromL[2*i+1]);
		printf("#%d from %d to %d: %11.3f\t %11.3f \n",i,rn,rank,GhostFromR[2*i],GhostFromR[2*i+1]);
	} 
	*/

	// force from local particles
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

	//froce from ghost particles, only need to look at ourgoing particles
	//[later] make a function call on ghost contribution
	for (j = 0; j < ngtl; j++) {
		int jidx = idGTL[j];
		for (k=0; k<dim; k++){
			rj[k] = GhostToL[j*dim+k]; //rj[k] = r_d[jidx*dim+k];
		}
		for (i = 0; i <ngfl ; i++){
			rel_len2 = 0.0;
			for (k = 0; k<dim; k++) {
				rel[k] = GhostFromL[i*dim+k]-rj[k];
				rel[k] -=  lxly[k]*nearbyint(rel[k]/lxly[k]);
				rel_len2 += pow(rel[k],2.0);
			}
			//Poor man's neighbor list
			if (rel_len2>Fcutoff2) continue;
			fij = 4.0*(6.0*pow(rel_len2,-4.0)-12.0*pow(rel_len2,-7.0));
			for (k = 0; k<dim; k++) {
				fijk = fij*rel[k];
				f_d[jidx*dim+k] += fijk;
			}
		}
	}
	for (j = 0; j < ngtr; j++) {
		int jidx = idGTR[j];
		for (k=0; k<dim; k++){
			rj[k] = GhostToR[j*dim+k]; //rj[k] = r_d[jidx*dim+k];
		}
		for (i = 0; i <ngfr; i++){
			rel_len2 = 0.0;
			for (k = 0; k<dim; k++) {
				rel[k] = GhostFromR[i*dim+k]-rj[k];
				rel[k] -=  lxly[k]*nearbyint(rel[k]/lxly[k]);
				rel_len2 += pow(rel[k],2.0);
			}
			//Poor man's neighbor list
			if (rel_len2>Fcutoff2) continue;
			fij = 4.0*(6.0*pow(rel_len2,-4.0)-12.0*pow(rel_len2,-7.0));
			for (k = 0; k<dim; k++) {
				fijk = fij*rel[k];
				f_d[jidx*dim+k] += fijk;
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
