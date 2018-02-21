/* This is an effort for estabulish a better MD code for Particle Based simulation course and for USPAS summer school*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mpi.h"
#include "constants.h"
#include "init.h"
#include "integrator.h"
#include "forcepotential.h"
#include "output.h"

int main(int argc,char *argv[]) {
	double R[N][dim]={{0.0}}, P[N][dim] = {{0.0}}, F[N][dim] = {{0.0}};
	int iter;
	double dt = 1.0e-3, realt = 0.0;
	int plotstride = 200;
	FILE *RPo,*Energyo,*To_warmup,*To_run;

	RPo = fopen("RandP.xyz","w+");
	Energyo = fopen("KE_PE_Etot.dat","w+");
	To_warmup = fopen("warmupT.dat","w+");
	To_run = fopen("runT.dat","w+");

	//Initialize random seed:
	srand(time(NULL));
	
	//=======MPI Initialization=======
	int 	rank,size,provided,i;
	MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//=======Initialization=======
	int ni[4] = {0}; // number of particles for each box, ni[0] for master
	double **ri = new double*[3];//outgoing particle R
	if (rank == 0){
		Init_R_2D_3box(R);
		Init_P_random(P);
		// counting number of particles in each cell and check
		int bid,boxid[N];
		for (i=0;i<N;i++){
			bid = (int)floor(R[i][0]/L)+1;
			boxid[i] = bid-1;
			ni[bid] += 1;
		}
		for (i=1;i<4;i++){
			if (ni[i] ==0){
				printf("no particle in box #%d \n",i);
				return 0;
			} else {
				printf("%d particles in box #%d \n",ni[i],i);
			}
		}
		//packing outgoing particles
		for (i=0;i<3;i++)	ri[i] = new double[ni[i+1]*2];
		int bidx[3]={0},idx;
		for (i=0;i<N;i++){
			bid = boxid[i];
			idx = 2*bidx[bid];
			ri[bid][idx] = R[i][0];
			ri[bid][idx+1] = R[i][1];
			bidx[bid] += 1;
		}
		//check
		for (i=0;i<3;i++){
			for (int j=0;j<3;j++){
				printf("%11.3f \t %11.3f\t",ri[i][2*j],ri[i][2*j+1]);
			}
			printf("\n");
		}
	}

	// send out particles
	int numlocal; //paticle number in local box
	// send out number of particles first
	MPI_Scatter(&ni[0],1,MPI_INT,&numlocal,1,MPI_INT,0,MPI_COMM_WORLD);
	
	// then send out the R and P
	double *R_local = new double[numlocal*2];
//	double *P_local = new double[numlocal*2];
	if (rank !=0){
		printf("This is box#%d, my numlocal is %d \n",rank,numlocal);
		// create array for R and P of local particles
		MPI_Recv(&R_local[0],2*numlocal,MPI_DOUBLE,0,rank,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		//recv check
		for (int j=0;j<3;j++){
			printf("#%d in %d: %11.3f \t %11.3f\n",j,rank,R_local[2*j],R_local[2*j+1]);
		}
	}else{ //rank == 0
		for (i=1;i<4;i++){
			MPI_Send(&ri[i-1][0],2*ni[i],MPI_DOUBLE,i,i,MPI_COMM_WORLD);
		}
	}


	//=======Warmup Run=======
/*
	for (iter = 1; iter<=(Itime+1); iter++) {
		double t_temp;
		VelocityVerlet(R,P,F,dt);
		t_temp = getT(P);
		//rescale(P);
		if (iter % 50 == 0){rescale_single(P);}
		fprintf(To_warmup,"%11.5f \n",t_temp);
	}
*/
	
	//=======Simulation=======
	ForceCalculation(R,F);
	for (iter = 1; iter<=(Ntime+1); iter++) {
		realt += dt;
		VelocityVerlet(R,P,F,dt);
		fprintf(To_run,"%11.5f \t %11.5f \n",realt,getT(P));	
		if ((iter % plotstride) == 0){
			output(R,P,realt,RPo,Energyo);
		//	printf("%11.5f \n",getT(P));	
		}
	}

	//=======MPI Finializaiton=======
	MPI_Finalize();

	fclose(RPo);
	fclose(Energyo);
	fclose(To_warmup);
	fclose(To_run);

	return 0;
}
