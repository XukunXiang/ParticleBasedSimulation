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
	MPI_Init_thread(&argc,&argv,MPI_THREAD_SINGLE,&provided);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//=======Initialization=======
	if (rank == 0){
		Init_R_2D_3box(R);
		Init_P_random(P);
		int ni[3] = {0},bid;
		// counting number of particles in each cell and check
		for (i=0;i<N;i++){
			bid = int(floor(R[i][0]/L));
			ni[bid] += 1;
		}
		for (i=0;i<3;i++){
			if (ni[i] ==0){
				printf("no particle in box #%d \n",i);
				return 0;
			} else {
				printf("%d particles in box #%d \n",ni[i],i);
			}
		}
	}
	// send out particles
	// send out number of particles first
	// then send out the R and P

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
