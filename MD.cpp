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
	//int plotstride = 2;
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
	double **pi = new double*[3];//outgoing particle P
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
		for (i=0;i<3;i++)	pi[i] = new double[ni[i+1]*2];
		int bidx[3]={0},idx;
		for (i=0;i<N;i++){
			bid = boxid[i];
			idx = 2*bidx[bid];
			ri[bid][idx] = R[i][0];
			ri[bid][idx+1] = R[i][1];
			pi[bid][idx] = P[i][0];
			pi[bid][idx+1] = P[i][1];
			bidx[bid] += 1;
		}
	}

	// send out particles
	int numlocal=0; //paticle number in local box
	// send out number of particles first
	MPI_Scatter(&ni[0],1,MPI_INT,&numlocal,1,MPI_INT,0,MPI_COMM_WORLD);
	
	// then send out the R and P
	double *R_local, *P_local;
	if (rank !=0){
		R_local = new double[numlocal*2];
		P_local = new double[numlocal*2];
		printf("This is box#%d, my numlocal is %d \n",rank,numlocal);
		// create array for R and P of local particles
		MPI_Recv(R_local,2*numlocal,MPI_DOUBLE,0,rank,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Recv(P_local,2*numlocal,MPI_DOUBLE,0,rank,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}else{ //rank == 0
		for (i=1;i<4;i++){
			MPI_Send(&ri[i-1][0],2*ni[i],MPI_DOUBLE,i,i,MPI_COMM_WORLD);
			MPI_Send(&pi[i-1][0],2*ni[i],MPI_DOUBLE,i,i,MPI_COMM_WORLD);
		}
	}
	// clean up
	if (rank ==0){
		for (i=0;i<3;i++){
			delete [] ri[i];
			delete [] pi[i];
		}
	}
	delete[] ri;
	delete[] pi;

	//=======Simulation=======
	double *F_local = new double[2*numlocal];
	MPI_Datatype stype;
	if (rank != 0)	ForceCalculation_MPI_3box(R_local,F_local,numlocal);
	for (iter = 1; iter<=(Ntime+1); iter++) {
		realt += dt;
		if (rank !=0) {
			printf("***before VV*****rank#%d: %7.3f  %d\n",rank,R_local[0],numlocal);
			VelocityVerlet_MPI_3box(R_local,P_local,F_local,dt,numlocal);
			printf("***after VV*****rank#%d: %7.3f  %d\n",rank,R_local[0],numlocal);
		}

//		fprintf(To_run,"%11.5f \t %11.5f \n",realt,getT(P));	
		if ((iter % plotstride) == 0){
			printf("iter:%d, numlocal = %d \n",iter,numlocal);
			//***gather data back to master
			double *Rarray = new double[N*2];
			double *Parray = new double[N*2];
			int *NumInBox,*displs,*rcounts;
			//******gather numlocal first
			if (rank==0) {
				NumInBox = new int[4];
				displs = new int[4];
				rcounts = new int[4];
			}
			MPI_Gather(&numlocal,1,MPI_INT,&NumInBox[0],1,MPI_INT,0,MPI_COMM_WORLD);
			if (rank==0) {
				int offset = 0;
				for (i=0;i<4;i++){
					displs[i] = offset;
					offset += 2*NumInBox[i];
					rcounts[i] = 2*NumInBox[i];
					printf("%d \t",NumInBox[i]);
				}
				printf("\n");
			}
			//******then gather data with shift
			MPI_Type_vector(numlocal,2,2,MPI_DOUBLE,&stype);
			MPI_Type_commit(&stype);
			MPI_Gatherv(R_local,1,stype,Rarray,rcounts,displs,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Gatherv(P_local,1,stype,Parray,rcounts,displs,MPI_DOUBLE,0,MPI_COMM_WORLD);
			//clean up
			if (rank==0) {
				delete [] NumInBox;
				delete [] displs;
				delete [] rcounts;
				delete [] Rarray;
				delete [] Parray;
//				output(Rarray,Parray,realt,RPo,Energyo);
			}
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
