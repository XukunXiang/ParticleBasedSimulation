/* This is an effort for estabulish a better MD code for Particle Based simulation course and for USPAS summer school*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "constants.h"
#include "init.h"
#include "integrator.h"
#include "forcepotential.h"
#include "output.h"

int main() {
	double R[N][dim]={{0.0}}, P[N][dim] = {{0.0}}, F[N][dim] = {{0.0}};
	int iter;
	double dt = 1.0e-3, realt = 0.0;
	int plotstride = 200;
	FILE *RPo,*Energyo,*To_warmup,*To_run;

	clock_t cpu_start, cpu_end;
	double cpu_dt;

	RPo = fopen("RandP.xyz","w+");
	Energyo = fopen("KE_PE_Etot.dat","w+");
	To_warmup = fopen("warmupT.dat","w+");
	To_run = fopen("runT.dat","w+");

	//initialize random seed:
	srand(time(NULL));
	
	//Initialization
//	Init_R(R);
	Init_R_FCC(R);
//	Init_P(P);
	Init_P_random(P);
///	output(R,P,realt,RPo,Energyo);

	//Warmup Run
	for (iter = 1; iter<=(Itime+1); iter++) {
		double t_temp;
		VelocityVerlet(R,P,F,dt);
		t_temp = getT(P);
//		rescale(P);
		if (iter % 50 == 0){rescale_single(P);}
		fprintf(To_warmup,"%11.5f \n",t_temp);
	}
	
	//Simulation
	ForceCalculation(R,F);
cpu_start = clock();
	for (iter = 1; iter<=(Ntime+1); iter++) {
		realt += dt;
		VelocityVerlet(R,P,F,dt);
		fprintf(To_run,"%11.5f \t %11.5f \n",realt,getT(P));	
		if ((iter % plotstride) == 0){
			output(R,P,realt,RPo,Energyo);
		//	printf("%11.5f \n",getT(P));	
		}
	}
cpu_end = clock();
cpu_dt = ((double)(cpu_end-cpu_start)/CLOCKS_PER_SEC)*1000.0;
printf("Time: %.3f ms\n",cpu_dt/Ntime);

	fclose(RPo);
	fclose(Energyo);
	fclose(To_warmup);
	fclose(To_run);

	return 0;
}
