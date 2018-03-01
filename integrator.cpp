#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "constants.h"
#include "integrator.h"
#include "forcepotential.h"
#include "mpi.h"
#include <stdio.h>

void VelocityVerlet(double R[N][dim],double P[N][dim],double F[N][dim],double dt){
	int i,j;
	double hdt = 0.5*dt,dtm = dt/m;
	//forwardVR
	for (i=0; i<N; i++){
		for (j=0; j<dim; j++){
			P[i][j] += F[i][j]*hdt;
			R[i][j] = R[i][j]+P[i][j]*dtm;
			R[i][j] = fmod(R[i][j],L);
		}
	}
	ForceCalculation(R,F);
	//forwardV
	for (i=0; i<N; i++){
		for (j=0; j<dim; j++){
			P[i][j] += F[i][j]*hdt;
		}
	}
}

//void VelocityVerlet_MPI_3box(double r_d[],double p_d[],double f_d[],double dt,int &numlocal){
void VelocityVerlet_MPI_3box(double *r_d,double *p_d,double *f_d,double dt,int &numlocal){
	int i,j;
	double hdt = 0.5*dt,dtm = dt/m;
 
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	int *mlist = new int[numlocal]; // migration status list
	int gl = 0, gr = 0, gs = 0;
	double left = (rank-1)*L,right = rank*L;
	double rt;
	//forwardVR
	int i2;
	for (i=0; i<numlocal; i++){
		i2 = i*2;
		//x direction
		p_d[i2] += f_d[i2]*hdt;
	//	rt = r_d[i2]+p_d[i2]*dtm;
		rt = r_d[i2]+1.0;
		if (rt < left) {
			r_d[i2] = fmod(rt+lxly[0],lxly[0]);
			mlist[i] = -1;
			gl += 1;
		}else if (rt > right) {
			r_d[i2] = fmod(rt+lxly[0],lxly[0]);
			mlist[i] = 1;
			gr += 1;
		}else {  
			r_d[i2] = rt;
			mlist[i] = 0;
			gs += 1;
		}
		//y direction 
		p_d[i2+1] += f_d[i2+1]*hdt;
		r_d[i2+1] = r_d[i2+1]+p_d[i2+1]*dtm;
		r_d[i2+1] = fmod(r_d[i2+1]+lxly[1],lxly[1]);//PBC in y
	}

	MPI_Status sta[4];
	MPI_Request rq[4];
	//send/recv number first
	//**since we are on rank=1,2,3 and size=4
	//left neighbor(with periodic BC) for [0,1,2] in size=3 should be (rank-1+size)%size
	//so in our case, ((rank-1)-1+(size-1))%(size-1)+1
	int cfl,cfr; // come from left and come from right
	int ln,rn; // left neighbor and right neighbor
	ln = (rank+size-3)%(size-1)+1;
	//(rank-1+1+size-1)%(size-1)+1
	rn = (rank+size-1)%(size-1)+1;
	//send/recv particle numbers
	//send/recv ghost particle number first, tag is 10*source+dest
	MPI_Irecv(&cfl,1,MPI_INT,ln,10*ln+rank,MPI_COMM_WORLD,&rq[0]);
	MPI_Irecv(&cfr,1,MPI_INT,rn,10*rn+rank,MPI_COMM_WORLD,&rq[1]);
	MPI_Isend(&gl,1,MPI_INT,ln,10*rank+ln,MPI_COMM_WORLD,&rq[2]);
	MPI_Isend(&gr,1,MPI_INT,rn,10*rank+rn,MPI_COMM_WORLD,&rq[3]);
	MPI_Waitall(4,rq,sta);

	int numlocal_new = gs + cfl + cfr;
	printf("rank%d: old=%d gl=%d gs=%d gr=%d cfl=%d cfr=%d new=%d\n",rank,numlocal,gl,gs,gr,cfl,cfr,numlocal_new);	
	double *newrd = new double[2*numlocal_new];
	double *newpd = new double[2*numlocal_new];
	double *recvrp = new double[4*(cfl+cfr)];
	double *pgl = new double[4*gl];
	double *pgr = new double[4*gr];
	//packing outgoing particles
	int sidx=0,lidx=0,ridx=0;
	for (i=0;i<numlocal;i++){
		int idx2 = i*2;
		if (mlist[i] == 0) {
			newrd[sidx] = r_d[idx2];
			newrd[sidx+1] = r_d[idx2+1];
			newpd[sidx] = p_d[idx2];
			newpd[sidx+1] = p_d[idx2+1];
			sidx += 2;
		} else if (mlist[i] == 1){
			pgr[ridx] = r_d[idx2];
			pgr[ridx+1] = r_d[idx2+1];
			pgr[ridx+2] = p_d[idx2];
			pgr[ridx+3] = p_d[idx2+1];
			ridx += 4;
		} else { // mlist[i]==-1
			pgl[lidx] = r_d[idx2];
			pgl[lidx+1] = r_d[idx2+1];
			pgl[lidx+2] = p_d[idx2];
			pgl[lidx+3] = p_d[idx2+1];
			lidx += 4;
		}
	}
	//send/recv particle RandP
	MPI_Irecv(&recvrp[0],4*cfl,MPI_DOUBLE,ln,10*ln+rank,MPI_COMM_WORLD,&rq[0]);
	MPI_Irecv(&recvrp[4*cfl],4*cfr,MPI_DOUBLE,rn,10*rn+rank,MPI_COMM_WORLD,&rq[1]);
	MPI_Isend(pgl,4*gl,MPI_DOUBLE,ln,10*rank+ln,MPI_COMM_WORLD,&rq[2]);
	MPI_Isend(pgr,4*gr,MPI_DOUBLE,rn,10*rank+rn,MPI_COMM_WORLD,&rq[3]);
	MPI_Waitall(4,rq,sta);

	//unpack
	for (i=0;i<(cfl+cfr);i++){
		int idx2 = 2*sidx,i4 = i*4;
		newrd[idx2] = recvrp[i4];
		newrd[idx2+1] = recvrp[i4+1];
		newpd[idx2] = recvrp[i4+2];
		newpd[idx2+1] = recvrp[i4+3];
		sidx += 1;
	}
	delete [] recvrp;
	delete [] pgl;
	delete [] pgr;

	//update new R and P
	delete [] r_d;
	delete [] p_d;
	delete [] f_d;
	r_d = newrd;
	p_d = newpd;
	for (i=0;i<numlocal;i++){
		printf("%d@%d: newr = %7.3f\n",i,rank,r_d[i*2]);
	}
	numlocal = numlocal_new;
	f_d = new double[2*numlocal];
	ForceCalculation_MPI_3box(r_d,f_d,numlocal);
	//forwardV
	for (i=0; i<numlocal; i++){
		for (j=0; j<dim; j++){
			p_d[i*dim+j] += f_d[i*dim+j]*hdt;
		}
	}
	printf("*****before VV exit%11.3f*********\n",r_d[0]);
}

