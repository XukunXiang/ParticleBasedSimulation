#ifndef __INTEGRATOR_H_INCLUDED__
#define	__INTEGRATOR_H_INCLUDED__
#include "constants.h"

void VelocityVerlet(double R[][dim],double P[][dim],double F[][dim],double dt);
void VelocityVerlet_MPI_3box(double *r_d,double *p_d,double *f_d,double dt,int &numlocal);
//void VelocityVerlet_MPI_3box(double r_d[],double p_d[],double f_d[],double dt,int &numlocal);

#endif
