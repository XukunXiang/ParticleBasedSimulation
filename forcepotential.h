#ifndef __FORCEPOTENTIAL_H_INCLUDED__
#define	__FORCEPOTENTIAL_H_INCLUDED__
#include "constants.h"

void ForceCalculation(double R[][dim],double F[][dim]);
void ForceCalculation_MPI_3box(double r_d[],double f_d[],int numlocal);
double getPE(double R[][dim]);

#endif
