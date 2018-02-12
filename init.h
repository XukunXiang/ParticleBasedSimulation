#ifndef __INIT_H_INCLUDED__
#define	__INIT_H_INCLUDED__
#include "constants.h"

double getrand();

void Init_R(double R[][dim]);
void Init_R_FCC(double R[][3]);

void Init_P(double P[][dim]);
void Init_P_random(double P[][dim]);

double rescale(double P[][dim]);
double rescale_single(double P[][dim]);
double getT(double P[][dim]);

#endif
