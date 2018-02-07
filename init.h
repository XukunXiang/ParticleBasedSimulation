#ifndef __INIT_H_INCLUDED__
#define	__INIT_H_INCLUDED__

double getrand();

void Init_R(double R[][3]);
void Init_R_FCC(double R[][3]);

void Init_P(double P[][3]);
void Init_P_random(double P[][3]);

double rescale(double P[][3]);
double rescale_single(double P[][3]);
double getT(double P[][3]);

#endif
