#ifndef MATC_H
#define MATC_H
// #include "problems.h"
/*
 Populate your function declarations here

*/

/*
*   Test function spherical function
*	x : decision vector of length n
*	y  : objective vector of length m
*	n : length of x
*	m : length of y
*/
extern void CL1(double* f, double* x);
extern void Deb512b(double*f, double* x);
extern void Deb521b(double*f, double* x);
extern void DTLZ1n2(double* f, double* x);
extern void DTLZ3n2(double* f, double* x);
extern void DTLZ5n2(double* f, double* x);
extern void Far1(double* f, double* x);
extern void Fonseca(double* f, double* x);
extern void I4(double* f, double* z);
extern void Jin1(double* f, double* x);
extern void Kursawe(double* f, double* x);
extern void L1ZDT4(double* f, double* x);
extern void L3ZDT1(double* f, double* x);
extern void L3ZDT2(double* f, double* x);
extern void L3ZDT3(double* f, double* x);
extern void L3ZDT4(double* f, double* x);
extern void L3ZDT6(double* f, double* x);
extern void VFM1(double* f, double* x);
extern void VU1(double* f, double* x);
extern void VU2(double* f, double* x);
extern void ZDT1(double* f, double* x);
extern void ZDT2(double* f, double* x);
extern void ZDT3(double* f, double* x);
extern void ZDT4(double* f, double* x);
extern void ZDT6(double* f, double* x);
extern void ZLT1(double* f, double* x);
extern void Deb512a(double* f, double* x);
extern void Deb521a(double *f, double *x);
extern void DTLZ1(double *f, double *x);
extern void DTLZ3(double *f, double *x);
extern void DTLZ5(double *f, double *x);
extern void ex005(double *f, double *x);
extern void FES3(double *f, double *x);
extern void I2(double *f, double *z);
extern void I3(double *f, double *z);
extern void IM1(double *f, double *x);
extern void Jin3(double *f, double *x);
extern void Jin4(double *f, double *x);
extern void L2ZDT1(double *f, double *x);
extern void L2ZDT2(double *f, double *x);
extern void L2ZDT4(double *f, double *x);
extern void L2ZDT6(double *f, double *x);
extern void lovison1(double *f, double *x);
extern void lovison3(double *f, double *x);
extern void lovison5(double *f, double *x);
extern void OKA1(double *f, double *x);
extern void OKA2(double *f, double *x);
extern void Sch1(double *f, double *x);
extern void SK1(double *f, double *x);
extern void SP1(double *f, double *x);
extern void SSFYY2(double *f, double *x);
extern void TKLY1(double *f, double *x);
extern void WFG6(double* f, double* zz);
extern void WFG7(double* f, double* zz);
extern void WFG8(double* f, double* zz);
extern void Deb53(double *f, double *x);
extern void BK1(double* f, double* x);
extern void Deb41(double* f, double* x);
extern void Deb512c(double*f, double* x);
extern void DG01(double*f, double* x);
extern void DTLZ2(double* f, double* x);
extern void DTLZ4(double* f, double* x);
extern void DTLZ6(double* f, double* x);
extern void FES1(double* f, double* x);
extern void I1(double* f, double* x);
extern void I5(double* f, double* z);
extern void Jin2(double* f, double* x);
extern void L2ZDT3(double* f, double* x);
extern void LE1(double* f, double* x);
extern void lovison2(double* f, double* x);
extern void lovison4(double* f, double* x);
extern void lovison6(double* f, double* x);
extern void MOP2(double* f, double* x);
extern void MOP6(double* f, double* x);
extern void QV1(double* f, double* x);
extern void SK2(double* f, double* x);
extern void SSFYY1(double* f, double* x);
extern void WFG3(double* f, double* x);
extern void MOP1(double* f, double* x);
extern void MOP3(double* f, double* x);
extern void MOP4(double* f, double* x);
extern void MOP5(double* f, double* x);
extern void MOP7(double* f, double* x);
extern void MLF1(double* f, double* x);
extern void MLF2(double* f, double* x);
extern void Deb513(double *f, double *x);
extern void DPAM1(double *f, double *x);
extern void DTLZ2n2(double *f, double *x);
extern void DTLZ4n2(double *f, double *x);
extern void DTLZ6n2(double *f, double *x);
extern void FES2(double *f, double *x);
extern void LRS1(double *f, double *x);
extern void MHHM1(double *f, double *x);
extern void MHHM2(double *f, double *x);
extern void WFG1(double* f, double* zz);
extern void WFG2(double* f, double* zz);
extern void WFG4(double* f, double* zz);
extern void WFG5(double* f, double* zz);
extern void IKK1(double *f, double *x);
extern void WFG9(double* f, double* zz);

#endif