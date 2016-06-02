

/*
* matc.c - MATLAB External Interfaces to 100mops
*
* Multiplies an input scalar (multiplier)
* times a 1xN matrix (inMatrix)
* and outputs a 1xN matrix (outMatrix)
*
* The calling syntax is:
*
*       outMatrix = arrayProduct(multiplier, inMatrix)
*
* This is a MEX file for MATLAB.
*/

#include "mex.h"
#include <string.h>
// #include "problems.h"
#include "matc.h"
#include <stdlib.h>
#include <math.h>
#define NUM_PROBLEMS 100
#define STR_LEN 1024
#define MAX_DIM 30

struct problemsInformation {
    size_t nVars;
    size_t nObjs;
    double lowerBounds[MAX_DIM];
    double upperBounds[MAX_DIM];
    char *name;
};

typedef struct problemsInformation problemsInfo;
typedef void (*bmobProblem)(double *, double *);

bmobProblem handles[NUM_PROBLEMS]={&CL1,&Deb512b,&Deb521b,&DTLZ1n2,&DTLZ3n2,&DTLZ5n2,&Far1,&Fonseca,&I4,&Jin1,&Jin3,&Kursawe,&L1ZDT4,&L3ZDT1,&L3ZDT2,&L3ZDT3,&L3ZDT4,&L3ZDT6,
  &VFM1,&VU1,&VU2,&ZDT1,&ZDT2,&ZDT3,&ZDT4,&ZDT6,&ZLT1,&Deb512a,&Deb521a,&Deb53,&DTLZ1,&DTLZ3,&DTLZ5,&ex005,&FES3,&I2,&I3,&IM1,&Jin4,&L2ZDT1,&L2ZDT2,&L2ZDT4,&L2ZDT6,&lovison1,
  &lovison3,&lovison5,&OKA1,&OKA2,&Sch1,&SK1,&SP1,&SSFYY2,&TKLY1,&WFG6,&WFG7,&WFG8,&BK1,&Deb41,&Deb512c,&DG01,&DTLZ2,&DTLZ4,&DTLZ6,&FES1,&I1,&I5,&L2ZDT3,&Jin2,&LE1,&lovison2,
  &lovison4,&lovison6,&MOP2,&MOP6,&QV1,&SK2,&SSFYY1,&WFG3,&MOP1,&MOP3,&MOP4,&MOP5,&MOP7,&MLF1,&MLF2,&Deb513,&DPAM1,&DTLZ2n2,&DTLZ4n2,&DTLZ6n2,&FES2,&LRS1,&MHHM1,&MHHM2,
  &WFG1,&WFG2,&WFG4,&WFG5,&WFG9,&IKK1}; /* for a better indexing of the problems rather than ifelse*/

// problemsInfo getInfo(char *problem, problemsInfo pInfo);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
  int nrhs, const mxArray *prhs[])
  {


    problemsInfo pInfo[NUM_PROBLEMS] =  {
      { 4, 2, {1.0, sqrt(2), sqrt(2), 1}, {3.0, 3.0, 3.0, 3.0}, "CL1"},
      { 2, 2, {0,0}, {1.0,1.0}, "Deb512b"},
      { 2, 2, {0,0}, {1.0,1.0}, "Deb521b"},
      { 2, 2, {0,0}, {1.0,1.0}, "DTLZ1n2"},
      { 2, 2, {0,0}, {1.0,1.0}, "DTLZ3n2"},
      { 2, 2, {0,0}, {1.0,1.0}, "DTLZ5n2"},
      { 2, 2, {-1.0,-1.0}, {1.0,1.0}, "Far1"},
      { 2, 2, {-4.0,-4.0}, {4.0,4.0}, "Fonseca"},
      { 8, 3, {0,0,0,0,0,0,0,0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,}, "I4"},
      { 2, 2, {0,0}, {1.0,1.0}, "Jin1"},
      { 2, 2, {0.0,0.0}, {1.0,1.0}, "Jin3"},
      { 3, 2, {-5.0,-5.0,-5.0}, {5.0,5.0,5.0}, "Kursawe"},
      { 10, 2, {0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0}, {1.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0}, "L1ZDT4"},
      { 30, 2, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "L3ZDT1"},
      { 30, 2, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "L3ZDT2"},
      { 30, 2, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "L3ZDT3"},
      { 30, 2, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "L3ZDT4"},
      { 10, 2, {0,0,0,0,0,0,0,0,0,0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "L3ZDT6"},
      { 2, 3, {-2.0,-2.0}, {2.0,2.0}, "VFM1"},
      { 2, 2, {-3.0,-3.0}, {3.0,3.0}, "VU1"},
      { 2, 2, {-3.0,-3.0}, {3.0,3.0}, "VU2"},
      { 30, 2, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "ZDT1"},
      { 30, 2, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "ZDT2"},
      { 30, 2, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "ZDT3"},
      { 10, 2, {0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0}, {1.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0}, "ZDT4"},
      { 10, 2, {0,0,0,0,0,0,0,0,0,0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "ZDT6"},
      { 10, 3, {-1000.0,-1000.0,-1000.0,-1000.0,-1000.0,-1000.0,-1000.0,-1000.0,-1000.0,-1000.0}, {1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0}, "ZLT1"},
      { 2, 2, {0.0,0.0}, {1.0,1.0}, "Deb512a"},
      { 2, 2, {0.0,0.0}, {1.0,1.0}, "Deb521a"},
      { 2, 2, {0.0,0.0}, {1.0,1.0}, "Deb53"},
      { 7, 3, {0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "DTLZ1"},
      { 12, 3, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "DTLZ3"},
      { 12, 3, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "DTLZ5"},
      { 2, 2, {-1.0,1.0}, {2.0,2.0}, "ex005"},
      { 10, 4, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "FES3"},
      { 8, 3, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "I2"},
      { 8, 3, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "I3"},
      { 2, 2, {1.0,1.0}, {4.0,2.0}, "IM1"},
      { 2, 2, {0.0,0.0}, {1.0,1.0}, "Jin4"},
      { 30, 2, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "L2ZDT1"},
      { 30, 2, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "L2ZDT2"},
      { 30, 2, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "L2ZDT4"},
      { 10, 2, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "L2ZDT6"},
      { 2, 2, {0.0,0.0}, {3.0,3.0}, "lovison1"},
      { 2, 2, {0.0,-4.0}, {6.0,4.0}, "lovison3"},
      { 3, 3, {-1.0,-1.0,-1.0}, {4.0,4.0,4.0}, "lovison5"},
      { 2, 2, {1.55291427062,-1.62620802141}, {7.62200523018,5.79555495773}, "OKA1"},
      { 3, 2, {-4*atan(1),-5.0, -5.0}, {4*atan(1),5, 5}, "OKA2"},
      { 1, 2, {0.0}, {5.0}, "Sch1"},
      { 1, 2, {-10.0}, {10.0}, "SK1"},
      { 2, 2, {-1.0,-1.0}, {5.0,5.0}, "SP1"},
      { 1, 2, {-100.0}, {100.0}, "SSFYY2"},
      { 4, 2, {0.1,0.0,0.0,0.0}, {1.0,1.0,1.0,1.0}, "TKLY1"},
      { 8, 3, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0}, "WFG6"},
      { 8, 3, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0}, "WFG7"},
      { 8, 3, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0}, "WFG8"},
      { 2, 2, {-5.0,-5.0}, {10.0,10.0}, "BK1"},
      { 2, 2, {0.1, 0.0}, {1.0,1.0}, "Deb41"},
      { 2, 2, {0.0,0.0}, {1.0,1.0}, "Deb512c"},
      { 1, 2, {-10.0}, {13.0}, "DG01"},
      { 12, 3, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "DTLZ2"},
      { 12, 3, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "DTLZ4"},
      { 22, 3, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "DTLZ6"},
      { 10, 2, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "FES1"},
      { 8, 3, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "I1"},
      { 8, 3, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "I5"},
      { 30, 2, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "L2ZDT3"},
      { 2, 2, {0.0,0.0}, {1.0,1.0}, "Jin2"},
      { 2, 2, {-5.0,-5.0}, {10.0,10.0}, "LE1"},
      { 2, 2, {-0.5,-0.5}, {0.0,0.5}, "lovison2"},
      { 2, 2, {0.0,-1.0}, {6.0,1.0}, "lovison4"},
      { 3, 3, {-1.0,-1.0, -1.0}, {4.0,4.0,4.0}, "lovison6"},
      { 4, 2, {-4.0,-4.0,-4.0,-4.0}, {4.0,4.0,4.0,4.0}, "MOP2"},
      { 2, 2, {0.0,0.0}, {1.0,1.0}, "MOP6"},
      { 10, 2, {-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12}, {5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12}, "QV1"},
      { 4, 2, {-10.0,-10.0,-10.0,-10.0}, {10.0,10.0,10.0,10.0}, "SK2"},
      { 2, 2, {-100.0,-100.0}, {100.0,100.0}, "SSFYY1"},
      { 8, 3, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0}, "WFG3"},
      { 1, 2, {-100000.0}, {100000.0}, "MOP1"},
      { 2, 2, {-4*atan(1),-4*atan(1)}, {4*atan(1),4*atan(1)}, "MOP3"},
      { 3, 2, {-5.0,-5.0,-5.0}, {5.0,5.0,5.0}, "MOP4"},
      { 2, 3, {-30.0,-30.0}, {30.0,30.0}, "MOP5"},
      { 2, 3, {-400.0,-400.0}, {400.0,400.0}, "MOP7"},
      { 1, 2, {0.0}, {20.0}, "MLF1"},
      { 2, 2, {-2.0,-2.0}, {2.0,2.0}, "MLF2"},
      { 2, 2, {0.0,0.0}, {1.0,1.0}, "Deb513"},
      { 10, 2, {-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3}, {0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3}, "DPAM1"},
      { 2, 2, {0.0,0.0}, {1.0,1.0}, "DTLZ2n2"},
      { 2, 2, {0.0,0.0}, {1.0,1.0}, "DTLZ4n2"},
      { 2, 2, {0.0,0.0}, {1.0,1.0}, "DTLZ6n2"},
      { 10, 3, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, "FES2"},
      { 2, 2, {-50.0,-50.0}, {50.0,50.0}, "LRS1"},
      { 1, 3, {0.0}, {1.0}, "MHHM1"},
      { 2, 3, {0.0, 0.0}, {1.0, 1.0}, "MHHM2"},
      { 8, 3, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0}, "WFG1"},
      { 8, 3, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0}, "WFG2"},
      { 8, 3, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0}, "WFG4"},
      { 8, 3, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0}, "WFG5"},
      { 8, 3, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0}, "WFG9"},
      { 2, 3, {-50.0,-50.0}, {50.0,50.0}, "IKK1"}
    };
    // char *problem;    /* input problem id */
    int problemID;
    double *x, *inp;          /* 1xn input vector */
    size_t n;           /* size of vector */
    int len, i;
    problemsInfo info;
    /* TO DO:
    We should have, besides functions to get objective values, other functions that return
    the bounds and the number of variables/objectives per problems so the algorithm knows
    what to send and what to receive.
    So for example n and m below should be known ahead by the algorithm accessing the problem (e.g., BK1)
    Here I am just assuming the algorithm knows ahead
    */
    size_t m;           /* TO DO: this should depend on the problem, size of output vector */
    double *y, *m_inp, *l_inp, *u_inp, *x_inp;      /* output vector */
    // checking arguments
    if(nrhs > 2 || nlhs != 5)
    {
      printf("\nSyntax:   [x, l, u, m, y] = matc(problem, x)");
      plhs[0]    = mxCreateDoubleMatrix(0 , 0 ,  mxREAL);
      return;
    }

    if(nrhs > 2) {
      mexErrMsgIdAndTxt("matc:nrhs",
      "Two inputs required.");
    }

    if(nlhs != 5) {
      mexErrMsgIdAndTxt("matc:nlhs",
      "Five output required.");
    }

    /* make sure the first input argument is string */
    /* input must be a string */
    if ( mxIsDouble(prhs[0]) != 1)
    mexErrMsgIdAndTxt( "MATLAB:revord:inputNotInteger","problem must be an integer.");


    /* Read Input */
    /* copy the string data from prhs[0] into a C string problem.    */
    // problem = mxArrayToString(prhs[0]);
    inp = mxGetPr(prhs[0]);
    problemID = (int)(inp[0]-1);

    // for(i = 0; i < NUM_PROBLEMS; i++){
    //     if(strcmp(pInfo[i].name, problem) == 0){
    //       info =  pInfo[i];
    //       break;
    //     }
    // }
    // if(info == NULL)
    //   return;
    info = pInfo[problemID];
    m = info.nObjs;

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL); /* x */
    plhs[1] = mxCreateDoubleMatrix(1,info.nVars,mxREAL); /* l */
    plhs[2] = mxCreateDoubleMatrix(1,info.nVars,mxREAL); /* u */
    plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL); /* m */
    plhs[4] = mxCreateDoubleMatrix(1,m,mxREAL); /* y */
    /* get a pointer to the real data in the output matrix */
    x_inp = mxGetPr(plhs[0]);
    l_inp = mxGetPr(plhs[1]);
    u_inp = mxGetPr(plhs[2]);
    m_inp = mxGetPr(plhs[3]);
    y = mxGetPr(plhs[4]);

    x_inp[0] = info.nVars;
    m_inp[0] = m;
    for(i = 0; i < info.nVars; i++){
        l_inp[i] = info.lowerBounds[i];
        u_inp[i] = info.upperBounds[i];
    }

    if(nrhs == 1){
        return;
    }
    
    if( !mxIsDouble(prhs[1]) ||
    mxIsComplex(prhs[1])) {
      mexErrMsgIdAndTxt("matc:notDouble",
      "x vector must be type double.");
    }
    /* create a pointer to the real data in the input matrix  */
    x = mxGetPr(prhs[1]);

    /* get dimensions of the input matrix */
    n = mxGetN(prhs[1]);
    
    /* call the computational routine */
    (*handles[problemID])(y,x);

    // if(strcmp("CL1",problem)==0){
    //   CL1(y, x);
    // }
    // else if(strcmp("Deb512b",problem)==0){
    //   Deb512b(y, x);
    // }
    // else if(strcmp("Deb521b",problem)==0){
    //   Deb521b(y, x);
    // }
    // else if(strcmp("DTLZ1n2",problem)==0){
    //   DTLZ1n2(y, x);
    // }
    // else if(strcmp("DTLZ3n2",problem)==0){
    //   DTLZ3n2(y, x);
    // }
    // else if(strcmp("DTLZ5n2",problem)==0){
    //   DTLZ5n2(y, x);
    // }
    // else if(strcmp("Far1",problem)==0){
    //   Far1(y, x);
    // }
    // else if(strcmp("Fonseca",problem)==0){
    //   Fonseca(y, x);
    // }
    // else if(strcmp("I4",problem)==0){
    //   I4(y, x);
    // }
    // else if(strcmp("Jin1",problem)==0){
    //   Jin1(y, x);
    // }
    // else if(strcmp("Jin3",problem)==0){
    //   Jin3(y, x);
    // }
    // else if(strcmp("Kursawe",problem)==0){
    //   Kursawe(y, x);
    // }
    // else if(strcmp("L1ZDT4",problem)==0){
    //   L1ZDT4(y, x);
    // }
    // else if(strcmp("L3ZDT1",problem)==0){
    //   L3ZDT1(y, x);
    // }
    // else if(strcmp("L3ZDT2",problem)==0){
    //   L3ZDT2(y, x);
    // }
    // else if(strcmp("L3ZDT3",problem)==0){
    //   L3ZDT3(y, x);
    // }
    // else if(strcmp("L3ZDT4",problem)==0){
    //   L3ZDT4(y, x);
    // }
    // else if(strcmp("L3ZDT6",problem)==0){
    //   L3ZDT6(y, x);
    // }
    // else if(strcmp("VFM1",problem)==0){
    //   VFM1(y, x);
    // }
    // else if(strcmp("VU1",problem)==0){
    //   VU1(y, x);
    // }
    // else if(strcmp("VU2",problem)==0){
    //   VU2(y, x);
    // }
    // else if(strcmp("ZDT1",problem)==0){
    //   ZDT1(y, x);
    // }
    // else if(strcmp("ZDT2",problem)==0){
    //   ZDT2(y, x);
    // }
    // else if(strcmp("ZDT3",problem)==0){
    //   ZDT3(y, x);
    // }
    // else if(strcmp("ZDT4",problem)==0){
    //   ZDT4(y, x);
    // }
    // else if(strcmp("ZDT6",problem)==0){
    //   ZDT6(y, x);
    // }
    // else if(strcmp("ZLT1",problem)==0){
    //   ZLT1(y, x);
    // }
    // else if(strcmp("BK1", problem) == 0){
    //   BK1(y, x);
    // }else if(strcmp("Deb41", problem) == 0){
    //   Deb41(y, x);
    // }else if(strcmp("Deb512c", problem) == 0){
    //   Deb512c(y, x);
    // }else if(strcmp("DG01", problem) == 0){
    //   DG01(y, x);
    // }else if(strcmp("FES1", problem) == 0){
    //   FES1(y, x);
    // }else if(strcmp("Jin2", problem) == 0){
    //   Jin2(y, x);
    // }else if(strcmp("LE1", problem) == 0){
    //   LE1(y, x);
    // }else if(strcmp("lovison2", problem) == 0){
    //   lovison2(y, x);
    // }else if(strcmp("lovison4", problem) == 0){
    //   lovison4(y, x);
    // }else if(strcmp("MOP2", problem) == 0){
    //   MOP2(y, x);
    // }else if(strcmp("MOP6", problem) == 0){
    //   MOP6(y, x);
    // }else if(strcmp("SK2", problem) == 0){
    //   SK2(y, x);
    // }else if(strcmp("MOP1", problem) == 0){
    //   MOP1(y, x);
    // }else if(strcmp("MOP3", problem) == 0){
    //   MOP3(y, x);
    // }else if(strcmp("MOP4", problem) == 0){
    //   MOP4(y, x);
    // }else if(strcmp("MOP5", problem) == 0){
    //   MOP5(y, x);
    // }else if(strcmp("MOP7", problem) == 0){
    //   MOP7(y, x);
    // }else if(strcmp("MLF1", problem) == 0){
    //   MLF1(y, x);
    // }else if(strcmp("MLF2", problem) == 0){
    //   MLF2(y, x);
    // }else if(strcmp("I1", problem) == 0){
    //   I1(y, x);
    // }else if(strcmp("WFG3", problem) == 0){
    //   WFG3(y, x);
    // }else if(strcmp("DTLZ2", problem) == 0){
    //   DTLZ2(y, x, n, m);
    // }else if(strcmp("DTLZ4", problem) == 0){
    //   DTLZ4(y, x, n, m);
    // }else if(strcmp("DTLZ6", problem) == 0){
    //   DTLZ6(y, x, n, m);
    // }else if(strcmp("lovison6", problem) == 0){
    //   lovison6(y, x);
    // }else if(strcmp("L2ZDT3", problem) == 0){
    //   L2ZDT3(y, x);
    // }else if(strcmp("I5", problem) == 0){
    //   I5(y, x);
    // }
    // else if (strcmp("Deb512a",problem)==0)
    // Deb512a(y,x);
    // else if (strcmp("Deb521a",problem)==0)
    // Deb521a(y,x);
    // else if (strcmp("DTLZ1",problem)==0)
    // {
    //   DTLZ1(y,x);
    // }
    // else if (strcmp("DTLZ3",problem)==0)
    // {
    //   DTLZ3(y,x);
    // }
    // else if (strcmp("DTLZ5",problem)==0)
    // {
    //   DTLZ5(y,x);
    // }
    // else if (strcmp("ex005",problem)==0)
    // {
    //   ex005(y,x);
    // }
    // else if (strcmp("FES3",problem)==0)
    // {
    //   FES3(y,x);
    // }
    // else if (strcmp("I2",problem)==0)
    // {
    //   I2(y,x);
    // }
    // else if (strcmp("I3",problem)==0)
    // {
    //   I3(y,x);
    // }
    // else if (strcmp("IM1",problem)==0)
    // {
    //   IM1(y,x);
    // }
    // else if (strcmp("Jin4",problem)==0)
    // {
    //   Jin4(y,x);
    // }
    // else if (strcmp("L2ZDT1",problem)==0)
    // {
    //   L2ZDT1(y,x);
    // }
    // else if (strcmp("L2ZDT2",problem)==0)
    // {
    //   L2ZDT2(y,x);
    // }
    // else if (strcmp("L2ZDT4",problem)==0)
    // {
    //   L2ZDT4(y,x);
    // }
    // else if (strcmp("L2ZDT6",problem)==0)
    // {
    //   L2ZDT6(y,x);
    // }
    // else if (strcmp("lovison1",problem)==0)
    // {
    //   lovison1(y,x);
    // }
    // else if (strcmp("lovison3",problem)==0)
    // {
    //   lovison3(y,x);
    // }
    // else if (strcmp("lovison5",problem)==0)
    // {
    //   lovison5(y,x);
    // }
    // else if (strcmp("OKA1",problem)==0)
    // {
    //   OKA1(y,x);
    // }
    // else if (strcmp("OKA2",problem)==0)
    // {
    //   OKA2(y,x);
    // }
    // else if (strcmp("QV1",problem)==0)
    // {
    //   QV1(y,x);
    // }
    // else if (strcmp("Sch1",problem)==0)
    // {
    //   Sch1(y,x);
    // }
    // else if (strcmp("SK1",problem)==0)
    // {
    //   SK1(y,x);
    // }
    // else if (strcmp("SP1",problem)==0)
    // {
    //   SP1(y,x);
    // }
    // else if (strcmp("SSFYY1",problem)==0)
    // {
    //   SSFYY1(y,x);
    // }
    // else if (strcmp("SSFYY2",problem)==0)
    // {
    //   SSFYY2(y,x);
    // }
    // else if (strcmp("TKLY1",problem)==0)
    // {
    //   TKLY1(y,x);
    // }
    // else if (strcmp("Deb53",problem)==0)
    // {
    //   Deb53(y,x);
    // }
    // else if (strcmp("WFG6",problem)==0)
    // {
    //   WFG6(y,x);
    // }
    // else if (strcmp("WFG7",problem)==0)
    // {
    //   WFG7(y,x);
    // }
    // else if (strcmp("WFG8",problem)==0)
    // {
    //   WFG8(y,x);
    // }
    // else if (strcmp("WFG9",problem)==0)
    // {
    //   WFG9(y,x);
    // }
    // else if (strcmp("Deb513",problem)==0)
    // {
    //   Deb513(y,x);
    // }
    // else if (strcmp("DPAM1",problem)==0)
    // {
    //   DPAM1(y,x);
    // }
    // else if (strcmp("DTLZ2n2",problem)==0)
    // {
    //   DTLZ2n2(y,x);
    // }
    // else if (strcmp("DTLZ4n2",problem)==0)
    // {
    //   DTLZ4n2(y,x);
    // }
    // else if (strcmp("DTLZ6n2",problem)==0)
    // {
    //   DTLZ6n2(y,x);
    // }
    // else if (strcmp("FES2",problem)==0)
    // {
    //   FES2(y,x);
    // }
    // else if (strcmp("LRS1",problem)==0)
    // {
    //   LRS1(y,x);
    // }
    // else if (strcmp("MHHM1",problem)==0)
    // {
    //   MHHM1(y,x);
    // }
    // else if (strcmp("MHHM2",problem)==0)
    // {
    //   MHHM2(y,x);
    // }
    // else if (strcmp("WFG1",problem)==0)
    // {
    //   WFG1(y,x);
    // }
    // else if (strcmp("WFG2",problem)==0)
    // {
    //   WFG2(y,x);
    // }
    // else if (strcmp("WFG4",problem)==0)
    // {
    //   WFG4(y,x);
    // }
    // else if (strcmp("WFG5",problem)==0)
    // {
    //   WFG5(y,x);
    // }
    // else if (strcmp("IKK1",problem)==0)
    // {
    //   IKK1(y,x);
    // }
    // else {
    //   printf("No such problem");
    //   plhs[4]    = mxCreateDoubleMatrix(0 , 0 ,  mxREAL);
    //   return;
    // }
    // len = sizeof(y) / sizeof(y[0]); 
    // for(i = len-1; i >= 0 ; i--){
    //   y[i + 4] = y[i];
    // }

    // y[0] = *info.nVars;
    // y[1] = *info.lowerBounds;
    // y[2] = *info.upperBounds;
    // y[3] = m;
    // return;


  }

// problemsInfo getInfo(char *problem, problemsInfo *pInfo){
//     int i;
//     for(i = 0; i < NUM_PROBLEMS; i++){
//     if(strcmp(pInfo[i].name, problem) == 0)
//       return pInfo[i];
//     }
//     printf("ERROR : Problem not found\n");
//     problemsInfo p;
//     return p;
// }

