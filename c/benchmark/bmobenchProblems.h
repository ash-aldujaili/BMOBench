
#ifndef _bmobProblems_H
#define _bmobProblems_H
#include <stdlib.h>

/*
Defines
*/

#define NUM_PROBLEMS 100
#define STR_LEN 1024
#define MAX_DIM 30
/*
Variables
*/
struct problemsInformation {
	size_t nVars;
	size_t nObjs;
	double lowerBounds[MAX_DIM];
	double upperBounds[MAX_DIM];
	char *name;
};

typedef struct problemsInformation problemsInfo;
problemsInfo pInfo[NUM_PROBLEMS];
size_t mopCountFEvals;

/*
Problems
*/

/* function calls to each of the 100 problems */
void p1(double* x, double* y);
void p2(double* x, double* y);
void p3(double* x, double* y);
void p4(double* x, double* y);
void p5(double* x, double* y);
void p6(double* x, double* y);
void p7(double* x, double* y);
void p8(double* x, double* y);
void p9(double* x, double* y);
void p10(double* x, double* y);
void p11(double* x, double* y);
void p12(double* x, double* y);
void p13(double* x, double* y);
void p14(double* x, double* y);
void p15(double* x, double* y);
void p16(double* x, double* y);
void p17(double* x, double* y);
void p18(double* x, double* y);
void p19(double* x, double* y);
void p20(double* x, double* y);
void p21(double* x, double* y);
void p22(double* x, double* y);
void p23(double* x, double* y);
void p24(double* x, double* y);
void p25(double* x, double* y);
void p26(double* x, double* y);
void p27(double* x, double* y);
void p28(double* x, double* y);
void p29(double* x, double* y);
void p30(double* x, double* y);
void p31(double* x, double* y);
void p32(double* x, double* y);
void p33(double* x, double* y);
void p34(double* x, double* y);
void p35(double* x, double* y);
void p36(double* x, double* y);
void p37(double* x, double* y);
void p38(double* x, double* y);
void p39(double* x, double* y);
void p40(double* x, double* y);
void p41(double* x, double* y);
void p42(double* x, double* y);
void p43(double* x, double* y);
void p44(double* x, double* y);
void p45(double* x, double* y);
void p46(double* x, double* y);
void p47(double* x, double* y);
void p48(double* x, double* y);
void p49(double* x, double* y);
void p50(double* x, double* y);
void p51(double* x, double* y);
void p52(double* x, double* y);
void p53(double* x, double* y);
void p54(double* x, double* y);
void p55(double* x, double* y);
void p56(double* x, double* y);
void p57(double* x, double* y);
void p58(double* x, double* y);
void p59(double* x, double* y);
void p60(double* x, double* y);
void p61(double* x, double* y);
void p62(double* x, double* y);
void p63(double* x, double* y);
void p64(double* x, double* y);
void p65(double* x, double* y);
void p66(double* x, double* y);
void p67(double* x, double* y);
void p68(double* x, double* y);
void p69(double* x, double* y);
void p70(double* x, double* y);
void p71(double* x, double* y);
void p72(double* x, double* y);
void p73(double* x, double* y);
void p74(double* x, double* y);
void p75(double* x, double* y);
void p76(double* x, double* y);
void p77(double* x, double* y);
void p78(double* x, double* y);
void p79(double* x, double* y);
void p80(double* x, double* y);
void p81(double* x, double* y);
void p82(double* x, double* y);
void p83(double* x, double* y);
void p84(double* x, double* y);
void p85(double* x, double* y);
void p86(double* x, double* y);
void p87(double* x, double* y);
void p88(double* x, double* y);
void p89(double* x, double* y);
void p90(double* x, double* y);
void p91(double* x, double* y);
void p92(double* x, double* y);
void p93(double* x, double* y);
void p94(double* x, double* y);
void p95(double* x, double* y);
void p96(double* x, double* y);
void p97(double* x, double* y);
void p98(double* x, double* y);
void p99(double* x, double* y);
void p100(double* x, double* y);



/* problem handlers */


typedef void (*bmobProblem)(double *, double *);
bmobProblem handles[NUM_PROBLEMS]; /* for a better indexing of the problems rather than ifelse*/

/*
Problem evaluation call
*/
void pEvaluate( double *x, double *y);



#endif
