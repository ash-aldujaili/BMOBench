#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "bmobenchProblems.h"
#include <stdbool.h>
//#include "global.h"
/* #include "paretofiltering.h"
#include "myarray.h"
#include "wrapbbob.h"
#include "myrandom.h" */
#include <unistd.h>
#include <windows.h>
// define


#define UNIFRAND (((double) rand())/ RAND_MAX)

// MODIRECT algorithm
void MORANDOM(size_t nVars, size_t nObjs, size_t maxFuncEvals, double *loBound, double* upBound) {


//  char isVerbose = 1;
	double *objVector;
	double *pointVector;

	// allocate
	objVector = malloc(sizeof(double) * nObjs );
	pointVector = malloc(sizeof(double) * nVars );

	int i;
	for(i=0; i<nObjs; ++i)	objVector[i] = 0;
	for(i=0; i<nVars; ++i)	pointVector[i] = 0;

 // sample and Evaluate
 for(size_t i=0; i < maxFuncEvals; i++){
   // sample
   for(size_t j=0; j < nVars; j++) pointVector[j] = UNIFRAND * (upBound[j]-loBound[j]) + loBound[j];
   // evaluate
   pEvaluate(pointVector, objVector);

	 //
	 //if (isVerbose){
		//  printf("y-vector:");
		//  for (size_t j=0; j < nObjs- 1; j ++) printf("%f,", objVector[j]);
		//  printf("%f).\n", objVector[nObjs - 1]);
		//  printf("x-vector:");
		//  for (size_t j=0; j < nVars- 1; j ++) printf("%f,", pointVector[j]);
		//  printf("%f).\n", pointVector[nVars -1]);

	 //}
 }

	// free memory:
	free(objVector);
	free(pointVector);

	//array1Dp_destruct(&pPoints);
	//array1Dp_destruct(&pObjs);

}
