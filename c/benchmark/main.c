#include <stdlib.h>
#include <stdio.h>
#include "globaldeclare.h"
#include "recorder.h"
#include "bmobenchProblems.h"



/* ===========================
include your algorithm files
=============================*/
#include "MORANDOM.h"
/* ===========================
For setting the number of runs and
the budget multiplier, refer to globaldeclare.h
=============================*/



/* To set the NUM_RUNS and BUDGET_MULTIPLIER, refer to globaldeclare.h */
int main(){
	// seed it
	srand(1111);
	// algorithm stuff
	char algName[STR_LEN] = "MORANDOM"; /* put your algorithm name here*/
	double *l; /* lower bounds */
	double *u; /* upper bounds */
	algoName = my_strdup(algName);
	// create the exp results folder if required
	creatFolder(outputDir);
	// loop over problems
	for(unsigned int p = 0; p < NUM_PROBLEMS ; p++){
		// get problems info
		nObjs = pInfo[p].nObjs;
		nVars = pInfo[p].nVars;
		problem_id = p;

		printf("%zu %s\n", problem_id, pInfo[p].name);

		maxFuncEvals = BUDGET_MULTIPLIER * nVars;
		l = (double*) malloc(sizeof(double) * nVars);
		u = (double*) malloc(sizeof(double) * nVars);
		for(unsigned i = 0; i < nVars; i++){
			l[i] = pInfo[p].lowerBounds[i];
			u[i] = pInfo[p].upperBounds[i];
		}
		/* loop over runs */
		for(unsigned int r = 0; r < NUM_RUNS ; r++){
			iRun = r + 1;
			/***************************************
			     put your algorithm call here
			***************************************/
			MORANDOM(nVars, nObjs, maxFuncEvals, l, u);
			// recorder stuff
			printArchive();
			freeArchive();
		}
		/* free stuff */
		free(l);
		free(u);
	}
	/* free algorithm name */
	free(algoName);
}
