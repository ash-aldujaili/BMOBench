#include <math.h>
#include "mex.h"
#define printf mexPrintf

/*
    returns the logical Pareto membership of the last point among a set of points.

    synopsis:  isnondominated = isnondominated(objMat)

    created by Abdullah Al-Dujail: ash.aldujaili@gmail.com
    
    adapted from Yi Cao: y.cao@cranfield.ac.uk
    
    for compiling type 

    mex isnondominated.c
    
*/


void isnondominated(bool *front, double * M, unsigned int row, unsigned int col);

void mexFunction( int nlhs, mxArray *plhs[] , int nrhs, const mxArray *prhs[] )
{
	bool *front;
	double * M;
	unsigned int row, col;
	const int  *dims;
    
	if(nrhs == 0 || nlhs > 1)
	{
	    printf("\nsynopsis:   front = paretofront(X)");
	    plhs[0]    = mxCreateDoubleMatrix(0 , 0 ,  mxREAL);
	    return;
	}
	
	M = mxGetPr(prhs[0]);
	dims = mxGetDimensions(prhs[0]);
	row = dims[0];
	col = dims[1];
	
	
	
	/* ----- output ----- */

	plhs[0]    = mxCreateLogicalMatrix (1 , 1);
	front = (bool *) mxGetPr(plhs[0]);
	
	
	/* main call */
	isnondominated(front,  M, row, col);
}

void isnondominated(bool * front, double * M, unsigned int row, unsigned int col)
{
    unsigned int i,j;
    unsigned int gCount, sCount;
    
    front[0] = true;
    for(i= 0; i < row - 1; i++){
        gCount = 0;
        sCount = 0;
        for(j = 0; j < col; j++){
            if (M[i+j*row] < M[(row-1)+j*row]) gCount++;
            else sCount++;
            /*printf("comparing %d column of %d element :%f vs %f\n",i,j, M[i+j*row] ,M[(row-1)+j*row]);*/
            if (sCount > 0 && gCount > 0) break;/*nondominated, skip*/
        }
        if (sCount > 0 && gCount > 0) continue;
        if (gCount > 0) {
            /*printf("hi");*/
            front[0]= false;
            break;
        }
    }
}