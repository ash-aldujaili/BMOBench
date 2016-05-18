// By : Abdullah Al-Dujaili adapted from http://www.ntu.edu.sg/home/epnsugan/
// to compute the IGD and GD quality indicator incrementally as well as single scalar
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#define max(x,y) (((x)>(y))? (x) : (y))
#define min(x,y) (((x)<(y))? (x) : (y))

/*
* Scalar IGD
* returns IGD(A,B)
*/
double igd(double* A, double *B, size_t numelA, size_t  numelB, size_t nObjs)
{
	int i, j, k;
	double d,min,dis;
	
	
	/* Step 2: calculate the IGD value for feasible points */
	dis=0.0;
	for( i=0; i<numelB; i++ )
	{
		min = 1.0E200;
		for( j=0; j< numelA; j++ )
		{
			d=0.0;
			for( k=0; k<nObjs; k++ )
				d += ( B[i*nObjs+k] - A[j*nObjs+k] ) * ( B[i*nObjs+k] - A[j*nObjs+k] );
			if( d < min ) min = d;
		}
		dis += sqrt(min);
	}
	return dis/(double)(numelB);
}

/*
* Scalar GD
* returns GD(A,B)
*/
double gd(double* A, double *B, size_t numelA, size_t  numelB, size_t nObjs)
{
	int i, j, k;
	double d,min,dis;
	
	dis=0.0;
	for( i=0; i<numelA; i++ )
	{
		min = 1.0E200;
		for( j=0; j< numelB; j++ )
		{
			d=0.0;
			for( k=0; k<nObjs; k++ )
				d += ( A[i*nObjs+k] - B[j*nObjs+k] ) * ( A[i*nObjs+k] - B[j*nObjs+k] );
			if( d < min ) min = d;
		}
		dis += sqrt(min);
	}
	return dis/(double)(numelA);
}

/*
* Incremental GD
* returns GD(A,B) over incremental subsets of A with respect to the whole set B,
* this function is faster but at the expense of additional space
*/
void incr_gd(double *gd_i, double* A, double *B, size_t numelA, size_t  numelB, size_t nObjs)
{
	int i, j, k;
	double d,min,dis;
	
	dis=0.0;
	for( i=0; i<numelA; i++ )
	{
		min = 1.0E200;
		for( j=0; j< numelB; j++ )
		{
			d=0.0;
			for( k=0; k<nObjs; k++ )
				d += ( A[i*nObjs+k] - B[j*nObjs+k] ) * ( A[i*nObjs+k] - B[j*nObjs+k] );
			if( d < min ) min = d;
		}
		dis += sqrt(min);
		gd_i[i] = dis / (double) (i+1);
	}
	// return void
	return;
}

/*
* Incremental IGD
* returns IGD(A,B) over incremental subsets of A with respect to the whole set B,
* this function is faster but at the expense of additional space
*/
void incr_igd(double *igd_i, double* A, double *B, size_t numelA, size_t  numelB, size_t nObjs)
{
	int i, j, k;
	double d,min,dis;

	dis=0.0;
	for( i=0; i<numelB; i++ )
	{
		min = 1.0E200;
		for( j=0; j< numelA; j++ )
		{
			d=0.0;
			for( k=0; k<nObjs; k++ )
				d += ( B[i*nObjs+k] - A[j*nObjs+k] ) * ( B[i*nObjs+k] - A[j*nObjs+k] );
			if( d < min ) min = d;
			igd_i[j] +=sqrt(min);
		}
	}
	// postprocessing code
	for( i=0; i<numelA; i++ )
		igd_i[i] = igd_i[i]/ (double) numelB;
	return;
}