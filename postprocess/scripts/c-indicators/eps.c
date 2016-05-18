// By : Abdullah Al-Dujaili
// Different functionality for computing the epsilon indicator (Ziztler et al. 2000)

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#define max(x,y) (((x)>(y))? (x) : (y))
#define min(x,y) (((x)<(y))? (x) : (y))

// returns I_{\epsilon+}(A,B) over incremental subsets of A with respect to the whole set B,
// this function is faster but at the expense of additional space
void fast_incr_eps(double* eps_i, double* A, double *B, size_t numelA, size_t  numelB, size_t nObjs, char isAdditive) {
    
      double *eps_j = (double*) malloc(numelA*sizeof(double));
      
      // pre-processing
      for (size_t j = 0; j < numelA; j++) 
		if (isAdditive)  
			eps_i[j] = 0.0;
		else
			eps_i[j] = -INFINITY;

      // iterate over the reference set
      for(size_t i = 0; i < numelB; i++) {
        for (size_t j = 0; j < numelA; j++) eps_j[j] = INFINITY;
    	  // iterate over the approximation set
      	for(size_t j = 0; j < numelA; j++) {
    	      double eps_k = 0.0;
    		    for(size_t k = 0; k < nObjs; k++) eps_k = max(eps_k, A[j * nObjs + k] - B[i * nObjs + k]);
    		    if (j > 0) 
    		      eps_j[j] = min(eps_j[j-1], eps_k);
    		    else 
    		      eps_j[j] = min(eps_j[j], eps_k);
    	  }
    	  // post-processing
        for (size_t j = 0; j < numelA; j++) eps_i[j] = max(eps_i[j], eps_j[j]);

      } 
      
            // free
      free(eps_j);
}


// returns I_{\epsilon+}(A,B): a single number
double eps(double* A, double *B, size_t numelA, size_t  numelB, size_t nObjs) {
    
    double eps_i = 0.0;
    // iterate over the reference set
    for(size_t i = 0; i < numelB; i++) {
    	double eps_j = INFINITY;
    	// iterate over the approximation set
    	for(size_t j = 0; j < numelA; j++) {
    	    double eps_k = 0.0;
    		for(size_t k = 0; k < nObjs; k++) eps_k = max(eps_k, A[j * nObjs + k] - B[i * nObjs + k]);
    		eps_j = min(eps_j, eps_k);
    	}
    	eps_i = max(eps_i, eps_j);
    } 
    return eps_i;
}


// returns I_{\epsilon+}(A,B) over incremental subsets of A with respect to the whole set B,
void incr_eps(double* eps_i, double* A, double *B, size_t numelA, size_t  numelB, size_t nObjs) {
    
    for (size_t l = 0; l < numelA; l++) {
      eps_i[l] = 0.0;
      // iterate over the reference set
      for(size_t i = 0; i < numelB; i++) {
    	  double eps_j = INFINITY;
    	  // iterate over the approximation set
      	for(size_t j = 0; j < l+1; j++) {
    	      double eps_k = 0.0;
    		  for(size_t k = 0; k < nObjs; k++) eps_k = max(eps_k, A[j * nObjs + k] - B[i * nObjs + k]);
    		  eps_j = min(eps_j, eps_k);
    	  }
    	  eps_i[l] = max(eps_i[l], eps_j);
      } 
    }
}