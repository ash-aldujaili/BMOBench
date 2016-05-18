// paretofront returns the logical Pareto membership of a set of points
// synopsis:  frontFlag = paretofront(objMat)
// This file has two functions:
// pf_selective: by Abdullah Al-Dujaili, allows to pass an array and process its elements selectively, i.e pass an array of flag that tells which element to compare among.
// pf_cao: Created by Yi Cao: y.cao@cranfield.ac.uk, does not allow for a selective processing of elements within the passed array

#include <stdio.h>
#include <stdlib.h>  // memory, e.g. malloc
#include <stdbool.h> // to use the bool datatype, required C99
#include <math.h>


void pf_selective(bool *frontFlag, double *obj, size_t nPoints, size_t nObjs) {
    //bool *checklist, colDominatedFlag;
    size_t nPosDiff;
    size_t nNegDiff;
	double delta = 0.0;
    //checklist = (bool*)malloc(nrow*sizeof(bool));
    // initialize
    //for (size_t p = 0; p < nPoints; p++) frontFlag[p] = true;
    
    for (size_t p = 0; p < nPoints - 1; p++) {
      if (!frontFlag[p]) continue; // skip this point if not of interest
      // loop over other points
      for (size_t q = p + 1; q < nPoints; q++) {
        // skip this point if it is dominated or of no interest
        if (!frontFlag[q]) continue;
        // otherwise compare their objective-differences
        nPosDiff = 0;
        nNegDiff = 0;
        for (size_t j = 0; j < nObjs; j++) {
			delta = obj[p * nObjs + j] - obj[q * nObjs + j];
			if (delta > 0) nNegDiff = nNegDiff + 1;
			else //if (delta < 0)commenting this to filter out identical items similar to cao's method
				nPosDiff = nPosDiff + 1;
		    // speed up step
		    if (nNegDiff > 0 && nPosDiff > 0) break; // incomparable
        }
        if (nNegDiff > 0 && nPosDiff > 0) continue; // incomparable
        else if (nNegDiff > 0) { // p is dominated
          frontFlag[p] = false;
          break; 
        }
        else if (nPosDiff > 0) frontFlag[q] = false; // q is dominated
        
      }
    }   
    //free(checklist); 
}


void pf_cao(bool *frontFlag, double *obj, unsigned nrow, unsigned ncol) {
    unsigned t, s, i, j, j1, j2;
    bool *checklist, colDominatedFlag;
    
    checklist = (bool*)malloc(nrow*sizeof(bool));
    
    for(t=0; t<nrow; t++)
        checklist[t] = true;
    for(s=0; s<nrow; s++) {
        t = s;
        if (!checklist[t])
            continue;
        checklist[t] = false;
        colDominatedFlag = true;
        for(i=t+1; i<nrow; i++) {
            if (!checklist[i])
                continue;
            checklist[i] = false;
            for (j=0,j1=i,j2=t; j<ncol; j++,j1+=nrow,j2+=nrow) {
                if (obj[j1] < obj[j2]) {
                    checklist[i] = true;
                    break;
                }
            }
            if (!checklist[i])
                continue;
            colDominatedFlag = false;
            for (j=0,j1=i,j2=t; j<ncol; j++,j1+=nrow,j2+=nrow) {
                if (obj[j1] > obj[j2]) {
                    colDominatedFlag = true;
                    break;
                }
            }
            if (!colDominatedFlag) { //swap active index continue checking
                frontFlag[t] = false;
                checklist[i] = false;
                colDominatedFlag = true;
                t = i;
            }
        }
        frontFlag[t] = colDominatedFlag;
        if (t>s) {
            for (i=s+1; i<t; i++) {
                if (!checklist[i])
                    continue;
                checklist[i] = false;
                for (j=0,j1=i,j2=t; j<ncol; j++,j1+=nrow,j2+=nrow) {
                    if (obj[j1] < obj[j2]) {
                        checklist[i] = true;
                        break;
                    }
                }
            }
        }
    }
    free(checklist); 
}