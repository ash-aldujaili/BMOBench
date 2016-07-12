/*From MOBBOB-GECCO, Brockhoff, Dimo, Thanh-Do Tran, and Nikolaus Hansen. "Benchmarking numerical multiobjective optimizers revisited." Proceedings of the 2015 Annual Conference on Genetic and Evolutionary Computation. ACM, 2015.*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h> // to use the bool datatype, required C99

#include "checkDominance.h"


/**
 * Unconstrained non-domination checking for minimization
 * Input: pointer to the objective vector of each solution
 * Ouput: 1 if point1 dominates point2
 *        0 if point1 and point2 are non-dominated
 *       -1 if point2 dominates point1
 *       -2 if point1 is identical to point1
 */
int checkDominance(double *point1Obj, double *point2Obj, size_t nObjs) {
    bool flag1 = false;
    bool flag2 = false;
    for (size_t i=0; i < nObjs; i++) {
        if (point1Obj[i] < point2Obj[i]) {
            flag1 = true;
        } else if (point1Obj[i] > point2Obj[i]) {
            flag2 = true;
        }
    }
    
    if (flag1 && !flag2) {
        return 1;
    } else if (!flag1 && flag2) {
        return -1;
    } else {  /* (flag1 && flag2) || (!flag1 && !flag2) */
        for (size_t i=0; i < nObjs; i++) {
            if (point1Obj[i] != point2Obj[i])
                return 0;
        }
        return -2; // the two points are identical
    }
}

