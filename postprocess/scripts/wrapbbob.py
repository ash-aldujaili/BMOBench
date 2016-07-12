#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# From MOBBOB-GECCO, Brockhoff, Dimo, Thanh-Do Tran, and Nikolaus Hansen. "Benchmarking numerical multiobjective optimizers revisited." Proceedings of the 2015 Annual Conference on Genetic and Evolutionary Computation. ACM, 2015.*/

from __future__ import division
import numpy as np

#import fgeneric
import bbobbenchmarks as bm


class MOBBOB:
    """Scalarizing multiobjective problem
    """
    def __init__(self, obj1_id, obj1_instance, obj2_id, obj2_instance, dimension):
        # Query the optimum from the benchmark
        # -------------------------------------
        self.obj1, self.opt1_obj1 = bm.instantiate(obj1_id, iinstance=obj1_instance)
        self.obj2, self.opt2_obj2 = bm.instantiate(obj2_id, iinstance=obj2_instance)
        # Activate the evaluation system to get the true optimum from the self object
        self.fcurrent = np.zeros(2) # store function values
        self.fcurrent[0] = self.obj1.evaluate(np.zeros((1, dimension)))
        self.fcurrent[1] = self.obj2.evaluate(np.zeros((1, dimension)))
        self.opt1_x = self.obj1.xopt
        self.opt2_x = self.obj2.xopt
        self.opt1_obj2 = self.obj2.evaluate(self.opt1_x)
        self.opt2_obj1 = self.obj1.evaluate(self.opt2_x)
    

def setInstances(obj1_id, obj1_instance, obj2_id, obj2_instance, dimension, axislim=[-5,5], verbal=True):
    ### obj1_instance = 1
    if obj2_id < obj1_id: print "obj2_id=%d < obj1_id=%d" % (obj2_id, obj1_id) # start from obj1_id onwards
    ### if obj2_id == obj1_id: obj2_instance = 5
    ### if obj2_id != obj1_id: obj2_instance = 5
    '''
    if obj2_id == obj1_id:
        obj2_instance = 1
        if obj1_id == 20:
            obj2_instance = 2
    if obj1_id==17 and obj1_instance==0 and obj2_id==18 and obj2_instance==0:
        obj2_instance = 1
    '''
    
    mop = MOBBOB(obj1_id, obj1_instance, obj2_id, obj2_instance, dimension)
    
    """
    # Query the optimum from the benchmark and try another instance if needed
    for t in xrange(10): # max 10 trials for picking obj2_instance such that obj1 != obj2
        obj2_instance = obj2_instance + t;
        mop = MOBBOB(obj1_id, obj1_instance, obj2_id, obj2_instance, dimension)
        # Two optima must be distinguished in search space AND
        # distinguished in obj space AND not vertically nor horizontally equivalent
        #if ( np.any(mop.opt1_x != mop.opt2_x) and ((mop.opt1_obj1 != mop.opt2_obj1) and (mop.opt1_obj2 != mop.opt2_obj2)) ):
        if ( (not np.array_equal(mop.opt1_x, mop.opt2_x)) and ((mop.opt1_obj1 != mop.opt2_obj1) and (mop.opt1_obj2 != mop.opt2_obj2)) ):
            break # safe to proceed as obj1 != obj2
        print "WARNING: F%d.%d and F%d.%d (%dD) have the same optimum in search/objective space. Instance changed." % \
                        (obj1_id, obj1_instance, obj2_id, obj2_instance, dimension)
        if t == 9: # failed on the last trial as well
            print "WARNING: Two objectives are the same across several instances!"
    """
    
    print "\n***********************************************************************"
    print "Working with F%d.%d and F%d.%d in [%d,%d]^%d space." % (obj1_id, obj1_instance, obj2_id, obj2_instance, axislim[0], axislim[1], dimension)
    print "opt1_obj = (%.2f, %.2f)" % (mop.opt1_obj1, mop.opt1_obj2)
    print "     at opt1_x =", mop.opt1_x
    print "opt2_obj = (%.2f, %.2f)" % (mop.opt2_obj1, mop.opt2_obj2)
    print "     at opt2_x =", mop.opt2_x
    print "***********************************************************************"
    
    return mop, obj1_instance, obj2_instance
    


if __name__ == "__main__":
    print "I am a module and should only be called from another program!" 
    # raw_input()
