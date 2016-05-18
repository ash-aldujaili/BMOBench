#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# This code was given by Thanh-Do Tran, Dimo Brockhoff, INRIA

from __future__ import division
import numpy as np


def compute(dataSet, refSet, method="additive"):
    """All objectives are subject to minimization.
    
    dataSet: Is a 2D numpy array (nData x nObj) containing the output of an MO run
             Should be nondominated points.
    
    refSet:  Is a 2D numpy array (nRef x nObj) containing the reference set
             with which the indicator is calculated.
    
    method:  Specifies the indicator type to use: "additive" or "multiplicative".
             Default is "additive" (i.e. the additive epsilon indicator).
    """
    
    if method == "additive":
        value = -np.inf
        additiveEpsilon = True
    elif method == "multiplicative":
        value = 0
        additiveEpsilon = False
    else:
        print "Wrong 'method' parameter!"
        quit()
    
    nData, nObj = dataSet.shape
    nRef, dimSet = refSet.shape
    if nObj < 1 or nObj != dimSet:
        print "Number of objectives mismatched."
        print "Dimensions of dataSet and refSet must be equal and > 1"
        quit()
    
    # For any interpretation of the epsilon indicator,
    # the outer loop must be over the points in refSet.
    for i in xrange(nRef):
        for j in xrange(nData):
            eps_j = epsilonDominance(dataSet[j,:], refSet[i,:], nObj, additiveEpsilon)
            if j == 0:
                eps_i = eps_j   # epsilon from dataSet to point i in refSet
            elif eps_i > eps_j: # takes the SMALLEST epsilon over all points of dataSet to point i in refSet
                eps_i = eps_j
        
        if i == 0:
            value = eps_i   # epsilon from dataSet to refSet
        elif value < eps_i: # takes the LARGEST epsilon over all points of refSet
            value = eps_i   # thereby when transformed by a factor ${value}, the whole refSet is weakly dominated by dataSet
    
    return value


def epsilonDominance(dataPoint, refPoint, nObj, additive=True):
    for k in xrange(nObj):
        if additive: # additive epsilon
            eps_k = dataPoint[k] - refPoint[k] # distance in this k-th objective
        else: # multiplicative
            if (refPoint[k] < 0 and dataPoint[k] > 0) or \
               (refPoint[k] > 0 and dataPoint[k] < 0) or \
               (refPoint[k] == 0 and dataPoint[k] == 0):
                print "dataPoint and refPoint have to be > 0"
                quit()
            eps_k = dataPoint[k] / refPoint[k]
        if k == 0:
            epsilon = eps_k   # epsilon from point j in dataSet to point i in refSet
        elif epsilon < eps_k: # takes the LARGEST over all the objectives
            epsilon = eps_k
    return epsilon


if __name__ == "__main__":
    print "I am a module to be called by another program!"
    raw_input()
