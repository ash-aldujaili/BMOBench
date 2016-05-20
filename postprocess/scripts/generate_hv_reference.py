#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division
import os
import sys # to access command-line arguments
import shutil # to delete a folder
import time
import numpy as np

#----------------------------------------------
import matplotlib
# Use a non-interactive backend such as Agg (for PNGs), PDF, SVG, or PS.
matplotlib.use('Agg') # make sure to call this before pyplot
#----------------------------------------------
import matplotlib.pyplot as plt
from matplotlib import cm
import csv

import bbobbenchmarks as bm
import wrapbbob as mobbob # set instances, combine bbob functions, get optima
import readformatedtxt # read my txt file into numpy array
#import readhvref # read hv_ref_all_DIMs.txt
# non-dominance operator (c-wrapper)
#import paretofront as pf
#from paretofront import paretofront
import pymoutils
from pymoutils import *
# MPI stuff
#from mpi4py import MPI
import ecdf


def generate_hv_reference():
    """
    ind: specifies the indicator used
    isMean : true to compute the normalized data profile over multiple runs else the best artificial profile
    data : a list of dicts, one dict per algorithm describing it.
    """
    #====================================#
    # Read info from BBOB
    #====================================#
    # DO NOT CHANGE THIS !!
    PROBLEM_LIST = ['BK1','ex005','Deb41','Deb512a','Deb512b','Deb512c','Deb513','Deb521a','Deb521b','Deb53','ZDT1','ZDT2','ZDT3','ZDT4','ZDT6','DTLZ1','DTLZ2','DTLZ3','DTLZ4','DTLZ5','DTLZ6','DTLZ1n2','DTLZ2n2','DTLZ3n2','DTLZ4n2','DTLZ5n2','DTLZ6n2','Kursawe','Fonseca','L1ZDT4','L2ZDT1','L2ZDT2','L2ZDT3','L2ZDT4','L2ZDT6','L3ZDT1','L3ZDT2','L3ZDT3','L3ZDT4','L3ZDT6','WFG1','WFG2','WFG3','WFG4','WFG5','WFG6','WFG7','WFG8','WFG9','I1','I2','I3','I4','I5','MOP1','MOP2','MOP3','MOP4','MOP5','MOP6','MOP7','DPAM1','DG01','Far1','FES1','FES2','FES3','IKK1','Jin1','Jin2','Jin3','Jin4','OKA1','OKA2','LRS1','IM1','LE1','MHHM1','MHHM2','MLF1','MLF2','QV1','Sch1','SP1','SSFYY1','SSFYY2','SK1','SK2','TKLY1','VU1','VU2','VFM1','ZLT1','CL1','lovison1','lovison2','lovison3','lovison4','lovison5','lovison6']
  
    FILENAME_REFERENCE_SET = '%s_ref_set.prtun'
    FILENAME_ROI = '%s.BOUND'
    normalizedMOP = True
    
    #======= TO BE EDITED =====================================================
    REF_TXT_DIR = "../../problems/reference-sets"
    ROI_TXT_DIR = "../../problems/roi" 	
    #==========================================================================
    
    #=============================================#
    #  Loop over MOPs #
    #=============================================#
    hv_ref = np.zeros((len(PROBLEM_LIST),1))
    for pidx, pid in enumerate(PROBLEM_LIST):

        # Set reference point for HV computation
        fname = FILENAME_ROI % (pid)
        roi_fname = ROI_TXT_DIR + "/" + fname
        roi = np.genfromtxt(roi_fname)
        nadir = roi[1,:]
        ideal = roi[0,:]
        
        normfactor = nadir - ideal # for normalizing the objective vector wrt ROI
        # read the reference set
        fname = FILENAME_REFERENCE_SET % (pid)
        ref_fname = REF_TXT_DIR + "/" + fname
        ref = np.genfromtxt(ref_fname)
        
        # normalize:
        if normalizedMOP:
            referencePoint = list(nadir/normfactor)
        if normalizedMOP:
            ref = ref/ normfactor[None,:] 
        
        hv_ref[pidx]= compute_pyhv(ref, referencePoint)


    np.savetxt("hv_reference.txt", hv_ref, fmt='%f' )
  

def main():
    generate_hv_reference()
        




if __name__ == "__main__":
    startTime = time.time()
	main()
    print "It took:", time.time() - startTime, "seconds."
    print "Done." 
    # print "Press Enter to continue ..." 
    # raw_input()
