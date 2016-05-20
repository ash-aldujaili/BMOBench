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


def ecdf_extracting(ind = 'hv', isMean = True, data = []):
    """
    ind: specifies the indicator used
    isMean : true to compute the normalized data profile over multiple runs else the best artificial profile
    data : a list of dicts, one dict per algorithm describing it.
    """
    #====================================#
    # Read info from BBOB
    #====================================#
    #FUNC_LIST = range(1,25) # from f1 to f24
    # DO NOT CHANGE THIS !!
    PROBLEM_LIST = ['BK1','ex005','Deb41','Deb512a','Deb512b','Deb512c','Deb513','Deb521a','Deb521b','Deb53','ZDT1','ZDT2','ZDT3','ZDT4','ZDT6','DTLZ1','DTLZ2','DTLZ3','DTLZ4','DTLZ5','DTLZ6','DTLZ1n2','DTLZ2n2','DTLZ3n2','DTLZ4n2','DTLZ5n2','DTLZ6n2','Kursawe','Fonseca','L1ZDT4','L2ZDT1','L2ZDT2','L2ZDT3','L2ZDT4','L2ZDT6','L3ZDT1','L3ZDT2','L3ZDT3','L3ZDT4','L3ZDT6','WFG1','WFG2','WFG3','WFG4','WFG5','WFG6','WFG7','WFG8','WFG9','I1','I2','I3','I4','I5','MOP1','MOP2','MOP3','MOP4','MOP5','MOP6','MOP7','DPAM1','DG01','Far1','FES1','FES2','FES3','IKK1','Jin1','Jin2','Jin3','Jin4','OKA1','OKA2','LRS1','IM1','LE1','MHHM1','MHHM2','MLF1','MLF2','QV1','Sch1','SP1','SSFYY1','SSFYY2','SK1','SK2','TKLY1','VU1','VU2','VFM1','ZLT1','CL1','lovison1','lovison2','lovison3','lovison4','lovison5','lovison6']
    DIM_LIST = [2,2,2,2,2,2,2,2,2,2,30,30,30,10,10,7,12,12,12,12,22,2,2,2,2,2,2,3,2,10,30,30,30,30,10,30,30,30,30,10,8,8,8,8,8,8,8,8,8,8,8,8,8,8,1,4,2,3,2,2,2,10,1,2,10,10,10,2,2,2,2,2,2,3,2,2,2,1,2,1,2,10,1,2,2,1,1,4,4,2,2,2,10,4,2,2,2,2,3,3]
    #DIM = dimension
    #xlim = [-5, 5]
    #INST_LIST = ((2,4), (3,5), (7,8), (9,10), (11,12))
    NUM_RUNS = 10
    NUM_TARGETS = 70
    FILENAME_REFERENCE_SET = '%s_ref_set.prtun'
    FILENAME_HV_REF = 'hv_reference.txt'
    FILENAME_ROI = '%s.BOUND'
   
   
    FILENAME_PROFILE_DATABASE = 'db_runprofiles_%s.pkl' % ind
    FILENAME_ECDF_DATABASE = "db_ecdfs_%s_maxfevalDiff_%druns_%dtargets.pkl" % (ind, NUM_RUNS, NUM_TARGETS)  
        
    readDatabase  = False  # set to False to avoid reading PKL and thus must process raw TXT data
    normalizedMOP = True
    
    #======= TO BE EDITED =====================================================
    #OBJ1_INSTANCE = 1
    #OBJ2_INSTANCE = 5
    INPUT_TXT_DIR = "../../EXP_RESULTS"
    REF_TXT_DIR = "../../problems/reference-sets"
    ROI_TXT_DIR = "../../problems/roi"
    
    if data == None:
      data = []
      data.append({
        'name': "MO-SOO",
        'algo': "MO-SOO",
        'dir' : INPUT_TXT_DIR,
        'nRun': 1,
        'nFev': 1e3,
      })
      data.append({
        'name': "MO-DIRECT",
        'algo': "MO-DIRECT",
        'dir' : INPUT_TXT_DIR,
        'nRun': 1,
        'nFev': 1e3,
      })   	
    #==========================================================================
    objcol = 1
    timecol = 0

    # read hv reference
    hv_ref_all_problems = np.genfromtxt(FILENAME_HV_REF) 
    #=============================================#
    #  Loop over MOPs #
    #=============================================#
    for pidx, pid in enumerate(PROBLEM_LIST):
        DIM = DIM_LIST[pidx]
        if readDatabase:
            # Read once to load the whole database for all combinations, avoid re-reading
            try: all_profiles
            except NameError:
                import cPickle as pickle
                f = open(FILENAME_PROFILE_DATABASE, 'rb')
                all_profiles = pickle.load(f)
                f.close()
            
        # Otherwise CREATE the database of run profiles
        else:
            #buffer = np.array([]) # help not dirty the next run!!!
            for ialg, alg in enumerate([d['algo'] for d in data]):
                print("Processing %s's data for %s" % (alg, pid))
                # Define an identifier of algorithm in the produced database
                alg_id = '%d: %s' % (ialg, alg)
                
                #for iinst, instance in enumerate(INST_LIST):
                #obj1_instance = instance[0]
                #obj2_instance = instance[1]
                
                if ind == 'hv':
                  hv_ref = hv_ref_all_problems[pidx]
                elif ind == 'eps' or ind == 'igd' or ind == 'gd':
                  fname = FILENAME_REFERENCE_SET % (pid)
                  ref_fname = REF_TXT_DIR + "/" + fname
                  ref = np.genfromtxt(ref_fname)
                
                # Working on the given combined BBOB functions
                #mop, _, _ = mobbob.setInstances(obj1_id, obj1_instance, obj2_id, obj2_instance, DIM, xlim, verbal=True)
                
                # Set reference point for HV computation
                fname = FILENAME_ROI % (pid)
                roi_fname = ROI_TXT_DIR + "/" + fname
                roi = np.genfromtxt(roi_fname)
                nadir = roi[1,:]
                ideal = roi[0,:]

                normfactor = nadir - ideal # for normalizing the objective vector wrt ROI
                if ind == 'hv':
                  referencePoint = list(nadir)
                  #idealPoint = list(ideal)
                  if normalizedMOP:
                      referencePoint = list(nadir/normfactor)
                      #idealPoint = list(ideal/normfactor)
                elif ind == 'eps' or ind == 'igd' or ind == 'gd': # add the two optimal two the reference point
                  if normalizedMOP:
                      ref = ref/ normfactor[None,:] #np.vstack((op1/ normfactor[None,:],op2/ normfactor[None,:], eps_ref))
                
                isHit = False
                for irun in xrange(data[ialg]['nRun']):
                    # Check existence of input TXT file
                    fname = "%s_%dD_%s_nfev%.1e_run%d.txt" % (
                             pid, DIM,
                             data[ialg]['name'], data[ialg]['nFev'], irun+1)
                    inputFile = data[ialg]['dir'] + "/" + fname
                    
                    if not os.path.isfile(inputFile):
                        print "No input file <%s>" % (fname)
                        continue
                    
                    # Read the data from text file
                    #print(fname)
                    buf = readformatedtxt.readtxt(filename=inputFile, ichunk='all')
                    
                    # In case the input text file has no entry!!!: a work around in case all
                    # are empy record the nadir point 
                    #if len(buf.shape) < 2: continue
                    if len(buf.shape) < 2:
                      print("Empty file")
                      buf = np.hstack((nadir,np.array([1]), np.array([0.00])))
                      buf = buf[None,]
                      
                    #isHit = True
                    # Normalizing all points by nadir coordinates
                    normObjVecs = buf[:,objcol:] / normfactor[None,:] if normalizedMOP else buf[:,objcol:]
                    #pfFlag, nPareto, hyperV = ParetoUpdateWithHV2obj.execute(normObjVecs, referencePoint) # commented to use on-the-fly HV
                    
                    # Which of the 2 below is used may produce slightly
                    # difference profile plots due to numerical issue
                    # concerning on-the-fly versus offline approaches!!!
                    #hv_diff = hv_ref[idx_obj1_id, idx_obj2_id] - hyperV
                    if ind == 'hv':
                      print("Computing the HV indicator")
                      hvfname = "%s_%dD_%s_nfev%.1e_run%d_hv.txt" % (
                             pid, DIM,
                             data[ialg]['name'], data[ialg]['nFev'], irun+1)
                      inputHVFile = data[ialg]['dir'] + "/" + hvfname
                      hv_profile = np.genfromtxt(inputHVFile)
                      #print len(hv_profile.shape)
                      if len(hv_profile.shape) == 1:
                          hv_profile = hv_profile[None,]
                      val = hv_ref - hv_profile[:,-1]# more precise
                    elif ind == 'eps':
                      print("Computing the Epsilon indicator")
                      val = compute_fast_incr_eps(normObjVecs, ref)
                    elif ind == 'igd':
                      print("Computing the IGD indicator")
                      val = compute_incr_igd(normObjVecs, ref)
                    elif ind == 'gd':
                      print("Computing the GD indicator")
                      val = compute_incr_gd(normObjVecs, ref)
                    # Take the time stamp, which is the last column
                    if ind == 'hv': # timestamp for hv depends on the hv file
                        timestamp = hv_profile[:,0].astype(int)
                    else:
                        timestamp = buf[:,timecol].astype(int)
                    print("Updating the profiles")
                    #-------------------------------------------------------
                    # Store profiles (timestamp and epsilon) to a big dictionary
                    # Init the database if not existing
                    try: all_profiles
                    except NameError:
                        all_profiles = {}
                    
                    # Init the dict for storing everything belonging to this obj1
                    try: all_profiles[pidx]
                    except KeyError:
                        #all_profiles.update({obj1_id: {}})
                        all_profiles[pidx] = {}
                    
                    # Init the sub-dict for storing everything belonging to this obj1 & obj2
                    #try: all_profiles[obj1_id][obj2_id]
                    #except KeyError:
                    #    #all_profiles[obj1_id].update({obj2_id: {}})
                    #    all_profiles[obj1_id][obj2_id] = {}
                    
                    # Init the sub-sub-dict for storing everything belonging to this obj1 & obj2 & alg
                    try: all_profiles[pidx][alg_id]
                    except KeyError:
                        #all_profiles[obj1_id][obj2_id].update({alg_id: []})
                        all_profiles[pidx][alg_id] = []
                    
                    # Init the data structure for storing the core data of
                    # an alg on a combination and add core data to the database
                    #if ind == 'hv':
                    #print("Appending a new hv-based profile")
                    # TODO address the dimension scaling with respect to each problem as LCM(DIMS) is not a suitable way (currently scaling down)
                    all_profiles[pidx][alg_id].append(            {'value'    : list(val),
                                                                   'timestamp': list(timestamp//DIM),#'timestamp': list(np.ceil(timestamp/DIM).astype(int)),#'timestamp': list(timestamp//DIM),
                                                                   'nfevs'    : int(data[ialg]['nFev'])})
                    # elif ind == 'eps':
                    #   print("Appending a new eps-based profile")
                    #   all_profiles[pidx][alg_id].append({'value'    : list(eps),
                    #                                                'timestamp': list(timestamp),
                    #                                                'nfevs'    : int(data[ialg]['nFev'] * DIM)})
                                                               
                    #print all_profiles
                    #-------------------------------------------------------
                    print "Finished run %d of %s on problem %s" % (irun+1, data[ialg]['algo'], pid)
            #END OF for ialg
                                
        #-------------------------------------------------------------------
        # Target values to use
        try: list_of_target_values
        except NameError:
            if ind == 'hv':
              # Make 70 targets, linearly on the log-scale from 10^-0.1 to 10^-7
              arr = np.linspace(start=-0.8, stop=-3, num=NUM_TARGETS, endpoint=True)
              list_of_target_values = 10**arr
            elif ind == 'eps':
              arr = np.linspace(start=-0.8, stop=-2, num=NUM_TARGETS, endpoint=True)
              list_of_target_values = 10**arr
            elif ind == 'igd' or ind == 'gd':
              arr = np.linspace(start=-0.8, stop=-3, num=NUM_TARGETS, endpoint=True)
              list_of_target_values = 10**arr          
              
        #-------------------------------------------------------------------
        
        # Extract raw ECDFs from run profiles of algorithms on a combination
        rawECDFs, percECDFs = ecdf.compute_multi_algos(all_profiles[pidx], list_of_target_values, isMean = isMean)
        #rawECDFs = percECDFs
        # Attach the configuration of the algorithm runs to the raw database
        for ialg, alg_id in enumerate(all_profiles[pidx].keys()):
            rawECDFs[alg_id].update({'info': data[ialg]})
            # In this script, we don't need to store all profiles (which is totally huge)
            # thus we can remove the ECDF-extracted profiles just to release memory
            all_profiles[pidx][alg_id] = []
        
        # Store the raw ECDFs data to a big dictionary as gallery of all 300 combinations
        try: all_ecdfs
        except NameError:
            all_ecdfs = {'targets': list_of_target_values}
        try: all_ecdfs[pidx]
        except KeyError:
            all_ecdfs.update({pidx: {}})
        
        all_ecdfs[pidx].update(rawECDFs)
        #-------------------------------------------------------------------
    #END of loop over problems
    
    # Save the data to a pickle file
    import cPickle as pickle
    #outputfile = FILENAME_ECDF_DATABASE
    f = open(FILENAME_ECDF_DATABASE, 'wb')
    pickle.dump(all_ecdfs, f)
    f.close()
    
    # Copy the database to INPUT_TXT_DIR for archival
    #shutil.copy2(outputfile, INPUT_TXT_DIR)



def main(ind = 'eps'):
    # MPI stuff TO DO LATER
    #nproc = MPI.COMM_WORLD.Get_size()
    #iproc = MPI.COMM_WORLD.Get_rank()
    
    #for idx, d in enumerate(dims):
    #  if idx % nproc == iproc:
    ecdf_extracting( ind = ind, isMean =  True)
        
    #MPI.COMM_WORLD.Barrier()




if __name__ == "__main__":
    startTime = time.time()
    if len(sys.argv) > 1:
      main(ind = sys.argv[1])
    else:
      print "syntax: python main_extract_ecdfy.py hv|eps|gd|igd"
    print "It took:", time.time() - startTime, "seconds."
    print "Done." 
    # print "Press Enter to continue ..." 
    # raw_input()
