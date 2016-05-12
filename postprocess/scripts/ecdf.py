#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division
import os
import sys # to access command-line arguments
import shutil # to delete a folder
import time
import numpy as np
from collections import defaultdict # for best profile


def compute_multi_algos(multi_algo_profiles, targets, isMean = True):
    """
    This function will apply the 'compute_single_algo' to the single_profile
    of each individual algorithm's data within the collective 'multi_algo_profiles'.
    The output includes a raw ECDF and a perc ECDF (with proportional data):
      - raw ECDF is the modest data that can be stored in a database.
      - perc ECDF adds 'proportion' that can be used to extract plottable data.
    isMean specifies whether the best data profile out of all the runs should be chosen or the normalized
      - isMean True computes the normalized
    """
    rawECDFs = {}
    percECDFs = {}
    for algo in multi_algo_profiles.keys():
        print algo
        rawECDFs[algo], percECDFs[algo] = compute_single_algo(multi_algo_profiles[algo], targets, isMean = isMean)
    return rawECDFs, percECDFs


def compute_single_algo(run_profiles, targets, rawECDFs=None, isMean = True):
    """
    + 'run_profiles': run profiles for a SINGLE ALGORITHM, e.g. hypervolume (or epsilon) indicator over time.
       Be a list of items, each is the data for 1 run. Each item is a dictionary with:
        - 'value'    : []
        - 'timestamp': []
        - 'nfevs'    : n
       The values for these 2 entries must be lists of equal length.
    + 'targets': a set of target values applied to the run profiles for extracting hitting times
    + 'rawECDFs': a dictionary, in which under key 'data' being a LIST of database records.
       Each record is a dictionary representing a target-hitting point, represented by:
        - ['at']   : the time, and
        - ['nhits']: the number of target being hit simultaneously at that time
    + isMean specifies whether the best data profile out of all the runs should be chosen or the normalized
      - isMean True computes the normalized. As the input profiles represent multiple runs over multiple 
        instances of the same problem, when isMean is set to False, the profile is computed by normalizing
        the best of each instance i.e. nprofiles = NUM_INST when isMean = false else nprofiles = NUM_INST *       NUM_RUNS

    """
    # If there is no input rawECDFs for updating,
    # initialize the DATA STRUCTURE for the output DATABASE
    if rawECDFs is None:
        rawECDFs = {}
        rawECDFs['nprofiles'] = 0  # typically the number of runs
        rawECDFs['ntargets']  = 0  # number of target values
        rawECDFs['nfevs']     = 0  # the run length
        rawECDFs['data']      = [] # main database is here: a list of dicts {'at': time, 'nhits': count}
        #rawECDFs['proportion'] = [] # with this entry, rawECDFs becomes percECDFs
        #rawECDFs['time']       = [] # with this entry, rawECDFs becomes percECDFs
    
    # Update information regarding the sizes of data
    NUM_INST = 1 # the report run_profiles are bundled for one problem i.e. NUM_INST * NUM_RUNS: here we aggregate over all of them either by mean or the best.
    NUM_RUNS = len(run_profiles) // NUM_INST
    #print(len(run_profiles))
    if isMean:
      rawECDFs['nprofiles'] += len(run_profiles)
    else: # normalize over one profile
      rawECDFs['nprofiles'] += NUM_INST
      
    rawECDFs['ntargets']   += len(targets)
    
    nfevs = max([profile['nfevs'] for profile in run_profiles])
    if nfevs > rawECDFs['nfevs']:
        rawECDFs['nfevs'] = nfevs
    
    isHit = False # a flag to indicate if any target has met to avoid the error of empty profile
    
    if isMean: # choose the normalized ecdf
      # Loop over the runs (items of 'run_profiles') to extract hitting times data
      for profile in run_profiles:
          itarget = 0  # reset the itarget for each run, i.e. restart from 1st target
          # Trace this run profile by going thru its points
          for time, value in zip(profile['timestamp'], profile['value']):
              # When a point is encounter to reach some target value:
              # "value <= targets[itarget]" is the major condition, but is preceded by
              # "itarget < len(targets)" to avoid index error when all targets are hit.
              if itarget < len(targets) and value <= targets[itarget]:
                  itarget_hit = itarget
                  count = 0
                  # Check with all targets from 'itarget' when there is a target hit
                  # to ensure the proper counting of multiple-target hitting
                  for target in targets[itarget_hit:]:
                      if value <= target:
                          count   += 1
                          itarget += 1
                      else:
                          break
                  isHit = True
                  rawECDFs['data'].append({'at'   : time,
                                         'nhits': count})
                  
           
    else: # choose the best profile:
          # combine all the ecdf in a tuple set (timestmap,value) joined over all the runs
          for inst_id in xrange(NUM_INST):
            time_value = []
            for profile in run_profiles[inst_id  * NUM_RUNS: (inst_id + 1) * NUM_RUNS ]:
                time_value = time_value + [x for x in zip(profile['timestamp'], profile['value'])]
            # join them 
            time_value_dict = defaultdict(list)
            for key, value in time_value:
                time_value_dict[key].append(value)
            # go over the targets
            itarget = 0  # reset the itarget for each run, i.e. restart from 1st target
            # Trace the best run of the algorithm at any timestamp
            #print(time_value_dict)
            for time, vals in time_value_dict.iteritems():
                value = min(vals)
                if itarget < len(targets) and value <= targets[itarget]:
                    itarget_hit = itarget
                    count = 0
                    # Check with all targets from 'itarget' when there is a target hit
                    # to ensure the proper counting of multiple-target hitting
                    for target in targets[itarget_hit:]:
                       if value <= target:
                            count   += 1 
                            itarget += 1
                       else:
                           break
                    isHit = True
                    rawECDFs['data'].append({'at'   : time,
                                         'nhits': count})
    
    
       
    # Sort the records in the list w.r.t. the key (or attribute) 'at' of each element (which is a dict)
    rawECDFs['data'] = sorted(rawECDFs['data'], key = lambda point: point['at'])
    
    # check if the data profile has zero targets
    if not isHit:
      rawECDFs['data'].append({'at'   : 0,'nhits': 0})
      rawECDFs['data'].append({'at'   : 1,'nhits': 0})
      
    #if (len(rawECDFs['data'])==0):
    #    print(len(rawECDFs['data']))
    # Add the first, the second (helps in semilogx), and the last elements
    if rawECDFs['data'][0]['at'] != 0:
        rawECDFs['data'].insert(0, {'at'   : 0,
                                    'nhits': 0})
    if rawECDFs['data'][1]['at'] != 1:
        rawECDFs['data'].insert(1, {'at'   : 1,
                                    'nhits': 0})
    if rawECDFs['data'][-1]['at'] < rawECDFs['nfevs']:
        rawECDFs['data'].append({'at'   : rawECDFs['nfevs'],
                                 'nhits': np.nan})
    
    # Compute the proportional data for the ECDF
    import copy
    percECDFs = copy.deepcopy(rawECDFs)
    induce_proportion_single_algo(percECDFs) # dict is mutable, thus passed by reference
    
    return rawECDFs, percECDFs


def mergeHittingTimes(listOfRawECDFs, itemsToMerge=None):
    """
    Aggregate rawECDFs array.
    If 'itemsToMerge' is not specified, simply merge all hitting times in the list
    """
    result = []
    if itemsToMerge is None:
        for rawECDFs in listOfRawECDFs:
            result = result + rawECDFs  # concatenation of 2 lists
    else:
        for i in itemsToMerge:
            result = result + listOfRawECDFs[i]  # concatenation of 2 lists
    
    return sorted(result, key = lambda point: point['at'])


def induce_proportion_multi_algos(percECDFs):
    """
    Dict is mutable, thus passed by reference
    This function will apply the 'induce_proportion_single_algo' to the single rawECDFs
    of each individual algorithm's data within the collective multi_algo rawECDFs.
    The output is a perc ECDF (with 'proportion' and 'time' entries):
    """
    for algo in percECDFs.keys():
        induce_proportion_single_algo(percECDFs[algo]) # dict is mutable, thus passed by reference


def induce_proportion_single_algo(singleECDF):
    """
    Dict is mutable, thus passed by reference
    Compute the y-axis data for the ECDF curve.
    The input is raw ECDFs, the output is perc ECDFs with additional keys:
      - 'proportion'
      - 'time'
    Called by:
        for ialgo, algo in enumerate([d['algo'] for d in data]):
            number_of_runs = data[ialgo]['nRun']
            loghit[algo] = {'info': data[ialgo]}
            loghit[algo].update(ecdf.compute_single_algo(run_profiles=hv_diff[:, :number_of_runs, ialgo], targets=target_set))
            
            ecdf_data.update(ecdf.induce_proportion_single_algo(ecdf_data, algo))
    or
        ecdf_data.update(ecdf.induce_proportion_single_algo(ecdf_data))
        # or: percECDFs = ecdf.induce_proportion_single_algo(singleECDF)
    """
    
    # For each algorithm's database, compute the increment (i.e. proportion for one step up)
    # and then induce the whole ECDF proportional profile according to the hitting times database
    increment = 1.0 / (singleECDF['ntargets'] * singleECDF['nprofiles'])
    
    # Pre-set for the 0th f-eval, which should always be always 0 (as no "performance"
    # is observed), before possibly hitting some target(s) at the 1st f-eval.
    singleECDF['proportion'] = [0]
    
    # The raw ECDFs is now perc ECDFs (i.e. with 2 new entries of keys 'proportion' and 'time')
    singleECDF['time'] = [item['at'] for item in singleECDF['data']]
    
    # Compute the cumulative proportion over time, starting from the
    # second element, i.e. k+1, as the first was already pre-set to 0
    for k in xrange(len(singleECDF['data']) - 1): # start from 2nd element, thus len-1
        singleECDF['proportion'].append( singleECDF['proportion'][k] + singleECDF['data'][k+1]['nhits'] * increment )


def extract_plot_data(percECDFs, startFevPlot, items=None):
    """
    Extract the plottable portion of ECDF data, x-limited by ${startFevPlot}.
    If ${items} is not given, all the ECDFs in ${percECDFs} are processed,
    otherwise, only the one(s) specified by ${items} will be processed.
    Called by:
        ecdf_coord = {}
        for ialgo in xrange(len(data)):
            ecdf_coord.update(ecdf.extract_plot_data(percECDFs, startFevPlot, ialgo))
    or simply:
        ecdf_coord = ecdf.extract_plot_data(percECDFs, startFevPlot)
    """
    # Induce 'keys' from the input 'items' for selecting the database of
    # the algorithm(s) that will be processed within the 'percECDFs'
    if items is None:
        keys = percECDFs.keys()
    elif type(items) == type(str()): # only 1 item is given as a string, not list of strings
        keys = [items]
    else:
        keys = items
    
    ecdf_coord = {}
    for key in keys:
        #ecdf_coord[key] = {'x': [item['at'] for item in percECDFs[key]['data']],
        #                   'y': percECDFs[key]['proportion']}
        ecdf_coord.update({key: {'x': percECDFs[key]['time'], # [item['at'] for item in percECDFs[key]['data']],
                                 'y': percECDFs[key]['proportion']}})
        y_pre = None
        # Trim the ecdf_coord array w.r.t. the x-limit by ${startFevPlot}
        while ecdf_coord[key]['x'][0] < startFevPlot:
            del ecdf_coord[key]['x'][0]
            y_pre = ecdf_coord[key]['y'][0] # store a copy of y before deleting for later use
            del ecdf_coord[key]['y'][0]
        
        # Add ${startFevPlot} as an element of 'x' and repeat a predecessor in 'y'
        if ecdf_coord[key]['x'][0] > startFevPlot:
            ecdf_coord[key]['x'].insert(0, startFevPlot)
            if y_pre == None:
                ecdf_coord[key]['y'].insert(0, y_pre)
            else:
                ecdf_coord[key]['y'].insert(0, 0)
    
    return ecdf_coord


def aggregate_ecdfs_2objs_bbob(ecdfs_db, obj1=None, obj2=None):
    """
    Aggregate ECDFs within ONE database indicated by the lists ${obj1}, ${obj2}.
    Return: a dictionary with similar structure to ecdfs_db[1][2] or so.
    """
    # According to the current bi-objective BBOB setup, f2 starts from f1 onwards,
    # we cannot do a complete cross-combination of two groups of functions. Thus
    # we generate a dictionary ${mop} to hold selected combinations as follows:
    #   e.g. obj1 = [1, 4, 2, 3, 6] # no need to be in order
    #        obj2 = [3, 2, 1, 7, 8, 6]
    #    --> mop = {1: {1: True, 2: True, 3: True, 6: True, 7: True, 8: True},
    #               2: {2: True, 3: True, 6: True, 7: True, 8: True},
    #               3: {3: True, 6: True, 7: True, 8: True},
    #               4: {6: True, 7: True, 8: True},
    #               6: {6: True, 7: True, 8: True}}
    if obj1 is None or obj2 is None: # use all possible combinations
        mop = {}
        for f1 in ecdfs_db.keys():
            if isinstance(f1, int): # since there is a key 'targets' out there!
                mop[f1] = {} # or: mop.update(f1: {})
                for f2 in ecdfs_db[f1].keys():
                    mop[f1][f2] = True # or: mop[f1].update(f2: True)
    else: # use those indicated
        mop = {}
        for f1 in obj1:
            mop[f1] = {}
            for f2 in obj2:
                if f2 < f1: continue # as values for f2 start from f1 onwards
                mop[f1][f2] = True
    
    # Initialize the output by copying an existing data structure
    f1 = mop.keys()[0]
    f2 = mop[f1].keys()[0]
    import copy
    result = copy.deepcopy(ecdfs_db[f1][f2]) # CAUTION: assignment does not make a copy!!!
    # here the 'targets' entry is dropped out.
    
    # Way to go: for each algorithm, do the aggregation over the specified groups of functions
    for algo in ecdfs_db[f1][f2].keys():
        number_of_profiles = 0
        all_hitting_times = []
        for f1 in mop.keys():
            for f2 in mop[f1].keys():
                #print(f1,f2)
                #if f1 == 9 and f2 == 9:
                #    print("Hi");
                # Check the existence of this combination:
                try: ecdfs_db[f1][f2]
                except:
                    print "Combination of %d and %d does not exist, I ignore it" % (f1,f2)
                    continue # simply ignore this combination
                # Deal with the entry 'ntargets'
                if result[algo]['ntargets'] != ecdfs_db[f1][f2][algo]['ntargets']:
                    print "You are aggregating ecdfs of different 'ntargets', which is not yet supported by this routine!"
                    quit()
                # Deal with the entry 'nfevs'
                if result[algo]['nfevs'] < ecdfs_db[f1][f2][algo]['nfevs']:
                    result[algo]['nfevs'] = ecdfs_db[f1][f2][algo]['nfevs']
                # Cumulate the entries 'nprofiles' and 'data'
                number_of_profiles += ecdfs_db[f1][f2][algo]['nprofiles']
                all_hitting_times = all_hitting_times + ecdfs_db[f1][f2][algo]['data']
        
        result[algo]['nprofiles'] = number_of_profiles
        result[algo]['data'] = sorted(all_hitting_times, key = lambda point: point['at'])
        
        # Remove the multiplicity of the added datum {'at'   : rawECDFs['nfevs'],
        #                                             'nhits': np.nan}
        # buggy condition does not remove all the nans
        # while np.isnan(result[algo]['data'][-1]['nhits']) and \
        #       np.isnan(result[algo]['data'][-2]['nhits']) and \
        #       result[algo]['data'][-2]['at'] == result[algo]['data'][-2]['at']:
        #     del result[algo]['data'][-1]
        # rectified code : filter the nans
        
        where_n_fevs = np.where(np.array([x['at'] for x in result[algo]['data']])==result[algo]['nfevs'])[0]
        nhits_at_fevs = [result[algo]['data'][id]['nhits'] for id in where_n_fevs]
        where_nan = np.where(np.isnan(nhits_at_fevs))[0]
        
        #print(where_nan)
        if len(where_nan) == 0:
            break
        elif len(where_nan) < len(nhits_at_fevs):# delete all nan
            
            for id in where_nan[::-1]:
                #print("Shape of data %d" % len(result[algo]['data']))
                #print("Idx removed %d" % where_n_fevs[id])
                del result[algo]['data'][where_n_fevs[id]]
        elif len(where_nan) == len(nhits_at_fevs): # keep one nan
            for id in where_nan[::-1]:
                del result[algo]['data'][where_n_fevs[id]]
        else:
            print("There are nans other than the end of the profile, this shouldnt happen !!")
            where_nan[0][0][0]


        
        # dict is mutable, thus passed by reference, no need of return
        # i.e. no more: result[algo].update(induce_proportion_single_algo(result[algo]))
        induce_proportion_single_algo(result[algo])
    
    return result


def aggregate_ecdfs_across_databases_2objs_bbob(list_of_ecdfs_dbs):
    """
    Aggregate the *common* ECDFs across multiple databases.
    Abort those ECDFs that do not exist over all input databases.
    Return: a dictionary with similar structure to ecdfs_db[1][2] or so.
    """
    # First, generate a dictionary ${mop} to hold *common* combinations as follows:
    #   e.g. obj1 = [1, 2, 3, 4, 6]
    #        obj2 = [1, 2, 3, 6, 7, 8]
    #    --> mop = {1: {1: True, 2: True, 3: True, 6: True, 7: True, 8: True},
    #               2: {2: True, 3: True, 6: True, 7: True, 8: True},
    #               3: {3: True, 6: True, 7: True, 8: True},
    #               4: {6: True, 7: True, 8: True},
    #               6: {6: True, 7: True, 8: True}}
    mop = {}
    for f1 in list_of_ecdfs_dbs[0].keys():
        # Skip this key if this is 'targets' key, not a funcID key
        if not isinstance(f1, int): continue
        # Ensure f1 exists in all databases
        abort = False
        for database in list_of_ecdfs_dbs[1:]:
            if not database.has_key(f1):  # if not (f1 in database.keys()):
                abort = True
                break
        if abort: continue
        # Now knowing f1 does exist over all databases
        mop[f1] = {}  # or: mop.update(f1: {})
        for f2 in list_of_ecdfs_dbs[0][f1].keys():
            # Ensure f2 exists in all databases
            abort = False
            for database in list_of_ecdfs_dbs[1:]:
                if not database[f1].has_key(f2):  # if not (f2 in database[f1].keys()):
                    abort = True
                    break
            if abort: continue
            # Now knowing f2 does exist over all databases
            mop[f1][f2] = True  # or: mop[f1].update(f2: True)
    
    # Initialize the output ecdfs database by copying an existing data structure
    import copy
    output = copy.deepcopy(list_of_ecdfs_dbs[0]) # CAUTION: assignment does not make a copy!!!
    # Delete funcID keys that are not common across databases, i.e. not in the map 'mop'
    # here the 'targets' entry gets deleted.
    for f1 in output.keys():
        if not mop.has_key(f1):
            del output[f1]
        else:
            for f2 in output[f1].keys():
                if not mop[f1].has_key(f2):
                    del output[f1][f2]
    
    for f1 in mop.keys():
        for f2 in mop[f1].keys():
            for algo in list_of_ecdfs_dbs[0][f1][f2].keys():
                number_of_profiles = 0
                all_hitting_times = []
                # Way to go: under each combination and for each algorithm,
                # do the aggregation over the given databases
                for database in list_of_ecdfs_dbs:
                    # Deal with the entry 'ntargets'
                    if output[f1][f2][algo]['ntargets'] != database[f1][f2][algo]['ntargets']:
                        print "You are aggregating ecdfs of different 'ntargets', which is not yet supported by this routine!"
                        quit()
                    # Deal with the entry 'nfevs'
                    if output[f1][f2][algo]['nfevs'] < database[f1][f2][algo]['nfevs']:
                        output[f1][f2][algo]['nfevs'] = database[f1][f2][algo]['nfevs']
                    # Cumulate the entries 'nprofiles' and 'data'
                    number_of_profiles += database[f1][f2][algo]['nprofiles']
                    all_hitting_times = all_hitting_times + database[f1][f2][algo]['data']
                
                output[f1][f2][algo]['nprofiles'] = number_of_profiles
                output[f1][f2][algo]['data'] = sorted(all_hitting_times, key = lambda point: point['at'])
                
                # Remove the multiplicity of the added {'at': rawECDFs['nfevs'], 'nhits': np.nan}
                while np.isnan(output[f1][f2][algo]['data'][-1]['nhits']) and \
                      np.isnan(output[f1][f2][algo]['data'][-2]['nhits']) and \
                      output[f1][f2][algo]['data'][-2]['at'] == output[f1][f2][algo]['data'][-2]['at']:
                    del output[f1][f2][algo]['data'][-1]
                
                # dict is mutable, thus passed by reference, no need of return
                # i.e. no more: output[f1][f2][algo].update(induce_proportion_single_algo(output[f1][f2][algo]))
                induce_proportion_single_algo(output[f1][f2][algo])
    
    return output



def aggregate_ecdfs_problems(ecdfs_db, problems=None):
    """
    Aggregate ECDFs within ONE database indicated by the lists problems.
    Return: a dictionary with similar structure to ecdfs_db or so.
    """
    # According to the current bi-objective BBOB setup, f2 starts from f1 onwards,
    # we cannot do a complete cross-combination of two groups of functions. Thus
    # we generate a dictionary ${mop} to hold selected combinations as follows:
    #   e.g. obj1 = [1, 4, 2, 3, 6] # no need to be in order
    #        obj2 = [3, 2, 1, 7, 8, 6]
    #    --> mop = {1: {1: True, 2: True, 3: True, 6: True, 7: True, 8: True},
    #               2: {2: True, 3: True, 6: True, 7: True, 8: True},
    #               3: {3: True, 6: True, 7: True, 8: True},
    #               4: {6: True, 7: True, 8: True},
    #               6: {6: True, 7: True, 8: True}}
    if problems is None: # use all possible combinations
        mop = {}
        for pidx in ecdfs_db.keys():
            if isinstance(pidx, int): # since there is a key 'targets' out there!
                mop[pidx] = True # or: mop.update(f1: {})
                #for f2 in ecdfs_db[f1].keys():
                #    mop[f1][f2] = True # or: mop[f1].update(f2: True)
    else: # use those indicated
        mop = {}
        for pidx in problems:
            mop[pidx] = True
    
    # Initialize the output by copying an existing data structure
    pid0 = mop.keys()[0]
    import copy
    result = copy.deepcopy(ecdfs_db[pid0]) # CAUTION: assignment does not make a copy!!!
    # here the 'targets' entry is dropped out.
    
    # Way to go: for each algorithm, do the aggregation over the specified groups of functions
    for algo in ecdfs_db[pid0].keys():
        number_of_profiles = 0
        all_hitting_times = []
        for pidx in mop.keys():
            #print(f1,f2)
            #if f1 == 9 and f2 == 9:
            #    print("Hi");
            # Check the existence of this combination:
            try: ecdfs_db[pidx]
            except:
                print "ECDF of %d does not exist, I ignore it" % (pidx)
                continue # simply ignore this combination
            # Deal with the entry 'ntargets'
            if result[algo]['ntargets'] != ecdfs_db[pidx][algo]['ntargets']:
                print "You are aggregating ecdfs of different 'ntargets', which is not yet supported by this routine!"
                quit()
            # Deal with the entry 'nfevs'
            if result[algo]['nfevs'] < ecdfs_db[pidx][algo]['nfevs']:
                result[algo]['nfevs'] = ecdfs_db[pidx][algo]['nfevs']
            # Cumulate the entries 'nprofiles' and 'data'
            number_of_profiles += ecdfs_db[pidx][algo]['nprofiles']
            all_hitting_times = all_hitting_times + ecdfs_db[pidx][algo]['data']
        
        result[algo]['nprofiles'] = number_of_profiles
        result[algo]['data'] = sorted(all_hitting_times, key = lambda point: point['at'])
        
        # Remove the multiplicity of the added datum {'at'   : rawECDFs['nfevs'],
        #                                             'nhits': np.nan}
        # buggy condition does not remove all the nans
        # while np.isnan(result[algo]['data'][-1]['nhits']) and \
        #       np.isnan(result[algo]['data'][-2]['nhits']) and \
        #       result[algo]['data'][-2]['at'] == result[algo]['data'][-2]['at']:
        #     del result[algo]['data'][-1]
        # rectified code : filter the nans
        
        where_n_fevs = np.where(np.array([x['at'] for x in result[algo]['data']])==result[algo]['nfevs'])[0]
        nhits_at_fevs = [result[algo]['data'][id]['nhits'] for id in where_n_fevs]
        where_nan = np.where(np.isnan(nhits_at_fevs))[0]
        
        #print(where_nan)
        if len(where_nan) == 0:
            pass
        elif len(where_nan) < len(nhits_at_fevs):# delete all nan
            
            for id in where_nan[::-1]:
                #print("Shape of data %d" % len(result[algo]['data']))
                #print("Idx removed %d" % where_n_fevs[id])
                del result[algo]['data'][where_n_fevs[id]]
        elif len(where_nan) == len(nhits_at_fevs): # keep one nan
            for id in where_nan[::-1]:
                del result[algo]['data'][where_n_fevs[id]]
        else:
            print("There are nans other than the end of the profile, this shouldnt happen !!")
            where_nan[0][0][0]
            
        
        where_n_fevs = np.where(np.array([x['at'] for x in result[algo]['data']])!=result[algo]['nfevs'])[0]
        nhits_at_fevs = [result[algo]['data'][id]['nhits'] for id in where_n_fevs]
        where_nan = np.where(np.isnan(nhits_at_fevs))[0]        
        
        
        if len(where_nan) == 0:
            pass
        else:
            for id in where_nan[::-1]:
                #print("Shape of data %d" % len(result[algo]['data']))
                #print("Idx removed %d" % where_n_fevs[id])
                del result[algo]['data'][where_n_fevs[id]]


        
        # dict is mutable, thus passed by reference, no need of return
        # i.e. no more: result[algo].update(induce_proportion_single_algo(result[algo]))
        induce_proportion_single_algo(result[algo])
    
    return result














def aggregate_ecdfs_across_dimensions_2objs_bbob(list_of_ecdfs_dbs, dims=[]):
    """
    Inputs are a list of ECDFs databases of different dimensions
           and a list of corresponding dimensions
    """
    # Take the least common multiple of the dimensions
    dim_lcm = least_common_multiple(*dims)
    
    # Make copies of input databases
    import copy
    databases = []
    for db in list_of_ecdfs_dbs:
        databases.append(copy.deepcopy(db))
    
    # Scale up the 'at' values and 'nfevs' value by a factor (dim_lcm / dim)
    for idb, database in enumerate(databases):
        for f1 in database.keys():
            # Skip this key if this is 'targets' key, not a funcID key
            if not isinstance(f1, int): continue
            # Now knowing f1 does exist over all databases
            for f2 in database[f1].keys():
                for algo in database[f1][f2]:
                    multiple = int(dim_lcm/dims[idb])
                    databases[idb][f1][f2][algo]['nfevs'] *= multiple
                    for i in xrange(len(database[f1][f2][algo]['data'])):
                        databases[idb][f1][f2][algo]['data'][i]['at'] *= multiple
    
    # Now ready for the aggregation across the databases
    agg_database = aggregate_ecdfs_across_databases_2objs_bbob(databases)
    agg_dim = dim_lcm
    
    return agg_database, agg_dim


#================= COMPUTE THE LEAST COMMON MULTIPLE ===========================
def gcd(a, b):
    """Return greatest common divisor using Euclid's Algorithm."""
    while b:
        a, b = b, a % b
    return a

def lcm(a, b):
    """Return lowest common multiple."""
    return a * b // gcd(a, b)

def least_common_multiple(*args):
    """Return lcm of args."""
    return reduce(lcm, args)

# Test cases:
#print least_common_multiple(2,3,4,5,6)     # this is correct
#print least_common_multiple([2,3,4,5,6])   # NOT a correct way to pass arguments
#print least_common_multiple(*[2,3,4,5,6])  # this the correct
#print least_common_multiple(*range(1, 20)) # this is correct
#===============================================================================



if __name__ == "__main__":
    print "I am a module and should only be called from another program!" 
    # raw_input()
