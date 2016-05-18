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
import texify # compile the main latex file to pdf

import ecdf
from ecdf import least_common_multiple
#import fancyplot
import fancyplot_bigFont as fancyplot

# MPI stuff
#from mpi4py import MPI

def composeLatex(cmd=None, file=None, stringlist=None):
    """Compose the content of a LaTeX file"""
    if cmd == "init":
        pfile = open(file, "w")
        # Write header code to LaTeX files
        pfile.write(r"""
\documentclass[a4paper]{article}
\usepackage[utf8x]{inputenc}
\usepackage{fullpage}
\usepackage[left=2mm,right=2mm,top=3mm,bottom=4mm,footskip=1mm,includefoot,heightrounded]{geometry}
\usepackage[usenames,svgnames,dvipsnames,table]{xcolor}
\usepackage{graphicx}
\usepackage{pdflscape}

\begin{document}
\begin{landscape}
\begin{center}

Intentionally left blank by a $\backslash${\tt clearpage}, otherwise the first table does not appear as expected.
\clearpage

""")
        return pfile
        
    elif cmd == "fini":
        # Write footer code to LaTeX files
        file.write("""
\end{center}
\end{landscape}
\end{document}
""")
        file.close()
        
    else:
        for string in stringlist:
            file.write("\includegraphics[width=7cm, trim = 0mm 0mm 0mm 0mm, clip]{%s}\n" % (string))
        file.write("\n" + r"\vspace{1mm}" + "\n\n")
        


def ecdf_processing(ind = 'eps'):
    #====================================#
    # Read info from BBOB
    #====================================#
    PROBLEM_LIST = ['BK1','ex005','Deb41','Deb512a','Deb512b','Deb512c','Deb513','Deb521a','Deb521b','Deb53','ZDT1','ZDT2','ZDT3','ZDT4','ZDT6','DTLZ1','DTLZ2','DTLZ3','DTLZ4','DTLZ5','DTLZ6','DTLZ1n2','DTLZ2n2','DTLZ3n2','DTLZ4n2','DTLZ5n2','DTLZ6n2','Kursawe','Fonseca','L1ZDT4','L2ZDT1','L2ZDT2','L2ZDT3','L2ZDT4','L2ZDT6','L3ZDT1','L3ZDT2','L3ZDT3','L3ZDT4','L3ZDT6','WFG1','WFG2','WFG3','WFG4','WFG5','WFG6','WFG7','WFG8','WFG9','I1','I2','I3','I4','I5','MOP1','MOP2','MOP3','MOP4','MOP5','MOP6','MOP7','DPAM1','DG01','Far1','FES1','FES2','FES3','IKK1','Jin1','Jin2','Jin3','Jin4','OKA1','OKA2','LRS1','IM1','LE1','MHHM1','MHHM2','MLF1','MLF2','QV1','Sch1','SP1','SSFYY1','SSFYY2','SK1','SK2','TKLY1','VU1','VU2','VFM1','ZLT1','CL1','lovison1','lovison2','lovison3','lovison4','lovison5','lovison6']
    DIM_LIST = [2,2,2,2,2,2,2,2,2,2,30,30,30,10,10,7,12,12,12,12,22,2,2,2,2,2,2,3,2,10,30,30,30,30,10,30,30,30,30,10,8,8,8,8,8,8,8,8,8,8,8,8,8,8,1,4,2,3,2,2,2,10,1,2,10,10,10,2,2,2,2,2,2,3,2,2,2,1,2,1,2,10,1,2,2,1,1,4,4,2,2,2,10,4,2,2,2,2,3,3]

    #xlim = [-5, 5]
    #INST_LIST = ((2,4), (3,5), (7,8), (9,10), (11,12))
    NUM_RUNS = 10
    NUM_TARGETS = 70
    
    #======= TO BE EDITED =====================================================
    FILENAME_ECDF_DATABASE = 'db_ecdfs_%s_maxfevalDiff_%druns_%dtargets.pkl'
    suffix = "_%druns_aggregate" % (NUM_RUNS)
    #==========================================================================

    
    fullRangeY    = True
    compileLatex  = True
    compileLogx   = True
    
    MAIN_TEX_FIGS = "plot_ecdf_%s%s.tex" % (ind,suffix) # used to control the compilation of LaTeX
    DIR_FIGS = "../" + "postproc" + "/" + "figs%s" % suffix
    DIR_TEX_FIGS = os.path.split(DIR_FIGS)[0] + "/" + "latex4figs%s" % suffix
    
    # Check the existence of output folders
    if not os.path.isdir(DIR_FIGS): os.makedirs(DIR_FIGS)
    if not os.path.isdir(DIR_TEX_FIGS): os.makedirs(DIR_TEX_FIGS)
    
    ftex_main_figs_pdf = composeLatex(cmd='init', file=DIR_TEX_FIGS + "/" + MAIN_TEX_FIGS)
    
    #-------------------------------------------------------------------
    # Read the database dictionary back from the pickle
    try: all_ecdfs
    except NameError:
        import pickle
        f = open(FILENAME_ECDF_DATABASE % (ind,NUM_RUNS, NUM_TARGETS), 'rb')
        all_ecdfs = pickle.load(f)
        f.close()
    #-------------------------------------------------------------------
    
    #g1func = g2func = range(1, 25)
    #g1func = range(1, 6)
    #g2func = range(6, 10)
    
    groups = [{'all functions'  : range(1, 25)},
              {'separable'      : range(1, 6)},
              {'moderate'       : range(6, 10)},
              {'ill-conditioned': range(10, 15)},
              {'multi-modal'    : range(15, 20)},
              {'w/s multi-modal': range(20, 25)}]
              
    #-- Problem category:

    L=set(['ex005','Deb41','Deb512a','Deb512b','Deb512c','Deb513','Deb521a','Deb521b','Deb53','DTLZ1n2','DTLZ2n2','DTLZ3n2','DTLZ4n2','DTLZ5n2','DTLZ6n2','Kursawe','Fonseca','MOP1','MOP2','MOP3','MOP4','MOP5','MOP6','MOP7','BK1','DG01','Far1','IKK1','Jin1','Jin2','Jin3','Jin4','OKA1','OKA2','LRS1','IM1','LE1','MHHM1','MHHM2','MLF1','MLF2','Sch1','SP1','SSFYY1','SSFYY2','SK1','SK2','TKLY1','VU1','VU2','VFM1','CL1','lovison1','lovison2','lovison3','lovison4','lovison5','lovison6'])
    H=set(['DPAM1','DTLZ1','DTLZ2','DTLZ3','DTLZ4','DTLZ5','DTLZ6','FES1','FES2','FES3','I1','I2','I3','I4','I5','L1ZDT4','L2ZDT1','L2ZDT2','L2ZDT3','L2ZDT4','L2ZDT6','L3ZDT1','L3ZDT2','L3ZDT3','L3ZDT4','L3ZDT6','QV1','WFG1','WFG2','WFG3','WFG4','WFG5','WFG6','WFG7','WFG8','WFG9','ZDT1','ZDT2','ZDT3','ZDT4','ZDT6','ZLT1'])

    S = set(['BK1','FES1','FES2','FES3','LE1','LRS1','MHHM2','MOP1','MOP2','MOP4','MOP6','QV1','SK1','SSFYY1','VFM1','VU1','VU2','WFG1','WFG4','WFG5','WFG7','ZDT1','ZDT2','ZDT3','ZDT4','ZDT6','ZLT1','I1','Fonseca'])
    NS = set(['DPAM1','Far1','MLF2','MOP5','SP1','WFG2','WFG3','WFG6','WFG8','WFG9','I2','I3','I4','I5'])

    U=set(['ex005','BK1','DTLZ2','DTLZ2n2','DTLZ4n2','DTLZ5n2','DTLZ6n2','DTLZ4','DTLZ5','DTLZ6','FES1','FES2','FES3','IKK1','IM1','LE1','MHHM1','MHHM2','LRS1','MOP1','MOP2','MOP7','SP1','SSFYY1','VFM1','VU1','VU2','WFG1','WFG3','WFG6','WFG7','WFG8','ZDT1','ZDT2','ZLT1','I1','I2','I3','I4','I5','Fonseca','Jin1','Jin2','Jin3','Jin4'])
    M = set(['DG01','DTLZ1','DTLZ3''DTLZ1n2','DTLZ3n2','Far1','MLF1','MLF2','QV1','SK1','WFG4','ZDT6'])
    All = set(PROBLEM_LIST)
    Mixed = All.difference(S.union(NS).intersection(U.union(M)))
    
    groups={'low-dimensioanl':[pid for pid, problem in enumerate(PROBLEM_LIST) if any(p in problem for p in list(L))],
            'high-dimensional':[pid for pid, problem in enumerate(PROBLEM_LIST) if any(p in problem for p in list(H))],
            'separable' :[pid for pid, problem in enumerate(PROBLEM_LIST) if any(p in problem for p in list(S))],
            'non-separable':[pid for pid, problem in enumerate(PROBLEM_LIST) if any(p in problem for p in list(NS))],
            'uni-modal':[pid for pid, problem in enumerate(PROBLEM_LIST) if any(p in problem for p in list(U))],
            'multi-modal':[pid for pid, problem in enumerate(PROBLEM_LIST) if any(p in problem for p in list(M))],
            'mixed':[pid for pid, problem in enumerate(PROBLEM_LIST) if any(p in problem for p in list(Mixed))],
            'all': np.arange(0,len(PROBLEM_LIST))}


              
    nGrps = len(groups)
    
    for gidx, gid in enumerate(groups):
        print gid
        dims = [DIM_LIST[pid] for pid in groups[gid]]
        DIM = 1#least_common_multiple(*dims)
        #agg_ecdfs = ecdf.aggregate_ecdfs_2objs_bbob(all_ecdfs, g1func, g2func)
        agg_ecdfs = ecdf.aggregate_ecdfs_problems(all_ecdfs, problems=groups[gid],dims= [DIM_LIST[id] for id in groups[gid]])
        #agg_ecdfs = ecdf.aggregate_ecdfs_2objs_bbob(all_ecdfs) # simply aggregates all ecdfs in the database
        
        print 'nprofiles:', agg_ecdfs[agg_ecdfs.keys()[0]]['nprofiles']
        print 'ndata:', len(agg_ecdfs[agg_ecdfs.keys()[0]]['data'])
        print '================================================'
        
        # Adjust the ecdf data such that we can plot from ${startFevPlot}
        startFevPlot = -1  # DIM - 1
        startFevLim  = 1
        ecdf_coord = ecdf.extract_plot_data(agg_ecdfs, startFevPlot)
        
        # Modify x's such that we can display the x-axis as "# f-evals / DIM"
        for algo in ecdf_coord.keys():
            ecdf_coord[algo]['x'] = [item for item in ecdf_coord[algo]['x']] #change for this
            #ecdf_coord[algo]['x'] = [item for item in ecdf_coord[algo]['x']]
        #-------------------------------------------------------------------
        
        # Plot commands start here
        pltObject = plottingConfig() # configure the plot
        #pltObject.legends = [key[3:] for key in ecdf_coord.keys()]
        pltObject.legends = [key[3:].replace('MO-CMA-ES','MO-CMA') for key in ecdf_coord.keys()]
        #pltObject.xlabel = "# f-evals / DIM"
        pltObject.xlabel = "# $\mathbf{f}$-evaluations / dimension"
        pltObject.ylabel = "Proportion of (function+target)"
        #pltObject.title = u"%s (F1) & %s (F2), %dD" % (g1name, g2name, DIM)
        pltObject.title = ""
        pltObject.note = ["%s" % (gid)]
        
        # Construct output filenames
        fname = "%s_%druns" % (gid, NUM_RUNS)
        fnameEcdf     = fname + "_ecdf_%s.pdf" % ind
        fnameEcdfLogx = fname + "_ecdf_%s_logx.pdf" % ind
        
        # Plot to the PDF files
        #pltObject.plotting(DIR_FIGS + "/" + fnameEcdf,     ecdf_coord, xstart=startFevLim/DIM, logx=False, yfull=fullRangeY)
        #if ig1 == 0:
        #  print("Plotting all")
        #pltObject.plotting(DIR_FIGS + "/" + fnameEcdfLogx, ecdf_coord, xstart=startFevLim/DIM, logx=True,  yfull=fullRangeY)
        pltObject.plotting(DIR_FIGS + "/" + fnameEcdfLogx, ecdf_coord, xstart=startFevLim, logx=True,  yfull=fullRangeY)
        #-------------------------------------------------------------------
        

        # MPI stuff
        #MPI.COMM_WORLD.Barrier()
        
        # Define a layout of PDFs and include them in TeX file
        try: name_buf
        except NameError:
            name_buf = []
        name_buf = name_buf + ([fnameEcdfLogx] if compileLogx else [fnameEcdf])
        
        if len(name_buf) == 4:
                composeLatex(cmd='add', file=ftex_main_figs_pdf,
                            stringlist=["../" + os.path.split(DIR_FIGS)[1] + "/" + name_buf[0],
                                        "../" + os.path.split(DIR_FIGS)[1] + "/" + name_buf[1],
                                        "../" + os.path.split(DIR_FIGS)[1] + "/" + name_buf[2],
                                        "../" + os.path.split(DIR_FIGS)[1] + "/" + name_buf[3]])
                name_buf = []
        #-------------------------------------------------------------------
    
    composeLatex(cmd='fini', file=ftex_main_figs_pdf)
    
    # Compile the main LaTeX file to PDF
    if compileLatex:
        texify.pdfLatex(directory=DIR_TEX_FIGS, filename=MAIN_TEX_FIGS, cleanup=[".log", ".aux"])


def ecdf_processing_allindicators(indicators = ['eps', 'gd','igd']):
    '''
        aggregate and plots ecdfs over all the problems and all the indicators into a single plot
    '''
    #====================================#
    # Read info from BBOB
    #====================================#
    PROBLEM_LIST = ['BK1','ex005','Deb41','Deb512a','Deb512b','Deb512c','Deb513','Deb521a','Deb521b','Deb53','ZDT1','ZDT2','ZDT3','ZDT4','ZDT6','DTLZ1','DTLZ2','DTLZ3','DTLZ4','DTLZ5','DTLZ6','DTLZ1n2','DTLZ2n2','DTLZ3n2','DTLZ4n2','DTLZ5n2','DTLZ6n2','Kursawe','Fonseca','L1ZDT4','L2ZDT1','L2ZDT2','L2ZDT3','L2ZDT4','L2ZDT6','L3ZDT1','L3ZDT2','L3ZDT3','L3ZDT4','L3ZDT6','WFG1','WFG2','WFG3','WFG4','WFG5','WFG6','WFG7','WFG8','WFG9','I1','I2','I3','I4','I5','MOP1','MOP2','MOP3','MOP4','MOP5','MOP6','MOP7','DPAM1','DG01','Far1','FES1','FES2','FES3','IKK1','Jin1','Jin2','Jin3','Jin4','OKA1','OKA2','LRS1','IM1','LE1','MHHM1','MHHM2','MLF1','MLF2','QV1','Sch1','SP1','SSFYY1','SSFYY2','SK1','SK2','TKLY1','VU1','VU2','VFM1','ZLT1','CL1','lovison1','lovison2','lovison3','lovison4','lovison5','lovison6']
    DIM_LIST = [2,2,2,2,2,2,2,2,2,2,30,30,30,10,10,7,12,12,12,12,22,2,2,2,2,2,2,3,2,10,30,30,30,30,10,30,30,30,30,10,8,8,8,8,8,8,8,8,8,8,8,8,8,8,1,4,2,3,2,2,2,10,1,2,10,10,10,2,2,2,2,2,2,3,2,2,2,1,2,1,2,10,1,2,2,1,1,4,4,2,2,2,10,4,2,2,2,2,3,3]

    #xlim = [-5, 5]
    #INST_LIST = ((2,4), (3,5), (7,8), (9,10), (11,12))
    NUM_RUNS = 10
    NUM_TARGETS = 70
    
    #======= TO BE EDITED =====================================================
    FILENAME_ECDF_DATABASE = 'db_ecdfs_%s_maxfevalDiff_%druns_%dtargets.pkl'
    suffix = "_%druns_aggregate" % (NUM_RUNS)
    #==========================================================================

    
    fullRangeY    = True
    compileLatex  = True
    compileLogx   = True
    
    MAIN_TEX_FIGS = "plot_ecdf_allind%s.tex" % (suffix) # used to control the compilation of LaTeX
    DIR_FIGS = "../" + "postproc" + "/" + "figs%s" % suffix
    DIR_TEX_FIGS = os.path.split(DIR_FIGS)[0] + "/" + "latex4figs%s" % suffix
    
    # Check the existence of output folders
    if not os.path.isdir(DIR_FIGS): os.makedirs(DIR_FIGS)
    if not os.path.isdir(DIR_TEX_FIGS): os.makedirs(DIR_TEX_FIGS)
    
    ftex_main_figs_pdf = composeLatex(cmd='init', file=DIR_TEX_FIGS + "/" + MAIN_TEX_FIGS)
    
    #-------------------------------------------------------------------
    # Read the database dictionary back from the pickle
    try: all_ecdfs
    except NameError:
        import pickle
        list_of_ecdfs_dbs = []
        for ind in indicators:
            f = open(FILENAME_ECDF_DATABASE % (ind,NUM_RUNS, NUM_TARGETS), 'rb')
            #temp = ecdf.aggregate_ecdfs_problems(pickle.load(f))
            list_of_ecdfs_dbs.append(ecdf.aggregate_ecdfs_problems(pickle.load(f), problems=np.arange(0,len(PROBLEM_LIST)), dims= DIM_LIST))
            f.close()
    #-------------------------------------------------------------------
    # TO DO: this is a work around by double processing the databased first across all the problems as above then across all indicators as shown below (processing one profile set per indicator hence the [0])
    # aggregate the data profiles:
    agg_ecdfs = ecdf.aggregate_ecdfs_across_indicators_problems(list_of_ecdfs_dbs, problems = [0])      

    # plot the single profile 
    DIM = 1#least_common_multiple(*dims)
    #agg_ecdfs = ecdf.aggregate_ecdfs_problems(all_ecdfs)
    
        
    print 'nprofiles:', agg_ecdfs[agg_ecdfs.keys()[0]]['nprofiles']
    print 'ndata:', len(agg_ecdfs[agg_ecdfs.keys()[0]]['data'])
    print '================================================'
    
    # Adjust the ecdf data such that we can plot from ${startFevPlot}
    startFevPlot = -1  # DIM - 1
    startFevLim  = 1
    ecdf_coord = ecdf.extract_plot_data(agg_ecdfs, startFevPlot)
    
    # Modify x's such that we can display the x-axis as "# f-evals / DIM"
    for algo in ecdf_coord.keys():
        ecdf_coord[algo]['x'] = [item for item in ecdf_coord[algo]['x']] #change for this
        #ecdf_coord[algo]['x'] = [item for item in ecdf_coord[algo]['x']]
    #-------------------------------------------------------------------
    
    # Plot commands start here
    pltObject = plottingConfig() # configure the plot
    #pltObject.legends = [key[3:] for key in ecdf_coord.keys()]
    pltObject.legends = [key[3:].replace('MO-CMA-ES','MO-CMA') for key in ecdf_coord.keys()]
    #pltObject.xlabel = "# f-evals / DIM"
    pltObject.xlabel = "# $\mathbf{f}$-evaluations"
    pltObject.ylabel = "Proportion of (function+target)"
    #pltObject.title = u"%s (F1) & %s (F2), %dD" % (g1name, g2name, DIM)
    pltObject.title = "All Indicators"
    pltObject.note = ["all"]
    
    # Construct output filenames
    fname = "all_%druns" % (NUM_RUNS)
    fnameEcdf     = fname + "_ecdf_allind.pdf"
    fnameEcdfLogx = fname + "_ecdf_allind_logx.pdf"
    
    # Plot to the PDF files
    #pltObject.plotting(DIR_FIGS + "/" + fnameEcdf,     ecdf_coord, xstart=startFevLim/DIM, logx=False, yfull=fullRangeY)
    #if ig1 == 0:
    #  print("Plotting all")
    #pltObject.plotting(DIR_FIGS + "/" + fnameEcdfLogx, ecdf_coord, xstart=startFevLim/DIM, logx=True,  yfull=fullRangeY)
    pltObject.plotting(DIR_FIGS + "/" + fnameEcdfLogx, ecdf_coord, xstart=startFevLim, logx=True,  yfull=fullRangeY)
    #-------------------------------------------------------------------
    

    # MPI stuff
    #MPI.COMM_WORLD.Barrier()
    
    # Define a layout of PDFs and include them in TeX file
    try: name_buf
    except NameError:
        name_buf = []
    name_buf = name_buf + ([fnameEcdfLogx] if compileLogx else [fnameEcdf])
    
    if len(name_buf) == 4:
            composeLatex(cmd='add', file=ftex_main_figs_pdf,
                        stringlist=["../" + os.path.split(DIR_FIGS)[1] + "/" + name_buf[0],
                                    "../" + os.path.split(DIR_FIGS)[1] + "/" + name_buf[1],
                                    "../" + os.path.split(DIR_FIGS)[1] + "/" + name_buf[2],
                                    "../" + os.path.split(DIR_FIGS)[1] + "/" + name_buf[3]])
            name_buf = []
    #-------------------------------------------------------------------
    
    composeLatex(cmd='fini', file=ftex_main_figs_pdf)
    
    # Compile the main LaTeX file to PDF
    if compileLatex:
        texify.pdfLatex(directory=DIR_TEX_FIGS, filename=MAIN_TEX_FIGS, cleanup=[".log", ".aux"])

def ecdf_processing_alldims(ind='eps'):
    #====================================#
    # Read info from BBOB
    #====================================#
    FUNC_LIST = range(1,25) # from f1 to f24
    dims = [2,3,5,10,20]
    DIM = least_common_multiple(*dims)
    xlim = [-5, 5]
    INST_LIST = ((2,4), (3,5), (7,8), (9,10), (11,12))
    NUM_RUNS = 10
    
    #======= TO BE EDITED =====================================================
    FILENAME_ECDF_DATABASE = 'db_ecdfs_%s_maxfevalDiff_%druns_%dtargets_%dD.pkl'
    suffix = "_%druns_alldims_aggregate" % (NUM_RUNS)
    #==========================================================================
    
    fullRangeY    = True
    compileLatex  = True
    compileLogx   = True
    
    MAIN_TEX_FIGS = "plot_ecdf_%s%s.tex" % (ind,suffix) # used to control the compilation of LaTeX
    DIR_FIGS = "../" + "postproc" + "/" + "figs%s" % suffix
    DIR_TEX_FIGS = os.path.split(DIR_FIGS)[0] + "/" + "latex4figs%s" % suffix
    
    # Check the existence of output folders
    if not os.path.isdir(DIR_FIGS): os.makedirs(DIR_FIGS)
    if not os.path.isdir(DIR_TEX_FIGS): os.makedirs(DIR_TEX_FIGS)
    
    ftex_main_figs_pdf = composeLatex(cmd='init', file=DIR_TEX_FIGS + "/" + MAIN_TEX_FIGS)
    
    #-------------------------------------------------------------------
    # Read the databases dictionary back from the pickle one by one
    try: list_of_ecdfs_dbs
    except NameError:
        list_of_ecdfs_dbs = []
        import pickle
        for d in dims:
          f = open(FILENAME_ECDF_DATABASE % (ind,NUM_RUNS, 70, d), 'rb')
          list_of_ecdfs_dbs.append(pickle.load(f))
          f.close()
    all_ecdfs, _ = ecdf.aggregate_ecdfs_across_dimensions_2objs_bbob(list_of_ecdfs_dbs, dims=[2,3,5,10,20])
    #-------------------------------------------------------------------
    
    #g1func = g2func = range(1, 25)
    #g1func = range(1, 6)
    #g2func = range(6, 10)
    
    groups = [{'all functions'  : range(1, 25)},
              {'separable'      : range(1, 6)},
              {'moderate'       : range(6, 10)},
              {'ill-conditioned': range(10, 15)},
              {'multi-modal'    : range(15, 20)},
              {'w/s multi-modal': range(20, 25)}]

    
    for ig1, dg1 in enumerate(groups):
        for ig2, dg2 in enumerate(groups):
            # The 'all functions' combined only with 'all functions'
            # to make an aggreation over all 300 combinations
            if ig1 == 0 and ig2 != 0: continue
            
            # Select a subset to work on
            #if ig1 != 2: continue
            #if ig2 != 4: continue
            
            # For group-group combinations
            # we have 300 'onwards' combinations, not the whole 24*24
            if ig2 < ig1: continue
            
            # Extract the function list of group and the group name
            # Note, [0] to extract the list inside the list of dict's values
            g1func = dg1.values()[0]
            g2func = dg2.values()[0]
            g1name = dg1.keys()[0]
            g2name = dg2.keys()[0]
            print 'g1func:', g1func
            print 'g2func:', g2func
            print '-----------------------'
            
            # STEPS TO PRODUCE PLOTTABLE DATA FOR AN *AGGREGATED* ECDF
            #-------------------------------------------------------------------
            # *Aggregate* the selected ECDFs (in fact the hitting times)
            # and compute proportional ECDF (i.e. y-axis) as well
            agg_ecdfs = ecdf.aggregate_ecdfs_2objs_bbob(all_ecdfs, g1func, g2func)
            #agg_ecdfs = ecdf.aggregate_ecdfs_2objs_bbob(all_ecdfs) # simply aggregates all ecdfs in the database
            
            print 'nprofiles:', agg_ecdfs[agg_ecdfs.keys()[0]]['nprofiles']
            print 'ndata:', len(agg_ecdfs[agg_ecdfs.keys()[0]]['data'])
            print '================================================'
            
            # Adjust the ecdf data such that we can plot from ${startFevPlot}
            startFevPlot = 10*DIM - 1  # DIM - 1
            startFevLim  = 10*DIM
            ecdf_coord = ecdf.extract_plot_data(agg_ecdfs, startFevPlot)
            
            # Modify x's such that we can display the x-axis as "# f-evals / DIM"
            for algo in ecdf_coord.keys():
                ecdf_coord[algo]['x'] = [item/DIM for item in ecdf_coord[algo]['x']]
            #-------------------------------------------------------------------
            
            # Plot commands start here
            pltObject = plottingConfig() # configure the plot
            #pltObject.legends = [key[3:] for key in ecdf_coord.keys()]
            pltObject.legends = [key[3:].replace('MO-CMA-ES','MO-CMA') for key in ecdf_coord.keys()]

            #pltObject.xlabel = "# f-evals / DIM"
            pltObject.xlabel = "# $\mathbf{f}$-evaluations / lcm(dimensions)"
            pltObject.ylabel = "Proportion of (function+target)"
            #pltObject.title = u"%s (F1) & %s (F2), %dD" % (g1name, g2name, DIM)
            pltObject.title = "dimension = $\{2,3,5,10,20\}$"
            pltObject.note = ["F1: $f_{%d} \\rightarrow f_{%d}$ (%s)" % (g1func[0], g1func[-1], g1name),
                              "F2: $f_{%d} \\rightarrow f_{%d}$ (%s)" % (g2func[0], g2func[-1], g2name)]
            
            # Construct output filenames
            fname = "f%d--f%d_vs_f%d--f%d_%dD_%dinsts_%druns" % (g1func[0], g1func[-1], g2func[0], g2func[-1], DIM, len(INST_LIST), NUM_RUNS)
            fnameEcdf     = fname + "_ecdf_%s.pdf" % ind
            fnameEcdfLogx = fname + "_ecdf_%s_logx.pdf" % ind
            
            # Plot to the PDF files
            #pltObject.plotting(DIR_FIGS + "/" + fnameEcdf,     ecdf_coord, xstart=startFevLim/DIM, logx=False, yfull=fullRangeY)
            pltObject.plotting(DIR_FIGS + "/" + fnameEcdfLogx, ecdf_coord, xstart=startFevLim/DIM, logx=True,  yfull=fullRangeY)
            #-------------------------------------------------------------------
            

            
            
            # Define a layout of PDFs and include them in TeX file
            try: name_buf
            except NameError:
                name_buf = []
            name_buf = name_buf + ([fnameEcdfLogx] if compileLogx else [fnameEcdf])
            if len(name_buf) == 4:
                composeLatex(cmd='add', file=ftex_main_figs_pdf,
                             stringlist=["../" + os.path.split(DIR_FIGS)[1] + "/" + name_buf[0],
                                         "../" + os.path.split(DIR_FIGS)[1] + "/" + name_buf[1],
                                         "../" + os.path.split(DIR_FIGS)[1] + "/" + name_buf[2],
                                         "../" + os.path.split(DIR_FIGS)[1] + "/" + name_buf[3]])
                name_buf = []
            #-------------------------------------------------------------------
    
    composeLatex(cmd='fini', file=ftex_main_figs_pdf)
    
    # Compile the main LaTeX file to PDF
    if compileLatex:
        texify.pdfLatex(directory=DIR_TEX_FIGS, filename=MAIN_TEX_FIGS, cleanup=[".log", ".aux"])


def plottingConfig():
    # Plotting commands start here
    pltObject = fancyplot.PlotWithFancyLegend(figsize=(7,5)) # changing this may effect all!!!
    pltObject.position = [0.112, 0.132, 0.635, 0.803] # [left, bottom, width, height]
    pltObject.lineWidth  = 2.7
    pltObject.markerSize = 10
    pltObject.nMarkers   = 5
    pltObject.xExtRatio    = 10/100.
    pltObject.xSegLenRatio = 60/100.  # 70/100.
    pltObject.yShrinkRatio = 30/100.
    pltObject.labelBottomTop = [0.85, 0.85] # a factor determines the the spread of the right legends
    # Newly added properties
    pltObject.topFontSize        = 21 # fontsize for plot title
    pltObject.bottomFontSize     = 22 # fontsize for x-label
    pltObject.leftFontSize       = 20 # fontsize for y-label
    pltObject.rightFontSize      = 17 # fontsize for legends on the right
    pltObject.annotationFontSize = 21 # fontsize for additional texts inside the plot
    pltObject.xTickLabelSize = 17 # fontsize for xtick labels (i.e. numbers on the x-axis)
    pltObject.yTickLabelSize = 17 # fontsize for ytick labels (i.e. numbers on the y-axis)
    pltObject.gridLineWeight = 0.5
    return pltObject


def main(ind = 'eps'):
    #dims = [2,3,5,10,20]
    #isHV = False
    #ecdf_processing_alldims()
    #ecdf_processing_alldims(isHV)
    # MPI stuff
    #nproc = MPI.COMM_WORLD.Get_size()
    #iproc = MPI.COMM_WORLD.Get_rank()
    
    #for idx, d in enumerate(dims):
    #  if idx % nproc == iproc:
    ecdf_processing(ind = ind)
    ecdf_processing_allindicators()    
    #MPI.COMM_WORLD.Barrier()


if __name__ == "__main__":
    startTime = time.time()
    if len(sys.argv) > 1:
      main(ind = sys.argv[1])
    else:
      print "syntax: python main_plot_aggregated_ecdf.py hv|eps|gd|igd"
    print "It took:", time.time() - startTime, "seconds."
    print "Done." 
    # print "Press Enter to continue ..." 
    # raw_input()
