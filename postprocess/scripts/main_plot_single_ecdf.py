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
    #FUNC_LIST = range(1,25) # from f1 to f24
    #DIM = dimension
    #xlim = [-5, 5]
    #INST_LIST = ((2,4), (3,5), (7,8), (9,10), (11,12))
    NUM_RUNS = 10
    NUM_TARGETS = 70
    #======= TO BE EDITED =====================================================      
    FILENAME_ECDF_DATABASE = 'db_ecdfs_%s_maxfevalDiff_%druns_%dtargets.pkl'
    suffix = "_%druns_single" % (NUM_RUNS)
    #==========================================================================
    
    fullRangeY    = True
    compileLatex  = True
    showLogOnly   = True
    
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
    
    #=============================================#
    #  Loop over MOPs and combine BBOB functions  #
    #=============================================#
    
    #=============================================#
    #  Loop over MOPs #
    #=============================================#
    for pidx, pid in enumerate(PROBLEM_LIST):
        DIM = DIM_LIST[pidx]
            
        print "Working on %s" % (pid)
        
        # STEPS TO PRODUCE PLOTTABLE DATA FOR A *SINGLE* ECDF
        #-------------------------------------------------------------------
        # Compute proportional data for ECDF (i.e. y-axis)
        # dict is mutable, thus passed by reference, no need of a return
        ecdf.induce_proportion_multi_algos(all_ecdfs[pidx])
        
        print 'nprofiles:', all_ecdfs[pidx][all_ecdfs[pidx].keys()[0]]['nprofiles']
        print 'ndata:', len(all_ecdfs[pidx][all_ecdfs[pidx].keys()[0]]['data'])
        print '================================================'
        
        # Adjust the ecdf data such that we can plot from ${startFevPlot}
        startFevPlot = 1 # DIM - 1
        startFevLim  = 1
        ecdf_coord = ecdf.extract_plot_data(all_ecdfs[pidx], startFevPlot)
        
        # Modify x's such that we can display the x-axis as "# f-evals / DIM"
        for algo in ecdf_coord.keys():
            ecdf_coord[algo]['x'] = [item for item in ecdf_coord[algo]['x']]
        #-------------------------------------------------------------------
        
        # Load algorithm run information
        info = [all_ecdfs[pidx][algo]['info'] for algo in all_ecdfs[pidx].keys()]
        
        # Plot commands start here
        pltObject = plottingConfig() # configure the plot
        #pltObject.legends = [key[3:] for key in ecdf_coord.keys()]
        pltObject.legends = [key[3:].replace('MO-CMA-ES','MO-CMA') for key in ecdf_coord.keys()]
        #pltObject.xlabel = "# f-evals / DIM"
        pltObject.xlabel = "# $\mathbf{f}$-evaluations / dimension"
        #pltObject.ylabel = "%d targets * %d runs" % (len(list_of_target_values), info[0]['nRun'])
        pltObject.ylabel = "Proportion of targets in %d trials" % (NUM_RUNS)
        #pltObject.title = ur"$f_{%d.%d}$ & $f_{%d.%d}$ (%dD)" % (obj1_id, obj1_instance, obj2_id, obj2_instance, DIM)
        #pltObject.title = ur"$f_{%d}$ & $f_{%d}$, %dD" % (obj1_id, obj2_id, DIM)
        pltObject.title = "dimension = %d" % DIM
        pltObject.note = ["%s" % pid]
        
        # Construct output filenames
        fname = "%s_%dD_%s_%druns" % \
                (pid, DIM, "3algos", info[0]['nRun'])
        fname = fname.replace(".", "-")
        fnameEcdf     = fname + "_ecdf_%s.pdf" % ind
        fnameEcdfLogx = fname + "_ecdf_%s_logx.pdf" % ind
        
        # Plot to the PDF files
        #pltObject.plotting(DIR_FIGS + "/" + fnameEcdf,     ecdf_coord, xstart=startFevLim/DIM, logx=False, yfull=fullRangeY)
        pltObject.plotting(DIR_FIGS + "/" + fnameEcdfLogx, ecdf_coord, xstart=startFevLim, logx=True,  yfull=fullRangeY)
        #-------------------------------------------------------------------
        
        # Define a layout of PDFs and include them in TeX file
        try: name_buf
        except NameError:
            name_buf = []
        name_buf = name_buf + ([fnameEcdfLogx] if showLogOnly else [fnameEcdf, fnameEcdfLogx])
        if len(name_buf) == 4: # four figs in a row
            if showLogOnly:
                composeLatex(cmd='add', file=ftex_main_figs_pdf,
                                stringlist=["../" + os.path.split(DIR_FIGS)[1] + "/" + name_buf[0],
                                            "../" + os.path.split(DIR_FIGS)[1] + "/" + name_buf[1],
                                            "../" + os.path.split(DIR_FIGS)[1] + "/" + name_buf[2],
                                            "../" + os.path.split(DIR_FIGS)[1] + "/" + name_buf[3]])
            else:
                composeLatex(cmd='add', file=ftex_main_figs_pdf,
                                stringlist=["../" + os.path.split(DIR_FIGS)[1] + "/" + name_buf[0],
                                            "../" + os.path.split(DIR_FIGS)[1] + "/" + name_buf[2],
                                            "../" + os.path.split(DIR_FIGS)[1] + "/" + name_buf[1],
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
    pltObject.annotationFontSize = 23 # fontsize for additional texts inside the plot
    pltObject.xTickLabelSize = 17 # fontsize for xtick labels (i.e. numbers on the x-axis)
    pltObject.yTickLabelSize = 17 # fontsize for ytick labels (i.e. numbers on the y-axis)
    pltObject.gridLineWeight = 0.5
    return pltObject


def main(ind= 'eps'):
    #dims = [2,3,5,10,20]
    #isHV = False
    #ecdf_processing_alldims()
    # MPI stuff
    #nproc = MPI.COMM_WORLD.Get_size()
    #iproc = MPI.COMM_WORLD.Get_rank()
    
    #for idx, d in enumerate(dims):
    #  if idx % nproc == iproc:
    ecdf_processing(ind = ind)
        
    #MPI.COMM_WORLD.Barrier()


if __name__ == "__main__":
    startTime = time.time()
    if len(sys.argv) > 1:
      main(ind = sys.argv[1])
    else:
      print "syntax: python main_plot_single_ecdf.py hv|eps|gd|igd"
    print "It took:", time.time() - startTime, "seconds."
    print "Done." 
    # print "Press Enter to continue ..." 
    # raw_input()
