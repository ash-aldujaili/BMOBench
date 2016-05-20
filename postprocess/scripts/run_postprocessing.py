#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division
import os
import sys # to access command-line arguments
import shutil # to delete a folder
import time

import main_extract_ecdf as extract_ecdf
import main_plot_aggregated_ecdf as plot_aggregated
import main_plot_single_ecdf as plot_single
def main():
    # ===============================================================
    # PUT YOUR SETTINGS HERE:
    data = []
    data.append({
        'name': "MODIRECT",
        'algo': "MO-DIRECT",
        'dir' : "../../EXP_RESULTS",
        'nRun': 1,
        'nFev': 1e3,
    })
    # use one or more of the following indicators: ['eps', 'gd', 'igd', 'hv']
    indicators = ['eps', 'gd', 'igd', 'hv']
    # ================================================================
    for ind in indicators:
        extract_ecdf.ecdf_extracting(ind = ind, isMean = False, data = data)
        plot_aggregated.main(ind = ind)
        plot_single.main(ind = ind)
    
    # Process all the indicators/ all problems in a single plot
    plot_aggregated.ecdf_processing_allindicators(indicators = indicators) 
    # move the aggregate data profiles to the paper template:
    NUM_RUNS = 10
    suffix = "_%druns_aggregate" % (NUM_RUNS)
    DIR_SRC = "../" + "postproc" + "/" + "figs%s" % suffix
    DIR_DST = "../" + "../" + "latex-template" + "/" + "performance-data"
    if os.path.exists(DIR_DST):
        shutil.rmtree(DIR_DST)
    shutil.move(DIR_SRC, DIR_DST)
if __name__ == "__main__":
    startTime = time.time()
    main()
    print "It took:", time.time() - startTime, "seconds."
    print "Done." 
    # print "Press Enter to continue ..." 
    # raw_input()