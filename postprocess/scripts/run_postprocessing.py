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
        'name': "MORANDOM",
        'algo': "MORANDOM",
        'dir' : "../../EXP_RESULTS",
        'nRun': 1,
        'nFev': 1e3,
    })
    isHV = False
    # ================================================================
    
    extract_ecdf.ecdf_extracting(isHV = isHV, isMean = True, data = data)
    plot_aggregated.main(isHV = isHV)
    plot_single.main(isHV = isHV)
    
    
if __name__ == "__main__":
    startTime = time.time()
    main()
    print "It took:", time.time() - startTime, "seconds."
    print "Done." 
    # print "Press Enter to continue ..." 
    # raw_input()