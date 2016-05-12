#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import numpy as np
import csv


def readtxt(filename, ichunk="all"):
    """Read the data from a formatted text file.
    A commented line must start with a #. Blank lines are permitted.
    A continuous chunk can be regarded as the final population of a run,
    different runs are separated by commented lines or blank lines.
    ichunk="all": the whole numerical data in the file are returned.
    ichunk=[idx]: only corresponding chunks (for those runs) are returned.
    """
    with open(filename) as f:
        buf = csv.reader(f, delimiter='\t', skipinitialspace=True)
        buf = list(buf)
        
        groupedBuf = []
        startNewChunk = True
        for line in buf:
            if line and (line[0][0] is not "#"): # if the i-th line is not an empty list, i.e. [], nor a commented line
                if startNewChunk is True:
                    chunk = []
                    startNewChunk = False
                chunk.append(line)
            elif startNewChunk is False:
                groupedBuf.append(chunk) # store the whole chunk as a list in the mother list groupedBuf
                startNewChunk = True
        if startNewChunk is False: # as the last chunk hasn't stored if there is no empty line at the end of file
            groupedBuf.append(chunk)
        
        if ichunk is "all": # Take ALL groups together
            buf = [item for sublist in groupedBuf for item in sublist]
        else: # Take only SOME of the groups
            igroup = ichunk # [1,3] # [0]
            buf = [item for sublist in [groupedBuf[i] for i in igroup] for item in sublist]
        
        # for i in reversed(xrange(len(buf))): # reversed check as buf is to be trimmed
        #     if (not buf[i]): # if the i-th line is an empty list, i.e. []
        #         del buf[i]
        #     elif buf[i][0][0] is "#": # if the i-th line is a commented line
        #         del buf[i]
                
        # USE THIS in case tabs and spaces for separating numbers are mixed
        for i in xrange(len(buf)): buf[i] = [s for string in buf[i] for s in string.split()]
        # buf = map(lambda row: filter(None, row), buf) # remove trailing tab (i.e. empty string in the list): fastest way?
        
        mydata = np.array(buf); mydata = mydata.astype(np.float)
    return mydata
    


if __name__ == "__main__":
    print "I am a module and should only be called from another program!" 
    # raw_input()
