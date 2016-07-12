#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# From MOBBOB-GECCO, Brockhoff, Dimo, Thanh-Do Tran, and Nikolaus Hansen. "Benchmarking numerical multiobjective optimizers revisited." Proceedings of the 2015 Annual Conference on Genetic and Evolutionary Computation. ACM, 2015.*/


def pdfLatex(directory, filename, cleanup=[]):
    """Control the OS to compile the main LaTeX file using 'pdflatex'
    and clean up junk files of types listed in the list 'cleanup'
    """
    import os
    workingDir = os.getcwd()
    os.chdir(directory)
    if os.path.isfile(filename):
        os.system("pdflatex %s" % (filename))
        # Cleaning up after compilation
        for ext in cleanup:
            f = filename[:-4] + ext
            os.remove(f)
            print "Cleaned up %s" % f
    else:
        print "\nThe required main TeX file '%s' not found. Skip compiling LaTeX." % (directory + "/" + filename)
    os.chdir(workingDir)
    

def epsLatex(directory, filename, cleanup=[]):
    """Control the OS to compile the main LaTeX file using: 'latex'-->'latex'-->'dvips'-->'ps2pdf'
    and clean up junk files of types listed in the list 'cleanup'
    """
    import os
    workingDir = os.getcwd()
    os.chdir(directory)
    if os.path.isfile(filename):
        os.system("latex %s" % (filename))
        os.system("latex %s" % (filename))
        os.system("dvips %s" % (filename[:-4] + '.dvi'))
        os.system("ps2pdf %s" % (filename[:-4] + '.ps'))
        # Cleaning up after compilation
        for ext in cleanup:
            f = filename[:-4] + ext
            os.remove(f)
            print "Cleaned up %s" % f
    else:
        print "\nThe required main TeX file '%s' not found. Skip compiling LaTeX." % (directory + "/" + filename)
    os.chdir(workingDir)
    

if __name__ == "__main__":
    print "I am a module and should only be called from another program!" 
    # raw_input()
