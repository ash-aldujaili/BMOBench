"""
A setup script for multi-objective c-2-python related functionalities

Abdullah Al-Dujaili, 2016

"""

import os
from os import system

from sys import platform as platform



# Define the compilation flags
if platform == "darwin": # mac os
	LDFLAGS = " -fPIC -std=c99 -dynamiclib "
else: # other systems tested so far
	LDFLAGS = " -fPIC -std=c99 -shared "



print "Cleaning files"
system("rm -vf *.so *.o *~")

print "Compiling libs"
system("gcc" + LDFLAGS + " c-indicators/pf.c -o libpf.so")
system("gcc" + LDFLAGS + " c-indicators/eps.c -o libeps.so")
system("gcc" + LDFLAGS + " c-indicators/hypervol.c -o libhv.so")
system("gcc" + LDFLAGS + " c-indicators/igd.c -o libgd.so")
