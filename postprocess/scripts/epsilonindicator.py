""" 
	wrapper for epsilonindicator.c that accepts a C double array as input
	using the numpy.ctypeslib to compute the epsilon indicator additive 
	of the approximation set with respect to the reference set

	Abdullah Al-Dujaili adapted from:
	http://www.scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html#id1
"""

import numpy as np
import numpy.ctypeslib as npct
from ctypes import c_int
from ctypes import c_double
from sys import platform as platform
import os
# input/outpu type for the epsilon indicator function
array_1d_double = npct.ndpointer(dtype=np.double, flags='C_CONTIGUOUS')


# load the library, using numpy mechanisms : use the name *.so
if platform == "darwin":
	libeps = npct.load_library("dlibepsilonindicator", ".") # classical 
	#libincreps = npct.load_library("libincrepsilonindicator.so", ".") # incremental 
	#libfeps = npct.load_library("lib_fast_eps.so", ".") # fast
else:
  libeps = npct.load_library("libepsilonindicator.so", ".") # classical 
	#libincreps = npct.load_library("libincrepsilonindicator", ".") # incremental 
	#libfeps = npct.load_library("lib_fast_eps", ".") # fast

# classical function
libeps.epsilonindicator.restype = c_double
libeps.epsilonindicator.argtypes = [array_1d_double, array_1d_double, c_int, c_int, c_int]

# incremental
libeps.incremental_epsilonindicator.restype = None
libeps.incremental_epsilonindicator.argtypes = [array_1d_double,array_1d_double, array_1d_double, c_int, c_int, c_int]

# fast (function name)
libeps.fast_eps_ind.restype = None
libeps.fast_eps_ind.argtypes = [array_1d_double,array_1d_double, array_1d_double, c_int, c_int, c_int]



def compute_eps(approximation_set, reference_set):
	"""
		returns the additive epsilon indicator value of the approximation set with respect to the reference set
	"""
	#assert approximation_set.shape[1] == reference_set.shape[1]
	return libeps.epsilonindicator(approximation_set, reference_set, approximation_set.shape[0], reference_set.shape[0], approximation_set.shape[1])


def compute_incr_eps(approximation_set, reference_set):
	"""
		returns a vector of the additive epsilon indicator values for incremental subsets of the approximation set with respect to the reference set.
		i.e. incr_epsilon[m] = epsilon(approximation_set[:m+1,:], reference_set)
	"""
	#assert approximation_set.shape[1] == reference_set.shape[1]
	incr_epsilon = np.array([0.0]* approximation_set.shape[0])
	libeps.incremental_epsilonindicator(incr_epsilon, np.ascontiguousarray(approximation_set), np.ascontiguousarray(reference_set), approximation_set.shape[0], reference_set.shape[0], approximation_set.shape[1])
	return incr_epsilon
  
  
def compute_fast_incr_eps(approximation_set, reference_set):
	"""
		returns a vector of the additive epsilon indicator values for incremental subsets of the approximation set with respect to the reference set.
		i.e. incr_epsilon[m] = epsilon(approximation_set[:m+1,:], reference_set) but in a much faster way at the cost of extra space
	"""
	#assert approximation_set.shape[1] == reference_set.shape[1]
	incr_epsilon = np.array([0.0]* approximation_set.shape[0])
	libeps.fast_eps_ind(incr_epsilon, np.ascontiguousarray(approximation_set), np.ascontiguousarray(reference_set), approximation_set.shape[0], reference_set.shape[0], approximation_set.shape[1])
	return incr_epsilon

