""" 
	wrapper for several c functions useful for multi-objective optimization computation
	Most of the functions deal accepts a C double array as input
	using the numpy.ctypeslib to compute the epsilon indicator additive 
	of the approximation set with respect to the reference set

	Abdullah Al-Dujaili adapted from:
	http://www.scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html#id1
"""

import numpy as np
import numpy.ctypeslib as npct
from ctypes import c_int
from ctypes import c_uint
from ctypes import c_double
from sys import platform as platform
import os
import hv
from hv import HyperVolume
##########################################################################
# define input/output type
array_1d_double = npct.ndpointer(dtype=np.double, flags='C_CONTIGUOUS')
array_1d_bool = npct.ndpointer(dtype=np.bool,  flags='C_CONTIGUOUS')
##########################################################################
# load the libraries, using numpy mechanisms : use the name *.so
if platform == "darwin":
	libeps = npct.load_library("libeps", ".")
	libpf  = npct.load_library("libpf",".")
	libhv  = npct.load_library("libhv",".")
	libgd  = npct.load_library("libgd",".")
else:
	libeps = npct.load_library("libeps.so", ".") 
	libpf  = npct.load_library("libpf.so",".")
	libhv  = npct.load_library("libhv.so",".")
	libgd  = npct.load_library("libgd.so",".")
##########################################################################	
# define the arg and res type for each of the functions:
# 1. epsilon library
libeps.eps.restype = c_double
libeps.eps.argtypes = [array_1d_double, array_1d_double, c_int, c_int, c_int]
libeps.incr_eps.restype = None
libeps.incr_eps.argtypes = [array_1d_double,array_1d_double, array_1d_double, c_int, c_int, c_int]
libeps.fast_incr_eps.restype = None
libeps.fast_incr_eps.argtypes = [array_1d_double,array_1d_double, array_1d_double, c_int, c_int, c_int]
# 2. pareto front library
libpf.pf_selective.restype = None
libpf.pf_selective.argtypes = [array_1d_bool, array_1d_double, c_int, c_int]
libpf.pf_cao.restype = None
libpf.pf_cao.argtypes = [array_1d_bool, array_1d_double, c_int, c_int]
# 3. hv library
libhv.incr_hv.restype = None
libhv.incr_hv.argtypes = [array_1d_double, array_1d_double, c_int, c_int]
# 4. GD library (GD,IGD)
libgd.gd.restype = c_double
libgd.gd.argtypes = [array_1d_double, array_1d_double, c_int, c_int, c_int]
libgd.igd.restype = c_double
libgd.igd.argtypes = [array_1d_double, array_1d_double, c_int, c_int, c_int]
libgd.incr_gd.restype = None
libgd.incr_gd.argtypes = [array_1d_double,array_1d_double, array_1d_double, c_int, c_int, c_int]
libgd.incr_igd.restype = None
libgd.incr_igd.argtypes = [array_1d_double,array_1d_double, array_1d_double, c_int, c_int, c_int]
##########################################################################
# Python wrappers:
def compute_pyhv(approximation_set, reference_point):
	"""
		returns the hypervolume of the approximation set with respect to the reference point based on Simon Wessing's code
	"""
	#referencePoint = [2, 2, 2]
	hv = HyperVolume(reference_point)
	#front = [[1,0,1], [0,1,0]]
	return hv.compute(approximation_set)

def compute_incr_hv_c(approximation_set, reference_point):
	"""
		returns the hypervolume of the approximation set with respect to the reference point based on Eckart Zitzler's code.
		currently not working
	"""
	#incr_hv_val = np.array([0.0]* approximation_set.shape[0])
	#normalized_approx_set = reference_point - approximation_set
	#libhv.incr_hv(incr_hv_val,np.ascontiguousarray(normalized_approx_set),approximation_set.shape[0], approximation_set.shape[1])
	#return incr_hv_val
	print("This is not working at the moment. It seems that: \n calling a c function, which calls other c functions, \n from python does not work out. \n The c source is available under src/hypervol.c\n Please email me if you work around it")
	
def compute_incr_hv(approximation_set, reference_point):
	"""
		returns a vector of the hypervolume indicator values for incremental subsets of the approximation set with respect to the reference point.
		i.e. incr_hv[m] = hv(approximation_set[:m+1,:], reference_point)
		This function is based on the python-based hypervolume code.
		TO DO: integrate it with the C-based implementation.
	"""	
	incr_hv_val = np.array([0.0]* approximation_set.shape[0])
	for i in xrange(len(incr_hv_val)):
		incr_hv_val[i]= compute_pyhv(approximation_set[:i+1], reference_point)
	return incr_hv_val
		
def compute_gd(approximation_set, reference_set):
	"""
		returns the generational distance indicator value of the approximation set with respect to the reference set
	"""
	#assert approximation_set.shape[1] == reference_set.shape[1]
	return libgd.gd(approximation_set, reference_set, approximation_set.shape[0], reference_set.shape[0], approximation_set.shape[1])


def compute_incr_gd(approximation_set, reference_set):
	"""
		returns a vector of the generational distance indicator values for incremental subsets of the approximation set with respect to the reference set.
		i.e. incr_gd[m] = gd(approximation_set[:m+1,:], reference_set)
	"""
	#assert approximation_set.shape[1] == reference_set.shape[1]
	incr_gd = np.array([0.0]* approximation_set.shape[0])
	libgd.incr_gd(incr_gd, np.ascontiguousarray(approximation_set), np.ascontiguousarray(reference_set), approximation_set.shape[0], reference_set.shape[0], approximation_set.shape[1])
	return incr_gd

def compute_igd(approximation_set, reference_set):
	"""
		returns the inverted generational distance indicator value of the approximation set with respect to the reference set
	"""
	#assert approximation_set.shape[1] == reference_set.shape[1]
	return libgd.igd(approximation_set, reference_set, approximation_set.shape[0], reference_set.shape[0], approximation_set.shape[1])


def compute_incr_igd(approximation_set, reference_set):
	"""
		returns a vector of the inverted generational distance indicator values for incremental subsets of the approximation set with respect to the reference set.
		i.e. incr_igd[m] = igd(approximation_set[:m+1,:], reference_set)
	"""
	#assert approximation_set.shape[1] == reference_set.shape[1]
	incr_igd = np.array([0.0]* approximation_set.shape[0])
	libgd.incr_igd(incr_igd, np.ascontiguousarray(approximation_set), np.ascontiguousarray(reference_set), approximation_set.shape[0], reference_set.shape[0], approximation_set.shape[1])
	return incr_igd

	
def compute_eps(approximation_set, reference_set):
	"""
		returns the additive epsilon indicator value of the approximation set with respect to the reference set
	"""
	#assert approximation_set.shape[1] == reference_set.shape[1]
	return libeps.eps(approximation_set, reference_set, approximation_set.shape[0], reference_set.shape[0], approximation_set.shape[1])


def compute_incr_eps(approximation_set, reference_set):
	"""
		returns a vector of the additive epsilon indicator values for incremental subsets of the approximation set with respect to the reference set.
		i.e. incr_epsilon[m] = epsilon(approximation_set[:m+1,:], reference_set)
	"""
	#assert approximation_set.shape[1] == reference_set.shape[1]
	incr_epsilon = np.array([0.0]* approximation_set.shape[0])
	libeps.incr_eps(incr_epsilon, np.ascontiguousarray(approximation_set), np.ascontiguousarray(reference_set), approximation_set.shape[0], reference_set.shape[0], approximation_set.shape[1])
	return incr_epsilon
  
  
def compute_fast_incr_eps(approximation_set, reference_set):
	"""
		returns a vector of the additive epsilon indicator values for incremental subsets of the approximation set with respect to the reference set.
		i.e. incr_epsilon[m] = epsilon(approximation_set[:m+1,:], reference_set) but in a much faster way at the cost of extra space
	"""
	#assert approximation_set.shape[1] == reference_set.shape[1]
	incr_epsilon = np.array([0.0]* approximation_set.shape[0])
	libeps.fast_incr_eps(incr_epsilon, np.ascontiguousarray(approximation_set), np.ascontiguousarray(reference_set), approximation_set.shape[0], reference_set.shape[0], approximation_set.shape[1])
	return incr_epsilon

def paretofront_cao(in_array):
	"""
	Returns the logical Pareto membership of a set of points.
	Takes a numpy array of nrows x ncols and returns a boolean vector of nrows x 1
	where 1 denotes a non-dominated vector, according to cao's method
	Yi Cao: y.cao@cranfield.ac.uk
	"""
	bool_array = np.array([False] * in_array.shape[0], dtype= np.bool)
	libpf.pf_cao(bool_array, np.ascontiguousarray(in_array.T), in_array.shape[0], in_array.shape[1])
	return bool_array


def paretofront(in_array):
	"""
	Returns the logical Pareto membership of a set of points
	Takes a numpy array of nrows x ncols and returns a boolean vector of nrows x 1
	where 1 denotes a non-dominated vector.
	"""
	bool_array = np.array([True] * in_array.shape[0], dtype= np.bool)
	libpf.pf_selective(bool_array, np.ascontiguousarray(in_array), in_array.shape[0], in_array.shape[1])
	return bool_array
##########################################################################
	