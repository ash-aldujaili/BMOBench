###############################################################################
#
#   As described by Huband et al. in "A review of multiobjective test problems
#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
#   Computing 10(5): 477-506, 2006.
#
#   Example Far1, see the previous cited paper for the original reference.
#
#   This file is part of a collection of problems developed for
#   derivative-free multiobjective optimization in
#   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
#   Direct Multisearch for Multiobjective Optimization, 2010.
#
#   Written by the authors in June 1, 2010.
#
###############################################################################

# x variable
var x{1..2} >=-1, <=1;

minimize f1:
	-2*exp(15*(-(x[1]-0.1)^2-x[2]^2))
    -exp(20*(-(x[1]-0.6)^2-(x[2]-0.6)^2))
    +exp(20*(-(x[1]+0.6)^2-(x[2]-0.6)^2))
    +exp(20*(-(x[1]-0.6)^2-(x[2]+0.6)^2))
    +exp(20*(-(x[1]+0.6)^2-(x[2]+0.6)^2));
minimize f2:
    2*exp(20*(-x[1]^2-x[2]^2))
    +exp(20*(-(x[1]-0.4)^2-(x[2]-0.6)^2))
    -exp(20*(-(x[1]+0.5)^2-(x[2]-0.7)^2))
    -exp(20*(-(x[1]-0.5)^2-(x[2]+0.7)^2))
    +exp(20*(-(x[1]+0.4)^2-(x[2]+0.8)^2));
