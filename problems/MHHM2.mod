###############################################################################
#
#   As described by Huband et al. in "A review of multiobjective test problems
#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
#   Computing 10(5): 477-506, 2006.
#
#   Example MHHM2, see the previous cited paper for the original reference.
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
var x{1..2} >=0, <=1;

minimize f1:
    (x[1]-0.8)^2+(x[2]-0.6)^2;
minimize f2:
    (x[1]-0.85)^2+(x[2]-0.7)^2;
minimize f3:
    (x[1]-0.9)^2+(x[2]-0.6)^2;

