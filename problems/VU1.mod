###############################################################################
#
#   As described by Huband et al. in "A review of multiobjective test problems
#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
#   Computing 10(5): 477-506, 2006.
#
#   Example VU1, see the previous cited paper for the original reference.
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
var x{1..2} >= -3, <=3;

minimize f1:
    1/(x[1]^2+x[2]^2+1);
minimize f2:
    x[1]^2+3*x[2]^2+1;


