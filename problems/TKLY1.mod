###############################################################################
#
#   As described by Huband et al. in "A review of multiobjective test problems
#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
#   Computing 10(5): 477-506, 2006.
#
#   Example TKLY1, see the previous cited paper for the original reference.
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
var x{1..4};

minimize f1:
    x[1];
minimize f2:
    (prod {i in 2..4} (2.0-exp(-((x[i]-0.1)/0.004)^2)-0.8*exp(-((x[i]-0.9)/0.4)^2)))/x[1];
    
subject to bound:
    0.1<=x[1]<=1;

subject to bounds{i in 2..4}:
    0<=x[i]<=1;

