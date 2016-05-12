###############################################################################
#
#   As described by Huband et al. in "A review of multiobjective test problems
#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
#   Computing 10(5): 477-506, 2006.
#
#   Example FES1, see the previous cited paper for the original reference.
#
#   In the above paper the number of variables was left undefined. 
#   We selected n=10 as default.
#
#   This file is part of a collection of problems developed for
#   derivative-free multiobjective optimization in
#   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
#   Direct Multisearch for Multiobjective Optimization, 2010.
#
#   Written by the authors in June 1, 2010.
#
###############################################################################

# parameter
param pi := 4*atan(1);
param n := 10;

# x variable
var x{1..n} >=0, <=1;

minimize f1:
    sum {i in 1..n} (abs(x[i]-exp((i/n)^2)/3)^0.5);
minimize f2:
    sum {i in 1..n} ((x[i]-0.5*cos(10*pi*i/n)-0.5)^2);
