###############################################################################
#
#   As described by Huband et al. in "A review of multiobjective test problems
#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
#   Computing 10(5): 477-506, 2006.
#
#   Example QV1, see the previous cited paper for the original reference.
#
#   In the original reference the number of variables was n=16. 
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

# parameters
param pi := 4*atan(1);
param n  := 10;

# x variable
var x{1..n} >=-5.12, <= 5.12;

minimize f1:
    (sum {i in 1..n} (x[i]^2-10*cos(2*pi*x[i])+10)/n)^0.25;
minimize f2:
    (sum {i in 1..n} ((x[i]-1.5)^2-10*cos(2*pi*(x[i]-1.5))+10)/n)^0.25;

