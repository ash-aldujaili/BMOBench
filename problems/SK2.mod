###############################################################################
#
#   As described by Huband et al. in "A review of multiobjective test problems
#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
#   Computing 10(5): 477-506, 2006.
#
#   Example SK2, see the previous cited paper for the original reference.
#
#   In the above paper/papers the variables bounds were not set.
#   We considered -10<=x[i]<=10, i=1,2,3,4.
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
var x{1..4} >=-10, <= 10;

#maximize f1:
#    -(x[1]-2)^2-(x[2]+3)^2-(x[3]-5)^2-(x[4]-4)^2+5;
#maximize f2:
#    (sum {i in 1..4} sin(x[i]))/(1+(sum {i in 1..4} x[i]^2)/100);

minimize f1:
    -(-(x[1]-2)^2-(x[2]+3)^2-(x[3]-5)^2-(x[4]-4)^2+5);
minimize f2:
    -((sum {i in 1..4} sin(x[i]))/(1+(sum {i in 1..4} x[i]^2)/100));