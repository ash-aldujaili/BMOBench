###############################################################################
#
#   As described by C.M. Fonseca and P.J. Fleming in "Multiobjective
#   Optimization and Multiple Constraint Handling with Evolutionary
#   Algorithms–Part I: A Unified Formulation", in IEEE Transactions 
#   on Systems, Man, and Cybernetics—Part A: Systems and Humans, 
#   vol. 28, no. 1, January 1998.
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
var x{1..2};


minimize f1:
    1-exp(-(x[1]-1)^2-(x[2]+1)^2);
minimize f2:
    1-exp(-(x[1]+1)^2-(x[2]-1)^2);


subject to bounds {i in 1..2}:
    -4.0 <= x[i] <= 4.0;
