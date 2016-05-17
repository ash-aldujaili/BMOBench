###############################################################################
#
#   As described by Huband et al. in "A review of multiobjective test problems
#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
#   Computing 10(5): 477-506, 2006.
#
#   Example MOP4, Van Valedhuizen's test suit.
#
#   This file is part of a collection of problems developed for
#   derivative-free multiobjective optimization in
#   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
#   Direct Multisearch for Multiobjective Optimization, 2010.
#
#   Written by the authors in June 1, 2010.
#
###############################################################################

# Define objective function variables
var x{1..3} >= -5, <=5;


# The objective functions
minimize fobj1:
    sum {i in 1..2} (-10*exp(-0.2*sqrt(x[i]^2+x[i+1]^2)));

minimize fobj2:
    sum {i in 1..3} (abs(x[i])^0.8+5*sin(x[i]^3));

