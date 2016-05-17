###############################################################################
#
#   As described by Huband et al. in "A review of multiobjective test problems
#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
#   Computing 10(5): 477-506, 2006.
#
#   Example MOP6, Van Valedhuizen's test suit.
#
#   This file is part of a collection of problems developed for
#   derivative-free multiobjective optimization in
#   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
#   Direct Multisearch for Multiobjective Optimization, 2010.
#
#   Written by the authors in June 1, 2010.
#
###############################################################################

# Pi
param pi := 4*atan(1);

# Define objective function variables
var x{1..2} >= 0, <=1;


# The objective functions
minimize fobj1:
    x[1];

minimize fobj2:
    (1+10*x[2])*(1-(x[1]/(1+10*x[2]))^2-x[1]/(1+10*x[2])*sin(8*pi*x[1]));


