###############################################################################
#
#   As described by Huband et al. in "A review of multiobjective test problems
#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
#   Computing 10(5): 477-506, 2006.
#
#   Example MOP3, Van Valedhuizen's test suit.
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
param A1 := 0.5*sin(1)-2*cos(1)+sin(2)-1.5*cos(2);
param A2 := 1.5*sin(1)-cos(1)+2*sin(2)-0.5*cos(2);

# Define objective function variables
var x{1..2} >= -pi, <=pi;

# Additional variables
var B1 = 0.5*sin(x[1])-2*cos(x[1])+sin(x[2])-1.5*cos(x[2]);
var B2 = 1.5*sin(x[1])-cos(x[1])+2*sin(x[2])-0.5*cos(x[2]);


# The objective functions
#maximize fobj1:
#    -1-(A1-B1)^2-(A2-B2)^2;

#maximize fobj2:
#    -(x[1]+3)^2-(x[2]+1)^2;

minimize fobj1:
    -(-1-(A1-B1)^2-(A2-B2)^2);

minimize fobj2:
    -(-(x[1]+3)^2-(x[2]+1)^2);