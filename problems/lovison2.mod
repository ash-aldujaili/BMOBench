###############################################################################
#
#   As described by A. Lovison in "A synthetic approach to multiobjective
#   optimization", arxiv Item: http://arxiv.org/abs/1002.0093.
#
#   Example 2.
#
#   In the above paper/papers the variables bounds were not set.
#   We considered -0.5<=x[1]<=0 and -0.5<=x[2]<=0.5.
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

#maximize f1:
#	-x[2];
#maximize f2:
#    (x[2]-x[1]^3)/(x[1]+1);
	
minimize f1:
	x[2];
minimize f2:
    -(x[2]-x[1]^3)/(x[1]+1);
    
subject to xbound1:
    -0.5 <= x[1] <= 0;
    
subject to xbound2:
    -0.5 <= x[2] <= 0.5;

