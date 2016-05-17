###############################################################################
#
#   As described by A. Lovison in "A synthetic approach to multiobjective
#   optimization", arxiv Item: http://arxiv.org/abs/1002.0093.
#
#   Example 5.
#
#   In the above paper/papers the variables bounds were not set.
#   We considered -1<=x[i]<=4, i=1,2,3.
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
param n := 3;           # number of variables
param m := 3;           # number of function used in the objectives
param C{1..n,1..m};     #:=Uniform(-1,1); # distinct non collinear points
param alpha{1..m,1..n}; #:=Uniform(0,1) # negative definite matrix diagonal
param beta{1..m};       #:=Uniform(-1,1);
param gamma{1..m};      #:=Uniform(-1,1);
param pi := 4*atan(1);

# x variable
var x{1..n} >= -1, <= 4;
var f{j in 1..m} = sum{i in 1..n} (-alpha[j,i]*(x[i]-C[i,j])^2);

#maximize u1:
#	f[1];
#maximize u2:
#    f[2]+beta[2]*sin(pi*(x[1]+x[2])/gamma[2]);
#maximize u3:
#    f[3]+beta[3]*cos(pi*(x[1]-x[2])/gamma[3]);
	
minimize u1:
	-f[1];
minimize u2:
    -(f[2]+beta[2]*sin(pi*(x[1]+x[2])/gamma[2]));
minimize u3:
    -(f[3]+beta[3]*cos(pi*(x[1]-x[2])/gamma[3]));	

data;

# Uniformly random generated data. Fixed in the numerical results
param C :=
1 1    0.218418
1 2   -0.620254
1 3    0.843784
2 1    0.914311
2 2   -0.788548
2 3    0.428212
3 1    0.103064
3 2   -0.47373
3 3   -0.300792;

param alpha :=
1 1   0.407247
1 2   0.665212
1 3   0.575807
2 1   0.942022
2 2   0.363525
2 3   0.00308876
3 1   0.755598
3 2   0.450103
3 3   0.170122;

param beta :=
1  0.575496
2  0.675617
3  0.180332;

param gamma :=
1  -0.593814
2  -0.492722
3   0.0646786;
