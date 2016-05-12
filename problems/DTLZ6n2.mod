###############################################################################
#
#   As described by K. Deb, L. Thiele, M. Laumanns and E. Zitzler in "Scalable
#   Multi-Objective Optimization Test Problems", Congress on Evolutionary 
#   Computation (CEC’2002): 825-830, 2002.
#
#   Example DTLZ6 with M=2 and n=2.
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
param M := 2, >=2;  # Number of objectives
param n := 2, >= M; # Number of variables
param k := n-M+1;
param pi := 4*atan(1);

# x variable
var x{1..n};

# g(x)
var gx = 1 + 9/k * sum {i in M..n} (x[i]);

# functions
var ffM = (1+gx) * (M - sum{i in 1..M-1} (x[i]/(1+gx)*(1+sin(3*pi*x[i]))));

minimize f {i in 1..M-1}:
    x[i];
minimize fM:
    ffM;


subject to bounds {i in 1..n}:
    0.0 <= x[i] <= 1.0;
