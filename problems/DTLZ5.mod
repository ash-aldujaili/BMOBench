###############################################################################
#
#   As described by K. Deb, L. Thiele, M. Laumanns and E. Zitzler in "Scalable
#   Multi-Objective Optimization Test Problems", Congress on Evolutionary 
#   Computation (CEC’2002): 825-830, 2002.
#
#   Example DTLZ5.
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
param M := 3, >=2;  # Number of objectives
param n := 12, >= M; # Number of variables
param k := n-M+1;
param pi := 4*atan(1);

# x variable
var x{1..n};

# g(x)
var gx = sum {i in M..n} (x[i]^0.1);

# theta
var theta {i in 2..n} = (pi/2)*(1+2*gx*x[i])/(2*(1+gx));

# functions
var ff1 = (1+gx)*(cos(0.5*pi*x[1]))*prod {j in 2..M-1} (cos(theta[j]));
var ff {i in 2..M-1} = (1+gx) * (cos(0.5*pi*x[1])) * (prod {j in 2..M-i} (cos(theta[j])))*(sin(theta[M-i+1]));
var ffM = (1+gx) * sin(0.5*pi*x[1]);

minimize f1:
    ff1;
minimize f {i in 2..M-1}:
	ff[i];
minimize fM:
    ffM;


subject to bounds {i in 1..n}:
    0.0 <= x[i] <= 1.0;
