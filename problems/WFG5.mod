###############################################################################
#
#   As described by Huband et al. in "A review of multiobjective test problems
#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
#   Computing 10(5): 477-506, 2006.
#           
#   Example WFG5.
#
#   This file is part of a collection of problems developed for
#   derivative-free multiobjective optimization in
#   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
#   Direct Multisearch for Multiobjective Optimization, 2010.
#
#   Written by the authors in June 1, 2010.
#
###############################################################################

param M := 3;
param k := 4;
param l := 4;
param n := k+l;

param pi := 4*atan(1);
param pi2:= 2*atan(1);

param S {m in 1..M} := 2*m;

# neq WFG3
param A {i in 1..M-1} := 1;


# problem variables
param zmax {i in 1..n} := 2*i;
var z{i in 1..n} >=0, <= zmax[i];

# transform z into [0,1] set
var y{i in 1..n} = z[i]/zmax[i];

# first level mapping
param AA := 0.35;
param BB := 0.001;
param CC := 0.05;
var t1{i in 1..n} = 1+(abs(y[i]-AA)-BB)*((floor(y[i]-AA+BB)*(1-CC+(AA-BB)/BB))/(AA-BB)+(floor(AA+BB-y[i])*(1-CC+(1-AA-BB)/BB))/(1-AA-BB)+1/BB);

# second level mapping
param w{i in 1..n} := 1;
var t2{i in 1..M} = if i<=M-1 then (sum {j in ((i-1)*k/(M-1)+1)..(i*k/(M-1))} (w[j]*t1[j]))/(sum {j in ((i-1)*k/(M-1)+1)..(i*k/(M-1))} w[j])
    else (sum {j in k+1..n} (w[j]*t1[j]))/(sum {j in k+1..n} w[j]);


# Define objective function variables
var x{i in 1..M} = if i<=M-1 then max(t2[M],A[i])*(t2[i]-0.5)+0.5
    else t2[M];

# Define objective function function h
var h{m in 1..M} = if m==1 then prod {i in 1..M-1} sin(x[i]*pi2)
    else if m<=M-1 then (prod {i in 1..M-m} sin(x[i]*pi2))*cos(x[M-m+1]*pi2)
        else cos(x[1]*pi2);


# The objective functions
minimize fobj {m in 1..M}:
    x[M]+S[m]*h[m];


