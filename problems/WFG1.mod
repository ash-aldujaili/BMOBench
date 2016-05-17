###############################################################################
#
#   As described by Huband et al. in "A review of multiobjective test problems
#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
#   Computing 10(5): 477-506, 2006.
#           
#   Example WFG1.
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
var t1{i in 1..n} = if i <= k then y[i] 
      else abs(y[i]-0.35)/abs(floor(0.35-y[i])+0.35);


# second level mapping
param AA := 0.8;
param BB := 0.75;
param CC := 0.85;
var t2{i in 1..n} = if i<=k then t1[i]
    else AA+ min(0,floor(t1[i]-BB))*(AA*(BB-t1[i]))/BB-min(0,floor(CC-t1[i]))*(1-AA)*(t1[i]-CC)/(1-CC);

# third level mapping
param AAA := 0.02;
var t3{i in 1..n} = t2[i]^AAA;

# forth level mapping
param w{i in 1..n} := 2*i;
var t4{i in 1..M} = if i<=M-1 then (sum {j in ((i-1)*k/(M-1)+1)..(i*k/(M-1))} (w[j]*t3[j]))/(sum {j in ((i-1)*k/(M-1)+1)..(i*k/(M-1))} w[j])
    else (sum {j in k+1..n} (w[j]*t3[j]))/(sum {j in k+1..n} w[j]);



# Define objective function variables
var x{i in 1..M} = if i<=M-1 then max(t4[M],A[i])*(t4[i]-0.5)+0.5
    else t4[M];

# Define objective function function h
param alpha :=1;
param AAAA :=5;
var h{m in 1..M} = if m==1 then prod {i in 1..M-1} (1-cos(x[i]*pi2))
    else if m<=M-1 then (prod {i in 1..M-m} (1-cos(x[i]*pi2)))*(1-sin(x[M-m+1]*pi2))
        else (1-x[1]-(cos(2*AAAA*pi*x[1]+pi2))/(2*AAAA*pi))^alpha;


# The objective functions
minimize fobj {m in 1..M}:
    x[M]+S[m]*h[m];

