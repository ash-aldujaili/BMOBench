###############################################################################
#
#   As described by Huband et al. in "A review of multiobjective test problems
#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
#   Computing 10(5): 477-506, 2006.
#           
#   Example WFG3.
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

# WFG3
param A {i in 1..M-1} := if i<=1 then 1 else 0;


# problem variables
param zmax {i in 1..n} := 2*i;
var z{i in 1..n} >=0, <= zmax[i];

# transform z into [0,1] set
var y{i in 1..n} = z[i]/zmax[i];

# first level mapping
var t1{i in 1..n} = if i <= k then y[i]
      else abs(y[i]-0.35)/abs(floor(0.35-y[i])+0.35);

# second level mapping
param AA := 2;
var t2{i in 1..k+l/2} = if i<=k then t1[i]
    else sum {ii in (k+2*(i-k)-1)..(k+2*(i-k))} (t1[ii]+sum {jj in 0..AA-2} abs(t1[ii]-t1[(k+2*(i-k)-1)+((ii+jj-(k+2*(i-k)-1)+1) mod (k+2*(i-k)-(k+2*(i-k)-1)+1))]))/((k+2*(i-k)-(k+2*(i-k)-1)+1)/AA*ceil(AA/2)*(1+2*AA-2*ceil(AA/2)));

# third level mapping
param w{i in 1..n} := 1;
var t3{i in 1..M} = if i<=M-1 then (sum {j in ((i-1)*k/(M-1)+1)..(i*k/(M-1))} (w[j]*t2[j]))/(sum {j in ((i-1)*k/(M-1)+1)..(i*k/(M-1))} w[j])
    else (sum {j in k+1..k+l/2} (w[j]*t2[j]))/(sum {j in k+1..k+l/2} w[j]);


# Define objective function variables
var x{i in 1..M} = if i<=M-1 then max(t3[M],A[i])*(t3[i]-0.5)+0.5
    else t3[M];

# Define objective function function h
var h{m in 1..M} = if m==1 then prod {i in 1..M-1} x[i]
    else if m<=M-1 then (prod {i in 1..M-m} x[i])*(1-x[M-m+1])
        else 1-x[1];


# The objective functions
minimize fobj {m in 1..M}:
    x[M]+S[m]*h[m];


