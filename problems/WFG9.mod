###############################################################################
#
#   As described by Huband et al. in "A review of multiobjective test problems
#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
#   Computing 10(5): 477-506, 2006.
#           
#   Example WFG9.
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
param w{1..n}, default 1.0;
param AA := 0.98/49.98;
param BB := 0.02;
param CC := 50;
var r_sum{i in 1..n-1} = (sum {j in i+1..n} (w[j]*y[j]))/(sum {j in i+1..n} w[j]);
var t1{i in 1..n} = if i<=n-1 then y[i]^(BB+(CC-BB)*(AA-(1-2*r_sum[i])*abs(floor(0.5-r_sum[i])+AA)))
    else y[i];

# second level mapping
param AAA := 0.35;
param BBB := 0.001;
param CCC := 0.05;
param AAAA := 30;
param BBBB := 95;
param CCCC := 0.35;
var t2{i in 1..n} = if i<=k then 1+(abs(t1[i]-AAA)-BBB)*((floor(t1[i]-AAA+BBB)*(1-CCC+(AAA-BBB)/BBB))/(AAA-BBB)+(floor(AAA+BBB-t1[i])*(1-CCC+(1-AAA-BBB)/BBB))/(1-AAA-BBB)+1/BBB)
    else (1+cos((4*AAAA+2)*pi*(0.5-abs(t1[i]-CCCC)/(2*(floor(CCCC-t1[i])+CCCC))))+4*BBBB*(abs(t1[i]-CCCC)/(2*(floor(CCCC-t1[i])+CCCC)))^2)/(BBBB+2);

# third level mapping
var t3{i in 1..M} = if i<=M-1 then sum {ii in ((i-1)*k/(M-1)+1)..(i*k/(M-1))} (t2[ii]+sum {jj in 0..(k/(M-1)-2)} abs(t2[ii]-t2[((i-1)*k/(M-1)+1)+((ii+jj-((i-1)*k/(M-1)+1)+1) mod ((i*k/(M-1))-((i-1)*k/(M-1)+1)+1))]))/(((i*k/(M-1))-((i-1)*k/(M-1)+1)+1)/(k/(M-1))*ceil(k/(M-1)/2)*(1+2*k/(M-1)-2*ceil(k/(M-1)/2)))
    else sum {ii in k+1..n} (t2[ii]+sum {jj in 0..(l-2)} abs(t2[ii]-t2[k+1+((ii+jj-(k+1)+1) mod (n-k))]))/(((n-k)/l)*ceil(l/2)*(1+2*l-2*ceil(l/2)));

# Define objective function variables
var x{i in 1..M} = if i<=M-1 then max(t3[M],A[i])*(t3[i]-0.5)+0.5
    else t3[M];

# Define objective function function h
var h{m in 1..M} = if m==1 then prod {i in 1..M-1} sin(x[i]*pi2)
    else if m<=M-1 then (prod {i in 1..M-m} sin(x[i]*pi2))*cos(x[M-m+1]*pi2)
        else cos(x[1]*pi2);


# The objective functions
minimize fobj {m in 1..M}:
    x[M]+S[m]*h[m];


