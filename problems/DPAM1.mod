###############################################################################
#
#   As described by Huband et al. in "A review of multiobjective test problems
#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
#   Computing 10(5): 477-506, 2006.
#
#   Example DPAM1, see the previous cited paper for the original reference.
#
#   In the above paper the number of variables was left undefined. 
#   We selected n=10 as default.
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
param n := 10;
param A{1..n,1..n}; # := Uniform(-1,1); # Fixed

# x variable
var x{1..n} >=-0.3, <=0.3;

# y variable
var y{i in 1..n} = sum {j in 1..n} A[i,j] * x[j];

# g(x)
var gx = 1 + 10*(n-1) + sum {i in 2..n} (y[i]^2-10*cos(4*pi*y[i]));

minimize f1:
	y[1];
minimize f2:
    gx*exp(-y[1]/gx);

data;


# Random generated matrix, but fixed for the numerical results
param A:        1            2            3           4            5           6
 :=
1     0.218418    -0.620254     0.843784    0.914311    -0.788548    0.428212
2     0.330423     0.151614     0.884043   -0.272951    -0.993822    0.511197
3     0.180332    -0.593814    -0.492722    0.0646786   -0.666503   -0.945716
4    -0.0265389   -0.920133     0.308861   -0.0437502   -0.374203    0.207359
5    -0.88565     -0.375906    -0.708948   -0.37902      0.576578    0.0194674
6     0.238261    -0.1596      -0.827302    0.669248     0.494475    0.691715
7    -0.218632    -0.865161    -0.715997    0.220772     0.692356    0.646453
8    -0.207987    -0.865931     0.613732   -0.525712    -0.995728    0.389633
9     0.60624      0.0951648   -0.160446   -0.394585    -0.167581    0.0679849
10    0.404396     0.449996     0.162711    0.294454    -0.563345   -0.114993

:        7            8             9            10         :=
1     0.103064    -0.47373     -0.300792     -0.185507
2    -0.0997948   -0.659756     0.575496      0.675617
3    -0.334582     0.611894     0.281032      0.508749
4    -0.219433     0.914104     0.184408      0.520599
5    -0.470262     0.572576     0.351245     -0.480477
6    -0.198585     0.0492812    0.959669      0.884086
7    -0.401724     0.615443    -0.0601957    -0.748176
8    -0.064173     0.662131    -0.707048     -0.340423
9     0.449799     0.733505    -0.00918638    0.00446808
10    0.549589    -0.775141     0.677726      0.610715;
