###############################################################################
#
#   As described by K. Deb and A. Sinha and S. Kukkonen in "Multi-objective
#   test problems, linkages, and evolutionary methodologies", GECCO'06}: 
#   Proceedings of the 8th Annual Conference on Genetic and Evolutionary 
#   Computation, 1141-1148, 2006.
#
#   Example T4, with linkage L1.
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
param m := 10, >=2;
param pi := 4*atan(1);
param A{1..m,1..m};

# x variable
var x{1..m};

# y variable
var y{i in 1..m} = sum {j in 1..m} A[i,j] * x[j];

# functions
var ff1 = y[1]^2;

# g(x)
var gx = 1 + 10*(m-1) + sum {i in 2..m} (y[i]^2-10*cos(4*pi*y[i]));

var h = 1-sqrt(ff1/gx);

minimize f1:
	ff1;
minimize f2:
    gx*h;

subject to bounds1:
	0.0 <= x[1] <= 1.0;

subject to bounds {i in 2..m}:
	-5.0 <= x[i] <= 5.0;

data;

param A :        1            2           3            4           5            6
 :=
1     1     0           0            0           0            0
2     0     0.884043   -0.272951    -0.993822    0.511197    -0.0997948
3     0    -0.492722    0.0646786   -0.666503   -0.945716    -0.334582
4     0     0.308861   -0.0437502   -0.374203    0.207359    -0.219433
5     0    -0.708948   -0.37902      0.576578    0.0194674   -0.470262
6     0    -0.827302    0.669248     0.494475    0.691715    -0.198585
7     0    -0.715997    0.220772     0.692356    0.646453    -0.401724
8     0     0.613732   -0.525712    -0.995728    0.389633    -0.064173
9     0    -0.160446   -0.394585    -0.167581    0.0679849    0.449799
10    0     0.162711    0.294454    -0.563345   -0.114993     0.549589

:        7             8             9            10        :=
1     0            0             0             0
2    -0.659756     0.575496      0.675617      0.180332
3     0.611894     0.281032      0.508749     -0.0265389
4     0.914104     0.184408      0.520599     -0.88565
5     0.572576     0.351245     -0.480477      0.238261
6     0.0492812    0.959669      0.884086     -0.218632
7     0.615443    -0.0601957    -0.748176     -0.207987
8     0.662131    -0.707048     -0.340423      0.60624
9     0.733505    -0.00918638    0.00446808    0.404396
10   -0.775141     0.677726      0.610715      0.0850755;
