###############################################################################
#
#   As described by F.Y. Cheng and X.S. Li, "Generalized center method for
#   multiobjective engineering optimization", Engineering Optimization,31:5,
#   641-661, 1999.
#
#   Example 2, four bar truss.
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
param F :=  10;
param E := 2*10^5;
param L := 200;
param sigma := 10;

# x variable
var x{1..4};


minimize f1:
    L*(2*x[1]+sqrt(2)*x[2]+sqrt(x[3])+x[4]);
minimize f2:
	F*L/E*(2/x[1]+(2*sqrt(2))/x[2]-(2*sqrt(2))/x[3]+2/x[4]);


subject to bound1:
    F/sigma <= x[1] <= 3*F/sigma;
subject to bound2:
    sqrt(2)*F/sigma <= x[2] <= 3*F/sigma;
subject to bound3:
    sqrt(2)*F/sigma <= x[3] <= 3*F/sigma;
subject to bound4:
    F/sigma <= x[4] <= 3*F/sigma;
