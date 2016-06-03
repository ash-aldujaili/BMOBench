#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define printf mexPrintf

/*
 *
 *  As described by F.Y. Cheng and X.S. Li, "Generalized center method for
 *  multiobjective engineering optimization", Engineering Optimization,31:5,
 *  641-661, 1999.
 *
 *  Example 2, four bar truss.
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void CL1(double* f, double* x){
    double F = 10;
    double E = 2*100000;
    double L = 200;
    double sigma = 10;

    double sqrt2 = sqrt(2);

    if( (x[0] >= (F/sigma) && x[0] <= (3*F)/sigma) && (x[1] >= ((sqrt2*F)/sigma) && x[1] <= ((3*F)/sigma)) && (x[2] >= ((sqrt2*F)/sigma) && x[2] <= (3*F/sigma)) && (x[3] >= (F/sigma) && x[3] <= ((3*F)/sigma)));
    else{
        printf("Invalid input range\n");
        return;
    }

    f[0] = L*(2*x[0]+sqrt(2)*x[1]+sqrt(x[2])+x[3]);
    f[1] = F*L/E*(2/x[0]+(2*sqrt(2))/x[1]-(2*sqrt(2))/x[2]+2/x[3]);
    return;
}

/*
 *
 *  As described by K. Deb in "Multi-objective Genetic Algorithms: Problem
 *  Difficulties and Construction of Test Problems", Evolutionary Computation
 *  7(3): 205-230, 1999.
 *
 *  Example 5.1.2 (Non-convex Pareto-optimal Front).
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void Deb512b(double* f, double* x){
    double beta = 1;
    double alpha = 4;

    double ff1 = 4*x[0];
    double gx, h;

    if( (x[0]<0 || x[0]>1 || x[1]<0 || x[1]>1) ) {
        printf("Invalid input range\n");
        return;
    }


    gx = (x[1]<=0.4) ? (4 - 3*exp(-1*pow((x[1]-0.2)/0.02, 2))) : (4 - 2*exp(-1*pow((x[1]-0.7)/0.2, 2)));
    h = (ff1<=beta*gx) ? (1-pow(ff1/(beta*gx),alpha)) : 0;

    f[0] = ff1;
    f[1] = gx*h;

    return;
}

/*
 *
 *  As described by K. Deb in "Multi-objective Genetic Algorithms: Problem
 *  Difficulties and Construction of Test Problems", Evolutionary Computation
 *  7(3): 205-230, 1999.
 *
 *  Example 5.2.1 (Biased Search Space).
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void Deb521b(double* f, double* x){

    double gamma = 1;

    double ff1 = x[0];
    double gx, h;

    if( (x[0]<0 || x[0]>1 || x[1]<0 || x[1]>1) ) {
        printf("Invalid input range\n");
        return;
    }

    gx = 1 + pow(x[1], gamma);
    h = 1 - pow(ff1/gx, 2);

    f[0] = ff1;
    f[1] = gx*h;

    return;
}

/*
 *
 *  As described by K. Deb, L. Thiele, M. Laumanns and E. Zitzler in "Scalable
 *  Multi-Objective Optimization Test Problems", Congress on Evolutionary
 *  Computation (CEC�2002): 825-830, 2002.
 *
 *  Example DTLZ1 with M=2 and n=2.
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void DTLZ1n2(double* f, double* x){

    int M = 2; // Number of Objectives; >=2;
    int n = 2; // Number of variables; >=M;
    int k = n-M+1;
    double pi = 4*atan(1);

    int i = 0;
    double sum = 0.0;
    double gx;

    double prod = 1.0;
    double ff1;

    for(i=0; i<n; ++i) {
        if( x[i]<0 || x[i]>1 ) {
            printf("Invalid input range\n");
            return;
        }
    }

    for(i=M-1; i<n; ++i)
        sum += pow(x[i]-0.5, 2) - cos(20*pi*(x[i]-0.5));

    gx = 100*(k+sum);

    for(i=0; i<M-1; ++i)
        prod *= x[i];

    ff1 = 0.5*(1+gx)*prod;

    f[0] = ff1;

    /*for(i=1; i<M; ++i){
      prod = 1;
      for(j=0; j<M-i; ++j){
        prod *= x[j];
      }
      f[i] = 0.5*(1+gx)*prod*(1-x[M-i-1]);
    } *//*Let's just leave this for now as long as M=2 */

    prod = 1;
    for(i=0; i<M-2; ++i)
      prod*=x[i];
    f[M-1] = 0.5*(1+gx)*prod*(1-x[0]);

    return;
}

/*
 *
 *  As described by K. Deb, L. Thiele, M. Laumanns and E. Zitzler in "Scalable
 *  Multi-Objective Optimization Test Problems", Congress on Evolutionary
 *  Computation (CEC�2002): 825-830, 2002.
 *
 *  Example DTLZ3 with M=2 and n=2.
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void DTLZ3n2(double* f, double* x){

    int M = 2; // Number of Objectives; >=2;
    int n = 2; // Number of variables; >=M;
    int k = n-M+1;
    double pi = 4*atan(1);

    int i = 0;
    double sum = 0.0;
    double gx, ff1;
    double prod = 1.0;

    for(i=0; i<n; ++i) {
        if( x[i]<0 || x[i]>1 ) {
            printf("Invalid input range\n");
            return;
        }
    }

    for(i=M-1; i<n; ++i)
        sum += pow(x[i]-0.5, 2) - cos(20*pi*(x[i]-0.5));

    gx = 100*(k+sum);

    for(i=0; i<M-1; ++i)
        prod *= cos(0.5*pi*x[i]);

    ff1 = (1+gx)*prod;

    f[0] = ff1;

    /*for(i=1; i<M; ++i){
      prod = 1;
      for(j=0; j<M-i-1; ++j){
        prod *= cos(0.5*pi*x[j]);
      }
      f[i] = (1+gx)*prod*(sin(0.5*pi*x[M-i]));
    } */ /*Let's just leave this for now as long as M=2 */

    prod = 1;
    for(i=0; i<M-2; ++i)
      prod*=cos(0.5*pi*x[i]);
    f[M-1] = (1+gx)*prod*(sin(0.5*pi*x[0]));

    return;
}

/*
 *
 *  As described by K. Deb, L. Thiele, M. Laumanns and E. Zitzler in "Scalable
 *  Multi-Objective Optimization Test Problems", Congress on Evolutionary
 *  Computation (CEC�2002): 825-830, 2002.
 *
 *  Example DTLZ5 with M=2 and n=2.
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void DTLZ5n2(double* f, double* x){

    int M = 2; // Number of Objectives; >=2;
    int n = 2; // Number of variables; >=M;
    double pi = 4*atan(1);

    double gx = 0.0;
    double* theta = (double*)malloc(sizeof(double)*(n));
    double prod = 1.0;
    double ff1;

    int i = 0;
    int j = 0;

    for(i=0; i<n; ++i) {
        if( x[i]<0 || x[i]>1 ) {
            printf("Invalid input range\n");
            return;
        }
    }

    for(i=M-1; i<n; ++i)
        gx += pow(x[i], 0.1);

    for(i=1; i<n; ++i)
      theta[i] = (pi/2)*(1+2*gx*x[i])/(2*(1+gx));

    for(i=1; i<M-1; ++i)
        prod *= cos(theta[i]);

    ff1 = (1+gx)*cos(0.5*pi*x[0])*prod;

    f[0] = ff1;

    if(M-1>1){
      for(i=1; i<M-1; ++i){
        prod = 1;
        for(j=1; j<M-i-1; ++j){
          prod *= cos(theta[j]);
        }
        f[i] = (1+gx)*cos(0.5*pi*x[0])*prod*(sin(theta[M-i]));
      }
    }

    f[M-1] = (1+gx)*sin(0.5*pi*x[0]);

    return;
}

/*
 *
 *  As described by Huband et al. in "A review of multiobjective test problems
 *  and a scalable test problem toolkit", IEEE Transactions on Evolutionary
 *  Computing 10(5): 477-506, 2006.
 *
 *  Example Far1, see the previous cited paper for the original reference.
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void Far1(double* f, double* x){

    int i;

    for(i=0; i<2; ++i) {
        if( x[i]<-1 || x[i]>1 ) {
            printf("Invalid input range\n");
            return;
        }
    }

    f[0] = -2*exp(15*(-1*pow(x[0]-0.1, 2) - pow(x[1], 2)))
              - exp(20*(-1*pow(x[0]-0.6, 2) - pow(x[1]-0.6, 2)))
              + exp(20*(-1*pow(x[0]+0.6, 2) - pow(x[1]-0.6, 2)))
              + exp(20*(-1*pow(x[0]-0.6, 2) - pow(x[1]+0.6, 2)))
              + exp(20*(-1*pow(x[0]+0.6, 2) - pow(x[1]+0.6, 2)));
    f[1] = 2*exp(20*(-1*pow(x[0], 2) - pow(x[1], 2)))
              + exp(20*(-1*pow(x[0]-0.4, 2) - pow(x[1]-0.6, 2)))
              - exp(20*(-1*pow(x[0]+0.5, 2) - pow(x[1]-0.7, 2)))
              - exp(20*(-1*pow(x[0]-0.5, 2) - pow(x[1]+0.7, 2)))
              + exp(20*(-1*pow(x[0]+0.4, 2) - pow(x[1]+0.8, 2)));;

    return;
}

/*
 *
 *  As described by C.M. Fonseca and P.J. Fleming in "Multiobjective
 *  Optimization and Multiple Constraint Handling with Evolutionary
 *  Algorithms�Part I: A Unified Formulation", in IEEE Transactions
 *  on Systems, Man, and Cybernetics�Part A: Systems and Humans,
 *  vol. 28, no. 1, January 1998.
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void Fonseca(double* f, double* x){

    int i;

    for(i=0; i<2; ++i) {
        if( x[i]<-4 || x[i]>4 ) {
            printf("Invalid input range\n");
            return;
        }
    }

    f[0] = 1 - exp(-1*pow(x[0]-1, 2) - pow(x[1]+1, 2));
    f[1] = 1 - exp(-1*pow(x[0]+1, 2) - pow(x[1]-1, 2));

    return;
}

/*
 *
 *  As described by Huband et al. in "A Scalable Multi-objective Test Problem
 *  Toolkit", in C. A. Coello Coello et al. (Eds.): EMO 2005, LNCS 3410,
 *  pp. 280�295, 2005, Springer-Verlag Berlin Heidelberg 2005.
 *
 *  Example I4.
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void I4(double* f, double* zz){

  /*Note: Input variable here is zz and not x */

    int M = 3;
    double k = 4;
    double l = 4;
    int n = k+l;

    double pi2 = 2*atan(1);

    double *S, *A;
    double *zmax, *x, *y, *t1, *t2, *t3, *h, *z;
    int i, j, ii, jj;
    double sum1, sum2;

    z = (double*)malloc(sizeof(double)*(n+1));
    for(i=n; i>=1; --i)
      z[i] = zz[i-1];

    S = (double*)malloc(sizeof(double)*(M+1));
    A = (double*)malloc(sizeof(double)*(M));

    for(i=1; i<=M-1; ++i){
      S[i] = 1;
      A[i] = 1;
    }
    S[M] = 1;

    zmax = (double*)malloc(sizeof(double)*(n+1));
    y = (double*)malloc(sizeof(double)*(n+1));
    t1 = (double*)malloc(sizeof(double)*(n+1));
    t2 = (double*)malloc(sizeof(double)*(n+1));
    t3 = (double*)malloc(sizeof(double)*(M+1));
    x = (double*)malloc(sizeof(double)*(M+1));
    h = (double*)malloc(sizeof(double)*(M+1));

    for(i=1; i<=n; ++i){
      zmax[i] = 1;
      y[i] = z[i]/zmax[i];
      t1[i] = y[i];
      if(i<=k)
        t2[i] = t1[i];
      else
        t2[i] = ((fabs(t1[i]-0.35)/fabs(floor(0.35-t1[i])+0.35)));
    }

    for(i=1; i<=n; ++i) {
        if( z[i]<0 || z[i]>zmax[i] ) {
            printf("%lf, Invalid input range\n", z[i]);
            return;
        }
    }

    for(i=1; i<=M; ++i){
      if(i<=(M-1)){
        sum1 = 0;

        for(ii=((i-1)*k/(M-1)+1); ii<=(i)*k/(M-1); ++ii){
          sum2 = 0;
          for(jj=0; jj<=(k/(M-1)-2); ++jj){
            sum2 += fabs(t2[ii] - t2[(int)(((i-1)*k/(M-1)+1)+((((int)(ii+jj-((i-1)*k/(M-1)+1)+1))) % (((int)(((i)*k/(M-1))-((i-1)*k/(M-1)+1)+1)))))]);
          }
          sum1 += (t2[ii] + sum2);
        }
        t3[i] = sum1/(((i*k/(M-1))-((i-1)*k/(M-1)+1)+1)/(k/(M-1))*ceil(k/(M-1)/2)*(1+2*k/(M-1)-2*ceil(k/(M-1)/2)));
      }
      else{
        sum1 = 0;
        sum2 = 0;
        for(ii=k+1; ii<=n; ++ii){
          sum2 = 0;
          for(jj=0; jj<=l-2; ++jj){
            sum2 += fabs(t2[ii] - t2[(int)(k+1+(((int)(ii+jj-(k+1)+1)) % ((int)(n-k))))]);
          }
          sum1 += t2[ii] + sum2;
        }
        t3[i] = sum1/(((n-k)/l)*ceil(l/2)*(1+2*l-2*ceil(l/2)));
      }
    }
    for(i=1; i<=M; ++i){
      if(i<=M-1){
        if(t3[M]>A[i])
          x[i] = t3[M]*(t3[i]-0.5)+0.5;
        else
          x[i] = A[i]*(t3[i]-0.5)+0.5;
      }
      else
        x[i] = t3[M];
    }

    for(i=0; i<=M; ++i){
      if(i==1){
        h[i] = 1;
        for(j=1; j<=M-1; ++j){
          h[i] *= sin(x[j]*pi2);
        }
      }
      else if(i<=M-1){
        h[i] = 1;
        for(j=1; j<=M-i; ++j)
          h[i] *= sin(x[j]*pi2);
        h[i] *= cos(x[M-i+1]*pi2);
      }
      else{
        h[i] = cos(x[1]*pi2);
      }
    }
    // f = (double*)malloc(sizeof(double)*M);
    for(i=1; i<=M; ++i){
      f[i-1] = x[M]+S[i]*h[i];
    }

    return;
}

/*
 *
 *  As described by Y. Jin, M. Olhofer and B. Sendhoff. "Dynamic weighted
 *  aggregation for evolutionary multi-objective optimization: Why does it
 *  work and how?", in Proceedings of Genetic and Evolutionary Computation
 *  Conference, pp.1042-1049, San Francisco, USA, 2001.
 *
 *  Test function 1, F1.
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void Jin1(double* f, double* x){

    int n = 2;
    int i;

    for(i=0; i<n; ++i) {
        if( x[i]<0 || x[i]>1 ) {
            printf("Invalid input range\n");
            return;
        }
    }

    f[0] = 0;
    f[1] = 0;

    for(i=0; i<n; ++i){
      f[0] += x[i]*x[i];
      f[1] += (x[i]-2)*(x[i]-2);
    }
    f[0] = f[0]/n;
    f[1] = f[1]/n;

    return;
}

/*
 *
 *  As described by F. Kursawe in "A variant of evolution strategies for
 *  vector optimization", in H. P. Schwefel and R. M�nner, editors, Parallel
 *  Problem Solving from Nature, 1st Workshop, PPSN I, volume 496 of Lecture
 *  Notes in Computer Science, pages 193-197, Berlin, Germany, Oct 1991,
 *  Springer-Verlag.
 *
 *  In the above paper the variables bounds were not set.
 *  We considered -5.0 <= x[i] <= 5.0, i=1,2,3.
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void Kursawe(double* f, double* x){

    int n = 3;
    int i;

    for(i=0; i<n; ++i) {
        if( x[i]<-5 || x[i]>5 ) {
            printf("Invalid input range\n");
            return;
        }
    }

    f[0] = 0;
    f[1] = 0;
    for(i=0; i<n-1; ++i){
      f[0] += (-10*exp(-0.2*sqrt(x[i]*x[i] + x[i+1]*x[i+1])));
      f[1] += (pow(fabs(x[i]), 0.8) + 5*sin(x[i])*sin(x[i])*sin(x[i]));
    }
    f[1] += (pow(fabs(x[n-1]), 0.8) + 5*sin(x[n-1])*sin(x[n-1])*sin(x[n-1]));

    return;
}

/*
 *
 *  As described by K. Deb and A. Sinha and S. Kukkonen in "Multi-objective
 *  test problems, linkages, and evolutionary methodologies", GECCO'06}:
 *  Proceedings of the 8th Annual Conference on Genetic and Evolutionary
 *  Computation, 1141-1148, 2006.
 *
 *  Example T4, with linkage L1.
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void L1ZDT4(double* f, double* x){
  int m = 10;
  double *y;
  double ff1, gx, h;
  double pi = 4*atan(1);

  int i, j;
  double A[100]={1,0,0,0,0,0,0,0,0,0,0,0.884043,-0.272951,-0.993822,0.511197,-0.0997948,-0.659756,0.575496,0.675617,0.180332,0,-0.492722,0.0646786,-0.666503,-0.945716,-0.334582,0.611894,0.281032,0.508749,-0.0265389,0,0.308861,-0.0437502,-0.374203,0.207359,-0.219433,0.914104,0.184408,0.520599,-0.88565,0,-0.708948,-0.37902,0.576578,0.0194674,-0.470262,0.572576,0.351245,-0.480477,0.238261,0,-0.827302,0.669248,0.494475,0.691715,-0.198585,0.0492812,0.959669,0.884086,-0.218632,0,-0.715997,0.220772,0.692356,0.646453,-0.401724,0.615443,-0.0601957,-0.748176,-0.207987,0,0.613732,-0.525712,-0.995728,0.389633,-0.064173,0.662131,-0.707048,-0.340423,0.60624,0,-0.160446,-0.394585,-0.167581,0.0679849,0.449799,0.733505,-0.00918638,0.00446808,0.404396,0,0.162711,0.294454,-0.563345,-0.114993,0.549589,-0.775141,0.677726,0.610715,0.0850755,};
  /*
  FILE *myFile;
  double **A = (double **)malloc(sizeof(double*)*m);
  myFile = fopen("./data/L1ZDT4.txt", "r");

  for(i = 0; i < m; i++){
    A[i] = (double *)malloc(sizeof(double)*m);
  }

	for (i=0; i<m; i++)
  {
    for (j=0; j<m; j++)
    {
      if(!	fscanf(myFile, "%lf\n", &A[i][j]))
	     break;
    }
  }
	fclose(myFile);
	*/
  if(x[0]<0 || x[0]>1){
    printf("Invalid input range\n");
    return;
  }

  for(i = 1; i < m; i++){
    if(x[i] >= -5 && x[i] <= 5)
      continue;
    else{
      printf("Invalid input range\n");
      return;
    }
  }

  y = (double*)malloc(sizeof(double)*m);
  for(i = 0; i < m; i++){
    y[i] = 0;
    for(j = 0; j < m; j++){
        y[i] = y[i] + A[i*m+j]*x[j];
    }
  }

  ff1 = y[0]*y[0];
  gx = 0;

  for(i = 1; i < m; i++){
    gx = gx + (y[i]*y[i] - 10*cos(4*pi*y[i]) );
  }
  gx = gx + 10*(m-1);
  gx = gx + 1;

  h = 1 - sqrt(ff1/gx);

  f[0] = ff1;
  f[1] = gx*h;

  return;
}

/*
 *
 *  As described by K. Deb and A. Sinha and S. Kukkonen in "Multi-objective
 *  test problems, linkages, and evolutionary methodologies", GECCO'06}:
 *  Proceedings of the 8th Annual Conference on Genetic and Evolutionary
 *  Computation, 1141-1148, 2006.
 *
 *  Example T1, with linkage L3.
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void L3ZDT1(double* f, double* x){
  int m = 30;
  double *y;
  double ff1, gx, h;

  double M[900] = {0.218418,-0.620254,0.843784,0.914311,-0.788548,0.428212,0.103064,-0.47373,-0.300792,-0.185507,0.330423,0.151614,0.884043,-0.272951,-0.993822,0.511197,-0.0997948,-0.659756,0.575496,0.675617,0.180332,-0.593814,-0.492722,0.0646786,-0.666503,-0.945716,-0.334582,0.611894,0.281032,0.508749,-0.0265389,-0.920133,0.308861,-0.0437502,-0.374203,0.207359,-0.219433,0.914104,0.184408,0.520599,-0.88565,-0.375906,-0.708948,-0.37902,0.576578,0.0194674,-0.470262,0.572576,0.351245,-0.480477,0.238261,-0.1596,-0.827302,0.669248,0.494475,0.691715,-0.198585,0.0492812,0.959669,0.884086,-0.218632,-0.865161,-0.715997,0.220772,0.692356,0.646453,-0.401724,0.615443,-0.0601957,-0.748176,-0.207987,-0.865931,0.613732,-0.525712,-0.995728,0.389633,-0.064173,0.662131,-0.707048,-0.340423,0.60624,0.0951648,-0.160446,-0.394585,-0.167581,0.0679849,0.449799,0.733505,-0.00918638,0.00446808,0.404396,0.449996,0.162711,0.294454,-0.563345,-0.114993,0.549589,-0.775141,0.677726,0.610715,0.0850755,0.0419388,-0.323614,-0.973719,-0.680238,-0.270873,-0.209617,0.968436,0.908798,0.975851,-0.994918,-0.0621977,0.628171,0.761228,0.34372,-0.792042,-0.144765,-0.965748,0.0133606,-0.0260565,-0.742377,0.426391,0.408202,0.633885,-0.0351053,-0.723444,-0.577654,0.0276004,0.0712472,-0.622791,0.155451,0.442717,-0.792786,0.925785,0.670266,-0.865566,-0.638281,0.333094,0.477628,0.47261,0.23151,0.82132,-0.589803,0.796275,0.57713,0.101149,0.970191,0.532821,0.814769,-0.0687269,0.712758,-0.191812,-0.390938,0.952828,0.921519,0.923094,0.93011,-0.945394,-0.0934027,0.964123,-0.795609,-0.289563,0.614236,-0.670585,0.466877,0.144597,-0.206416,0.6937,-0.967958,-0.0951247,-0.942473,-0.610767,-0.655472,-0.0960986,-0.302779,-0.734976,-0.342188,-0.315861,-0.912834,0.24499,0.0969326,0.089775,-0.241157,0.0835558,-0.420236,-0.686633,-0.711276,-0.00325377,0.435196,-0.710002,0.00283691,-0.168757,-0.134045,-0.655235,0.172361,0.998291,0.376291,-0.962215,-0.363174,-0.88777,-0.519929,-0.560554,-0.984415,0.601529,-0.984103,-0.228237,-0.578066,0.307023,0.606123,0.959635,0.00225943,0.0101814,0.441456,0.0633629,0.406631,-0.0100638,-0.177972,-0.491075,0.537035,-0.924987,-0.699424,0.742285,0.0181443,0.718971,-0.0308272,0.086931,0.524476,0.956457,0.143024,0.616481,0.217909,-0.128427,-0.262427,-0.938208,-0.52479,0.12919,0.721925,0.766492,0.470845,-0.0976466,0.507807,0.804148,0.963269,0.357128,-0.832565,-0.312441,0.327779,0.184745,0.246139,-0.936814,-0.931734,-0.0327827,0.319293,0.044473,-0.641645,0.596118,-0.293934,-0.63373,0.409658,0.759892,-0.257078,0.939616,-0.227661,0.115754,0.10964,-0.240557,0.66842,0.855535,-0.451536,0.264961,-0.61366,-0.204783,-0.842476,-0.249524,-0.0985226,0.0671501,-0.527707,-0.509489,-0.883254,0.14851,-0.906465,0.496238,-0.853211,-0.779234,-0.979515,0.827175,0.228969,-0.402829,-0.970118,0.762559,0.506495,0.460303,0.897304,0.686003,0.739986,0.15731,0.281697,-0.922955,-0.780824,0.449716,0.125225,0.487649,0.147046,0.679639,0.593707,-0.311828,-0.797099,-0.35815,0.95808,0.907244,0.772426,0.720574,-0.873217,0.371431,-0.826029,0.942716,0.70609,-0.658158,-0.782185,-0.806743,-0.627986,-0.405551,-0.258495,-0.796524,0.222498,0.087545,-0.0917108,-0.62542,-0.110256,0.0417576,0.24476,0.941339,-0.613783,0.402772,0.300775,-0.820314,-0.894233,-0.405896,0.0735439,0.486645,-0.394355,0.125097,-0.316386,-0.701215,-0.845742,0.2065,-0.413743,0.406725,-0.423813,-0.941255,-0.558804,0.312326,0.345314,0.319143,-0.644653,-0.0408415,0.176461,0.740113,0.470737,-0.914927,-0.591523,-0.606614,-0.181873,0.692975,0.50208,-0.536704,0.359652,0.839082,0.56817,-0.0776788,-0.00332785,0.459538,-0.518313,-0.270738,-0.629958,-0.755084,-0.721573,0.431107,-0.221877,0.32543,0.163743,0.0759916,0.695064,-0.656856,0.074163,0.264319,-0.73174,0.731548,-0.489341,0.678946,0.0271269,0.804879,-0.402973,0.800373,0.760082,-0.878364,0.176801,-0.548932,-0.225601,-0.164912,-0.208143,0.7768,-0.542743,-0.156021,0.671736,0.878648,-0.419588,-0.0752896,0.0299447,-0.494459,-0.72415,0.35978,-0.32646,-0.96605,0.0127605,0.563174,-0.814853,-0.949609,-0.526794,-0.801902,-0.753397,0.617418,0.689874,0.983384,0.668786,0.0304653,-0.625221,-0.13318,0.827343,-0.101358,-0.999522,-0.0525574,-0.458319,0.587409,-0.334639,0.0759643,0.0255827,0.128944,0.17317,-0.284309,0.287161,-0.550725,-0.433083,-0.242821,0.878879,0.691699,-0.660499,0.389985,0.599856,-0.711442,-0.798697,-0.244945,-0.942649,0.402856,-0.494672,0.439941,-0.88216,0.170196,0.650734,-0.0982391,-0.468732,0.342133,-0.838071,-0.832362,0.658177,-0.565361,0.149473,0.69331,-0.491848,0.74916,0.526025,-0.155339,0.0998096,0.468761,0.324649,0.128488,0.544144,-0.495222,0.965229,-0.79314,-0.545421,-0.500243,0.154371,0.170017,-0.259108,-0.868862,-0.50731,-0.848317,0.835712,0.616391,-0.442608,-0.158,0.313451,0.703748,-0.755984,-0.249443,0.491564,0.985068,0.678644,0.808324,0.81975,-0.435823,-0.839855,0.00282368,-0.569165,0.0884339,-0.222144,0.499412,-0.565198,0.64824,0.956914,-0.0620912,0.634479,0.928617,0.464664,0.377022,0.63047,-0.198619,-0.576153,0.565373,-0.524245,-0.187299,-0.614524,0.429316,-0.491171,0.399495,-0.333898,-0.646636,-0.0189709,-0.339605,-0.798791,0.0494081,0.367012,0.852545,0.43557,0.150039,-0.0454542,0.604861,-0.598288,-0.500696,0.249008,0.370711,-0.633174,-0.0121906,0.42006,0.169373,-0.975542,-0.0297852,0.80481,0.638317,-0.670967,0.935792,-0.35605,0.175773,0.878601,-0.275168,-0.932517,-0.372497,-0.0732907,-0.185493,-0.357004,0.314786,-0.229239,0.530256,-0.51327,0.44187,0.940309,-0.240334,-0.0276121,0.74383,-0.630329,-0.763376,0.62538,0.818945,0.891598,0.680494,0.471868,-0.769787,-0.878099,-0.973724,0.354362,-0.1792,-0.225034,-0.44548,0.598865,0.544005,-0.478962,0.327193,-0.525784,0.903179,-0.899248,0.156514,0.154329,0.499808,-0.836327,-0.802627,0.378082,-0.112673,-0.47926,-0.3355,-0.699445,0.237731,-0.324597,-0.800406,-0.42585,-0.710739,-0.144068,-0.828545,-0.800912,0.184654,-0.63675,-0.16696,0.240427,-0.513443,0.812664,0.744943,0.970612,0.00172899,-0.726378,-0.0985012,0.224232,0.16495,0.560077,-0.813112,0.112894,-0.0955366,0.0187107,0.913887,0.123076,0.550338,0.400334,-0.367816,0.198455,-0.983183,0.278976,0.714817,0.307911,0.812861,-0.403497,-0.784382,-0.161823,-0.120835,0.323172,0.583739,0.732924,-0.220603,-0.594121,0.935093,-0.216736,0.659318,-0.750417,-0.284773,-0.271496,0.491731,-0.712174,-0.763681,0.0781023,0.951666,0.734031,0.826912,0.57488,-0.361951,-0.0739728,0.91438,-0.391653,0.0193887,0.412634,-0.169813,0.471794,0.660792,-0.350906,-0.612644,0.347876,0.112573,-0.501126,0.456761,-0.109004,0.289352,-0.566504,0.585042,0.584934,0.923676,0.895312,-0.161036,-0.995895,0.0853141,-0.583368,-0.157612,0.234119,0.875043,0.430805,0.706102,0.423887,0.296828,-0.265607,0.338806,-0.15829,0.642516,0.355126,0.174447,-0.975015,0.869905,-0.145285,-0.484002,-0.475966,-0.67704,0.996452,-0.0685748,-0.851985,0.416498,0.791047,-0.211323,-0.302819,0.640735,-0.317908,-0.116586,-0.896382,-0.817317,-0.948837,-0.597427,0.975863,-0.971371,-0.124115,0.4339,-0.254671,0.298068,-0.349803,-0.73185,0.488106,-0.0495073,0.253969,0.168116,0.148772,0.889593,-0.512213,-0.165437,0.666773,-0.976304,-0.170024,0.905794,0.473908,-0.855725,-0.0413591,-0.508661,0.443453,0.842925,-0.144503,0.936699,-0.443935,-0.182996,0.803564,0.960386,-0.0323329,0.638181,-0.895684,-0.360502,0.0646258,-0.202449,-0.717228,0.970489,0.404608,-0.0861868,-0.879417,-0.866462,-0.938336,-0.799685,0.213464,-0.932344,-0.668236,0.751366,-0.22712,-0.407783,0.657463,0.0970092,-0.579109,-0.868866,-0.504041,0.926483,0.169978,-0.00842563,-0.530324,0.282745,0.0255867,0.287686,0.410417,-0.766576,-0.536566,-0.628288,0.69665,0.820713,-0.506162,-0.404114,0.640099,-0.956123,-0.576586,0.435502,-0.470676,-0.367062,-0.831765,-0.294942,0.518991,0.922338,0.337886,-0.67474,-0.725667,0.916684,0.39175,0.759081,0.496979,-0.200691,0.0417966,-0.687391,0.438773,0.287357,0.316636,-0.262311,-0.0755541,-0.442313,0.621378,0.670105,0.060982,0.944162,0.643442,-0.750684,-0.639973,0.217424,0.592823,0.929094,-0.239135,-0.41628,0.570893,-0.0798988,-0.917135,-0.749545,-0.982047,0.0626998,-0.977963,0.660401,0.470569,-0.0528868,-0.00138645,0.931065,-0.748519,0.304188,-0.266153,0.672524,-0.105179,-0.874749,-0.154355,-0.774656,-0.69654,0.433098,0.615897,-0.387919,-0.429779,0.650202,0.122306,-0.237727,0.626817,-0.227929,0.405916,0.483328,0.282047,-0.262206,0.784123,0.83125,-0.662272,0.702768,0.875814,-0.701221,0.553793,0.471795,0.769147,0.059668,-0.841617,-0.191179,-0.972471,-0.825361,0.779826,-0.917201,0.43272,0.10301,0.358771,0.793448,-0.0379954,-0.870112,0.600442,-0.990603,0.549151,0.512146,-0.795843,0.490091,0.372046,-0.549437,0.0964285,0.753047,-0.86284,-0.589688,0.178612,-0.720358,};

  int i, j;
  /*
  FILE *myFile;
  double **M = (double **)malloc(sizeof(double*)*m);
  myFile = fopen("./data/L3ZDT1.txt", "r");

  for(i = 0; i < m; i++){
    M[i] = (double *)malloc(sizeof(double)*m);
  }

  for (i=0; i<m; i++)
  {
    for (j=0; j<m; j++)
    {
      if(!	fscanf(myFile, "%lf\n", &M[i][j])){
	       printf("Burp\n");
         break;
      }
    }
  }
  fclose(myFile);
*/
  for(i = 0; i < m; i++){
    if(x[i] >= 0 && x[i] <= 1)
      continue;
    else{
      printf("Invalid input range\n");
      return;
    }
  }

  y = (double*)malloc(sizeof(double)*m);
  for(i = 0; i < m; i++){
    y[i] = 0;
    for(j = 0; j < m; j++){
      y[i] = y[i] + ( M[i*m+j] * (x[j]*x[j]) );
    }
  }

  ff1 = y[0]*y[0];
  gx = 0;

  for(i = 1; i < m; i++){
    gx = gx + y[i]*y[i];
  }
  gx = 1+(9.0/(m-1))*gx;
  h = 1 - sqrt(ff1/gx);
  f[0] = ff1;
  f[1] = gx*h;

  return;
}

/*
 *
 *  As described by K. Deb and A. Sinha and S. Kukkonen in "Multi-objective
 *  test problems, linkages, and evolutionary methodologies", GECCO'06}:
 *  Proceedings of the 8th Annual Conference on Genetic and Evolutionary
 *  Computation, 1141-1148, 2006.
 *
 *  Example T2, with linkage L3.
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void L3ZDT2(double* f, double* x){
  int m = 30;
  double ff1, gx, h;
  double *y;

  double M[900] = {0.218418,-0.620254,0.843784,0.914311,-0.788548,0.428212,0.103064,-0.47373,-0.300792,-0.185507,0.330423,0.151614,0.884043,-0.272951,-0.993822,0.511197,-0.0997948,-0.659756,0.575496,0.675617,0.180332,-0.593814,-0.492722,0.0646786,-0.666503,-0.945716,-0.334582,0.611894,0.281032,0.508749,-0.0265389,-0.920133,0.308861,-0.0437502,-0.374203,0.207359,-0.219433,0.914104,0.184408,0.520599,-0.88565,-0.375906,-0.708948,-0.37902,0.576578,0.0194674,-0.470262,0.572576,0.351245,-0.480477,0.238261,-0.1596,-0.827302,0.669248,0.494475,0.691715,-0.198585,0.0492812,0.959669,0.884086,-0.218632,-0.865161,-0.715997,0.220772,0.692356,0.646453,-0.401724,0.615443,-0.0601957,-0.748176,-0.207987,-0.865931,0.613732,-0.525712,-0.995728,0.389633,-0.064173,0.662131,-0.707048,-0.340423,0.60624,0.0951648,-0.160446,-0.394585,-0.167581,0.0679849,0.449799,0.733505,-0.00918638,0.00446808,0.404396,0.449996,0.162711,0.294454,-0.563345,-0.114993,0.549589,-0.775141,0.677726,0.610715,0.0850755,0.0419388,-0.323614,-0.973719,-0.680238,-0.270873,-0.209617,0.968436,0.908798,0.975851,-0.994918,-0.0621977,0.628171,0.761228,0.34372,-0.792042,-0.144765,-0.965748,0.0133606,-0.0260565,-0.742377,0.426391,0.408202,0.633885,-0.0351053,-0.723444,-0.577654,0.0276004,0.0712472,-0.622791,0.155451,0.442717,-0.792786,0.925785,0.670266,-0.865566,-0.638281,0.333094,0.477628,0.47261,0.23151,0.82132,-0.589803,0.796275,0.57713,0.101149,0.970191,0.532821,0.814769,-0.0687269,0.712758,-0.191812,-0.390938,0.952828,0.921519,0.923094,0.93011,-0.945394,-0.0934027,0.964123,-0.795609,-0.289563,0.614236,-0.670585,0.466877,0.144597,-0.206416,0.6937,-0.967958,-0.0951247,-0.942473,-0.610767,-0.655472,-0.0960986,-0.302779,-0.734976,-0.342188,-0.315861,-0.912834,0.24499,0.0969326,0.089775,-0.241157,0.0835558,-0.420236,-0.686633,-0.711276,-0.00325377,0.435196,-0.710002,0.00283691,-0.168757,-0.134045,-0.655235,0.172361,0.998291,0.376291,-0.962215,-0.363174,-0.88777,-0.519929,-0.560554,-0.984415,0.601529,-0.984103,-0.228237,-0.578066,0.307023,0.606123,0.959635,0.00225943,0.0101814,0.441456,0.0633629,0.406631,-0.0100638,-0.177972,-0.491075,0.537035,-0.924987,-0.699424,0.742285,0.0181443,0.718971,-0.0308272,0.086931,0.524476,0.956457,0.143024,0.616481,0.217909,-0.128427,-0.262427,-0.938208,-0.52479,0.12919,0.721925,0.766492,0.470845,-0.0976466,0.507807,0.804148,0.963269,0.357128,-0.832565,-0.312441,0.327779,0.184745,0.246139,-0.936814,-0.931734,-0.0327827,0.319293,0.044473,-0.641645,0.596118,-0.293934,-0.63373,0.409658,0.759892,-0.257078,0.939616,-0.227661,0.115754,0.10964,-0.240557,0.66842,0.855535,-0.451536,0.264961,-0.61366,-0.204783,-0.842476,-0.249524,-0.0985226,0.0671501,-0.527707,-0.509489,-0.883254,0.14851,-0.906465,0.496238,-0.853211,-0.779234,-0.979515,0.827175,0.228969,-0.402829,-0.970118,0.762559,0.506495,0.460303,0.897304,0.686003,0.739986,0.15731,0.281697,-0.922955,-0.780824,0.449716,0.125225,0.487649,0.147046,0.679639,0.593707,-0.311828,-0.797099,-0.35815,0.95808,0.907244,0.772426,0.720574,-0.873217,0.371431,-0.826029,0.942716,0.70609,-0.658158,-0.782185,-0.806743,-0.627986,-0.405551,-0.258495,-0.796524,0.222498,0.087545,-0.0917108,-0.62542,-0.110256,0.0417576,0.24476,0.941339,-0.613783,0.402772,0.300775,-0.820314,-0.894233,-0.405896,0.0735439,0.486645,-0.394355,0.125097,-0.316386,-0.701215,-0.845742,0.2065,-0.413743,0.406725,-0.423813,-0.941255,-0.558804,0.312326,0.345314,0.319143,-0.644653,-0.0408415,0.176461,0.740113,0.470737,-0.914927,-0.591523,-0.606614,-0.181873,0.692975,0.50208,-0.536704,0.359652,0.839082,0.56817,-0.0776788,-0.00332785,0.459538,-0.518313,-0.270738,-0.629958,-0.755084,-0.721573,0.431107,-0.221877,0.32543,0.163743,0.0759916,0.695064,-0.656856,0.074163,0.264319,-0.73174,0.731548,-0.489341,0.678946,0.0271269,0.804879,-0.402973,0.800373,0.760082,-0.878364,0.176801,-0.548932,-0.225601,-0.164912,-0.208143,0.7768,-0.542743,-0.156021,0.671736,0.878648,-0.419588,-0.0752896,0.0299447,-0.494459,-0.72415,0.35978,-0.32646,-0.96605,0.0127605,0.563174,-0.814853,-0.949609,-0.526794,-0.801902,-0.753397,0.617418,0.689874,0.983384,0.668786,0.0304653,-0.625221,-0.13318,0.827343,-0.101358,-0.999522,-0.0525574,-0.458319,0.587409,-0.334639,0.0759643,0.0255827,0.128944,0.17317,-0.284309,0.287161,-0.550725,-0.433083,-0.242821,0.878879,0.691699,-0.660499,0.389985,0.599856,-0.711442,-0.798697,-0.244945,-0.942649,0.402856,-0.494672,0.439941,-0.88216,0.170196,0.650734,-0.0982391,-0.468732,0.342133,-0.838071,-0.832362,0.658177,-0.565361,0.149473,0.69331,-0.491848,0.74916,0.526025,-0.155339,0.0998096,0.468761,0.324649,0.128488,0.544144,-0.495222,0.965229,-0.79314,-0.545421,-0.500243,0.154371,0.170017,-0.259108,-0.868862,-0.50731,-0.848317,0.835712,0.616391,-0.442608,-0.158,0.313451,0.703748,-0.755984,-0.249443,0.491564,0.985068,0.678644,0.808324,0.81975,-0.435823,-0.839855,0.00282368,-0.569165,0.0884339,-0.222144,0.499412,-0.565198,0.64824,0.956914,-0.0620912,0.634479,0.928617,0.464664,0.377022,0.63047,-0.198619,-0.576153,0.565373,-0.524245,-0.187299,-0.614524,0.429316,-0.491171,0.399495,-0.333898,-0.646636,-0.0189709,-0.339605,-0.798791,0.0494081,0.367012,0.852545,0.43557,0.150039,-0.0454542,0.604861,-0.598288,-0.500696,0.249008,0.370711,-0.633174,-0.0121906,0.42006,0.169373,-0.975542,-0.0297852,0.80481,0.638317,-0.670967,0.935792,-0.35605,0.175773,0.878601,-0.275168,-0.932517,-0.372497,-0.0732907,-0.185493,-0.357004,0.314786,-0.229239,0.530256,-0.51327,0.44187,0.940309,-0.240334,-0.0276121,0.74383,-0.630329,-0.763376,0.62538,0.818945,0.891598,0.680494,0.471868,-0.769787,-0.878099,-0.973724,0.354362,-0.1792,-0.225034,-0.44548,0.598865,0.544005,-0.478962,0.327193,-0.525784,0.903179,-0.899248,0.156514,0.154329,0.499808,-0.836327,-0.802627,0.378082,-0.112673,-0.47926,-0.3355,-0.699445,0.237731,-0.324597,-0.800406,-0.42585,-0.710739,-0.144068,-0.828545,-0.800912,0.184654,-0.63675,-0.16696,0.240427,-0.513443,0.812664,0.744943,0.970612,0.00172899,-0.726378,-0.0985012,0.224232,0.16495,0.560077,-0.813112,0.112894,-0.0955366,0.0187107,0.913887,0.123076,0.550338,0.400334,-0.367816,0.198455,-0.983183,0.278976,0.714817,0.307911,0.812861,-0.403497,-0.784382,-0.161823,-0.120835,0.323172,0.583739,0.732924,-0.220603,-0.594121,0.935093,-0.216736,0.659318,-0.750417,-0.284773,-0.271496,0.491731,-0.712174,-0.763681,0.0781023,0.951666,0.734031,0.826912,0.57488,-0.361951,-0.0739728,0.91438,-0.391653,0.0193887,0.412634,-0.169813,0.471794,0.660792,-0.350906,-0.612644,0.347876,0.112573,-0.501126,0.456761,-0.109004,0.289352,-0.566504,0.585042,0.584934,0.923676,0.895312,-0.161036,-0.995895,0.0853141,-0.583368,-0.157612,0.234119,0.875043,0.430805,0.706102,0.423887,0.296828,-0.265607,0.338806,-0.15829,0.642516,0.355126,0.174447,-0.975015,0.869905,-0.145285,-0.484002,-0.475966,-0.67704,0.996452,-0.0685748,-0.851985,0.416498,0.791047,-0.211323,-0.302819,0.640735,-0.317908,-0.116586,-0.896382,-0.817317,-0.948837,-0.597427,0.975863,-0.971371,-0.124115,0.4339,-0.254671,0.298068,-0.349803,-0.73185,0.488106,-0.0495073,0.253969,0.168116,0.148772,0.889593,-0.512213,-0.165437,0.666773,-0.976304,-0.170024,0.905794,0.473908,-0.855725,-0.0413591,-0.508661,0.443453,0.842925,-0.144503,0.936699,-0.443935,-0.182996,0.803564,0.960386,-0.0323329,0.638181,-0.895684,-0.360502,0.0646258,-0.202449,-0.717228,0.970489,0.404608,-0.0861868,-0.879417,-0.866462,-0.938336,-0.799685,0.213464,-0.932344,-0.668236,0.751366,-0.22712,-0.407783,0.657463,0.0970092,-0.579109,-0.868866,-0.504041,0.926483,0.169978,-0.00842563,-0.530324,0.282745,0.0255867,0.287686,0.410417,-0.766576,-0.536566,-0.628288,0.69665,0.820713,-0.506162,-0.404114,0.640099,-0.956123,-0.576586,0.435502,-0.470676,-0.367062,-0.831765,-0.294942,0.518991,0.922338,0.337886,-0.67474,-0.725667,0.916684,0.39175,0.759081,0.496979,-0.200691,0.0417966,-0.687391,0.438773,0.287357,0.316636,-0.262311,-0.0755541,-0.442313,0.621378,0.670105,0.060982,0.944162,0.643442,-0.750684,-0.639973,0.217424,0.592823,0.929094,-0.239135,-0.41628,0.570893,-0.0798988,-0.917135,-0.749545,-0.982047,0.0626998,-0.977963,0.660401,0.470569,-0.0528868,-0.00138645,0.931065,-0.748519,0.304188,-0.266153,0.672524,-0.105179,-0.874749,-0.154355,-0.774656,-0.69654,0.433098,0.615897,-0.387919,-0.429779,0.650202,0.122306,-0.237727,0.626817,-0.227929,0.405916,0.483328,0.282047,-0.262206,0.784123,0.83125,-0.662272,0.702768,0.875814,-0.701221,0.553793,0.471795,0.769147,0.059668,-0.841617,-0.191179,-0.972471,-0.825361,0.779826,-0.917201,0.43272,0.10301,0.358771,0.793448,-0.0379954,-0.870112,0.600442,-0.990603,0.549151,0.512146,-0.795843,0.490091,0.372046,-0.549437,0.0964285,0.753047,-0.86284,-0.589688,0.178612,-0.720358,};
  int i, j;
  // FILE *myFile;
  // for(i = 0; i < m; i++){
  //   M[i] = (double *)malloc(sizeof(double)*m);
  // }
  // myFile = fopen("./data/L3ZDT1.txt", "r");

  // for (i=0; i<m; i++)
  // {
  //   for (j=0; j<m; j++)
  //   {
  //     if(!	fscanf(myFile, "%lf\n", &M[i][j]))
   //     break;
  //   }
  // }
  // fclose(myFile);

  for(i = 0; i < m; i++){
    if(x[i] >= 0 && x[i] <= 1)
      continue;
    else{
      printf("Invalid input range\n");
      return;
    }
  }

  y = (double*)malloc(sizeof(double)*m);
  for(i = 0; i < m; i++){
    y[i] = 0;
    for(j = 0; j < m; j++){
      y[i] = y[i] + ( M[i*m+j] * (x[j]*x[j]) );
    }
  }

  ff1 = y[0]*y[0];
  gx = 0;

  for(i = 1; i < m; i++){
    gx = gx + y[i]*y[i];
  }
  gx = gx * (9.0/(m-1));
  gx = gx + 1;

  h = 1 - ( (ff1/gx)*(ff1/gx) );

  f[0] = ff1;
  f[1] = gx*h;

  free(y);

  return;
}

/*
 *
 *  As described by K. Deb and A. Sinha and S. Kukkonen in "Multi-objective
 *  test problems, linkages, and evolutionary methodologies", GECCO'06}:
 *  Proceedings of the 8th Annual Conference on Genetic and Evolutionary
 *  Computation, 1141-1148, 2006.
 *
 *  Example T3, with linkage L3.
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void L3ZDT3(double* f, double* x){
  int m = 30;

  double ff1, gx, h;
  double *y;

  double pi = 4*atan(1);

  int i, j;
  double M[900] = {0.218418,-0.620254,0.843784,0.914311,-0.788548,0.428212,0.103064,-0.47373,-0.300792,-0.185507,0.330423,0.151614,0.884043,-0.272951,-0.993822,0.511197,-0.0997948,-0.659756,0.575496,0.675617,0.180332,-0.593814,-0.492722,0.0646786,-0.666503,-0.945716,-0.334582,0.611894,0.281032,0.508749,-0.0265389,-0.920133,0.308861,-0.0437502,-0.374203,0.207359,-0.219433,0.914104,0.184408,0.520599,-0.88565,-0.375906,-0.708948,-0.37902,0.576578,0.0194674,-0.470262,0.572576,0.351245,-0.480477,0.238261,-0.1596,-0.827302,0.669248,0.494475,0.691715,-0.198585,0.0492812,0.959669,0.884086,-0.218632,-0.865161,-0.715997,0.220772,0.692356,0.646453,-0.401724,0.615443,-0.0601957,-0.748176,-0.207987,-0.865931,0.613732,-0.525712,-0.995728,0.389633,-0.064173,0.662131,-0.707048,-0.340423,0.60624,0.0951648,-0.160446,-0.394585,-0.167581,0.0679849,0.449799,0.733505,-0.00918638,0.00446808,0.404396,0.449996,0.162711,0.294454,-0.563345,-0.114993,0.549589,-0.775141,0.677726,0.610715,0.0850755,0.0419388,-0.323614,-0.973719,-0.680238,-0.270873,-0.209617,0.968436,0.908798,0.975851,-0.994918,-0.0621977,0.628171,0.761228,0.34372,-0.792042,-0.144765,-0.965748,0.0133606,-0.0260565,-0.742377,0.426391,0.408202,0.633885,-0.0351053,-0.723444,-0.577654,0.0276004,0.0712472,-0.622791,0.155451,0.442717,-0.792786,0.925785,0.670266,-0.865566,-0.638281,0.333094,0.477628,0.47261,0.23151,0.82132,-0.589803,0.796275,0.57713,0.101149,0.970191,0.532821,0.814769,-0.0687269,0.712758,-0.191812,-0.390938,0.952828,0.921519,0.923094,0.93011,-0.945394,-0.0934027,0.964123,-0.795609,-0.289563,0.614236,-0.670585,0.466877,0.144597,-0.206416,0.6937,-0.967958,-0.0951247,-0.942473,-0.610767,-0.655472,-0.0960986,-0.302779,-0.734976,-0.342188,-0.315861,-0.912834,0.24499,0.0969326,0.089775,-0.241157,0.0835558,-0.420236,-0.686633,-0.711276,-0.00325377,0.435196,-0.710002,0.00283691,-0.168757,-0.134045,-0.655235,0.172361,0.998291,0.376291,-0.962215,-0.363174,-0.88777,-0.519929,-0.560554,-0.984415,0.601529,-0.984103,-0.228237,-0.578066,0.307023,0.606123,0.959635,0.00225943,0.0101814,0.441456,0.0633629,0.406631,-0.0100638,-0.177972,-0.491075,0.537035,-0.924987,-0.699424,0.742285,0.0181443,0.718971,-0.0308272,0.086931,0.524476,0.956457,0.143024,0.616481,0.217909,-0.128427,-0.262427,-0.938208,-0.52479,0.12919,0.721925,0.766492,0.470845,-0.0976466,0.507807,0.804148,0.963269,0.357128,-0.832565,-0.312441,0.327779,0.184745,0.246139,-0.936814,-0.931734,-0.0327827,0.319293,0.044473,-0.641645,0.596118,-0.293934,-0.63373,0.409658,0.759892,-0.257078,0.939616,-0.227661,0.115754,0.10964,-0.240557,0.66842,0.855535,-0.451536,0.264961,-0.61366,-0.204783,-0.842476,-0.249524,-0.0985226,0.0671501,-0.527707,-0.509489,-0.883254,0.14851,-0.906465,0.496238,-0.853211,-0.779234,-0.979515,0.827175,0.228969,-0.402829,-0.970118,0.762559,0.506495,0.460303,0.897304,0.686003,0.739986,0.15731,0.281697,-0.922955,-0.780824,0.449716,0.125225,0.487649,0.147046,0.679639,0.593707,-0.311828,-0.797099,-0.35815,0.95808,0.907244,0.772426,0.720574,-0.873217,0.371431,-0.826029,0.942716,0.70609,-0.658158,-0.782185,-0.806743,-0.627986,-0.405551,-0.258495,-0.796524,0.222498,0.087545,-0.0917108,-0.62542,-0.110256,0.0417576,0.24476,0.941339,-0.613783,0.402772,0.300775,-0.820314,-0.894233,-0.405896,0.0735439,0.486645,-0.394355,0.125097,-0.316386,-0.701215,-0.845742,0.2065,-0.413743,0.406725,-0.423813,-0.941255,-0.558804,0.312326,0.345314,0.319143,-0.644653,-0.0408415,0.176461,0.740113,0.470737,-0.914927,-0.591523,-0.606614,-0.181873,0.692975,0.50208,-0.536704,0.359652,0.839082,0.56817,-0.0776788,-0.00332785,0.459538,-0.518313,-0.270738,-0.629958,-0.755084,-0.721573,0.431107,-0.221877,0.32543,0.163743,0.0759916,0.695064,-0.656856,0.074163,0.264319,-0.73174,0.731548,-0.489341,0.678946,0.0271269,0.804879,-0.402973,0.800373,0.760082,-0.878364,0.176801,-0.548932,-0.225601,-0.164912,-0.208143,0.7768,-0.542743,-0.156021,0.671736,0.878648,-0.419588,-0.0752896,0.0299447,-0.494459,-0.72415,0.35978,-0.32646,-0.96605,0.0127605,0.563174,-0.814853,-0.949609,-0.526794,-0.801902,-0.753397,0.617418,0.689874,0.983384,0.668786,0.0304653,-0.625221,-0.13318,0.827343,-0.101358,-0.999522,-0.0525574,-0.458319,0.587409,-0.334639,0.0759643,0.0255827,0.128944,0.17317,-0.284309,0.287161,-0.550725,-0.433083,-0.242821,0.878879,0.691699,-0.660499,0.389985,0.599856,-0.711442,-0.798697,-0.244945,-0.942649,0.402856,-0.494672,0.439941,-0.88216,0.170196,0.650734,-0.0982391,-0.468732,0.342133,-0.838071,-0.832362,0.658177,-0.565361,0.149473,0.69331,-0.491848,0.74916,0.526025,-0.155339,0.0998096,0.468761,0.324649,0.128488,0.544144,-0.495222,0.965229,-0.79314,-0.545421,-0.500243,0.154371,0.170017,-0.259108,-0.868862,-0.50731,-0.848317,0.835712,0.616391,-0.442608,-0.158,0.313451,0.703748,-0.755984,-0.249443,0.491564,0.985068,0.678644,0.808324,0.81975,-0.435823,-0.839855,0.00282368,-0.569165,0.0884339,-0.222144,0.499412,-0.565198,0.64824,0.956914,-0.0620912,0.634479,0.928617,0.464664,0.377022,0.63047,-0.198619,-0.576153,0.565373,-0.524245,-0.187299,-0.614524,0.429316,-0.491171,0.399495,-0.333898,-0.646636,-0.0189709,-0.339605,-0.798791,0.0494081,0.367012,0.852545,0.43557,0.150039,-0.0454542,0.604861,-0.598288,-0.500696,0.249008,0.370711,-0.633174,-0.0121906,0.42006,0.169373,-0.975542,-0.0297852,0.80481,0.638317,-0.670967,0.935792,-0.35605,0.175773,0.878601,-0.275168,-0.932517,-0.372497,-0.0732907,-0.185493,-0.357004,0.314786,-0.229239,0.530256,-0.51327,0.44187,0.940309,-0.240334,-0.0276121,0.74383,-0.630329,-0.763376,0.62538,0.818945,0.891598,0.680494,0.471868,-0.769787,-0.878099,-0.973724,0.354362,-0.1792,-0.225034,-0.44548,0.598865,0.544005,-0.478962,0.327193,-0.525784,0.903179,-0.899248,0.156514,0.154329,0.499808,-0.836327,-0.802627,0.378082,-0.112673,-0.47926,-0.3355,-0.699445,0.237731,-0.324597,-0.800406,-0.42585,-0.710739,-0.144068,-0.828545,-0.800912,0.184654,-0.63675,-0.16696,0.240427,-0.513443,0.812664,0.744943,0.970612,0.00172899,-0.726378,-0.0985012,0.224232,0.16495,0.560077,-0.813112,0.112894,-0.0955366,0.0187107,0.913887,0.123076,0.550338,0.400334,-0.367816,0.198455,-0.983183,0.278976,0.714817,0.307911,0.812861,-0.403497,-0.784382,-0.161823,-0.120835,0.323172,0.583739,0.732924,-0.220603,-0.594121,0.935093,-0.216736,0.659318,-0.750417,-0.284773,-0.271496,0.491731,-0.712174,-0.763681,0.0781023,0.951666,0.734031,0.826912,0.57488,-0.361951,-0.0739728,0.91438,-0.391653,0.0193887,0.412634,-0.169813,0.471794,0.660792,-0.350906,-0.612644,0.347876,0.112573,-0.501126,0.456761,-0.109004,0.289352,-0.566504,0.585042,0.584934,0.923676,0.895312,-0.161036,-0.995895,0.0853141,-0.583368,-0.157612,0.234119,0.875043,0.430805,0.706102,0.423887,0.296828,-0.265607,0.338806,-0.15829,0.642516,0.355126,0.174447,-0.975015,0.869905,-0.145285,-0.484002,-0.475966,-0.67704,0.996452,-0.0685748,-0.851985,0.416498,0.791047,-0.211323,-0.302819,0.640735,-0.317908,-0.116586,-0.896382,-0.817317,-0.948837,-0.597427,0.975863,-0.971371,-0.124115,0.4339,-0.254671,0.298068,-0.349803,-0.73185,0.488106,-0.0495073,0.253969,0.168116,0.148772,0.889593,-0.512213,-0.165437,0.666773,-0.976304,-0.170024,0.905794,0.473908,-0.855725,-0.0413591,-0.508661,0.443453,0.842925,-0.144503,0.936699,-0.443935,-0.182996,0.803564,0.960386,-0.0323329,0.638181,-0.895684,-0.360502,0.0646258,-0.202449,-0.717228,0.970489,0.404608,-0.0861868,-0.879417,-0.866462,-0.938336,-0.799685,0.213464,-0.932344,-0.668236,0.751366,-0.22712,-0.407783,0.657463,0.0970092,-0.579109,-0.868866,-0.504041,0.926483,0.169978,-0.00842563,-0.530324,0.282745,0.0255867,0.287686,0.410417,-0.766576,-0.536566,-0.628288,0.69665,0.820713,-0.506162,-0.404114,0.640099,-0.956123,-0.576586,0.435502,-0.470676,-0.367062,-0.831765,-0.294942,0.518991,0.922338,0.337886,-0.67474,-0.725667,0.916684,0.39175,0.759081,0.496979,-0.200691,0.0417966,-0.687391,0.438773,0.287357,0.316636,-0.262311,-0.0755541,-0.442313,0.621378,0.670105,0.060982,0.944162,0.643442,-0.750684,-0.639973,0.217424,0.592823,0.929094,-0.239135,-0.41628,0.570893,-0.0798988,-0.917135,-0.749545,-0.982047,0.0626998,-0.977963,0.660401,0.470569,-0.0528868,-0.00138645,0.931065,-0.748519,0.304188,-0.266153,0.672524,-0.105179,-0.874749,-0.154355,-0.774656,-0.69654,0.433098,0.615897,-0.387919,-0.429779,0.650202,0.122306,-0.237727,0.626817,-0.227929,0.405916,0.483328,0.282047,-0.262206,0.784123,0.83125,-0.662272,0.702768,0.875814,-0.701221,0.553793,0.471795,0.769147,0.059668,-0.841617,-0.191179,-0.972471,-0.825361,0.779826,-0.917201,0.43272,0.10301,0.358771,0.793448,-0.0379954,-0.870112,0.600442,-0.990603,0.549151,0.512146,-0.795843,0.490091,0.372046,-0.549437,0.0964285,0.753047,-0.86284,-0.589688,0.178612,-0.720358,};
  // FILE *myFile;
  // for(i = 0; i < m; i++){
  //   M[i] = (double *)malloc(sizeof(double)*m);
  // }
  // myFile = fopen("./data/L3ZDT1.txt", "r");

  // for (i=0; i<m; i++)
  // {
  //   for (j=0; j<m; j++)
  //   {
  //     if(!	fscanf(myFile, "%lf\n", &M[i][j]))
	 //     break;
  //   }
  // }
  // fclose(myFile);

  for(i = 0; i < m; i++){
    if(x[i] >= 0 && x[i] <= 1)
      continue;
    else{
      printf("Invalid input range\n");
      return;
    }
  }

  y = (double*)malloc(sizeof(double)*m);
  for(i = 0; i < m; i++){
    y[i] = 0;
    for(j = 0; j < m; j++){
      y[i] = y[i] + ( M[i*m+j] * (x[j]*x[j]) );
    }
  }

  ff1 = y[0]*y[0];
  gx = 0;

  for(i = 1; i < m; i++){
    gx = gx + y[i]*y[i];
  }
  gx = gx * (9.0/(m-1));
  gx = gx + 1;

  h = 1 - sqrt(ff1/gx)-(ff1/gx)*sin(10*pi*ff1);

  f[0] = ff1;
  f[1] = gx*h;

  free(y);

  return;
}

/*
 *
 *  As described by K. Deb and A. Sinha and S. Kukkonen in "Multi-objective
 *  test problems, linkages, and evolutionary methodologies", GECCO'06}:
 *  Proceedings of the 8th Annual Conference on Genetic and Evolutionary
 *  Computation, 1141-1148, 2006.
 *
 *  Example T4, with linkage L3.
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void L3ZDT4(double* f, double* x){
  int m = 30;
  double ff1, gx, h;
  double *y;

  double pi = 4*atan(1);

  double M[900] = {0.218418,-0.620254,0.843784,0.914311,-0.788548,0.428212,0.103064,-0.47373,-0.300792,-0.185507,0.330423,0.151614,0.884043,-0.272951,-0.993822,0.511197,-0.0997948,-0.659756,0.575496,0.675617,0.180332,-0.593814,-0.492722,0.0646786,-0.666503,-0.945716,-0.334582,0.611894,0.281032,0.508749,-0.0265389,-0.920133,0.308861,-0.0437502,-0.374203,0.207359,-0.219433,0.914104,0.184408,0.520599,-0.88565,-0.375906,-0.708948,-0.37902,0.576578,0.0194674,-0.470262,0.572576,0.351245,-0.480477,0.238261,-0.1596,-0.827302,0.669248,0.494475,0.691715,-0.198585,0.0492812,0.959669,0.884086,-0.218632,-0.865161,-0.715997,0.220772,0.692356,0.646453,-0.401724,0.615443,-0.0601957,-0.748176,-0.207987,-0.865931,0.613732,-0.525712,-0.995728,0.389633,-0.064173,0.662131,-0.707048,-0.340423,0.60624,0.0951648,-0.160446,-0.394585,-0.167581,0.0679849,0.449799,0.733505,-0.00918638,0.00446808,0.404396,0.449996,0.162711,0.294454,-0.563345,-0.114993,0.549589,-0.775141,0.677726,0.610715,0.0850755,0.0419388,-0.323614,-0.973719,-0.680238,-0.270873,-0.209617,0.968436,0.908798,0.975851,-0.994918,-0.0621977,0.628171,0.761228,0.34372,-0.792042,-0.144765,-0.965748,0.0133606,-0.0260565,-0.742377,0.426391,0.408202,0.633885,-0.0351053,-0.723444,-0.577654,0.0276004,0.0712472,-0.622791,0.155451,0.442717,-0.792786,0.925785,0.670266,-0.865566,-0.638281,0.333094,0.477628,0.47261,0.23151,0.82132,-0.589803,0.796275,0.57713,0.101149,0.970191,0.532821,0.814769,-0.0687269,0.712758,-0.191812,-0.390938,0.952828,0.921519,0.923094,0.93011,-0.945394,-0.0934027,0.964123,-0.795609,-0.289563,0.614236,-0.670585,0.466877,0.144597,-0.206416,0.6937,-0.967958,-0.0951247,-0.942473,-0.610767,-0.655472,-0.0960986,-0.302779,-0.734976,-0.342188,-0.315861,-0.912834,0.24499,0.0969326,0.089775,-0.241157,0.0835558,-0.420236,-0.686633,-0.711276,-0.00325377,0.435196,-0.710002,0.00283691,-0.168757,-0.134045,-0.655235,0.172361,0.998291,0.376291,-0.962215,-0.363174,-0.88777,-0.519929,-0.560554,-0.984415,0.601529,-0.984103,-0.228237,-0.578066,0.307023,0.606123,0.959635,0.00225943,0.0101814,0.441456,0.0633629,0.406631,-0.0100638,-0.177972,-0.491075,0.537035,-0.924987,-0.699424,0.742285,0.0181443,0.718971,-0.0308272,0.086931,0.524476,0.956457,0.143024,0.616481,0.217909,-0.128427,-0.262427,-0.938208,-0.52479,0.12919,0.721925,0.766492,0.470845,-0.0976466,0.507807,0.804148,0.963269,0.357128,-0.832565,-0.312441,0.327779,0.184745,0.246139,-0.936814,-0.931734,-0.0327827,0.319293,0.044473,-0.641645,0.596118,-0.293934,-0.63373,0.409658,0.759892,-0.257078,0.939616,-0.227661,0.115754,0.10964,-0.240557,0.66842,0.855535,-0.451536,0.264961,-0.61366,-0.204783,-0.842476,-0.249524,-0.0985226,0.0671501,-0.527707,-0.509489,-0.883254,0.14851,-0.906465,0.496238,-0.853211,-0.779234,-0.979515,0.827175,0.228969,-0.402829,-0.970118,0.762559,0.506495,0.460303,0.897304,0.686003,0.739986,0.15731,0.281697,-0.922955,-0.780824,0.449716,0.125225,0.487649,0.147046,0.679639,0.593707,-0.311828,-0.797099,-0.35815,0.95808,0.907244,0.772426,0.720574,-0.873217,0.371431,-0.826029,0.942716,0.70609,-0.658158,-0.782185,-0.806743,-0.627986,-0.405551,-0.258495,-0.796524,0.222498,0.087545,-0.0917108,-0.62542,-0.110256,0.0417576,0.24476,0.941339,-0.613783,0.402772,0.300775,-0.820314,-0.894233,-0.405896,0.0735439,0.486645,-0.394355,0.125097,-0.316386,-0.701215,-0.845742,0.2065,-0.413743,0.406725,-0.423813,-0.941255,-0.558804,0.312326,0.345314,0.319143,-0.644653,-0.0408415,0.176461,0.740113,0.470737,-0.914927,-0.591523,-0.606614,-0.181873,0.692975,0.50208,-0.536704,0.359652,0.839082,0.56817,-0.0776788,-0.00332785,0.459538,-0.518313,-0.270738,-0.629958,-0.755084,-0.721573,0.431107,-0.221877,0.32543,0.163743,0.0759916,0.695064,-0.656856,0.074163,0.264319,-0.73174,0.731548,-0.489341,0.678946,0.0271269,0.804879,-0.402973,0.800373,0.760082,-0.878364,0.176801,-0.548932,-0.225601,-0.164912,-0.208143,0.7768,-0.542743,-0.156021,0.671736,0.878648,-0.419588,-0.0752896,0.0299447,-0.494459,-0.72415,0.35978,-0.32646,-0.96605,0.0127605,0.563174,-0.814853,-0.949609,-0.526794,-0.801902,-0.753397,0.617418,0.689874,0.983384,0.668786,0.0304653,-0.625221,-0.13318,0.827343,-0.101358,-0.999522,-0.0525574,-0.458319,0.587409,-0.334639,0.0759643,0.0255827,0.128944,0.17317,-0.284309,0.287161,-0.550725,-0.433083,-0.242821,0.878879,0.691699,-0.660499,0.389985,0.599856,-0.711442,-0.798697,-0.244945,-0.942649,0.402856,-0.494672,0.439941,-0.88216,0.170196,0.650734,-0.0982391,-0.468732,0.342133,-0.838071,-0.832362,0.658177,-0.565361,0.149473,0.69331,-0.491848,0.74916,0.526025,-0.155339,0.0998096,0.468761,0.324649,0.128488,0.544144,-0.495222,0.965229,-0.79314,-0.545421,-0.500243,0.154371,0.170017,-0.259108,-0.868862,-0.50731,-0.848317,0.835712,0.616391,-0.442608,-0.158,0.313451,0.703748,-0.755984,-0.249443,0.491564,0.985068,0.678644,0.808324,0.81975,-0.435823,-0.839855,0.00282368,-0.569165,0.0884339,-0.222144,0.499412,-0.565198,0.64824,0.956914,-0.0620912,0.634479,0.928617,0.464664,0.377022,0.63047,-0.198619,-0.576153,0.565373,-0.524245,-0.187299,-0.614524,0.429316,-0.491171,0.399495,-0.333898,-0.646636,-0.0189709,-0.339605,-0.798791,0.0494081,0.367012,0.852545,0.43557,0.150039,-0.0454542,0.604861,-0.598288,-0.500696,0.249008,0.370711,-0.633174,-0.0121906,0.42006,0.169373,-0.975542,-0.0297852,0.80481,0.638317,-0.670967,0.935792,-0.35605,0.175773,0.878601,-0.275168,-0.932517,-0.372497,-0.0732907,-0.185493,-0.357004,0.314786,-0.229239,0.530256,-0.51327,0.44187,0.940309,-0.240334,-0.0276121,0.74383,-0.630329,-0.763376,0.62538,0.818945,0.891598,0.680494,0.471868,-0.769787,-0.878099,-0.973724,0.354362,-0.1792,-0.225034,-0.44548,0.598865,0.544005,-0.478962,0.327193,-0.525784,0.903179,-0.899248,0.156514,0.154329,0.499808,-0.836327,-0.802627,0.378082,-0.112673,-0.47926,-0.3355,-0.699445,0.237731,-0.324597,-0.800406,-0.42585,-0.710739,-0.144068,-0.828545,-0.800912,0.184654,-0.63675,-0.16696,0.240427,-0.513443,0.812664,0.744943,0.970612,0.00172899,-0.726378,-0.0985012,0.224232,0.16495,0.560077,-0.813112,0.112894,-0.0955366,0.0187107,0.913887,0.123076,0.550338,0.400334,-0.367816,0.198455,-0.983183,0.278976,0.714817,0.307911,0.812861,-0.403497,-0.784382,-0.161823,-0.120835,0.323172,0.583739,0.732924,-0.220603,-0.594121,0.935093,-0.216736,0.659318,-0.750417,-0.284773,-0.271496,0.491731,-0.712174,-0.763681,0.0781023,0.951666,0.734031,0.826912,0.57488,-0.361951,-0.0739728,0.91438,-0.391653,0.0193887,0.412634,-0.169813,0.471794,0.660792,-0.350906,-0.612644,0.347876,0.112573,-0.501126,0.456761,-0.109004,0.289352,-0.566504,0.585042,0.584934,0.923676,0.895312,-0.161036,-0.995895,0.0853141,-0.583368,-0.157612,0.234119,0.875043,0.430805,0.706102,0.423887,0.296828,-0.265607,0.338806,-0.15829,0.642516,0.355126,0.174447,-0.975015,0.869905,-0.145285,-0.484002,-0.475966,-0.67704,0.996452,-0.0685748,-0.851985,0.416498,0.791047,-0.211323,-0.302819,0.640735,-0.317908,-0.116586,-0.896382,-0.817317,-0.948837,-0.597427,0.975863,-0.971371,-0.124115,0.4339,-0.254671,0.298068,-0.349803,-0.73185,0.488106,-0.0495073,0.253969,0.168116,0.148772,0.889593,-0.512213,-0.165437,0.666773,-0.976304,-0.170024,0.905794,0.473908,-0.855725,-0.0413591,-0.508661,0.443453,0.842925,-0.144503,0.936699,-0.443935,-0.182996,0.803564,0.960386,-0.0323329,0.638181,-0.895684,-0.360502,0.0646258,-0.202449,-0.717228,0.970489,0.404608,-0.0861868,-0.879417,-0.866462,-0.938336,-0.799685,0.213464,-0.932344,-0.668236,0.751366,-0.22712,-0.407783,0.657463,0.0970092,-0.579109,-0.868866,-0.504041,0.926483,0.169978,-0.00842563,-0.530324,0.282745,0.0255867,0.287686,0.410417,-0.766576,-0.536566,-0.628288,0.69665,0.820713,-0.506162,-0.404114,0.640099,-0.956123,-0.576586,0.435502,-0.470676,-0.367062,-0.831765,-0.294942,0.518991,0.922338,0.337886,-0.67474,-0.725667,0.916684,0.39175,0.759081,0.496979,-0.200691,0.0417966,-0.687391,0.438773,0.287357,0.316636,-0.262311,-0.0755541,-0.442313,0.621378,0.670105,0.060982,0.944162,0.643442,-0.750684,-0.639973,0.217424,0.592823,0.929094,-0.239135,-0.41628,0.570893,-0.0798988,-0.917135,-0.749545,-0.982047,0.0626998,-0.977963,0.660401,0.470569,-0.0528868,-0.00138645,0.931065,-0.748519,0.304188,-0.266153,0.672524,-0.105179,-0.874749,-0.154355,-0.774656,-0.69654,0.433098,0.615897,-0.387919,-0.429779,0.650202,0.122306,-0.237727,0.626817,-0.227929,0.405916,0.483328,0.282047,-0.262206,0.784123,0.83125,-0.662272,0.702768,0.875814,-0.701221,0.553793,0.471795,0.769147,0.059668,-0.841617,-0.191179,-0.972471,-0.825361,0.779826,-0.917201,0.43272,0.10301,0.358771,0.793448,-0.0379954,-0.870112,0.600442,-0.990603,0.549151,0.512146,-0.795843,0.490091,0.372046,-0.549437,0.0964285,0.753047,-0.86284,-0.589688,0.178612,-0.720358,};
  int i, j;
  // FILE *myFile;
  // for(i = 0; i < m; i++){
  //   M[i] = (double *)malloc(sizeof(double)*m);
  // }
  // myFile = fopen("./data/L3ZDT1.txt", "r");
  // for (i=0; i<m; i++)
  // {
  //   for (j=0; j<m; j++)
  //   {
  //     if(!	fscanf(myFile, "%lf\n", &M[i][j]))
	 //     break;
  //   }
  // }
  // fclose(myFile);

  for(i = 0; i < m; i++){
    if(x[i] >= 0 && x[i] <= 1)
      continue;
    else{
      printf("Invalid input range\n");
      return;
    }
  }

  y = (double*)malloc(sizeof(double)*m);
  for(i = 0; i < m; i++){
    y[i] = 0;
    for(j = 0; j < m; j++){
      y[i] = y[i] + ( M[i*m+j] * (x[j]*x[j]) );
    }
  }

  ff1 = y[0]*y[0];
  gx = 0;

  for(i = 1; i < m; i++){
    gx += y[i]*y[i] - 10*cos(4*pi*y[i]);
  }
  gx = gx + (10.0*(m-1));
  gx = gx + 1;

  h = 1 - sqrt(ff1/gx);

  f[0] = ff1;
  f[1] = gx*h;

  free(y);

  return;
}

/*
 *
 *  As described by K. Deb and A. Sinha and S. Kukkonen in "Multi-objective
 *  test problems, linkages, and evolutionary methodologies", GECCO'06}:
 *  Proceedings of the 8th Annual Conference on Genetic and Evolutionary
 *  Computation, 1141-1148, 2006.
 *
 *  Example T6, with linkage L3.
 *
 *  In the above paper the number of variables was set to 30.
 *  We selected n=10 according to the dimension of problem ZDT6.
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void L3ZDT6(double* f, double* x){
  int m = 10;

  double *y;
  double ff1, gx, h;

  int i, j;
  double M[900] = {0.218418,-0.620254,0.843784,0.914311,-0.788548,0.428212,0.103064,-0.47373,-0.300792,-0.185507,0.330423,0.151614,0.884043,-0.272951,-0.993822,0.511197,-0.0997948,-0.659756,0.575496,0.675617,0.180332,-0.593814,-0.492722,0.0646786,-0.666503,-0.945716,-0.334582,0.611894,0.281032,0.508749,-0.0265389,-0.920133,0.308861,-0.0437502,-0.374203,0.207359,-0.219433,0.914104,0.184408,0.520599,-0.88565,-0.375906,-0.708948,-0.37902,0.576578,0.0194674,-0.470262,0.572576,0.351245,-0.480477,0.238261,-0.1596,-0.827302,0.669248,0.494475,0.691715,-0.198585,0.0492812,0.959669,0.884086,-0.218632,-0.865161,-0.715997,0.220772,0.692356,0.646453,-0.401724,0.615443,-0.0601957,-0.748176,-0.207987,-0.865931,0.613732,-0.525712,-0.995728,0.389633,-0.064173,0.662131,-0.707048,-0.340423,0.60624,0.0951648,-0.160446,-0.394585,-0.167581,0.0679849,0.449799,0.733505,-0.00918638,0.00446808,0.404396,0.449996,0.162711,0.294454,-0.563345,-0.114993,0.549589,-0.775141,0.677726,0.610715,0.0850755,0.0419388,-0.323614,-0.973719,-0.680238,-0.270873,-0.209617,0.968436,0.908798,0.975851,-0.994918,-0.0621977,0.628171,0.761228,0.34372,-0.792042,-0.144765,-0.965748,0.0133606,-0.0260565,-0.742377,0.426391,0.408202,0.633885,-0.0351053,-0.723444,-0.577654,0.0276004,0.0712472,-0.622791,0.155451,0.442717,-0.792786,0.925785,0.670266,-0.865566,-0.638281,0.333094,0.477628,0.47261,0.23151,0.82132,-0.589803,0.796275,0.57713,0.101149,0.970191,0.532821,0.814769,-0.0687269,0.712758,-0.191812,-0.390938,0.952828,0.921519,0.923094,0.93011,-0.945394,-0.0934027,0.964123,-0.795609,-0.289563,0.614236,-0.670585,0.466877,0.144597,-0.206416,0.6937,-0.967958,-0.0951247,-0.942473,-0.610767,-0.655472,-0.0960986,-0.302779,-0.734976,-0.342188,-0.315861,-0.912834,0.24499,0.0969326,0.089775,-0.241157,0.0835558,-0.420236,-0.686633,-0.711276,-0.00325377,0.435196,-0.710002,0.00283691,-0.168757,-0.134045,-0.655235,0.172361,0.998291,0.376291,-0.962215,-0.363174,-0.88777,-0.519929,-0.560554,-0.984415,0.601529,-0.984103,-0.228237,-0.578066,0.307023,0.606123,0.959635,0.00225943,0.0101814,0.441456,0.0633629,0.406631,-0.0100638,-0.177972,-0.491075,0.537035,-0.924987,-0.699424,0.742285,0.0181443,0.718971,-0.0308272,0.086931,0.524476,0.956457,0.143024,0.616481,0.217909,-0.128427,-0.262427,-0.938208,-0.52479,0.12919,0.721925,0.766492,0.470845,-0.0976466,0.507807,0.804148,0.963269,0.357128,-0.832565,-0.312441,0.327779,0.184745,0.246139,-0.936814,-0.931734,-0.0327827,0.319293,0.044473,-0.641645,0.596118,-0.293934,-0.63373,0.409658,0.759892,-0.257078,0.939616,-0.227661,0.115754,0.10964,-0.240557,0.66842,0.855535,-0.451536,0.264961,-0.61366,-0.204783,-0.842476,-0.249524,-0.0985226,0.0671501,-0.527707,-0.509489,-0.883254,0.14851,-0.906465,0.496238,-0.853211,-0.779234,-0.979515,0.827175,0.228969,-0.402829,-0.970118,0.762559,0.506495,0.460303,0.897304,0.686003,0.739986,0.15731,0.281697,-0.922955,-0.780824,0.449716,0.125225,0.487649,0.147046,0.679639,0.593707,-0.311828,-0.797099,-0.35815,0.95808,0.907244,0.772426,0.720574,-0.873217,0.371431,-0.826029,0.942716,0.70609,-0.658158,-0.782185,-0.806743,-0.627986,-0.405551,-0.258495,-0.796524,0.222498,0.087545,-0.0917108,-0.62542,-0.110256,0.0417576,0.24476,0.941339,-0.613783,0.402772,0.300775,-0.820314,-0.894233,-0.405896,0.0735439,0.486645,-0.394355,0.125097,-0.316386,-0.701215,-0.845742,0.2065,-0.413743,0.406725,-0.423813,-0.941255,-0.558804,0.312326,0.345314,0.319143,-0.644653,-0.0408415,0.176461,0.740113,0.470737,-0.914927,-0.591523,-0.606614,-0.181873,0.692975,0.50208,-0.536704,0.359652,0.839082,0.56817,-0.0776788,-0.00332785,0.459538,-0.518313,-0.270738,-0.629958,-0.755084,-0.721573,0.431107,-0.221877,0.32543,0.163743,0.0759916,0.695064,-0.656856,0.074163,0.264319,-0.73174,0.731548,-0.489341,0.678946,0.0271269,0.804879,-0.402973,0.800373,0.760082,-0.878364,0.176801,-0.548932,-0.225601,-0.164912,-0.208143,0.7768,-0.542743,-0.156021,0.671736,0.878648,-0.419588,-0.0752896,0.0299447,-0.494459,-0.72415,0.35978,-0.32646,-0.96605,0.0127605,0.563174,-0.814853,-0.949609,-0.526794,-0.801902,-0.753397,0.617418,0.689874,0.983384,0.668786,0.0304653,-0.625221,-0.13318,0.827343,-0.101358,-0.999522,-0.0525574,-0.458319,0.587409,-0.334639,0.0759643,0.0255827,0.128944,0.17317,-0.284309,0.287161,-0.550725,-0.433083,-0.242821,0.878879,0.691699,-0.660499,0.389985,0.599856,-0.711442,-0.798697,-0.244945,-0.942649,0.402856,-0.494672,0.439941,-0.88216,0.170196,0.650734,-0.0982391,-0.468732,0.342133,-0.838071,-0.832362,0.658177,-0.565361,0.149473,0.69331,-0.491848,0.74916,0.526025,-0.155339,0.0998096,0.468761,0.324649,0.128488,0.544144,-0.495222,0.965229,-0.79314,-0.545421,-0.500243,0.154371,0.170017,-0.259108,-0.868862,-0.50731,-0.848317,0.835712,0.616391,-0.442608,-0.158,0.313451,0.703748,-0.755984,-0.249443,0.491564,0.985068,0.678644,0.808324,0.81975,-0.435823,-0.839855,0.00282368,-0.569165,0.0884339,-0.222144,0.499412,-0.565198,0.64824,0.956914,-0.0620912,0.634479,0.928617,0.464664,0.377022,0.63047,-0.198619,-0.576153,0.565373,-0.524245,-0.187299,-0.614524,0.429316,-0.491171,0.399495,-0.333898,-0.646636,-0.0189709,-0.339605,-0.798791,0.0494081,0.367012,0.852545,0.43557,0.150039,-0.0454542,0.604861,-0.598288,-0.500696,0.249008,0.370711,-0.633174,-0.0121906,0.42006,0.169373,-0.975542,-0.0297852,0.80481,0.638317,-0.670967,0.935792,-0.35605,0.175773,0.878601,-0.275168,-0.932517,-0.372497,-0.0732907,-0.185493,-0.357004,0.314786,-0.229239,0.530256,-0.51327,0.44187,0.940309,-0.240334,-0.0276121,0.74383,-0.630329,-0.763376,0.62538,0.818945,0.891598,0.680494,0.471868,-0.769787,-0.878099,-0.973724,0.354362,-0.1792,-0.225034,-0.44548,0.598865,0.544005,-0.478962,0.327193,-0.525784,0.903179,-0.899248,0.156514,0.154329,0.499808,-0.836327,-0.802627,0.378082,-0.112673,-0.47926,-0.3355,-0.699445,0.237731,-0.324597,-0.800406,-0.42585,-0.710739,-0.144068,-0.828545,-0.800912,0.184654,-0.63675,-0.16696,0.240427,-0.513443,0.812664,0.744943,0.970612,0.00172899,-0.726378,-0.0985012,0.224232,0.16495,0.560077,-0.813112,0.112894,-0.0955366,0.0187107,0.913887,0.123076,0.550338,0.400334,-0.367816,0.198455,-0.983183,0.278976,0.714817,0.307911,0.812861,-0.403497,-0.784382,-0.161823,-0.120835,0.323172,0.583739,0.732924,-0.220603,-0.594121,0.935093,-0.216736,0.659318,-0.750417,-0.284773,-0.271496,0.491731,-0.712174,-0.763681,0.0781023,0.951666,0.734031,0.826912,0.57488,-0.361951,-0.0739728,0.91438,-0.391653,0.0193887,0.412634,-0.169813,0.471794,0.660792,-0.350906,-0.612644,0.347876,0.112573,-0.501126,0.456761,-0.109004,0.289352,-0.566504,0.585042,0.584934,0.923676,0.895312,-0.161036,-0.995895,0.0853141,-0.583368,-0.157612,0.234119,0.875043,0.430805,0.706102,0.423887,0.296828,-0.265607,0.338806,-0.15829,0.642516,0.355126,0.174447,-0.975015,0.869905,-0.145285,-0.484002,-0.475966,-0.67704,0.996452,-0.0685748,-0.851985,0.416498,0.791047,-0.211323,-0.302819,0.640735,-0.317908,-0.116586,-0.896382,-0.817317,-0.948837,-0.597427,0.975863,-0.971371,-0.124115,0.4339,-0.254671,0.298068,-0.349803,-0.73185,0.488106,-0.0495073,0.253969,0.168116,0.148772,0.889593,-0.512213,-0.165437,0.666773,-0.976304,-0.170024,0.905794,0.473908,-0.855725,-0.0413591,-0.508661,0.443453,0.842925,-0.144503,0.936699,-0.443935,-0.182996,0.803564,0.960386,-0.0323329,0.638181,-0.895684,-0.360502,0.0646258,-0.202449,-0.717228,0.970489,0.404608,-0.0861868,-0.879417,-0.866462,-0.938336,-0.799685,0.213464,-0.932344,-0.668236,0.751366,-0.22712,-0.407783,0.657463,0.0970092,-0.579109,-0.868866,-0.504041,0.926483,0.169978,-0.00842563,-0.530324,0.282745,0.0255867,0.287686,0.410417,-0.766576,-0.536566,-0.628288,0.69665,0.820713,-0.506162,-0.404114,0.640099,-0.956123,-0.576586,0.435502,-0.470676,-0.367062,-0.831765,-0.294942,0.518991,0.922338,0.337886,-0.67474,-0.725667,0.916684,0.39175,0.759081,0.496979,-0.200691,0.0417966,-0.687391,0.438773,0.287357,0.316636,-0.262311,-0.0755541,-0.442313,0.621378,0.670105,0.060982,0.944162,0.643442,-0.750684,-0.639973,0.217424,0.592823,0.929094,-0.239135,-0.41628,0.570893,-0.0798988,-0.917135,-0.749545,-0.982047,0.0626998,-0.977963,0.660401,0.470569,-0.0528868,-0.00138645,0.931065,-0.748519,0.304188,-0.266153,0.672524,-0.105179,-0.874749,-0.154355,-0.774656,-0.69654,0.433098,0.615897,-0.387919,-0.429779,0.650202,0.122306,-0.237727,0.626817,-0.227929,0.405916,0.483328,0.282047,-0.262206,0.784123,0.83125,-0.662272,0.702768,0.875814,-0.701221,0.553793,0.471795,0.769147,0.059668,-0.841617,-0.191179,-0.972471,-0.825361,0.779826,-0.917201,0.43272,0.10301,0.358771,0.793448,-0.0379954,-0.870112,0.600442,-0.990603,0.549151,0.512146,-0.795843,0.490091,0.372046,-0.549437,0.0964285,0.753047,-0.86284,-0.589688,0.178612,-0.720358,};
  // FILE *myFile;

  // for(i = 0; i < m; i++){
  //   M[i] = (double *)malloc(sizeof(double)*m);
  // }
  // myFile = fopen("./data/L3ZDT1.txt", "r");
  // for (i=0; i<m; i++)
  // {
  //   for (j=0; j<m; j++)
  //   {
  //     if(!	fscanf(myFile, "%lf\n", &M[i][j]))
	 //     break;
  //   }
  // }
  // fclose(myFile);

  for(i = 0; i < m; i++){
    if(x[i] >= 0 && x[i] <= 1)
      continue;
    else{
      printf("Invalid input range\n");
      return;
    }
  }

  y = (double*)malloc(sizeof(double)*m);
  for(i = 0; i < m; i++){
    y[i] = 0;
    for(j = 0; j < m; j++){
      y[i] = y[i] + ( M[i*m+j] * (x[j]*x[j]) );
    }
  }

  ff1 = y[0]*y[0];
  gx = 0;

  for(i = 1; i < m; i++){
    gx += y[i]*y[i];
  }
  gx = 9.0*pow(gx/(m-1), 0.25);
  gx = gx + 1;

  h = 1 - (ff1/gx)*(ff1/gx);

  f[0] = ff1;
  f[1] = gx*h;

  free(y);

  return;
}

/*
 *
 *  As described by Huband et al. in "A review of multiobjective test problems
 *  and a scalable test problem toolkit", IEEE Transactions on Evolutionary
 *  Computing 10(5): 477-506, 2006.
 *
 *  Example VFM1, see the previous cited paper for the original reference.
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void VFM1(double* f, double* x){

    int i;

    for(i=0; i<2; ++i) {
        if( x[i]<-2 || x[i]>2 ) {
            printf("Invalid input range\n");
            return;
        }
    }

    f[0] = x[0]*x[0] + (x[1]-1)*(x[1]-1);
    f[1] = x[0]*x[0] + (x[1]+1)*(x[1]+1) + 1;
    f[2] = (x[0]-1)*(x[0]-1) + (x[1])*(x[1]) + 2;

    return;
}

/*
 *
 *  As described by Huband et al. in "A review of multiobjective test problems
 *  and a scalable test problem toolkit", IEEE Transactions on Evolutionary
 *  Computing 10(5): 477-506, 2006.
 *
 *  Example VU1, see the previous cited paper for the original reference.
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void VU1(double* f, double* x){

    int i;

    for(i=0; i<2; ++i) {
        if( x[i]<-3 || x[i]>3 ) {
            printf("Invalid input range\n");
            return;
        }
    }

    f[0] = 1/(x[0]*x[0] + x[1]*x[1] + 1);
    f[1] = x[0]*x[0] + 3*x[1]*x[1] + 1;

    return;
}

/*
 *
 *  As described by Huband et al. in "A review of multiobjective test problems
 *  and a scalable test problem toolkit", IEEE Transactions on Evolutionary
 *  Computing 10(5): 477-506, 2006.
 *
 *  Example VU2, see the previous cited paper for the original reference.
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void VU2(double* f, double* x){

    int i;

    for(i=0; i<2; ++i) {
        if( x[i]<-3 || x[i]>3 ) {
            printf("Invalid input range\n");
            return;
        }
    }

    f[1] = x[0] + x[1] + 1;
    f[0] = x[0]*x[0] + 2*x[1] - 1;

    return;
}

/*
 *
 *  As described by E. Zitzler, K. Deb, and L. Thiele in "Comparison of
 *  Multiobjective Evolutionary Algorithms: Empirical Results", Evolutionary
 *  Computation 8(2): 173-195, 2000.
 *
 *  Example T1.
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void ZDT1(double* f, double* x){
    int m = 30;
    double ff1, sum, gx, h;
    int i;

    for(i=0; i<m; ++i) {
        if( x[i]<0 || x[i]>1 ) {
            printf("Invalid input range\n");
            return;
        }
    }

    ff1 = x[0];
    sum = 0;
    for(i=1; i<m; ++i)
      sum += x[i];

    gx = 1 + 9.0/(m-1) * sum;
    h = 1 - sqrt(ff1/gx);

    f[0] = ff1;
    f[1] = gx*h;

    return;
}

/*
 *
 *  As described by E. Zitzler, K. Deb, and L. Thiele in "Comparison of
 *  Multiobjective Evolutionary Algorithms: Empirical Results", Evolutionary
 *  Computation 8(2): 173-195, 2000.
 *
 *  Example T2.
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void ZDT2(double* f, double* x){
    int m = 30;
    double ff1, sum, gx, h;
    int i;

    for(i=0; i<m; ++i) {
        if( x[i]<0 || x[i]>1 ) {
            printf("Invalid input range\n");
            return;
        }
    }

    ff1 = x[0];
    sum = 0;
    for(i=1; i<m; ++i)
      sum += x[i];

    gx = 1 + 9.0/(m-1) * sum;
    h = 1 - (ff1/gx)*(ff1/gx);

    f[0] = ff1;
    f[1] = gx*h;

    return;
}

/*
 *
 *  As described by E. Zitzler, K. Deb, and L. Thiele in "Comparison of
 *  Multiobjective Evolutionary Algorithms: Empirical Results", Evolutionary
 *  Computation 8(2): 173-195, 2000.
 *
 *  Example T3.
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void ZDT3(double* f, double* x){
    int m = 30;
    double ff1, pi, sum, gx, h;
    int i;

    for(i=0; i<m; ++i) {
        if( x[i]<0 || x[i]>1 ) {
            printf("Invalid input range\n");
            return;
        }
    }

    ff1 = x[0];
    pi = 4*atan(1);
    sum = 0;
    for(i=1; i<m; ++i)
      sum += x[i];

    gx = 1 + 9.0/(m-1) * sum;
    h = 1 - sqrt(ff1/gx) - (ff1/gx)*sin(10*pi*ff1);

    f[0] = ff1;
    f[1] = gx*h;

    return;
}

/*
 *
 *  As described by E. Zitzler, K. Deb, and L. Thiele in "Comparison of
 *  Multiobjective Evolutionary Algorithms: Empirical Results", Evolutionary
 *  Computation 8(2): 173-195, 2000.
 *
 *  Example T4.
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void ZDT4(double* f, double* x){
    int m = 10;
    double ff1, pi, sum, gx, h;
    int i;

    if( x[0]<0 || x[0]>1 ) {
        printf("Invalid input range\n");
        return;
    }

    for(i=1; i<m; ++i) {
        if( x[i]<-5 || x[i]>5 ) {
            printf("Invalid input range\n");
            return;
        }
    }

    ff1 = x[0];
    pi = 4*atan(1);
    sum = 0;
    for(i=1; i<m; ++i)
      sum += (x[i]*x[i] - 10*cos(4*pi*x[i]));

    gx = 1 + 10*(m-1) + sum;
    h = 1 - sqrt(ff1/gx);

    f[0] = ff1;
    f[1] = gx*h;

    return;
}

/*
 *
 *  As described by E. Zitzler, K. Deb, and L. Thiele in "Comparison of
 *  Multiobjective Evolutionary Algorithms: Empirical Results", Evolutionary
 *  Computation 8(2): 173-195, 2000.
 *
 *  Example T6.
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void ZDT6(double* f, double* x){
    int m = 10;
    double ff1, pi, sum, gx, h;
    int i;

    for(i=0; i<m; ++i) {
        if( x[i]<0 || x[i]>1 ) {
            printf("Invalid input range\n");
            return;
        }
    }
    pi = 4*atan(1);
    ff1 = 1 - exp(-4*x[0])*pow(sin(6*pi*x[0]), 6);

    sum = 0;

    for(i=1; i<m; ++i)
      sum += x[i];

    gx = 1 + 9*pow(sum/(m-1), 0.25);
    h = 1 - (ff1/gx)*(ff1/gx);

    f[0] = ff1;
    f[1] = gx*h;

    return;
}

/*
 *
 *  As described by Huband et al. in "A review of multiobjective test problems
 *  and a scalable test problem toolkit", IEEE Transactions on Evolutionary
 *  Computing 10(5): 477-506, 2006.
 *
 *  Example ZLT1, see the previous cited paper for the original reference.
 *
 *  In the above paper the number of variables was set equal to 100.
 *  We selected n=10 as default.
 *
 *  This file is part of a collection of problems developed for
 *  the BMOBench platform and based on
 *  derivative-free multiobjective optimization in
 *  A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
 *
 */
void ZLT1(double* f, double* x){
    int n = 10;
    int m = 3;
    double sum;
    int i, j;

    for(i=0; i<n; ++i) {
        if( x[i]<-1000 || x[i]>1000 ) {
            printf("Invalid input range\n");
            return;
        }
    }

    for(i=0; i<m; ++i){
      sum = 0;
      for(j=0; j<n; ++j)
        sum = sum + ((i!=j) ? (x[j]*x[j]) : (0));
      f[i] = pow(x[i]-1, 2) + sum;
    }

    return;
}


/*
*
*
*   As described by K. Deb in "Multi-objective Genetic Algorithms: Problem
*   Difficulties and Construction of Test Problems", Evolutionary Computation
*   7(3): 205-230, 1999.
*
*   Example 5.1.2 (Convex Pareto-optimal Front).
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void Deb512a(double* f, double* x){

	//Parameters
	double beta ;
	double alpha ;
	double ff1 ;
	double gx;
	double h;

	beta = 1;
	alpha = 0.25;
	ff1 = 4*x[0];
	h = 0;

	if(x[0] < 0.0 || x[0] >1.0 || x[1] < 0.0 || x[1] > 1.0 ){
		printf("Invalid input range\n");
		return;
	}

	if(x[1]<=0.4)
		gx = 4 - (3*exp(-1*pow((x[1] - 0.2)/0.02 , 2)));
	else
		gx = 4-2*exp(-1*pow((x[1]-0.7)/0.2 , 2));

	if(ff1 <= beta*gx)
		h = 1-pow(ff1/(beta*gx), alpha);

	f[0] =  ff1;
	f[1] = gx*h;

	return;

}


/*
*
*
*   As described by K. Deb in "Multi-objective Genetic Algorithms: Problem
*   Difficulties and Construction of Test Problems", Evolutionary Computation
*   7(3): 205-230, 1999.
*
*   Example 5.2.1 (Biased Search Space).
*
*   This file is part of a collection of problems developed for
*   derivative-free multiobjective optimization in
*   the BMOBench platform and based on
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void Deb521a(double *f, double *x) {

	double gamma;
	double ff1;
	double gx;
	double temp;
	double h;

	gamma = 0.25;
	ff1 = x[0];
	gx = 1 + pow(x[1], gamma);
	temp = ff1/gx;
	h = 1 - pow(temp, 2);

	if(x[0] < 0.0 || x[0] >1.0 || x[1] < 0.0 || x[1] > 1.0 ){
		printf("Invalid input range\n");
		return;
	}

	f[0] = ff1;
	f[1] = gx*h;

	return;
}


/*
*
*
*   As described by K. Deb, L. Thiele, M. Laumanns and E. Zitzler in "Scalable
*   Multi-Objective Optimization Test Problems", Congress on Evolutionary
*   Computation (CEC’2002): 825-830, 2002.
*
*   Example DTLZ1.
*
*   This file is part of a collection of problems developed for
*   derivative-free multiobjective optimization in
*   the BMOBench platform and based on
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void DTLZ1(double *f, double *x) {
	int M;
	int n;
	double k;
	double pi;
	int b;
	int i;
	int j;
	double sum;
	double prod;
	double gx;
	double ff1;
	double *ff;

	M = 3;
	n = 7;
	k = n - M + 1;
	pi = 4 * atan(1);
	b = 0;
	i = 0;
	j = 0;
	sum = 0;
	prod = 1;
	ff = (double*)malloc((M-1)*sizeof(double));

	for(i=0;i<n;i++){
		if(x[i]<0.0 || x[i]>1.0){
			printf("Invalid input range\n");
			return;
		}
	}

	if(M<2){
		printf("Invalid value of M");
	}

	if(n<M){
		printf("Invalid value of N");
	}

	for(i = M-1 ; i <= n-1 ; i++){
		sum = sum + (pow((x[i]-0.5), 2)-cos(20*pi*(x[i]-0.5)));
	}

	for(j = 0 ; j <= M-2 ; j++){
		prod = prod*x[j];
	}

	gx = 100.0*(k + sum);

	ff1 = 0.5*(1+gx)*prod;

	i = 0;
	j=0;
	for (i = 1; i <= M-1; i++)
	{
		prod = 1;
		for (j = 0; j < M - i - 1; j++)
		{
			prod = prod*x[j];
		}
		ff[i] = 0.5*(1 + gx)*prod*(1 - x[M - i - 1]);
	}

	f[0] = ff1;

	for(i=1;i<M;i++)
	{
		f[i] = ff[i];
	}
	return;
}


/*
*
*
*   As described by K. Deb, L. Thiele, M. Laumanns and E. Zitzler in "Scalable
*   Multi-Objective Optimization Test Problems", Congress on Evolutionary
*   Computation (CEC’2002): 825-830, 2002.
*
*   Example DTLZ3.
*
*   This file is part of a collection of problems developed for
*   derivative-free multiobjective optimization in
*   the BMOBench platform and based on
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void DTLZ3(double *f, double *x) {

	int M;
	int n ;
	double k;
	double pi;
	int b;
	int i;
	int j;
	double temp;
	double sum;
	double prod ;
	double gx;
	double ff1;
	double *ff;


	M = 3;
	n = 12;
	k = n - M + 1;
	pi = 4 * atan(1);
	b = 0;
	i = 0;
	j = 0;
	sum = 0;
	prod = 1;
 	ff = (double*)malloc((M-1)*sizeof(double));

	if(M<3){
		printf("Invalid value of M");
	}

	if(n<M){
		printf("Invalid value of N");
	}




	for(b=0;b<n;b++){
		if(x[i]<0.0 || x[i]>1.0){
			printf("Invalid input range\n");
			return;
		}
	}


	for(i = M-1 ; i <= n-1 ; i++){
		temp = pow((x[i]-0.5), 2)-cos(20*pi*(x[i]-0.5));
		sum = sum + temp;
	}

	for(j = 0 ; j <= M-2 ; j++){
		prod = prod*cos(0.5*pi*x[j]);
	}

	gx = 100*(k + sum);

	ff1 = (1+gx)*prod;

	i = 0;
	j=0;
	for (i = 1; i <= M-1; ++i)
	{
		prod = 1;
		for (j = 0; j < M - i - 1; ++j)
		{
			prod = prod*cos(0.5*pi*x[j]);
		}
		ff[i] = (1 + gx)*prod*(sin(0.5*pi*x[M - i - 1]));
	}

	f[0] = ff1;
	for(i=1;i<M;i++)
	{
		f[i] = ff[i];
	}

	return;
}


/*
*
*
*   As described by K. Deb, L. Thiele, M. Laumanns and E. Zitzler in "Scalable
*   Multi-Objective Optimization Test Problems", Congress on Evolutionary
*   Computation (CEC’2002): 825-830, 2002.
*
*   Example DTLZ5.
*
*   This file is part of a collection of problems developed for
*   derivative-free multiobjective optimization in
*   the BMOBench platform and based on
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void DTLZ5(double *f, double *x) {
  	int M;
	int n ;
	double pi;
	int b;
	int i;
	int j;
	double temp;
	double sum;
	double prod ;
	double gx;
	double *theta;
	double ff1;
	double *ff;
	double ffM;


	M = 3;
	n = 12;
	pi = 4 * atan(1);
	b = 0;
	i = 0;
	j = 0;
	sum = 0.0;
	prod = 1.0;
	theta = (double *)malloc((n)*sizeof(double));
 	ff = (double*)malloc((M)*sizeof(double));

	if(M<3){
		printf("Invalid value of M");
	}

	if(n<M){
		printf("Invalid value of N");
	}



	for(b=0;b<n;b++){
		if(x[b]<0.0 || x[b]>1.0){
			printf("Invalid input range\n");
			return;
		}
	}

	for(i = M-1 ; i <= n-1 ; i++){
		temp = pow(x[i], 0.1);
		sum = sum + temp;
	}

	gx = sum;


	for(i = 1;i<=n-1;i++){
		theta[i] = (pi/2.0)*(1.0 + (2.0*gx*x[i]))/(2.0*(1.0 + gx));
	}

	for(j = 1 ; j <= M-2 ; j++){
		prod = prod*cos(theta[j]);
	}



	ff1 = (1.0+gx)*prod*(cos(0.5*pi*x[0]));

	i = 0;
	j=0;
	for (i = 1; i <= M-2; ++i)
	{
		prod = 1.0;
		for (j = 1; j < M - i - 1; ++j)
		{
			prod = prod*cos(theta[j]);
		}
		ff[i] = (1 + gx)*cos(0.5*pi*x[0])*prod*(sin(theta[M - i -1]));
	}

	ffM = (1 + gx) * sin(0.5*pi*x[0]);

	f[0] = ff1;
	f[1] = ff[1];
	f[2] = ffM;

	return;
}


/*
* ex005.mod
* Original AMPL coding by Sven Leyffer, Argonne Natl. Lab.
*
* A simple multi-objective optimization problem (p. 281):
* C.-L. Hwang & A. S. Md. Masud, Multiple Objective
* Decision Making - Methods and Applications, No. 164 in
* Lecture Notes in Economics and Mathematical Systems,
* Springer, 1979.
*/
void ex005(double *f, double *x) {

	if(x[0] < -1 || x[0] > 2 || x[1] < 1 || x[1] > 2 ){
		printf("Invalid input range\n");
		return;
	}

	f[0] = pow(x[0], 2) - pow(x[1], 2);
	f[1] = x[0]/x[1];


	return;
}


/*
*
*
*   As described by Huband et al. in "A review of multiobjective test problems
*   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
*   Computing 10(5): 477-506, 2006.
*
*   Example FES3, see the previous cited paper for the original reference.
*
*   In the above paper the number of variables was left undefined.
*   We selected n=10 as default.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void FES3(double *f, double *x) {

	double n;
	double pi;
	int b;
	int i;
	double sum1;
	double sum2;
	double sum3;
	double sum4;
	double temp;

	n = 10;
	pi = 4 * atan(1);
	b = 0;

	for(b=0;b<n;b++){
		if(x[b]<0.0 || x[b]>1.0){
			printf("Invalid input range\n");
			return;
		}
	}

	sum1 = 0;
	for(i = 0; i <= n - 1; i++){
		temp = fabs(x[i] - exp(pow(((i+1)/n), 2))/3);
		sum1 = sum1 + pow(temp, 0.5);
	}

	sum2 = 0;
	for(i = 0; i <= n - 1; i++){
		temp = fabs(x[i] - pow(sin(i), 2)*pow(cos(i), 2));
		sum2 = sum2 + pow(temp, 0.5);
	}

	sum3 = 0;
	for(i = 0; i <= n - 1; i++){
		temp = fabs(x[i] - (0.25*cos(i)*cos(2*i)) - 0.5);
		sum3 = sum3 + pow(temp, 0.5);
	}

	sum4 = 0;
	for(i = 0; i <= n - 1; i++){
		temp = fabs(x[i] - (0.5*sin(1000*pi*(i+1)/n)) - 0.5);
		sum4 = sum4 + pow(temp, 2);
	}

	f[0] = sum1;
	f[1] = sum2;
	f[2] = sum3;
	f[3] = sum4;
	return;
}


/*
*
*
*   As described by Huband et al. in "A Scalable Multi-objective Test Problem
*   Toolkit", in C. A. Coello Coello et al. (Eds.): EMO 2005, LNCS 3410,
*   pp. 280–295, 2005, Springer-Verlag Berlin Heidelberg 2005.
*
*   Example I2.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void I2(double* f, double* zz){

    /*Note: Input variable here is z and not x*/

      int M;
      int k;
      int l;
      int n;

      double pi2;

      double *S, *A;
      double *zmax, *x, *y, *t1, *t2, *t3, *h, *w, *r_sum, *z;
      int i, j;
      double sum1, sum2, AA, BB, CC, temp;

      M = 3;
      k = 4;
      l = 4;
      n = k+l;
      pi2 = 2*atan(1);
      S = (double*)malloc(sizeof(double)*(M+1));
      A = (double*)malloc(sizeof(double)*(M));
      zmax = (double*)malloc(sizeof(double)*(n+1));
      y = (double*)malloc(sizeof(double)*(n+1));
      t1 = (double*)malloc(sizeof(double)*(n+1));
      t2 = (double*)malloc(sizeof(double)*(n+1));
      t3 = (double*)malloc(sizeof(double)*(M+1));
      x = (double*)malloc(sizeof(double)*(M+1));
      h = (double*)malloc(sizeof(double)*(M+1));
      w = (double *)malloc((n+1)*sizeof(double));
      r_sum = (double *)malloc((n+1)*sizeof(double));
      z = (double *)malloc(sizeof(double)*(n+1));

      AA = 0.98/49.98;
      BB = 0.02;
      CC = 50.0;

      for(i=n; i>=1; --i)
      {
        z[i] = zz[i-1];
    	}

      for(i=1; i<=M-1; ++i){
        S[i] = 1.0;
        A[i] = 1.0;
      }
      S[M] = 1.0;

      for(i=1; i<=n; ++i){
        zmax[i] = 1.0;
    	}

      for(i=1;i<=n;i++)
      {
        w[i] = 1.0;
      }

       for(i=1; i<=n; ++i){
        y[i] = (double)z[i]/zmax[i];
    	}

      for(i=1;i<=n-1;i++)
      {
        sum1=0.0;
        sum2=0.0;
        for(j = i+1;j<=n;j++)
        {
          sum1 = sum1 + w[j]*y[j];
          sum2 = sum2 + w[j];
        }
        r_sum[i] = sum1/sum2;
      }

      for(i=1; i<=n; ++i){
        if(i<n)
        {
          temp = BB + (CC-BB) * (AA-(1.0-2.0*(r_sum[i])) * (fabs(floor(0.5-(r_sum[i]))+AA)));
          t1[i] = pow(y[i],temp);
        }
        else
        {
          t1[i] = y[i];
        }

        if(i<=k)
        {
          t2[i] = t1[i];
      	}
        else
        {
          t2[i] = (fabs(t1[i]-0.35)/fabs(floor(0.35-t1[i])+0.35));
      	}
      }

      for(i=1; i<=n; ++i) {
          if( z[i]<0 || z[i]>zmax[i] ) {
              printf("%lf, Invalid input range\n", z[i]);
              return;
          }
      }

      for(i=1; i<=M; ++i){
        if(i<=(M-1)){
          sum1 = 0.0;
          sum2 = 0.0;
          for(j=((i-1)*k/(M-1)+1); j<=((i)*k/(M-1)); ++j){
            sum1 = sum1 + (w[j]*t2[j]);
            sum2 = sum2 + (w[j]);
          }
          t3[i] = sum1/sum2;
        }
        else{
          sum1 = 0.0;
          sum2 = 0.0;
          for(j=k+1; j<=n; ++j){
            sum1 = sum1 + (w[j]*t2[j]);
            sum2 = sum2 + (w[j]);
          }
          t3[i] = sum1/sum2;
        }
      }
      for(i=1; i<=M; ++i){
        if(i<=M-1){
          if(t3[M]>A[i])
          {
            x[i] = t3[M]*(t3[i]-0.5)+0.5;
        	}
          else
          {
            x[i] = A[i]*(t3[i]-0.5)+0.5;
        	}
        }
        else
          x[i] = t3[M];
      }

      for(i=1;i<=M;i++)
      {
      	if(i==1){
  	        h[i] = 1.0;
  	        for(j=1; j<=M-1; ++j)
  	        {
  	          h[i] = h[i] * (sin(x[j]*pi2));
  	      	}
        	}
        	else if(i<=M-1){
  	        h[i] = 1.0;
  	        for(j=1; j<=M-i; ++j)
  	        {
  	          h[i] = h[i] * sin(x[j]*pi2);
  	      	}
  	        h[i] = h[i] * cos(x[M-i+1]*pi2);
        	}
        	else{
          	h[i] = cos(x[1]*pi2);
        	}
      }

      // f = (double*)malloc(sizeof(double)*M);
      for(i=1; i<=M; ++i)
      {
        f[i-1] = x[M]+S[i]*h[i];
    	}

      return;
}


/*
*
*
*   As described by Huband et al. in "A Scalable Multi-objective Test Problem
*   Toolkit", in C. A. Coello Coello et al. (Eds.): EMO 2005, LNCS 3410,
*   pp. 280–295, 2005, Springer-Verlag Berlin Heidelberg 2005.
*
*   Example I3.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
*/
void I3(double* f, double* zz){

  /*Note: Input variable here is z and not x*/

    int M;
    int k;
    int l;
    int n;

    double pi2;

    double *S, *A;
    double *zmax, *x, *y, *t1, *t2, *t3, *h, *w, *r_sum, *z;
    int i, j;
    double sum1, sum2, AA, BB, CC, temp;

    M = 3;
    k = 4;
    l = 4;
    n = k+l;
    pi2 = 2*atan(1);
    S = (double*)malloc(sizeof(double)*(M+1));
    A = (double*)malloc(sizeof(double)*(M));
    zmax = (double*)malloc(sizeof(double)*(n+1));
    y = (double*)malloc(sizeof(double)*(n+1));
    t1 = (double*)malloc(sizeof(double)*(n+1));
    t2 = (double*)malloc(sizeof(double)*(n+1));
    t3 = (double*)malloc(sizeof(double)*(M+1));
    x = (double*)malloc(sizeof(double)*(M+1));
    h = (double*)malloc(sizeof(double)*(M+1));
    w = (double *)malloc((n+1)*sizeof(double));
    r_sum = (double *)malloc((n+1)*sizeof(double));
    z = (double *)malloc(sizeof(double)*(n+1));

    AA = 0.98/49.98;
    BB = 0.02;
    CC = 50.0;

    for(i=n; i>=1; --i)
    {
      z[i] = zz[i-1];
  	}

    for(i=1; i<=M-1; ++i){
      S[i] = 1.0;
      A[i] = 1.0;
    }
    S[M] = 1.0;

    for(i=1; i<=n; ++i){
      zmax[i] = 1.0;
  	}

    for(i=1;i<=n;i++)
    {
      w[i] = 1.0;
    }

     for(i=1; i<=n; ++i){
      y[i] = (double)z[i]/zmax[i];
  	}

    for(i=2;i<=n;i++)
    {
      sum1=0.0;
      sum2=0.0;
      for(j = 1;j<=i-1;j++)
      {
        sum1 = sum1 + w[j]*y[j];
        sum2 = sum2 + w[j];
      }
      r_sum[i] = sum1/sum2;
    }

    for(i=1; i<=n; ++i){
      if(i<2)
      {
        t1[i] = y[i];
      }
      else
      {
        temp = BB + (CC-BB) * (AA-(1.0-2.0*(r_sum[i])) * (fabs(floor(0.5-(r_sum[i]))+AA)));
        t1[i] = pow(y[i],temp);
      }

      if(i<=k)
      {
        t2[i] = t1[i];
    	}
      else
      {
        t2[i] = (fabs(t1[i]-0.35)/fabs(floor(0.35-t1[i])+0.35));
    	}
    }

    for(i=1; i<=n; ++i) {
        if( z[i]<0 || z[i]>zmax[i] ) {
            printf("%lf, Invalid input range\n", z[i]);
            return;
        }
    }

    for(i=1; i<=M; ++i){
      if(i<=(M-1)){
        sum1 = 0.0;
        sum2 = 0.0;
        for(j=((i-1)*k/(M-1)+1); j<=((i)*k/(M-1)); ++j){
          sum1 = sum1 + (w[j]*t2[j]);
          sum2 = sum2 + (w[j]);
        }
        t3[i] = sum1/sum2;
      }
      else{
        sum1 = 0.0;
        sum2 = 0.0;
        for(j=k+1; j<=n; ++j){
          sum1 = sum1 + (w[j]*t2[j]);
          sum2 = sum2 + (w[j]);
        }
        t3[i] = sum1/sum2;
      }
    }
    for(i=1; i<=M; ++i){
      if(i<=M-1){
        if(t3[M]>A[i])
        {
          x[i] = t3[M]*(t3[i]-0.5)+0.5;
      	}
        else
        {
          x[i] = A[i]*(t3[i]-0.5)+0.5;
      	}
      }
      else
        x[i] = t3[M];
    }

    for(i=1;i<=M;i++)
    {
    	if(i==1){
	        h[i] = 1.0;
	        for(j=1; j<=M-1; ++j)
	        {
	          h[i] = h[i] * (sin(x[j]*pi2));
	      	}
      	}
      	else if(i<=M-1){
	        h[i] = 1.0;
	        for(j=1; j<=M-i; ++j)
	        {
	          h[i] = h[i] * sin(x[j]*pi2);
	      	}
	        h[i] = h[i] * cos(x[M-i+1]*pi2);
      	}
      	else{
        	h[i] = cos(x[1]*pi2);
      	}
    }

    for(i=1; i<=M; ++i)
    {
      f[i-1] = x[M]+S[i]*h[i];
  	}

    return;
}


/*
*
*
*   As described by Huband et al. in "A review of multiobjective test problems
*   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
*   Computing 10(5): 477-506, 2006.
*
*   Example IM1, see the previous cited paper for the original reference.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void IM1(double *f, double *x) {

	if(x[0] < 1 || x[0] > 4 || x[1] < 1 || x[1] > 2 ){
		printf("Invalid input range\n");
		return;
	}

	f[0] = 2*sqrt(x[0]);
	f[1] = x[0]*(1 - x[1]) + 5;

	return;
}


/*
*
*
*   As described by Y. Jin, M. Olhofer and B. Sendhoff. "Dynamic weighted
*   aggregation for evolutionary multi-objective optimization: Why does it
*   work and how?", in Proceedings of Genetic and Evolutionary Computation
*   Conference, pp.1042-1049, San Francisco, USA, 2001.
*
*   Test function 4, F4.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void Jin4(double *f, double *x) {
  	double n;
	int b ;
	int i;
	double sum;
	double gx;

	n = 2;
	sum = 0;
	for(b=0;b<n;b++){
		if(x[b]<0.0 || x[b]>1){
			printf("Invalid input range\n");
			return;
		}
	}

	for(i = 1;i<=n-1;i++)
	{
		sum = sum + x[i];
	}

	gx = 1 + (9*sum)/(n - 1);

	f[1] = x[0];
	f[0] = gx*(1 - pow((x[0]/gx),0.25) - pow((x[0]/gx), 4));

	return;
}


/*
*
*
*   As described by Huband et al. in "A review of multiobjective test problems
*   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
*   Computing 10(5): 477-506, 2006.
*
*   Example WFG8.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void L2ZDT1(double *f, double *x) {

    double m =30;
  	double *y;
  	double sum;
  	int b;
  	int i, j;
  	double ff1;
  	double gx;
  	double h;

    double M[30][30]={{0.218418,-0.620254,0.843784,0.914311,-0.788548,0.428212,0.103064,-0.473730,-0.300792,-0.185507,0.330423,0.151614,0.884043,-0.272951,-0.993822,0.511197,-0.0997948,-0.659756,0.575496,0.675617,0.180332,-0.593814,-0.492722,0.0646786,-0.666503,-0.945716,-0.334582,0.611894,0.281032,0.508749},
    {-0.0265389,-0.920133,0.308861,-0.0437502,-0.374203,0.207359,-0.219433,0.914104,0.184408,0.520599,-0.885650,-0.375906,-0.708948,-0.379020,0.576578,0.0194674,-0.470262,0.572576,0.351245,-0.480477,0.238261,-0.159600,-0.827302,0.669248,0.494475,0.691715,-0.198585,0.0492812,0.959669,0.884086},
    {-0.218632,-0.865161,-0.715997,0.220772,0.692356,0.646453,-0.401724,0.615443,-0.0601957,-0.748176,-0.207987,-0.865931,0.613732,-0.525712,-0.995728,0.389633,-0.0641730,0.662131,-0.707048,-0.340423,0.606240,0.0951648,-0.160446,-0.394585,-0.167581,0.0679849,0.449799,0.733505,-0.00918638,0.00446808},
    {0.404396,0.449996,0.162711,0.294454,-0.563345,-0.114993,0.549589,-0.775141,0.677726,0.610715,0.0850755,0.0419388,-0.323614,-0.973719,-0.680238,-0.270873,-0.209617,0.968436,0.908798,0.975851,-0.994918,-0.0621977,0.628171,0.761228,0.343720,-0.792042,-0.144765,-0.965748,0.0133606,-0.0260565},
    {-0.742377,0.426391,0.408202,0.633885,-0.0351053,-0.723444,-0.577654,0.0276004,0.0712472,-0.622791,0.155451,0.442717,-0.792786,0.925785,0.670266,-0.865566,-0.638281,0.333094,0.477628,0.472610,0.231510,0.821320,-0.589803,0.796275,0.577130,0.101149,0.970191,0.532821,0.814769,-0.0687269},
    {0.712758,-0.191812,-0.390938,0.952828,0.921519,0.923094,0.930110,-0.945394,-0.0934027,0.964123,-0.795609,-0.289563,0.614236,-0.670585,0.466877,0.144597,-0.206416,0.693700,-0.967958,-0.0951247,-0.942473,-0.610767,-0.655472,-0.0960986,-0.302779,-0.734976,-0.342188,-0.315861,-0.912834,0.244990},
    {0.0969326,0.0897750,-0.241157,0.0835558,-0.420236,-0.686633,-0.711276,-0.00325377,0.435196,-0.710002,0.00283691,-0.168757,-0.134045,-0.655235,0.172361,0.998291,0.376291,-0.962215,-0.363174,-0.887770,-0.519929,-0.560554,-0.984415,0.601529,-0.984103,-0.228237,-0.578066,0.307023,0.606123,0.959635},
    {0.00225943,0.0101814,0.441456,0.0633629,0.406631,-0.0100638,-0.177972,-0.491075,0.537035,-0.924987,-0.699424,0.742285,0.0181443,0.718971,-0.0308272,0.0869310,0.524476,0.956457,0.143024,0.616481,0.217909,-0.128427,-0.262427,-0.938208,-0.524790,0.129190,0.721925,0.766492,0.470845,-0.0976466},
    {0.507807,0.804148,0.963269,0.357128,-0.832565,-0.312441,0.327779,0.184745,0.246139,-0.936814,-0.931734,-0.0327827,0.319293,0.0444730,-0.641645,0.596118,-0.293934,-0.633730,0.409658,0.759892,-0.257078,0.939616,-0.227661,0.115754,0.109640,-0.240557,0.668420,0.855535,-0.451536,0.264961},
    {-0.613660,-0.204783,-0.842476,-0.249524,-0.0985226,0.0671501,-0.527707,-0.509489,-0.883254,0.148510,-0.906465,0.496238,-0.853211,-0.779234,-0.979515,0.827175,0.228969,-0.402829,-0.970118,0.762559,0.506495,0.460303,0.897304,0.686003,0.739986,0.157310,0.281697,-0.922955,-0.780824,0.449716},
    {0.125225,0.487649,0.147046,0.679639,0.593707,-0.311828,-0.797099,-0.358150,0.958080,0.907244,0.772426,0.720574,-0.873217,0.371431,-0.826029,0.942716,0.706090,-0.658158,-0.782185,-0.806743,-0.627986,-0.405551,-0.258495,-0.796524,0.222498,0.0875450,-0.0917108,-0.625420,-0.110256,0.0417576},
    {0.244760,0.941339,-0.613783,0.402772,0.300775,-0.820314,-0.894233,-0.405896,0.0735439,0.486645,-0.394355,0.125097,-0.316386,-0.701215,-0.845742,0.206500,-0.413743,0.406725,-0.423813,-0.941255,-0.558804,0.312326,0.345314,0.319143,-0.644653,-0.0408415,0.176461,0.740113,0.470737,-0.914927},
    {-0.591523,-0.606614,-0.181873,0.692975,0.502080,-0.536704,0.359652,0.839082,0.568170,-0.0776788,-0.00332785,0.459538,-0.518313,-0.270738,-0.629958,-0.755084,-0.721573,0.431107,-0.221877,0.325430,0.163743,0.0759916,0.695064,-0.656856,0.0741630,0.264319,-0.731740,0.731548,-0.489341,0.678946},
    {0.0271269,0.804879,-0.402973,0.800373,0.760082,-0.878364,0.176801,-0.548932,-0.225601,-0.164912,-0.208143,0.776800,-0.542743,-0.156021,0.671736,0.878648,-0.419588,-0.0752896,0.0299447,-0.494459,-0.724150,0.359780,-0.326460,-0.966050,0.0127605,0.563174,-0.814853,-0.949609,-0.526794,-0.801902},
    {-0.753397,0.617418,0.689874,0.983384,0.668786,0.0304653,-0.625221,-0.133180,0.827343,-0.101358,-0.999522,-0.0525574,-0.458319,0.587409,-0.334639,0.0759643,0.0255827,0.128944,0.173170,-0.284309,0.287161,-0.550725,-0.433083,-0.242821,0.878879,0.691699,-0.660499,0.389985,0.599856,-0.711442},
    {-0.798697,-0.244945,-0.942649,0.402856,-0.494672,0.439941,-0.882160,0.170196,0.650734,-0.0982391,-0.468732,0.342133,-0.838071,-0.832362,0.658177,-0.565361,0.149473,0.693310,-0.491848,0.749160,0.526025,-0.155339,0.0998096,0.468761,0.324649,0.128488,0.544144,-0.495222,0.965229,-0.793140},
    {-0.545421,-0.500243,0.154371,0.170017,-0.259108,-0.868862,-0.507310,-0.848317,0.835712,0.616391,-0.442608,-0.158000,0.313451,0.703748,-0.755984,-0.249443,0.491564,0.985068,0.678644,0.808324,0.819750,-0.435823,-0.839855,0.00282368,-0.569165,0.0884339,-0.222144,0.499412,-0.565198,0.648240},
    {0.956914,-0.0620912,0.634479,0.928617,0.464664,0.377022,0.630470,-0.198619,-0.576153,0.565373,-0.524245,-0.187299,-0.614524,0.429316,-0.491171,0.399495,-0.333898,-0.646636,-0.0189709,-0.339605,-0.798791,0.0494081,0.367012,0.852545,0.435570,0.150039,-0.0454542,0.604861,-0.598288,-0.500696},
    {0.249008,0.370711,-0.633174,-0.0121906,0.420060,0.169373,-0.975542,-0.0297852,0.804810,0.638317,-0.670967,0.935792,-0.356050,0.175773,0.878601,-0.275168,-0.932517,-0.372497,-0.0732907,-0.185493,-0.357004,0.314786,-0.229239,0.530256,-0.513270,0.441870,0.940309,-0.240334,-0.0276121,0.743830},
    {-0.630329,-0.763376,0.625380,0.818945,0.891598,0.680494,0.471868,-0.769787,-0.878099,-0.973724,0.354362,-0.179200,-0.225034,-0.445480,0.598865,0.544005,-0.478962,0.327193,-0.525784,0.903179,-0.899248,0.156514,0.154329,0.499808,-0.836327,-0.802627,0.378082,-0.112673,-0.479260,-0.335500},
    {-0.699445,0.237731,-0.324597,-0.800406,-0.425850,-0.710739,-0.144068,-0.828545,-0.800912,0.184654,-0.636750,-0.166960,0.240427,-0.513443,0.812664,0.744943,0.970612,0.00172899,-0.726378,-0.0985012,0.224232,0.164950,0.560077,-0.813112,0.112894,-0.0955366,0.0187107,0.913887,0.123076,0.550338},
    {0.400334,-0.367816,0.198455,-0.983183,0.278976,0.714817,0.307911,0.812861,-0.403497,-0.784382,-0.161823,-0.120835,0.323172,0.583739,0.732924,-0.220603,-0.594121,0.935093,-0.216736,0.659318,-0.750417,-0.284773,-0.271496,0.491731,-0.712174,-0.763681,0.0781023,0.951666,0.734031,0.826912},
    {0.574880,-0.361951,-0.0739728,0.914380,-0.391653,0.0193887,0.412634,-0.169813,0.471794,0.660792,-0.350906,-0.612644,0.347876,0.112573,-0.501126,0.456761,-0.109004,0.289352,-0.566504,0.585042,0.584934,0.923676,0.895312,-0.161036,-0.995895,0.0853141,-0.583368,-0.157612,0.234119,0.875043},
    {0.430805,0.706102,0.423887,0.296828,-0.265607,0.338806,-0.158290,0.642516,0.355126,0.174447,-0.975015,0.869905,-0.145285,-0.484002,-0.475966,-0.677040,0.996452,-0.0685748,-0.851985,0.416498,0.791047,-0.211323,-0.302819,0.640735,-0.317908,-0.116586,-0.896382,-0.817317,-0.948837,-0.597427},
    {0.975863,-0.971371,-0.124115,0.433900,-0.254671,0.298068,-0.349803,-0.731850,0.488106,-0.0495073,0.253969,0.168116,0.148772,0.889593,-0.512213,-0.165437,0.666773,-0.976304,-0.170024,0.905794,0.473908,-0.855725,-0.0413591,-0.508661,0.443453,0.842925,-0.144503,0.936699,-0.443935,-0.182996},
    {0.803564,0.960386,-0.0323329,0.638181,-0.895684,-0.360502,0.0646258,-0.202449,-0.717228,0.970489,0.404608,-0.0861868,-0.879417,-0.866462,-0.938336,-0.799685,0.213464,-0.932344,-0.668236,0.751366,-0.227120,-0.407783,0.657463,0.0970092,-0.579109,-0.868866,-0.504041,0.926483,0.169978,-0.00842563},
    {-0.530324,0.282745,0.0255867,0.287686,0.410417,-0.766576,-0.536566,-0.628288,0.696650,0.820713,-0.506162,-0.404114,0.640099,-0.956123,-0.576586,0.435502,-0.470676,-0.367062,-0.831765,-0.294942,0.518991,0.922338,0.337886,-0.674740,-0.725667,0.916684,0.391750,0.759081,0.496979,-0.200691},
    {0.0417966,-0.687391,0.438773,0.287357,0.316636,-0.262311,-0.0755541,-0.442313,0.621378,0.670105,0.0609820,0.944162,0.643442,-0.750684,-0.639973,0.217424,0.592823,0.929094,-0.239135,-0.416280,0.570893,-0.0798988,-0.917135,-0.749545,-0.982047,0.0626998,-0.977963,0.660401,0.470569,-0.0528868},
    {-0.00138645,0.931065,-0.748519,0.304188,-0.266153,0.672524,-0.105179,-0.874749,-0.154355,-0.774656,-0.696540,0.433098,0.615897,-0.387919,-0.429779,0.650202,0.122306,-0.237727,0.626817,-0.227929,0.405916,0.483328,0.282047,-0.262206,0.784123,0.831250,-0.662272,0.702768,0.875814,-0.701221},
    {0.553793,0.471795,0.769147,0.0596680,-0.841617,-0.191179,-0.972471,-0.825361,0.779826,-0.917201,0.432720,0.103010,0.358771,0.793448,-0.0379954,-0.870112,0.600442,-0.990603,0.549151,0.512146,-0.795843,0.490091,0.372046,-0.549437,0.0964285,0.753047,-0.862840,-0.589688,0.178612,-0.720358}};

    /*
    	M = (double **)malloc(m * sizeof(double *));
        for (i=0; i<m; i++)
        {
             M[i] = (double *)malloc(m * sizeof(double));
        }
    	*/
      	y = (double *)malloc(m*sizeof(double));
        sum = 0;
        b = 0;
    	/*
        fp = fopen("./data/L2ZDT1.txt","r");
    */
    	if(m<2){
    		printf("Invalid value of m");
    	}

    	for(b=0;b<m;b++){
    		if(x[b]<0.0 || x[b]>1.0){
    			printf("Invalid input range\n");
    			return;
    		}
    	}

    /*
    	for(i=0;i<m;i++)
    	{
    		for(j=0;j<m;j++)
    		{
    			fscanf(fp, "%lf", &M[j][i]);
    		}
    	}
    */
    	for(i=0;i<m;i++)
    	{
    		sum = 0;
    		for(j=0;j<m;j++)
    		{
    			sum = sum + (M[i][j]*x[j]);
    		}
    		y[i] = sum;
    	}
    	sum = 0.0;
    	for(i=1;i<m;i++)
    	{
    		sum = sum + pow(y[i], 2);
    	}

    	ff1 = pow(y[0],2);
    	gx = 1 + 9/(m - 1) * sum;
    	h = 1 - sqrt(ff1/gx);

    	f[0] = ff1;
    	f[1] = gx*h;

    	//free(M);
    	free(y);
    	//fclose(fp);
    	return;
}


/*
*
*
*   As described by K. Deb and A. Sinha and S. Kukkonen in "Multi-objective
*   test problems, linkages, and evolutionary methodologies", GECCO'06}:
*   Proceedings of the 8th Annual Conference on Genetic and Evolutionary
*   Computation, 1141-1148, 2006.
*
*   Example T2, with linkage L2.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void L2ZDT2(double *f, double *x) {

  	double m = 30;
  	double *y;
  	double sum;
  	int b;
  	int i, j;
  	double ff1;
  	double gx;
  	double h;

    	double M[30][30]={{0.218418,-0.620254,0.843784,0.914311,-0.788548,0.428212,0.103064,-0.473730,-0.300792,-0.185507,0.330423,0.151614,0.884043,-0.272951,-0.993822,0.511197,-0.0997948,-0.659756,0.575496,0.675617,0.180332,-0.593814,-0.492722,0.0646786,-0.666503,-0.945716,-0.334582,0.611894,0.281032,0.508749},
    {-0.0265389,-0.920133,0.308861,-0.0437502,-0.374203,0.207359,-0.219433,0.914104,0.184408,0.520599,-0.885650,-0.375906,-0.708948,-0.379020,0.576578,0.0194674,-0.470262,0.572576,0.351245,-0.480477,0.238261,-0.159600,-0.827302,0.669248,0.494475,0.691715,-0.198585,0.0492812,0.959669,0.884086},
    {-0.218632,-0.865161,-0.715997,0.220772,0.692356,0.646453,-0.401724,0.615443,-0.0601957,-0.748176,-0.207987,-0.865931,0.613732,-0.525712,-0.995728,0.389633,-0.0641730,0.662131,-0.707048,-0.340423,0.606240,0.0951648,-0.160446,-0.394585,-0.167581,0.0679849,0.449799,0.733505,-0.00918638,0.00446808},
    {0.404396,0.449996,0.162711,0.294454,-0.563345,-0.114993,0.549589,-0.775141,0.677726,0.610715,0.0850755,0.0419388,-0.323614,-0.973719,-0.680238,-0.270873,-0.209617,0.968436,0.908798,0.975851,-0.994918,-0.0621977,0.628171,0.761228,0.343720,-0.792042,-0.144765,-0.965748,0.0133606,-0.0260565},
    {-0.742377,0.426391,0.408202,0.633885,-0.0351053,-0.723444,-0.577654,0.0276004,0.0712472,-0.622791,0.155451,0.442717,-0.792786,0.925785,0.670266,-0.865566,-0.638281,0.333094,0.477628,0.472610,0.231510,0.821320,-0.589803,0.796275,0.577130,0.101149,0.970191,0.532821,0.814769,-0.0687269},
    {0.712758,-0.191812,-0.390938,0.952828,0.921519,0.923094,0.930110,-0.945394,-0.0934027,0.964123,-0.795609,-0.289563,0.614236,-0.670585,0.466877,0.144597,-0.206416,0.693700,-0.967958,-0.0951247,-0.942473,-0.610767,-0.655472,-0.0960986,-0.302779,-0.734976,-0.342188,-0.315861,-0.912834,0.244990},
    {0.0969326,0.0897750,-0.241157,0.0835558,-0.420236,-0.686633,-0.711276,-0.00325377,0.435196,-0.710002,0.00283691,-0.168757,-0.134045,-0.655235,0.172361,0.998291,0.376291,-0.962215,-0.363174,-0.887770,-0.519929,-0.560554,-0.984415,0.601529,-0.984103,-0.228237,-0.578066,0.307023,0.606123,0.959635},
    {0.00225943,0.0101814,0.441456,0.0633629,0.406631,-0.0100638,-0.177972,-0.491075,0.537035,-0.924987,-0.699424,0.742285,0.0181443,0.718971,-0.0308272,0.0869310,0.524476,0.956457,0.143024,0.616481,0.217909,-0.128427,-0.262427,-0.938208,-0.524790,0.129190,0.721925,0.766492,0.470845,-0.0976466},
    {0.507807,0.804148,0.963269,0.357128,-0.832565,-0.312441,0.327779,0.184745,0.246139,-0.936814,-0.931734,-0.0327827,0.319293,0.0444730,-0.641645,0.596118,-0.293934,-0.633730,0.409658,0.759892,-0.257078,0.939616,-0.227661,0.115754,0.109640,-0.240557,0.668420,0.855535,-0.451536,0.264961},
    {-0.613660,-0.204783,-0.842476,-0.249524,-0.0985226,0.0671501,-0.527707,-0.509489,-0.883254,0.148510,-0.906465,0.496238,-0.853211,-0.779234,-0.979515,0.827175,0.228969,-0.402829,-0.970118,0.762559,0.506495,0.460303,0.897304,0.686003,0.739986,0.157310,0.281697,-0.922955,-0.780824,0.449716},
    {0.125225,0.487649,0.147046,0.679639,0.593707,-0.311828,-0.797099,-0.358150,0.958080,0.907244,0.772426,0.720574,-0.873217,0.371431,-0.826029,0.942716,0.706090,-0.658158,-0.782185,-0.806743,-0.627986,-0.405551,-0.258495,-0.796524,0.222498,0.0875450,-0.0917108,-0.625420,-0.110256,0.0417576},
    {0.244760,0.941339,-0.613783,0.402772,0.300775,-0.820314,-0.894233,-0.405896,0.0735439,0.486645,-0.394355,0.125097,-0.316386,-0.701215,-0.845742,0.206500,-0.413743,0.406725,-0.423813,-0.941255,-0.558804,0.312326,0.345314,0.319143,-0.644653,-0.0408415,0.176461,0.740113,0.470737,-0.914927},
    {-0.591523,-0.606614,-0.181873,0.692975,0.502080,-0.536704,0.359652,0.839082,0.568170,-0.0776788,-0.00332785,0.459538,-0.518313,-0.270738,-0.629958,-0.755084,-0.721573,0.431107,-0.221877,0.325430,0.163743,0.0759916,0.695064,-0.656856,0.0741630,0.264319,-0.731740,0.731548,-0.489341,0.678946},
    {0.0271269,0.804879,-0.402973,0.800373,0.760082,-0.878364,0.176801,-0.548932,-0.225601,-0.164912,-0.208143,0.776800,-0.542743,-0.156021,0.671736,0.878648,-0.419588,-0.0752896,0.0299447,-0.494459,-0.724150,0.359780,-0.326460,-0.966050,0.0127605,0.563174,-0.814853,-0.949609,-0.526794,-0.801902},
    {-0.753397,0.617418,0.689874,0.983384,0.668786,0.0304653,-0.625221,-0.133180,0.827343,-0.101358,-0.999522,-0.0525574,-0.458319,0.587409,-0.334639,0.0759643,0.0255827,0.128944,0.173170,-0.284309,0.287161,-0.550725,-0.433083,-0.242821,0.878879,0.691699,-0.660499,0.389985,0.599856,-0.711442},
    {-0.798697,-0.244945,-0.942649,0.402856,-0.494672,0.439941,-0.882160,0.170196,0.650734,-0.0982391,-0.468732,0.342133,-0.838071,-0.832362,0.658177,-0.565361,0.149473,0.693310,-0.491848,0.749160,0.526025,-0.155339,0.0998096,0.468761,0.324649,0.128488,0.544144,-0.495222,0.965229,-0.793140},
    {-0.545421,-0.500243,0.154371,0.170017,-0.259108,-0.868862,-0.507310,-0.848317,0.835712,0.616391,-0.442608,-0.158000,0.313451,0.703748,-0.755984,-0.249443,0.491564,0.985068,0.678644,0.808324,0.819750,-0.435823,-0.839855,0.00282368,-0.569165,0.0884339,-0.222144,0.499412,-0.565198,0.648240},
    {0.956914,-0.0620912,0.634479,0.928617,0.464664,0.377022,0.630470,-0.198619,-0.576153,0.565373,-0.524245,-0.187299,-0.614524,0.429316,-0.491171,0.399495,-0.333898,-0.646636,-0.0189709,-0.339605,-0.798791,0.0494081,0.367012,0.852545,0.435570,0.150039,-0.0454542,0.604861,-0.598288,-0.500696},
    {0.249008,0.370711,-0.633174,-0.0121906,0.420060,0.169373,-0.975542,-0.0297852,0.804810,0.638317,-0.670967,0.935792,-0.356050,0.175773,0.878601,-0.275168,-0.932517,-0.372497,-0.0732907,-0.185493,-0.357004,0.314786,-0.229239,0.530256,-0.513270,0.441870,0.940309,-0.240334,-0.0276121,0.743830},
    {-0.630329,-0.763376,0.625380,0.818945,0.891598,0.680494,0.471868,-0.769787,-0.878099,-0.973724,0.354362,-0.179200,-0.225034,-0.445480,0.598865,0.544005,-0.478962,0.327193,-0.525784,0.903179,-0.899248,0.156514,0.154329,0.499808,-0.836327,-0.802627,0.378082,-0.112673,-0.479260,-0.335500},
    {-0.699445,0.237731,-0.324597,-0.800406,-0.425850,-0.710739,-0.144068,-0.828545,-0.800912,0.184654,-0.636750,-0.166960,0.240427,-0.513443,0.812664,0.744943,0.970612,0.00172899,-0.726378,-0.0985012,0.224232,0.164950,0.560077,-0.813112,0.112894,-0.0955366,0.0187107,0.913887,0.123076,0.550338},
    {0.400334,-0.367816,0.198455,-0.983183,0.278976,0.714817,0.307911,0.812861,-0.403497,-0.784382,-0.161823,-0.120835,0.323172,0.583739,0.732924,-0.220603,-0.594121,0.935093,-0.216736,0.659318,-0.750417,-0.284773,-0.271496,0.491731,-0.712174,-0.763681,0.0781023,0.951666,0.734031,0.826912},
    {0.574880,-0.361951,-0.0739728,0.914380,-0.391653,0.0193887,0.412634,-0.169813,0.471794,0.660792,-0.350906,-0.612644,0.347876,0.112573,-0.501126,0.456761,-0.109004,0.289352,-0.566504,0.585042,0.584934,0.923676,0.895312,-0.161036,-0.995895,0.0853141,-0.583368,-0.157612,0.234119,0.875043},
    {0.430805,0.706102,0.423887,0.296828,-0.265607,0.338806,-0.158290,0.642516,0.355126,0.174447,-0.975015,0.869905,-0.145285,-0.484002,-0.475966,-0.677040,0.996452,-0.0685748,-0.851985,0.416498,0.791047,-0.211323,-0.302819,0.640735,-0.317908,-0.116586,-0.896382,-0.817317,-0.948837,-0.597427},
    {0.975863,-0.971371,-0.124115,0.433900,-0.254671,0.298068,-0.349803,-0.731850,0.488106,-0.0495073,0.253969,0.168116,0.148772,0.889593,-0.512213,-0.165437,0.666773,-0.976304,-0.170024,0.905794,0.473908,-0.855725,-0.0413591,-0.508661,0.443453,0.842925,-0.144503,0.936699,-0.443935,-0.182996},
    {0.803564,0.960386,-0.0323329,0.638181,-0.895684,-0.360502,0.0646258,-0.202449,-0.717228,0.970489,0.404608,-0.0861868,-0.879417,-0.866462,-0.938336,-0.799685,0.213464,-0.932344,-0.668236,0.751366,-0.227120,-0.407783,0.657463,0.0970092,-0.579109,-0.868866,-0.504041,0.926483,0.169978,-0.00842563},
    {-0.530324,0.282745,0.0255867,0.287686,0.410417,-0.766576,-0.536566,-0.628288,0.696650,0.820713,-0.506162,-0.404114,0.640099,-0.956123,-0.576586,0.435502,-0.470676,-0.367062,-0.831765,-0.294942,0.518991,0.922338,0.337886,-0.674740,-0.725667,0.916684,0.391750,0.759081,0.496979,-0.200691},
    {0.0417966,-0.687391,0.438773,0.287357,0.316636,-0.262311,-0.0755541,-0.442313,0.621378,0.670105,0.0609820,0.944162,0.643442,-0.750684,-0.639973,0.217424,0.592823,0.929094,-0.239135,-0.416280,0.570893,-0.0798988,-0.917135,-0.749545,-0.982047,0.0626998,-0.977963,0.660401,0.470569,-0.0528868},
    {-0.00138645,0.931065,-0.748519,0.304188,-0.266153,0.672524,-0.105179,-0.874749,-0.154355,-0.774656,-0.696540,0.433098,0.615897,-0.387919,-0.429779,0.650202,0.122306,-0.237727,0.626817,-0.227929,0.405916,0.483328,0.282047,-0.262206,0.784123,0.831250,-0.662272,0.702768,0.875814,-0.701221},
    {0.553793,0.471795,0.769147,0.0596680,-0.841617,-0.191179,-0.972471,-0.825361,0.779826,-0.917201,0.432720,0.103010,0.358771,0.793448,-0.0379954,-0.870112,0.600442,-0.990603,0.549151,0.512146,-0.795843,0.490091,0.372046,-0.549437,0.0964285,0.753047,-0.862840,-0.589688,0.178612,-0.720358}};

      // for (i=0; i<m; i++)
      // {
      //      M[i] = (double *)malloc(m * sizeof(double));
      // }
    	y = (double *)malloc(m*sizeof(double));
      sum = 0;
      b = 0;
      // fp = fopen("./data/L2ZDT1.txt","r");

  	if(m<2){
  		printf("Invalid value of m");
  	}

  	for(b=0;b<m;b++){
  		if(x[b]<0.0 || x[b]>1.0){
  			printf("Invalid input range\n");
  			return;
  		}
  	}


  	// for(i=0;i<m;i++)
  	// {
  	// 	for(j=0;j<m;j++)
  	// 	{
  	// 		fscanf(fp, "%lf", &M[j][i]);
  	// 	}
  	// }

  	for(i=0;i<m;i++)
  	{
  		sum = 0;
  		for(j=0;j<m;j++)
  		{
  			sum = sum + (M[i][j]*x[j]);
  		}
  		y[i] = sum;
  	}
  	sum = 0.0;
  	for(i=1;i<m;i++)
  	{
  		sum = sum + pow(y[i], 2);
  	}

  	ff1 = pow(y[0],2);
  	gx = 1 + 9/(m - 1) * sum;
  	h = 1 - pow((ff1/gx),2);

  	f[0] = ff1;
  	f[1] = gx*h;

    free(y);
    return;
}


/*
*
*
*   As described by K. Deb and A. Sinha and S. Kukkonen in "Multi-objective
*   test problems, linkages, and evolutionary methodologies", GECCO'06}:
*   Proceedings of the 8th Annual Conference on Genetic and Evolutionary
*   Computation, 1141-1148, 2006.
*
*   Example T4, with linkage L2.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void L2ZDT4(double *f, double *x) {

  	double m = 30;
  	double *y;
  	double sum;
  	double temp;
  	int b;
  	int i, j;
  	double ff1;
  	double gx;
  	double h;
  	double pi;
    	double M[30][30]={{0.218418,-0.620254,0.843784,0.914311,-0.788548,0.428212,0.103064,-0.473730,-0.300792,-0.185507,0.330423,0.151614,0.884043,-0.272951,-0.993822,0.511197,-0.0997948,-0.659756,0.575496,0.675617,0.180332,-0.593814,-0.492722,0.0646786,-0.666503,-0.945716,-0.334582,0.611894,0.281032,0.508749},
    {-0.0265389,-0.920133,0.308861,-0.0437502,-0.374203,0.207359,-0.219433,0.914104,0.184408,0.520599,-0.885650,-0.375906,-0.708948,-0.379020,0.576578,0.0194674,-0.470262,0.572576,0.351245,-0.480477,0.238261,-0.159600,-0.827302,0.669248,0.494475,0.691715,-0.198585,0.0492812,0.959669,0.884086},
    {-0.218632,-0.865161,-0.715997,0.220772,0.692356,0.646453,-0.401724,0.615443,-0.0601957,-0.748176,-0.207987,-0.865931,0.613732,-0.525712,-0.995728,0.389633,-0.0641730,0.662131,-0.707048,-0.340423,0.606240,0.0951648,-0.160446,-0.394585,-0.167581,0.0679849,0.449799,0.733505,-0.00918638,0.00446808},
    {0.404396,0.449996,0.162711,0.294454,-0.563345,-0.114993,0.549589,-0.775141,0.677726,0.610715,0.0850755,0.0419388,-0.323614,-0.973719,-0.680238,-0.270873,-0.209617,0.968436,0.908798,0.975851,-0.994918,-0.0621977,0.628171,0.761228,0.343720,-0.792042,-0.144765,-0.965748,0.0133606,-0.0260565},
    {-0.742377,0.426391,0.408202,0.633885,-0.0351053,-0.723444,-0.577654,0.0276004,0.0712472,-0.622791,0.155451,0.442717,-0.792786,0.925785,0.670266,-0.865566,-0.638281,0.333094,0.477628,0.472610,0.231510,0.821320,-0.589803,0.796275,0.577130,0.101149,0.970191,0.532821,0.814769,-0.0687269},
    {0.712758,-0.191812,-0.390938,0.952828,0.921519,0.923094,0.930110,-0.945394,-0.0934027,0.964123,-0.795609,-0.289563,0.614236,-0.670585,0.466877,0.144597,-0.206416,0.693700,-0.967958,-0.0951247,-0.942473,-0.610767,-0.655472,-0.0960986,-0.302779,-0.734976,-0.342188,-0.315861,-0.912834,0.244990},
    {0.0969326,0.0897750,-0.241157,0.0835558,-0.420236,-0.686633,-0.711276,-0.00325377,0.435196,-0.710002,0.00283691,-0.168757,-0.134045,-0.655235,0.172361,0.998291,0.376291,-0.962215,-0.363174,-0.887770,-0.519929,-0.560554,-0.984415,0.601529,-0.984103,-0.228237,-0.578066,0.307023,0.606123,0.959635},
    {0.00225943,0.0101814,0.441456,0.0633629,0.406631,-0.0100638,-0.177972,-0.491075,0.537035,-0.924987,-0.699424,0.742285,0.0181443,0.718971,-0.0308272,0.0869310,0.524476,0.956457,0.143024,0.616481,0.217909,-0.128427,-0.262427,-0.938208,-0.524790,0.129190,0.721925,0.766492,0.470845,-0.0976466},
    {0.507807,0.804148,0.963269,0.357128,-0.832565,-0.312441,0.327779,0.184745,0.246139,-0.936814,-0.931734,-0.0327827,0.319293,0.0444730,-0.641645,0.596118,-0.293934,-0.633730,0.409658,0.759892,-0.257078,0.939616,-0.227661,0.115754,0.109640,-0.240557,0.668420,0.855535,-0.451536,0.264961},
    {-0.613660,-0.204783,-0.842476,-0.249524,-0.0985226,0.0671501,-0.527707,-0.509489,-0.883254,0.148510,-0.906465,0.496238,-0.853211,-0.779234,-0.979515,0.827175,0.228969,-0.402829,-0.970118,0.762559,0.506495,0.460303,0.897304,0.686003,0.739986,0.157310,0.281697,-0.922955,-0.780824,0.449716},
    {0.125225,0.487649,0.147046,0.679639,0.593707,-0.311828,-0.797099,-0.358150,0.958080,0.907244,0.772426,0.720574,-0.873217,0.371431,-0.826029,0.942716,0.706090,-0.658158,-0.782185,-0.806743,-0.627986,-0.405551,-0.258495,-0.796524,0.222498,0.0875450,-0.0917108,-0.625420,-0.110256,0.0417576},
    {0.244760,0.941339,-0.613783,0.402772,0.300775,-0.820314,-0.894233,-0.405896,0.0735439,0.486645,-0.394355,0.125097,-0.316386,-0.701215,-0.845742,0.206500,-0.413743,0.406725,-0.423813,-0.941255,-0.558804,0.312326,0.345314,0.319143,-0.644653,-0.0408415,0.176461,0.740113,0.470737,-0.914927},
    {-0.591523,-0.606614,-0.181873,0.692975,0.502080,-0.536704,0.359652,0.839082,0.568170,-0.0776788,-0.00332785,0.459538,-0.518313,-0.270738,-0.629958,-0.755084,-0.721573,0.431107,-0.221877,0.325430,0.163743,0.0759916,0.695064,-0.656856,0.0741630,0.264319,-0.731740,0.731548,-0.489341,0.678946},
    {0.0271269,0.804879,-0.402973,0.800373,0.760082,-0.878364,0.176801,-0.548932,-0.225601,-0.164912,-0.208143,0.776800,-0.542743,-0.156021,0.671736,0.878648,-0.419588,-0.0752896,0.0299447,-0.494459,-0.724150,0.359780,-0.326460,-0.966050,0.0127605,0.563174,-0.814853,-0.949609,-0.526794,-0.801902},
    {-0.753397,0.617418,0.689874,0.983384,0.668786,0.0304653,-0.625221,-0.133180,0.827343,-0.101358,-0.999522,-0.0525574,-0.458319,0.587409,-0.334639,0.0759643,0.0255827,0.128944,0.173170,-0.284309,0.287161,-0.550725,-0.433083,-0.242821,0.878879,0.691699,-0.660499,0.389985,0.599856,-0.711442},
    {-0.798697,-0.244945,-0.942649,0.402856,-0.494672,0.439941,-0.882160,0.170196,0.650734,-0.0982391,-0.468732,0.342133,-0.838071,-0.832362,0.658177,-0.565361,0.149473,0.693310,-0.491848,0.749160,0.526025,-0.155339,0.0998096,0.468761,0.324649,0.128488,0.544144,-0.495222,0.965229,-0.793140},
    {-0.545421,-0.500243,0.154371,0.170017,-0.259108,-0.868862,-0.507310,-0.848317,0.835712,0.616391,-0.442608,-0.158000,0.313451,0.703748,-0.755984,-0.249443,0.491564,0.985068,0.678644,0.808324,0.819750,-0.435823,-0.839855,0.00282368,-0.569165,0.0884339,-0.222144,0.499412,-0.565198,0.648240},
    {0.956914,-0.0620912,0.634479,0.928617,0.464664,0.377022,0.630470,-0.198619,-0.576153,0.565373,-0.524245,-0.187299,-0.614524,0.429316,-0.491171,0.399495,-0.333898,-0.646636,-0.0189709,-0.339605,-0.798791,0.0494081,0.367012,0.852545,0.435570,0.150039,-0.0454542,0.604861,-0.598288,-0.500696},
    {0.249008,0.370711,-0.633174,-0.0121906,0.420060,0.169373,-0.975542,-0.0297852,0.804810,0.638317,-0.670967,0.935792,-0.356050,0.175773,0.878601,-0.275168,-0.932517,-0.372497,-0.0732907,-0.185493,-0.357004,0.314786,-0.229239,0.530256,-0.513270,0.441870,0.940309,-0.240334,-0.0276121,0.743830},
    {-0.630329,-0.763376,0.625380,0.818945,0.891598,0.680494,0.471868,-0.769787,-0.878099,-0.973724,0.354362,-0.179200,-0.225034,-0.445480,0.598865,0.544005,-0.478962,0.327193,-0.525784,0.903179,-0.899248,0.156514,0.154329,0.499808,-0.836327,-0.802627,0.378082,-0.112673,-0.479260,-0.335500},
    {-0.699445,0.237731,-0.324597,-0.800406,-0.425850,-0.710739,-0.144068,-0.828545,-0.800912,0.184654,-0.636750,-0.166960,0.240427,-0.513443,0.812664,0.744943,0.970612,0.00172899,-0.726378,-0.0985012,0.224232,0.164950,0.560077,-0.813112,0.112894,-0.0955366,0.0187107,0.913887,0.123076,0.550338},
    {0.400334,-0.367816,0.198455,-0.983183,0.278976,0.714817,0.307911,0.812861,-0.403497,-0.784382,-0.161823,-0.120835,0.323172,0.583739,0.732924,-0.220603,-0.594121,0.935093,-0.216736,0.659318,-0.750417,-0.284773,-0.271496,0.491731,-0.712174,-0.763681,0.0781023,0.951666,0.734031,0.826912},
    {0.574880,-0.361951,-0.0739728,0.914380,-0.391653,0.0193887,0.412634,-0.169813,0.471794,0.660792,-0.350906,-0.612644,0.347876,0.112573,-0.501126,0.456761,-0.109004,0.289352,-0.566504,0.585042,0.584934,0.923676,0.895312,-0.161036,-0.995895,0.0853141,-0.583368,-0.157612,0.234119,0.875043},
    {0.430805,0.706102,0.423887,0.296828,-0.265607,0.338806,-0.158290,0.642516,0.355126,0.174447,-0.975015,0.869905,-0.145285,-0.484002,-0.475966,-0.677040,0.996452,-0.0685748,-0.851985,0.416498,0.791047,-0.211323,-0.302819,0.640735,-0.317908,-0.116586,-0.896382,-0.817317,-0.948837,-0.597427},
    {0.975863,-0.971371,-0.124115,0.433900,-0.254671,0.298068,-0.349803,-0.731850,0.488106,-0.0495073,0.253969,0.168116,0.148772,0.889593,-0.512213,-0.165437,0.666773,-0.976304,-0.170024,0.905794,0.473908,-0.855725,-0.0413591,-0.508661,0.443453,0.842925,-0.144503,0.936699,-0.443935,-0.182996},
    {0.803564,0.960386,-0.0323329,0.638181,-0.895684,-0.360502,0.0646258,-0.202449,-0.717228,0.970489,0.404608,-0.0861868,-0.879417,-0.866462,-0.938336,-0.799685,0.213464,-0.932344,-0.668236,0.751366,-0.227120,-0.407783,0.657463,0.0970092,-0.579109,-0.868866,-0.504041,0.926483,0.169978,-0.00842563},
    {-0.530324,0.282745,0.0255867,0.287686,0.410417,-0.766576,-0.536566,-0.628288,0.696650,0.820713,-0.506162,-0.404114,0.640099,-0.956123,-0.576586,0.435502,-0.470676,-0.367062,-0.831765,-0.294942,0.518991,0.922338,0.337886,-0.674740,-0.725667,0.916684,0.391750,0.759081,0.496979,-0.200691},
    {0.0417966,-0.687391,0.438773,0.287357,0.316636,-0.262311,-0.0755541,-0.442313,0.621378,0.670105,0.0609820,0.944162,0.643442,-0.750684,-0.639973,0.217424,0.592823,0.929094,-0.239135,-0.416280,0.570893,-0.0798988,-0.917135,-0.749545,-0.982047,0.0626998,-0.977963,0.660401,0.470569,-0.0528868},
    {-0.00138645,0.931065,-0.748519,0.304188,-0.266153,0.672524,-0.105179,-0.874749,-0.154355,-0.774656,-0.696540,0.433098,0.615897,-0.387919,-0.429779,0.650202,0.122306,-0.237727,0.626817,-0.227929,0.405916,0.483328,0.282047,-0.262206,0.784123,0.831250,-0.662272,0.702768,0.875814,-0.701221},
    {0.553793,0.471795,0.769147,0.0596680,-0.841617,-0.191179,-0.972471,-0.825361,0.779826,-0.917201,0.432720,0.103010,0.358771,0.793448,-0.0379954,-0.870112,0.600442,-0.990603,0.549151,0.512146,-0.795843,0.490091,0.372046,-0.549437,0.0964285,0.753047,-0.862840,-0.589688,0.178612,-0.720358}};

      // for (i=0; i<m; i++)
      // {
      //      M[i] = (double *)malloc(m * sizeof(double));
      // }
    	y = (double *)malloc(m*sizeof(double));
      sum = 0;
      b = 0;
      pi = 4 * atan(1);
      // fp = fopen("./data/L2ZDT1.txt","r");

  	if(m<2){
  		printf("Invalid value of m");
  	}

  	for(b=0;b<m;b++){
  		if(x[b]<0.0 || x[b]>1.0){
  			printf("Invalid input range\n");
  			return;
  		}
  	}


  	// for(i=0;i<m;i++)
  	// {
  	// 	for(j=0;j<m;j++)
  	// 	{
  	// 		fscanf(fp, "%lf", &M[j][i]);
  	// 	}
  	// }

  	for(i=0;i<m;i++)
  	{
  		sum = 0;
  		for(j=0;j<m;j++)
  		{
  			sum = sum + (M[i][j]*x[j]);
  		}
  		y[i] = sum;
  	}
  	sum = 0.0;
  	for(i=1;i<m;i++)
  	{
  		temp = pow(y[i], 2) - 10*cos(4*pi*y[i]);
  		sum = sum + temp;
  	}

  	ff1 = pow(y[0],2);
  	gx = 1 + 10*(m - 1) + sum;
  	h = 1 - sqrt(ff1/gx);

  	f[0] = ff1;
  	f[1] = gx*h;

    free(y);
  	return;
}


/*
*
*
*   As described by K. Deb and A. Sinha and S. Kukkonen in "Multi-objective
*   test problems, linkages, and evolutionary methodologies", GECCO'06}:
*   Proceedings of the 8th Annual Conference on Genetic and Evolutionary
*   Computation, 1141-1148, 2006.
*
*   Example T6, with linkage L2.
*
*   In the above paper the number of variables was set to 30.
*   We selected n=10 according to the dimension of problem ZDT6.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void L2ZDT6(double *f, double *x) {

    double m =10;
  	double *y;
  	double sum;
  	int b;
  	int i, j;
  	double ff1;
  	double gx;
  	double h;
    	//FILE * fp;
  	/* depends on m*/
  	double M[10][10] = {{0.218418,-0.620254,0.843784,0.914311,-0.788548,0.428212,0.103064,-0.473730,-0.300792,-0.185507},
  {0.330423,0.151614,0.884043,-0.272951,-0.993822,0.511197,-0.0997948,-0.659756,0.575496,0.675617},
  {0.180332,-0.593814,-0.492722,0.0646786,-0.666503,-0.945716,-0.334582,0.611894,0.281032,0.508749},
  {-0.0265389,-0.920133,0.308861,-0.0437502,-0.374203,0.207359,-0.219433,0.914104,0.184408,0.520599},
  {-0.885650,-0.375906,-0.708948,-0.379020,0.576578,0.0194674,-0.470262,0.572576,0.351245,-0.480477},
  {0.238261,-0.159600,-0.827302,0.669248,0.494475,0.691715,-0.198585,0.0492812,0.959669,0.884086},
  {-0.218632,-0.865161,-0.715997,0.220772,0.692356,0.646453,-0.401724,0.615443,-0.0601957,-0.748176},
  {-0.207987,-0.865931,0.613732,-0.525712,-0.995728,0.389633,-0.0641730,0.662131,-0.707048,-0.340423},
  {0.606240,0.0951648,-0.160446,-0.394585,-0.167581,0.0679849,0.449799,0.733505,-0.00918638,0.00446808},
  {0.404396,0.449996,0.162711,0.294454,-0.563345,-0.114993,0.549589,-0.775141,0.677726,0.610715}};
  /*
  	M = (double **)malloc(m * sizeof(double *));
      for (i=0; i<m; i++)
      {
           M[i] = (double *)malloc(m * sizeof(double));
      }
  	*/
    	y = (double *)malloc(m*sizeof(double));
      sum = 0;
      b = 0;
  	/*
      fp = fopen("./data/L2ZDT6.txt","r");
  	*/
  	if(m<2){
  		printf("Invalid value of m");
  	}

  	for(b=0;b<m;b++){
  		if(x[b]<0.0 || x[b]>1.0){
  			printf("Invalid input range\n");
  			return;
  		}
  	}

  	/*
  	for(i=0;i<m;i++)
  	{
  		for(j=0;j<m;j++)
  		{
  			fscanf(fp, "%lf", &M[j][i]);
  		}
  	}
  	*/
  	for(i=0;i<m;i++)
  	{
  		sum = 0;
  		for(j=0;j<m;j++)
  		{
  			sum = sum + (M[i][j]*x[j]);
  		}
  		y[i] = sum;
  	}
  	sum = 0.0;
  	for(i=1;i<m;i++)
  	{
  		sum = sum + pow(y[i], 2);
  	}

  	ff1 = pow(y[0],2);
  	gx = 1 + 9*pow((sum/(m - 1)), 0.25);
  	h = 1 - pow((ff1/gx),2);

  	f[0] = ff1;
  	f[1] = gx*h;

  	//free(M);
  	free(y);
  	//fclose(fp);

  	return;
}


/*
*
*
*   As described by A. Lovison in "A synthetic approach to multiobjective
*   optimization", arxiv Item: http://arxiv.org/abs/1002.0093.
*
*   Example 1.
*
*   In the above paper/papers the variables bounds were not set.
*   We considered 0<=x[i]<=3, i=1,2.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void lovison1(double *f, double *x) {

	if(x[0] < 0 || x[0] > 3 || x[1] < 0 || x[1] > 3 ){
		printf("Invalid input range\n");
		return;
	}

	f[0] = -1*(-1.05*pow((x[0]), 2) - 0.98*(pow((x[1]), 2)));
	f[1] = -1*(-0.99*pow((x[0] - 3), 2) - 1.03*pow((x[1] - 2.5), 2));


	return;
}


/*
*
*
*   As described by A. Lovison in "A synthetic approach to multiobjective
*   optimization", arxiv Item: http://arxiv.org/abs/1002.0093.
*
*   Example 3.
*
*   In the above paper/papers the variables bounds were not set.
*   We considered 0<=x[1]<=6 and -4<=x[2]<=4.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void lovison3(double *f, double *x) {

	if(x[0] < 0 || x[0] > 6 || x[1] < -4 || x[1] > 4 ){
		printf("Invalid input range\n");
		return;
	}

	f[0] = -1*((-1*pow((x[0]), 2)) - (pow((x[1]), 2)));
	f[1] = -1*(-1*pow((x[0] - 6), 2) + pow((x[1] + 0.3), 2));


	return;
}


/*
*
*
*   As described by A. Lovison in "A synthetic approach to multiobjective
*   optimization", arxiv Item: http://arxiv.org/abs/1002.0093.
*
*   Example 5.
*
*   In the above paper/papers the variables bounds were not set.
*   We considered -1<=x[i]<=4, i=1,2,3.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void lovison5(double *f, double *x) {

  	//parameters
	 int n = 3;
	 int m = 3;
	//double **C;
	//double **alpha;
	//double *beta;
	//double *gamma ;
	double pi = 4 * atan(1);
	int i,j;
	// temp;
	double *func;
	//FILE * fp;
	double sum;
	double temp1;

    double C[3][3] = {{0.218418,-0.620254,0.843784},
    {0.914311,-0.788548,0.428212},
    {0.103064,-0.47373,-0.300792}};

    double alpha[3][3] = {{0.407247,0.665212,0.575807},
    {0.942022,0.363525,0.00308876},
    {0.755598,0.450103,0.170122}};

    double beta[3] = {0.575496,0.675617,0.180332};
    double gamma[3] = {-0.593814,-0.492722,0.0646786};

	sum = 0;
	//beta = (double *)malloc((m)*sizeof(double));
	//gamma = (double *)malloc((m)*sizeof(double));
	//C = (double **)malloc(n * sizeof(double *));
	func = (double *)malloc(m*sizeof(double));

    /*
    for (i=0; i<n; i++)
    {
         C[i] = (double *)malloc(m * sizeof(double));
    }

    alpha = (double **)malloc(m * sizeof(double *));
    for (i=0; i<m; i++)
    {
         alpha[i] = (double *)malloc(n * sizeof(double));
    }

	fp = fopen("./data/lovison5.txt","r");

	for(i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
		{
			fscanf(fp, "%lf", &C[i][j]);
		}
	}
	fscanf(fp,"%c", &temp);


	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			fscanf(fp,"%lf", &alpha[i][j]);
		}
	}
	fscanf(fp,"%c", &temp);

	for(i=0;i<m;i++)
	{
		fscanf(fp,"%lf",&beta[i]);
	}
	fscanf(fp,"%c", &temp);

	for(i=0;i<m;i++)
	{
		fscanf(fp,"%lf",&gamma[i]);
	}

    */


	for(i=0;i<=n-1;i++){
		if(x[i]<-1 || x[i]>4){
			printf("Invalid input range\n");
			return;
		}
	}

	for(j = 0;j<=m-1;j++)
	{
		sum = 0;
		for(i=0;i<=n-1;i++)
		{
			temp1 = pow((x[i] - C[i][j]), 2);
			sum = sum + (temp1*alpha[j][i]*(-1));
		}
		func[j] = sum;
	}

	f[0] = (-1)*func[0];
	f[1] = (-1)*(func[1] + beta[1]*sin(pi*(x[0] + x[1])/gamma[1]));
	f[2] = (-1)*(func[2] + beta[2]*cos(pi*(x[0] - x[1])/gamma[2]));
	//free(C);
	//free(alpha);
	//free(beta);
	//free(gamma);
	free(func);
	//fclose(fp);
	return;
}


/*
*
*
*   As described by T. Okabe, Y. Jin, M. Olhofer, and B. Sendhoff. "On test
*   functions for evolutionary multi-objective optimization.", Parallel
*   Problem Solving from Nature, VIII, LNCS 3242, Springer, pp.792-802,
*   September 2004.
*
*   Test function OKA1.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void OKA1(double *f, double *x) {
	double pi;
	double y[2];

	double b1_a;
	double b1_b;
	double b2_a;
	double b2_b;

	pi = 4 * atan(1);
	b1_a = 6*sin(pi/12);
	b1_b = 6*sin(pi/12) + 2*pi*cos(pi/12);
	b2_a = -2*pi*sin(pi/12);
	b2_b = 6*cos(pi/12);

	if(x[0]<b1_a || x[0]>b1_b || x[1]<b2_a || x[1]>b2_b){
		printf("Invalid input range\n");
        printf("L1:%1.15f\n",b1_a);
		printf("U1:%1.15f\n",b1_b);
        printf("L2:%1.15f\n",b2_a);
		printf("U2:%1.15f\n",b2_b);
		printf("X1:%1.15f\n",x[0]);
		printf("X2:%1.15f\n",x[1]);
		return;
	}

	y[0] = cos(pi/12)*x[0] - sin(pi/12)*x[1];
	y[1] = sin(pi/12)*x[0] + cos(pi/12)*x[1];

	f[0] = y[0];
	f[1] = sqrt(2*pi) - sqrt(fabs(y[0])) + 2*pow(fabs(y[1] - 3*cos(y[0]) - 3),(1.0/3.0));

	return;
}


/*
*
*
*   As described by T. Okabe, Y. Jin, M. Olhofer, and B. Sendhoff. "On test
*   functions for evolutionary multi-objective optimization.", Parallel
*   Problem Solving from Nature, VIII, LNCS 3242, Springer, pp.792-802,
*   September 2004.
*
*   Test function OKA2.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void OKA2(double *f, double *x) {
	double pi = 4 * atan(1);
	double b1_a = -1*pi;
	double b1_b = pi;
	double b2_a = -5;
	double b2_b = 5;

	if(x[0]<b1_a || x[0]>b1_b || x[1]<b2_a || x[1]>b2_b || x[2]<b2_a || x[2]>b2_b){
		printf("Invalid input range\n");
		return;
	}

	f[1] = x[0];
	f[0] = 1 - pow((x[0] + pi), 2)/(4*pow(pi, 2)) + pow((fabs(x[1] - 5*cos(x[0]))) , 1.0/3.0) + pow((fabs(x[2] - 5*sin(x[0]))) , 1.0/3.0);

	return;
}


/*
*
*
*   As described by Huband et al. in "A review of multiobjective test problems
*   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
*   Computing 10(5): 477-506, 2006.
*
*   Example Sch1, see the previous cited paper for the original reference.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void Sch1(double *f, double *x) {

	if(*x<0 || *x>5)
	{
		printf("Invalid input range!\n");
		return;
	}

	if((*x)<=1)
	{
		f[0] = -1*(*x);
	}
	else if((*x)<=3)
	{
		f[0] = -2 + (*x);
	}
	else if((*x)<=4)
	{
		f[0] = 4 - (*x);
	}
	else
	{
		f[0] = -4 + (*x);
	}

	f[1] = pow(((*x)-5),2);

	return;
}


/*
*
*
*   As described by Huband et al. in "A review of multiobjective test problems
*   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
*   Computing 10(5): 477-506, 2006.
*
*   Example SK1, see the previous cited paper for the original reference.
*   Function f2 differs in the original and in the cited references. The herein
*   codification follows the original reference.
*
*   In the above paper/papers the variables bounds were not set.
*   We considered -10<=x<=10.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void SK1(double *f, double *x) {
	if(x[0]<-10 || x[0]>10)
	{
		printf("Invalid input range!\n");
		return;
	}

	f[0] = -1*(-1*(pow((x[0]),4)) - 3*pow((x[0]),3) + 10*pow((x[0]),2) + 10*(x[0]) + 10);
	f[1] = -1*(0.5*(pow((x[0]),4)) + 2*pow((x[0]),3) + 10*pow((x[0]),2) - 10*(x[0]) + 5);

	return;
}


/*
*
*
*   As described by Huband et al. in "A review of multiobjective test problems
*   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
*   Computing 10(5): 477-506, 2006.
*
*   Example SP1, see the previous cited paper for the original reference.
*
*   In the above paper the variables bounds were not set.
*   We considered -1<=x[i]<=5, i=1,2.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void SP1(double *f, double *x) {
	int i;
	for(i=0;i<2;i++)
	{
		if(x[i]<-1 || x[i]>5)
		{
			printf("Invalid input range!\n");
			return;
		}
	}

	f[0] = pow((x[0] - 1),2) + pow((x[0] - x[1]),2);
	f[1] = pow((x[1] - 3),2) + pow((x[0] - x[1]),2);

	return;
}


/*
*
*
*   As described by Huband et al. in "A review of multiobjective test problems
*   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
*   Computing 10(5): 477-506, 2006.
*
*   Example SSFYY2, see the previous cited paper for the original reference.
*
*   This file is part of a collection of problems developed for
    the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void SSFYY2(double *f, double *x) {

  	double pi;
  	pi = 4 * atan(1);

	if(*x<(-100) || *x>100)
	{
		printf("Invalid input range!\n");
		return;
	}

	f[0] = 10 + pow(*x,2) - 10*cos((*x)*pi/2);
	f[1] = pow(((*x)-4),2);

	return;
}


/*
*
*
*   As described by Huband et al. in "A review of multiobjective test problems
*   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
*   Computing 10(5): 477-506, 2006.
*
*   Example TKLY1, see the previous cited paper for the original reference.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void TKLY1(double *f, double *x) {
	int i;
	double prod;
	double temp;

	prod = 1;
	temp = 0;

	if(x[0]<0.1 || x[0]>1)
	{
		printf("Invalid input range!\n");
		return;
	}

	for(i=1;i<4;i++)
	{
		if(x[i]<0 || x[i]>1)
		{
			printf("Invalid input range!\n");
			return;
		}
	}
	for(i=1;i<=3;i++)
	{
		temp = 2.0 - exp(-1*(pow(((x[i] - 0.1)/0.004),2))) - 0.8*exp(-1*pow(((x[i] - 0.9)/0.4),2));
		prod = prod*temp;
	}

	f[1] = x[0];
	f[0] = prod/x[0];


	return;
}


/*
*
*
*   As described by K. Deb in "Multi-objective Genetic Algorithms: Problem
*   Difficulties and Construction of Test Problems", Evolutionary Computation
*   7(3): 205-230, 1999.
*
*   Example 5.1.3 (Non-uniformly Represented Pareto-optimal Front).
*
*   In the above paper the variables bounds were not set.
*   We considered 0.0 <= x[i] <= 1.0, i=1,2. The function f1 corresponds
*   to equation (16) and gx to equation (10).
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void Deb53(double *f, double *x){
	int n;
	double beta;
	double alpha;
	double pi;
	double ff1;
	double gx;
	double h;
  int i;

	n = 2;
	beta = 1.0;
	alpha = 4.0;
	pi = 4*atan(1);

	ff1 = 1-exp(-4*x[0])*pow(sin(5*pi*x[0]), 4);

  for(i=0; i<n; ++i){
  	if (x[i]>1.0 || x[i]<0.0)
  	{
  		printf("Invalid input range.\n");
  		return;
  	}
  }

	if (x[1]<=0.4)
	{
		gx = 4-3*exp(-1*pow(((x[1]-0.2)/0.02), 2));
	}
	else
	{
		gx = 4-2*exp(-1*pow(((x[1]-0.7)/0.2), 2));
	}

	if (ff1<beta*gx)
	{
		h = 1-pow((ff1/(beta*gx)), alpha);
	}
	else
	{
		h = 0;
	}

	f[0] = ff1;
	f[1] = gx*h;
	return;
}


/*
*
*
*   As described by Huband et al. in "A review of multiobjective test problems
*   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
*   Computing 10(5): 477-506, 2006.
*
*   Example WFG6.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void WFG6(double* f, double* zz){
    int M;
    int k;
    int l;
    int n;
    int i, j;
    double pi2;
    double *S;
    double *A;
    double *zmax;
    double *y;
    double *t1;
    double *w;
    double *t2;
    double *x;
    double *h;
	double *z;
    int ii, jj;
    double sum1, sum2;

    M = 3;
    k = 4;
    l = 4;
    n = k+l;
    pi2 = 2.0*atan(1);
	z = (double *)malloc(sizeof(double)*(n+1));
    for(i=1; i<=n; i++)
    {
        z[i] = zz[i-1];
    }

    S = (double*)malloc(sizeof(double)*(M+1));
    for(i = 1; i <= M; i++)
    {
        S[i] = 2.0*(double)i;
    }

    A = (double*)malloc(sizeof(double)*(M));
    for(i = 1; i <= M-1; i++)
    {
        A[i] = 1.0;
    }

    zmax = (double*)malloc(sizeof(double)*(n+1));
    for(i = 1; i <= n; i++)
    {
        zmax[i] = 2.0*(double)i;
    }

    for(i = 1; i <= n; i++)
    {
        if(z[i] >= 0 && z[i] <= zmax[i])
        {
            continue;
        }
        else
        {
            printf("Invalid input range\n");
            return;
        }
    }

    y = (double*)malloc(sizeof(double)*(n+1));
    for(i = 1; i <= n; i++)
    {
        y[i] = z[i]/zmax[i];
    }

    t1 = (double*)malloc(sizeof(double)*(n+1));
    for(i = 1; i <= n; i++)
    {
        if (i<=k)
        {
            t1[i] = y[i];
        }
        else
        {
            t1[i] = fabs(y[i]-0.35)/fabs(floor(0.35-y[i])+0.35);
        }
    }

    w = (double*)malloc(sizeof(double)*(n+1));
    for(i = 1; i <= n; i++)
    {
        w[i] = 1.0;
    }

    t2 = (double*)malloc(sizeof(double)*(M+1));
    for(i = 1; i <= M; i++)
    {
        if(i <= M-1)
        {
            sum1 = 0.0;
            for(ii = ((i-1)*k/(M-1)+1); ii <=(i*k/(M-1)); ii++)
            {
                sum2 = 0.0;
                for(jj = 0; jj <= (k/(M-1)-2); jj++)
                {
                    sum2 = sum2 + (double)(fabs(t1[ii]-t1[(int)(((double)((i-1.0)*k)/(double)(M-1.0)+1.0)+fmod((ii+jj-((double)((i-1.0)*k)/(double)(M-1.0)+1.0)+1.0), (((double)(i*k)/(double)(M-1.0))-((double)((i-1.0)*k)/(double)(M-1.0)+1.0)+1.0)))]));
                }
                sum1 = sum1 + (double)(t1[ii]+sum2)/(double)((double)(((double)(i*k)/(double)(M-1.0))-((double)((i-1.0)*k)/(double)(M-1.0)+1.0)+1.0)/(double)(k/(M-1.0))*ceil((double)k/(double)(M-1.0)/(double)2.0)*(1.0+(double)(2.0*k)/(double)(M-1.0)-2.0*ceil((double)(k)/(double)(M-1.0)/(double)2.0)));
            }
            t2[i] = sum1;
        }
        else
        {
            sum1 = 0.0;
            for(ii = k+1; ii <= n; ii++)
            {
                sum2 = 0.0;
                for(jj = 0; jj <= l-2; jj++)
                {
                    sum2 = sum2 + (double)(fabs(t1[ii]-t1[(int)(k+1.0+fmod((ii+jj-(k+1.0)+1.0), (double)(n-k)))]));
                }
                sum1 = sum1 + (double)(t1[ii]+sum2)/(double)(((double)(n-k)/(double)l)*ceil((double)l/(double)2.0)*(1.0+2.0*(double)l-2.0*(double)(ceil((double)l/2.0))));
            }
            t2[i] = sum1;
        }
    }

    x = (double*)malloc(sizeof(double)*(M+1));
    for(i = 1; i<=M; i++)
    {
        if(i<=M-1)
        {
        	if(t2[M]>A[i])
        	{
            	x[i] = t2[M]*(t2[i]-0.5)+0.5;
            }
            else
            {
            	x[i] = A[i]*(t2[i]-0.5)+0.5;
            }
        }
        else
        {
            x[i] = t2[M];
        }
    }

    h = (double*)malloc(sizeof(double)*(M+1));
    for(i = 1; i <= M; i++)
    {
        if(i == 1)
        {
            h[i] = 1.0;
			for (j=1; j<=M-1; j++)
			{
            	h[i] = h[i]*(sin(x[j]*pi2));
			}
        }
        else if(i<=(M-1))
        {
			h[i] = cos(x[M-i+1]*pi2);
            for(j = 1; j <= M-i; j++)
            {
                h[i] = h[i]*(sin(x[j]*pi2));
            }
        }
        else
        {
            h[i] = cos(x[1]*pi2);
        }
    }

    for(i = 0; i < M; i++)
    {
        f[i] = x[M]+S[i+1]*h[i+1];
    }

    return;
}


/*
*
*
*   As described by Huband et al. in "A review of multiobjective test problems
*   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
*   Computing 10(5): 477-506, 2006.
*
*   Example WFG7.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void WFG7(double* f, double* zz){

  /*Note: Input variable here is z and not x*/

    int M;
    int k;
    int l;
    int n;

    double pi2;

    double *S, *A;
    double *zmax, *x, *y, *t1, *t2, *t3, *h, *w, *r_sum, *z;
    int i, j;
    double sum1, sum2, AA, BB, CC, temp;

    M = 3;
    k = 4;
    l = 4;
    n = k+l;
    pi2 = 2*atan(1);
    S = (double*)malloc(sizeof(double)*(M+1));
    A = (double*)malloc(sizeof(double)*(M));
    zmax = (double*)malloc(sizeof(double)*(n+1));
    y = (double*)malloc(sizeof(double)*(n+1));
    t1 = (double*)malloc(sizeof(double)*(n+1));
    t2 = (double*)malloc(sizeof(double)*(n+1));
    t3 = (double*)malloc(sizeof(double)*(M+1));
    x = (double*)malloc(sizeof(double)*(M+1));
    h = (double*)malloc(sizeof(double)*(M+1));
    w = (double *)malloc((n+1)*sizeof(double));
    r_sum = (double *)malloc((n+1)*sizeof(double));
    z = (double *)malloc(sizeof(double)*(n+1));

    AA = 0.98/49.98;
    BB = 0.02;
    CC = 50.0;

    for(i=n; i>=1; --i)
    {
      z[i] = zz[i-1];
    }

    for(i=1; i<=M-1; ++i){
      S[i] = 2.0*(double)i;
      A[i] = 1.0;
    }
    S[M] = 2.0*(double)M;

    for(i=1; i<=n; ++i){
      zmax[i] = 2.0*(double)i;
    }

    for(i=1;i<=n;i++)
    {
      w[i] = 1.0;
    }

     for(i=1; i<=n; ++i){
      y[i] = (double)z[i]/zmax[i];
    }

    for(i=1;i<=k;i++)
    {
      sum1=0.0;
      sum2=0.0;
      for(j = i+1;j<=n;j++)
      {
        sum1 = sum1 + w[j]*y[j];
        sum2 = sum2 + w[j];
      }
      r_sum[i] = sum1/sum2;
    }

    for(i=1; i<=n; ++i){
      if(i<=k)
      {
        temp = BB + (CC-BB) * (AA-(1.0-2.0*(r_sum[i])) * (fabs(floor(0.5-(r_sum[i]))+AA)));
        t1[i] = pow(y[i],temp);
      }
      else
      {
        t1[i] = y[i];
      }

      if(i<=k)
      {
        t2[i] = t1[i];
      }
      else
      {
        t2[i] = (fabs(t1[i]-0.35)/fabs(floor(0.35-t1[i])+0.35));
      }
    }

    for(i=1; i<=n; ++i) {
        if( z[i]<0 || z[i]>zmax[i] ) {
            printf("%lf, Invalid input range\n", z[i]);
            return;
        }
    }

    for(i=1; i<=M; ++i){
      if(i<=(M-1)){
        sum1 = 0.0;
        sum2 = 0.0;
        for(j=((i-1)*k/(M-1)+1); j<=((i)*k/(M-1)); ++j){
          sum1 = sum1 + (w[j]*t2[j]);
          sum2 = sum2 + (w[j]);
        }
        t3[i] = sum1/sum2;
      }
      else{
        sum1 = 0.0;
        sum2 = 0.0;
        for(j=k+1; j<=n; ++j){
          sum1 = sum1 + (w[j]*t2[j]);
          sum2 = sum2 + (w[j]);
        }
        t3[i] = sum1/sum2;
      }
    }
    for(i=1; i<=M; ++i){
      if(i<=M-1){
        if(t3[M]>A[i])
        {
          x[i] = t3[M]*(t3[i]-0.5)+0.5;
        }
        else
        {
          x[i] = A[i]*(t3[i]-0.5)+0.5;
        }
      }
      else
        x[i] = t3[M];
    }

    for(i=1;i<=M;i++)
    {
      if(i==1){
          h[i] = 1.0;
          for(j=1; j<=M-1; ++j)
          {
            h[i] = h[i] * (sin(x[j]*pi2));
          }
        }
        else if(i<=M-1){
          h[i] = 1.0;
          for(j=1; j<=M-i; ++j)
          {
            h[i] = h[i] * sin(x[j]*pi2);
          }
          h[i] = h[i] * cos(x[M-i+1]*pi2);
        }
        else{
          h[i] = cos(x[1]*pi2);
        }
    }

    for(i=1; i<=M; ++i)
    {
      f[i-1] = x[M]+S[i]*h[i];
    }

    return;
}


/*
*
*
*   As described by Huband et al. in "A review of multiobjective test problems
*   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
*   Computing 10(5): 477-506, 2006.
*
*   Example WFG8.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
*
 */
void WFG8(double* f, double* zz){

  /*Note: Input variable here is z and not x*/

    int M;
    int k;
    int l;
    int n;

    double pi2;

    double *S, *A;
    double *zmax, *x, *y, *t1, *t2, *t3, *h, *w, *r_sum, *z;
    int i, j;
    double sum1, sum2, AA, BB, CC, temp;

    M = 3;
    k = 4;
    l = 4;
    n = k+l;
    pi2 = 2*atan(1);
    S = (double*)malloc(sizeof(double)*(M+1));
    A = (double*)malloc(sizeof(double)*(M));
    zmax = (double*)malloc(sizeof(double)*(n+1));
    y = (double*)malloc(sizeof(double)*(n+1));
    t1 = (double*)malloc(sizeof(double)*(n+1));
    t2 = (double*)malloc(sizeof(double)*(n+1));
    t3 = (double*)malloc(sizeof(double)*(M+1));
    x = (double*)malloc(sizeof(double)*(M+1));
    h = (double*)malloc(sizeof(double)*(M+1));
    w = (double *)malloc((n+1)*sizeof(double));
    r_sum = (double *)malloc((n+1)*sizeof(double));
    z = (double *)malloc(sizeof(double)*(n+1));

    AA = 0.98/49.98;
    BB = 0.02;
    CC = 50.0;

    for(i=n; i>=1; --i)
    {
      z[i] = zz[i-1];
    }

    for(i=1; i<=M-1; ++i){
      S[i] = 2.0*(double)i;
      A[i] = 1.0;
    }
    S[M] = 2.0*(double)M;

    for(i=1; i<=n; ++i){
      zmax[i] = 2.0*(double)i;
    }

    for(i=1;i<=n;i++)
    {
      w[i] = 1.0;
    }

     for(i=1; i<=n; ++i){
      y[i] = (double)z[i]/zmax[i];
    }

    for(i=k+1;i<=n;i++)
    {
      sum1=0.0;
      sum2=0.0;
      for(j = 1;j<=i-1;j++)
      {
        sum1 = sum1 + w[j]*y[j];
        sum2 = sum2 + w[j];
      }
      r_sum[i] = sum1/sum2;
    }

    for(i=1; i<=n; ++i){
      if(i<=k)
      {
        t1[i] = y[i];
      }
      else
      {
        temp = BB + (CC-BB) * (AA-(1.0-2.0*(r_sum[i])) * (fabs(floor(0.5-(r_sum[i]))+AA)));
        t1[i] = pow(y[i],temp);
      }

      if(i<=k)
      {
        t2[i] = t1[i];
      }
      else
      {
        t2[i] = (fabs(t1[i]-0.35)/fabs(floor(0.35-t1[i])+0.35));
      }
    }

    for(i=1; i<=n; ++i) {
        if( z[i]<0 || z[i]>zmax[i] ) {
            printf("%lf, Invalid input range\n", z[i]);
            return;
        }
    }

    for(i=1; i<=M; ++i){
      if(i<=(M-1)){
        sum1 = 0.0;
        sum2 = 0.0;
        for(j=((i-1)*k/(M-1)+1); j<=((i)*k/(M-1)); ++j){
          sum1 = sum1 + (w[j]*t2[j]);
          sum2 = sum2 + (w[j]);
        }
        t3[i] = sum1/sum2;
      }
      else{
        sum1 = 0.0;
        sum2 = 0.0;
        for(j=k+1; j<=n; ++j){
          sum1 = sum1 + (w[j]*t2[j]);
          sum2 = sum2 + (w[j]);
        }
        t3[i] = sum1/sum2;
      }
    }
    for(i=1; i<=M; ++i){
      if(i<=M-1){
        if(t3[M]>A[i])
        {
          x[i] = t3[M]*(t3[i]-0.5)+0.5;
        }
        else
        {
          x[i] = A[i]*(t3[i]-0.5)+0.5;
        }
      }
      else
        x[i] = t3[M];
    }

    for(i=1;i<=M;i++)
    {
      if(i==1){
          h[i] = 1.0;
          for(j=1; j<=M-1; ++j)
          {
            h[i] = h[i] * (sin(x[j]*pi2));
          }
        }
        else if(i<=M-1){
          h[i] = 1.0;
          for(j=1; j<=M-i; ++j)
          {
            h[i] = h[i] * sin(x[j]*pi2);
          }
          h[i] = h[i] * cos(x[M-i+1]*pi2);
        }
        else{
          h[i] = cos(x[1]*pi2);
        }
    }

    for(i=1; i<=M; ++i)
    {
      f[i-1] = x[M]+S[i]*h[i];
    }

    return;
}


/*
 *
 *   As described by Huband et al. in "A review of multiobjective test problems
 *   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
 *   Computing 10(5): 477-506, 2006.
 *
 *   Example WFG3.
 *
 *   This file is part of a collection of problems developed for
 *   the BMOBench platform and based on
 *   derivative-free multiobjective optimization in
 *   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *   Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *   Written by the authors in June 1, 2010.
 *
 */

void WFG3(double* f, double* z){
	int M = 3;
	int k = 4;
	int l = 4;
	int n = k + l;
	double AA = 2;

	int m, i, ii, jj, j;

	double *S, *A, *zmax, *y, *t1, *t2, *t3, *w, *x, *h;

	double sum1, sum2;

	S = (double*)malloc(sizeof(double)*(M+1));
	for(m = 1; m <= M; m++)
		S[m] = 2*m;

	A = (double*)malloc(sizeof(double)*(M));
	for(i = 1; i <= M-1; i++)
		if( i <= 1 )
			A[i] = 1;
		else
			A[i] = 0;

	zmax = (double*)malloc(sizeof(double)*(n+1));
	for(i = 1; i <= n; i++)
		zmax[i] = 2*i;

	for(i = 1; i <= n; i++){
		if(z[i-1] >= 0 && z[i-1] <= zmax[i]);
		else{
			printf("Invalid input range\n");
			return;
		}
	}

	y = (double*)malloc(sizeof(double)*(n+1));
	for(i = 1; i <= n; i++)
		y[i] = z[i-1] / zmax[i];

	t1 = (double*)malloc(sizeof(double)*(n+1));
	for(i = 1; i <= n; i++){
		if(i <= k)
			t1[i] = y[i];
		else{
			t1[i] = fabs(y[i]-0.35)/fabs(floor(0.35-y[i])+0.35);
		}
	}

	t2 = (double*)malloc(sizeof(double)*((k+l/2)+1));
	for(i = 1; i <= k+l/2; i++){
		if(i <= k)
			t2[i] = t1[i];
		else{
			sum1 = 0;
			for(ii = (k+2*(i-k)-1); ii <= (k+2*(i-k)); ii++){
				sum2 = 0;
				for(jj = 0; jj <= AA-2; jj++)
					sum2 = sum2 + fabs(t1[ii]-t1[(k+2*(i-k)-1)+((ii+jj-(k+2*(i-k)-1)+1) % (k+2*(i-k)-(k+2*(i-k)-1)+1))]);
				sum1 = sum1 + t1[ii] + sum2;
			}
			t2[i] = sum1 / ((k+2*(i-k)-(k+2*(i-k)-1)+1)/AA*ceil(AA/2)*(1+2*AA-2*ceil(AA/2)));
		}
	}

	w = (double*)malloc(sizeof(double)*(n+1));
	for(i = 1; i <= n; i++)
		w[i] = 1;

	t3 = (double*)malloc(sizeof(double)*(M+1));
	for(i = 1; i <= M; i++){
		if(i <= M-1){
			sum1 = 0;
			for(j = ((i-1)*k/(M-1)+1); j <= (i*k/(M-1)); j++)
				sum1 = sum1 + (w[j]*t2[j]);
			sum2 = 0;
			for(j = ((i-1)*k/(M-1)+1); j <= (i*k/(M-1)); j++)
				sum2 = sum2 + w[j];
			t3[i] = sum1 / sum2;
		}
		else{
			sum1 = 0;
			for(j = k+1; j <= k+l/2; j++)
				sum1 = sum1 + (w[j]*t2[j]);
			sum2 = 0;
			for(j = k+1; j <= k+l/2; j++)
				sum2 = sum2 + w[j];
			t3[i] = sum1 / sum2;
		}
	}

	x = (double*)malloc(sizeof(double)*(M+1));
	for(i = 1; i <= M; i++){
		if(i <= M-1)
			x[i] = ((t3[M]>A[i])?t3[M]:A[i])*(t3[i]-0.5)+0.5;
		else
			x[i] = t3[M];
	}

	h = (double*)malloc(sizeof(double)*(M+1));
	for(m = 1; m <= M; m++){
		if(m == 1){
			h[m] = 1;
			for(i = 1; i <= M-1; i++)
				h[m] = h[m] * x[i];
		}else if(m <= M-1){
			h[m] = 1;
			for(i = 1; i <= M-m; i++)
				h[m] = h[m] * x[i];
			h[m] = h[m] * (1-x[M-m+1]);
		}else{
			h[m] = 1-x[1];
		}
	}

	for(m = 1; m <= M; m++)
		f[m-1] = x[M]+S[m]*h[m];

	return;
}

/*
 *
 *   As described by Huband et al. in "A review of multiobjective test problems
 *   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
 *   Computing 10(5): 477-506, 2006.
 *
 *   Example MOP1, Van Valedhuizen's test suit.
 *
 *   This file is part of a collection of problems developed for
 *   the BMOBench platform and based on
 *   derivative-free multiobjective optimization in
 *   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *   Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *   Written by the authors in June 1, 2010.
 *
 */

void MOP1(double* f, double* x){
  if(x[0] >= -100000 && x[0] <= 100000);
  else{
    printf("Invalid input range\n");
    return;
  }

  f[0] = x[0]*x[0];
  f[1] = (x[0]-2)*(x[0]-2);

  return;
}

/*
 *
 *   As described by Huband et al. in "A review of multiobjective test problems
 *   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
 *   Computing 10(5): 477-506, 2006.
 *
 *   Example MOP3, Van Valedhuizen's test suit.
 *
 *   This file is part of a collection of problems developed for
 *   the BMOBench platform and based on
 *   derivative-free multiobjective optimization in
 *   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *   Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *   Written by the authors in June 1, 2010.
 *
 */

void MOP3(double* f, double* x){
  double pi = 4 * atan(1);
  double A1 = 0.5*sin(1)-2*cos(1)+sin(2)-1.5*cos(2);
  double A2 = 1.5*sin(1)-cos(1)+2*sin(2)-0.5*cos(2);
  double B1 = 0.5*sin(x[0])-2*cos(x[0])+sin(x[1])-1.5*cos(x[1]);
  double B2 = 1.5*sin(x[0])-cos(x[0])+2*sin(x[1])-0.5*cos(x[1]);
  double n = 2;
  int i;

  for(i = 0; i < n; i++){
    if(x[i] >= -pi && x[i] <= pi)
      continue;
    else{
      printf("Invalid input range\n");
      return;
    }
  }

  f[0] = -(-1-pow((A1-B1), 2)-pow((A2-B2), 2));
  f[1] = -(-pow((x[0]+3), 2)-pow((x[1]+1), 2));

  return;
}

/*
 *
 *   As described by Huband et al. in "A review of multiobjective test problems
 *   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
 *   Computing 10(5): 477-506, 2006.
 *
 *   Example MOP4, Van Valedhuizen's test suit.
 *
 *   This file is part of a collection of problems developed for
 *   the BMOBench platform and based on
 *   derivative-free multiobjective optimization in
 *   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *   Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *   Written by the authors in June 1, 2010.
 *
 */

void MOP4(double* f, double* x){
  double n = 3;
  int i;

  for(i = 0; i < n; i++){
    if(x[i] >= -5 && x[i] <= 5)
      continue;
    else{
      printf("Invalid input range\n");
      return;
    }
  }

  f[0] = 0;
  for(i = 0; i < n-1; i++){
    f[0] = f[0] + (-10*exp(-0.2*sqrt(x[i]*x[i]+x[i+1]*x[i+1])));
  }

  f[1] = 0;
  for(i = 0; i < n; i++){
    f[1] = f[1] + (pow(fabs(x[i]), 0.8) + 5*sin(pow(x[i], 3)));
  }

  return;
}

/*
 *
 *   As described by Huband et al. in "A review of multiobjective test problems
 *   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
 *   Computing 10(5): 477-506, 2006.
 *
 *   Example MOP5, Van Valedhuizen's test suit.
 *
 *   This file is part of a collection of problems developed for
 *   the BMOBench platform and based on
 *   derivative-free multiobjective optimization in
 *   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *   Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *   Written by the authors in June 1, 2010.
 *
 */

void MOP5(double* f, double* x){
  double n = 2;
  int i;

  for(i = 0; i < n; i++){
    if(x[i] >= -30 && x[i] <= 30)
      continue;
    else{
      printf("Invalid input range\n");
      return;
    }
  }

  f[0] = 0.5*(x[0]*x[0]+x[1]*x[1])+sin(x[0]*x[0]+x[1]*x[1]);
  f[1] = pow((3*x[0]-2*x[1]+4), 2)/8+pow((x[0]-x[1]+1), 2)/27+15;
  f[2] = 1/(x[0]*x[0]+x[1]*x[1]+1)-1.1*exp(-x[0]*x[0]-x[1]*x[1]);

  return;
}

/*
 *
 *   As described by Huband et al. in "A review of multiobjective test problems
 *   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
 *   Computing 10(5): 477-506, 2006.
 *
 *   Example MOP7, Van Valedhuizen's test suit.
 *
 *   This file is part of a collection of problems developed for
 *   the BMOBench platform and based on
 *   derivative-free multiobjective optimization in
 *   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *   Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *   Written by the authors in June 1, 2010.
 *
 */

void MOP7(double* f, double* x){
  double n = 2;
  int i;

  for(i = 0; i < n; i++){
    if(x[i] >= -400 && x[i] <= 400)
      continue;
    else{
      printf("Invalid input range\n");
      return;
    }
  }

  f[0] = pow((x[0]-2), 2)/2+pow((x[1]+1), 2)/13+3;
  f[1] = pow((x[0]+x[1]-3), 2)/36+pow((-x[0]+x[1]+2), 2)/8-17;
  f[2] = pow((x[0]+2*x[1]-1), 2)/175+pow((-x[0]+2*x[1]), 2)/17-13;

  return;
}

/*
*
*   As described by Huband et al. in "A review of multiobjective test problems
*   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
*   Computing 10(5): 477-506, 2006.
*
*   Example MLF1, see the previous cited paper for the original reference.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
 */

void MLF1(double* f, double* x){
  double n = 1;
  int i;

  for(i = 0; i < n; i++){
    if(x[i] >= 0 && x[i] <= 20)
      continue;
    else{
      printf("Invalid input range\n");
      return;
    }
  }

  f[0] = (1+x[0]/20)*sin(x[0]);
  f[1] = (1+x[0]/20)*cos(x[0]);

  return;
}

/*
*
*   As described by Huband et al. in "A review of multiobjective test problems
*   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
*   Computing 10(5): 477-506, 2006.
*
*   Example MLF2, See the previous cited paper for the original reference.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
 */

void MLF2(double* f, double* x){
  double n = 2;

  int i;
  double *y;

  for(i = 0; i < n; i++){
    if(x[i] >= -2 && x[i] <= 2)
      continue;
    else{
      printf("Invalid input range\n");
      return;
    }
  }

  y = (double*)malloc(sizeof(double)*n);
  for(i = 0; i < n; i++)
    y[i] = 2*x[i];


  f[0] = -(5-(pow((x[0]*x[0]+x[1]-11), 2)+pow((x[0]+x[1]*x[1]-7), 2))/200);
  f[1] = -(5-(pow((y[0]*y[0]+y[1]-11), 2)+pow((y[0]+y[1]*y[1]-7), 2))/200);

  return;
}

/*
*
*   As described by Huband et al. in "A review of multiobjective test problems
*   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
*   Computing 10(5): 477-506, 2006.
*
*   Example SSFYY1, see the previous cited paper for the original reference.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
 */

void SSFYY1(double* f, double* x){
  double n = 2;
  int i;

  for(i = 0; i < n; i++){
    if(x[i] >= -100 && x[i] <= 100)
      continue;
    else{
      printf("Invalid input range\n");
      return;
    }
  }

  f[0] = x[0]*x[0] + x[1]*x[1];
  f[1] = (x[0]-1)*(x[0]-1) + (x[1]-2)*(x[1]-2);

  return;
}

/*
*
*   As described by Huband et al. in "A review of multiobjective test problems
*   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
*   Computing 10(5): 477-506, 2006.
*
*   Example SK2, see the previous cited paper for the original reference.
*
*   In the above paper/papers the variables bounds were not set.
*   We considered -10<=x[i]<=10, i=1,2,3,4.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
 */

void SK2(double* f, double* x){
  double n = 4;
  int i;
  double numerator = 0;
  double denominator = 0;

  for(i = 0; i < n; i++){
    if(x[i] >= -10 && x[i] <= 10)
      continue;
    else{
      printf("Invalid input range\n");
      return;
    }
  }

  f[0] = -(-pow((x[0]-2), 2)-pow((x[1]+3), 2)-pow((x[2]-5), 2)-pow((x[3]-4), 2)+5);
  f[1] = 0;
  for(i = 0; i < n; i++){
    numerator = numerator + sin(x[i]);
  }
  numerator = -1 * numerator;

  for(i = 0; i < n; i++){
    denominator = denominator + (x[i]*x[i])/100;
  }
  denominator = 1 + denominator;

  f[1] = numerator / denominator;

  return;
}

/*
*
*   As described by Huband et al. in "A review of multiobjective test problems
*   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
*   Computing 10(5): 477-506, 2006.
*
*   Example QV1, see the previous cited paper for the original reference.
*
*   In the original reference the number of variables was n=16.
*   We selected n=10 as default.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
 */

void QV1(double* f, double* x){
  double n = 10;
  int i;
  double pi = 4 * atan(1);

  for(i = 0; i < n; i++){
    if(x[i] >= -5.12 && x[i] <= 5.12)
      continue;
    else{
      printf("Invalid input range\n");
      return;
    }
  }

  f[0] = 0;
  f[1] = 0;

  for(i = 0; i < n; i++){
    f[0] = f[0] + (x[i]*x[i]-10*cos(2*pi*x[i])+10)/n;
  }
  f[0] = pow(f[0], 0.25);

  for(i = 0; i < n; i++){
    f[1] = f[1] + ((x[i] - 1.5)*(x[i] - 1.5)-10*cos(2*pi*(x[i] - 1.5))+10)/n;
  }
  f[1] = pow(f[1], 0.25);
  return;
}

/*
*
*   As described by Huband et al. in "A review of multiobjective test problems
*   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
*   Computing 10(5): 477-506, 2006.
*
*   Example MOP6, Van Valedhuizen's test suit.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
 */

void MOP6(double* f, double* x){
  double n = 2;
  int i;
  double pi = 4 * atan(1);

  for(i = 0; i < n; i++){
    if(x[i] >= 0 && x[i] <= 1)
      continue;
    else{
      printf("Invalid input range\n");
      return;
    }
  }

  f[1] = x[0];
  f[0] = (1+10*x[1])*(1-pow((x[0]/(1+10*x[1])), 2)-x[0]/(1+10*x[1])*sin(8*pi*x[0]));

  return;
}

/*
*   As described by Huband et al. in "A review of multiobjective test problems
*   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
*   Computing 10(5): 477-506, 2006.
*
*   Example MOP2, Van Valedhuizen's test suit.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
 */

void MOP2(double* f, double* x){
  double n = 4;
  int i;

  for(i = 0; i < n; i++){
    if(x[i] >= -4 && x[i] <= 4)
      continue;
    else{
      printf("Invalid input range\n");
      return;
    }
  }
  f[0] = 0;
  f[1] = 0;

  for(i = 0; i < n; i++){
    f[0] = f[0] + (x[i] - 1/sqrt(n))*(x[i] - 1/sqrt(n));
  }

  f[0] = exp(-f[0]);
  f[0] = 1 - f[0];

  for(i = 0; i < n; i++){
    f[1] = f[1] + (x[i] + 1/sqrt(n))*(x[i] + 1/sqrt(n));
  }

  f[1] = exp(-f[1]);
  f[1] = 1 - f[1];

  return;
}

/*
*
*   As described by A. Lovison in "A synthetic approach to multiobjective
*   optimization", arxiv Item: http://arxiv.org/abs/1002.0093.
*
*   Example 5.
*
*   In the above paper/papers the variables bounds were not set.
*   We considered -1<=x[i]<=4, i=1,2,3.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
 */

void lovison6(double* f, double* x){
  int n = 3;
  int m = 4;
  //FILE *myFile;
  //double **alpha = malloc(sizeof(double*)*m);
  //double **C = malloc(sizeof(double*)*n);
  int j;
  int i;
  //double *beta;
  //double *gamma;
  double *func;
  double pi = 4*atan(1);








  double C[3][4] = {{0.218418,-0.620254,0.843784,0.914311},
    {-0.788548,0.428212,0.103064,-0.47373},
    {-0.300792,-0.185507,0.330423,0.151614}};


  double alpha[4][3] = {{0.942022,0.363525,0.00308876},
                        {0.755598,0.450103,0.170122},
                        {0.787748,0.837808,0.590166},
                        {0.203093,0.253639,0.532339}};


   double beta[4] = {-0.666503,-0.945716,-0.334582,0.611894};
   double gamma[4] = {0.281032,0.508749,-0.0265389,-0.920133};


  for(i = 0; i < n; i++){
    if(x[i] >= -1 && x[i] <= 4)
      continue;
    else{
      printf("Invalid input range\n");
      return;
    }
  }
/*
  for(i = 0; i < n; i++)
    C[i] = malloc(sizeof(double*)*m);

  myFile = (FILE*)fopen("./data/lovison6_Avalues_C.txt", "r");

  for (i=0; i<n; i++)
  {
    for (j=0; j<m; j++)
    {
      fscanf(myFile, "%lf", &C[i][j]);
    }
  }
  fclose(myFile);

  for(i = 0; i < m; i++)
    alpha[i] = malloc(sizeof(double*)*n);


  myFile = (FILE*)fopen("./data/lovison6_Avalues_alpha.txt", "r");

  for (i=0; i<m; i++)
  {
    for (j=0; j<n; j++)
    {
      fscanf(myFile, "%lf", &alpha[i][j]);
    }
  }
  fclose(myFile);

  beta = (double*)malloc(sizeof(double)*m);

  myFile = (FILE*)fopen("./data/lovison6_Avalues_beta.txt", "r");
  for (i=0; i<m; i++)
  {
      fscanf(myFile, "%lf", &beta[i]);
  }
  fclose(myFile);

  gamma = (double*)malloc(sizeof(double)*m);
  myFile = (FILE*)fopen("./data/lovison6_Avalues_gamma.txt", "r");

  for (i=0; i<m; i++)
  {
      fscanf(myFile, "%lf", &gamma[i]);
  }
  fclose(myFile);

 */
  func = (double*)malloc(sizeof(double)*m);
  for(j = 0; j < m; j++){
    func[j] = 0;
    for(i = 0; i < n; i++){
      func[j] = func[j] + (-alpha[j][i]*pow(x[i] - C[i][j], 2));
    }
  }

  f[0] = -(func[0] + beta[0]*exp(func[3]/gamma[0]));
  f[1] = -(func[1] + beta[1]*sin(pi*(x[0]+x[1])/gamma[1]));
  f[2] = -(func[2] + beta[2]*cos(pi*(x[0]-x[1])/gamma[2]));

  return;
}

/*
*
*   As described by A. Lovison in "A synthetic approach to multiobjective
*   optimization", arxiv Item: http://arxiv.org/abs/1002.0093.
*
*   Example 4.
*
*   In the above paper/papers the variables bounds were not set.
*   We considered 0<=x[1]<=6 and -1<=x[2]<=1.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
 */

void lovison4(double* f, double* x){
  if((x[0] >= 0 && x[0] <= 6) && (x[1] >= -1 && x[1] <= 1));
  else{
    printf("Invalid input range\n");
    return;
  }

  f[0] = -(-x[0]*x[0]-x[1]*x[1]-4*(exp(-pow((x[0]+2), 2)-x[1]*x[1])+exp(-pow((x[0]-2), 2)-x[1]*x[1])));
  f[1] = pow((x[0]-6), 2)+pow((x[1]+0.5), 2);

  return;
}

/*
*
*   As described by A. Lovison in "A synthetic approach to multiobjective
*   optimization", arxiv Item: http://arxiv.org/abs/1002.0093.
*
*   Example 2.
*
*   In the above paper/papers the variables bounds were not set.
*   We considered -0.5<=x[1]<=0 and -0.5<=x[2]<=0.5.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
 */

void lovison2(double* f, double* x){
  if((x[0] >= -0.5 && x[0] <= 0) && (x[1] >= -0.5 && x[1] <= 0.5));
  else{
    printf("Invalid input range\n");
    return;
  }

  f[1] = x[1];
  f[0] = -(x[1] - pow(x[0], 3)) / (x[0] + 1);

  return;
}

/*
*
*   As described by Huband et al. in "A review of multiobjective test problems
*   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
*   Computing 10(5): 477-506, 2006.
*
*   Example LE1, see the previous cited paper for the original reference.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
 */

void LE1(double* f, double* x){
  int i;
  int n = 2;
  for(i = 0; i < n; i++){
    if(x[i] >= -5 && x[i] <= 10)
      continue;
    else{
      printf("Invalid input range\n");
      return;
    }
  }

  f[0] = pow((x[0]*x[0]+x[1]*x[1]), 0.125);
  f[1] = pow((pow((x[0]-0.5), 2)+pow((x[1]-0.5), 2)), 0.25);

  return;
}

/*
*
*   As described by K. Deb and A. Sinha and S. Kukkonen in "Multi-objective
*   test problems, linkages, and evolutionary methodologies", GECCO'06}:
*   Proceedings of the 8th Annual Conference on Genetic and Evolutionary
*   Computation, 1141-1148, 2006.
*
*   Example T3, with linkage L2.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
 */

void L2ZDT3(double* f, double* x){
  int m = 30;
  double ff1, gx, h;
  double *y;
  double pi = 4 * atan(1);

  int i, j;
  //FILE *myFile;
  //double **M = (double **)malloc(sizeof(double*)*m);
  double M[900] = {0.218418,-0.620254,0.843784,0.914311,-0.788548,0.428212,0.103064,-0.47373,-0.300792,-0.185507,0.330423,0.151614,0.884043,-0.272951,-0.993822,0.511197,-0.0997948,-0.659756,0.575496,0.675617,0.180332,-0.593814,-0.492722,0.0646786,-0.666503,-0.945716,-0.334582,0.611894,0.281032,0.508749,-0.0265389,-0.920133,0.308861,-0.0437502,-0.374203,0.207359,-0.219433,0.914104,0.184408,0.520599,-0.88565,-0.375906,-0.708948,-0.37902,0.576578,0.0194674,-0.470262,0.572576,0.351245,-0.480477,0.238261,-0.1596,-0.827302,0.669248,0.494475,0.691715,-0.198585,0.0492812,0.959669,0.884086,-0.218632,-0.865161,-0.715997,0.220772,0.692356,0.646453,-0.401724,0.615443,-0.0601957,-0.748176,-0.207987,-0.865931,0.613732,-0.525712,-0.995728,0.389633,-0.064173,0.662131,-0.707048,-0.340423,0.60624,0.0951648,-0.160446,-0.394585,-0.167581,0.0679849,0.449799,0.733505,-0.00918638,0.00446808,0.404396,0.449996,0.162711,0.294454,-0.563345,-0.114993,0.549589,-0.775141,0.677726,0.610715,0.0850755,0.0419388,-0.323614,-0.973719,-0.680238,-0.270873,-0.209617,0.968436,0.908798,0.975851,-0.994918,-0.0621977,0.628171,0.761228,0.34372,-0.792042,-0.144765,-0.965748,0.0133606,-0.0260565,-0.742377,0.426391,0.408202,0.633885,-0.0351053,-0.723444,-0.577654,0.0276004,0.0712472,-0.622791,0.155451,0.442717,-0.792786,0.925785,0.670266,-0.865566,-0.638281,0.333094,0.477628,0.47261,0.23151,0.82132,-0.589803,0.796275,0.57713,0.101149,0.970191,0.532821,0.814769,-0.0687269,0.712758,-0.191812,-0.390938,0.952828,0.921519,0.923094,0.93011,-0.945394,-0.0934027,0.964123,-0.795609,-0.289563,0.614236,-0.670585,0.466877,0.144597,-0.206416,0.6937,-0.967958,-0.0951247,-0.942473,-0.610767,-0.655472,-0.0960986,-0.302779,-0.734976,-0.342188,-0.315861,-0.912834,0.24499,0.0969326,0.089775,-0.241157,0.0835558,-0.420236,-0.686633,-0.711276,-0.00325377,0.435196,-0.710002,0.00283691,-0.168757,-0.134045,-0.655235,0.172361,0.998291,0.376291,-0.962215,-0.363174,-0.88777,-0.519929,-0.560554,-0.984415,0.601529,-0.984103,-0.228237,-0.578066,0.307023,0.606123,0.959635,0.00225943,0.0101814,0.441456,0.0633629,0.406631,-0.0100638,-0.177972,-0.491075,0.537035,-0.924987,-0.699424,0.742285,0.0181443,0.718971,-0.0308272,0.086931,0.524476,0.956457,0.143024,0.616481,0.217909,-0.128427,-0.262427,-0.938208,-0.52479,0.12919,0.721925,0.766492,0.470845,-0.0976466,0.507807,0.804148,0.963269,0.357128,-0.832565,-0.312441,0.327779,0.184745,0.246139,-0.936814,-0.931734,-0.0327827,0.319293,0.044473,-0.641645,0.596118,-0.293934,-0.63373,0.409658,0.759892,-0.257078,0.939616,-0.227661,0.115754,0.10964,-0.240557,0.66842,0.855535,-0.451536,0.264961,-0.61366,-0.204783,-0.842476,-0.249524,-0.0985226,0.0671501,-0.527707,-0.509489,-0.883254,0.14851,-0.906465,0.496238,-0.853211,-0.779234,-0.979515,0.827175,0.228969,-0.402829,-0.970118,0.762559,0.506495,0.460303,0.897304,0.686003,0.739986,0.15731,0.281697,-0.922955,-0.780824,0.449716,0.125225,0.487649,0.147046,0.679639,0.593707,-0.311828,-0.797099,-0.35815,0.95808,0.907244,0.772426,0.720574,-0.873217,0.371431,-0.826029,0.942716,0.70609,-0.658158,-0.782185,-0.806743,-0.627986,-0.405551,-0.258495,-0.796524,0.222498,0.087545,-0.0917108,-0.62542,-0.110256,0.0417576,0.24476,0.941339,-0.613783,0.402772,0.300775,-0.820314,-0.894233,-0.405896,0.0735439,0.486645,-0.394355,0.125097,-0.316386,-0.701215,-0.845742,0.2065,-0.413743,0.406725,-0.423813,-0.941255,-0.558804,0.312326,0.345314,0.319143,-0.644653,-0.0408415,0.176461,0.740113,0.470737,-0.914927,-0.591523,-0.606614,-0.181873,0.692975,0.50208,-0.536704,0.359652,0.839082,0.56817,-0.0776788,-0.00332785,0.459538,-0.518313,-0.270738,-0.629958,-0.755084,-0.721573,0.431107,-0.221877,0.32543,0.163743,0.0759916,0.695064,-0.656856,0.074163,0.264319,-0.73174,0.731548,-0.489341,0.678946,0.0271269,0.804879,-0.402973,0.800373,0.760082,-0.878364,0.176801,-0.548932,-0.225601,-0.164912,-0.208143,0.7768,-0.542743,-0.156021,0.671736,0.878648,-0.419588,-0.0752896,0.0299447,-0.494459,-0.72415,0.35978,-0.32646,-0.96605,0.0127605,0.563174,-0.814853,-0.949609,-0.526794,-0.801902,-0.753397,0.617418,0.689874,0.983384,0.668786,0.0304653,-0.625221,-0.13318,0.827343,-0.101358,-0.999522,-0.0525574,-0.458319,0.587409,-0.334639,0.0759643,0.0255827,0.128944,0.17317,-0.284309,0.287161,-0.550725,-0.433083,-0.242821,0.878879,0.691699,-0.660499,0.389985,0.599856,-0.711442,-0.798697,-0.244945,-0.942649,0.402856,-0.494672,0.439941,-0.88216,0.170196,0.650734,-0.0982391,-0.468732,0.342133,-0.838071,-0.832362,0.658177,-0.565361,0.149473,0.69331,-0.491848,0.74916,0.526025,-0.155339,0.0998096,0.468761,0.324649,0.128488,0.544144,-0.495222,0.965229,-0.79314,-0.545421,-0.500243,0.154371,0.170017,-0.259108,-0.868862,-0.50731,-0.848317,0.835712,0.616391,-0.442608,-0.158,0.313451,0.703748,-0.755984,-0.249443,0.491564,0.985068,0.678644,0.808324,0.81975,-0.435823,-0.839855,0.00282368,-0.569165,0.0884339,-0.222144,0.499412,-0.565198,0.64824,0.956914,-0.0620912,0.634479,0.928617,0.464664,0.377022,0.63047,-0.198619,-0.576153,0.565373,-0.524245,-0.187299,-0.614524,0.429316,-0.491171,0.399495,-0.333898,-0.646636,-0.0189709,-0.339605,-0.798791,0.0494081,0.367012,0.852545,0.43557,0.150039,-0.0454542,0.604861,-0.598288,-0.500696,0.249008,0.370711,-0.633174,-0.0121906,0.42006,0.169373,-0.975542,-0.0297852,0.80481,0.638317,-0.670967,0.935792,-0.35605,0.175773,0.878601,-0.275168,-0.932517,-0.372497,-0.0732907,-0.185493,-0.357004,0.314786,-0.229239,0.530256,-0.51327,0.44187,0.940309,-0.240334,-0.0276121,0.74383,-0.630329,-0.763376,0.62538,0.818945,0.891598,0.680494,0.471868,-0.769787,-0.878099,-0.973724,0.354362,-0.1792,-0.225034,-0.44548,0.598865,0.544005,-0.478962,0.327193,-0.525784,0.903179,-0.899248,0.156514,0.154329,0.499808,-0.836327,-0.802627,0.378082,-0.112673,-0.47926,-0.3355,-0.699445,0.237731,-0.324597,-0.800406,-0.42585,-0.710739,-0.144068,-0.828545,-0.800912,0.184654,-0.63675,-0.16696,0.240427,-0.513443,0.812664,0.744943,0.970612,0.00172899,-0.726378,-0.0985012,0.224232,0.16495,0.560077,-0.813112,0.112894,-0.0955366,0.0187107,0.913887,0.123076,0.550338,0.400334,-0.367816,0.198455,-0.983183,0.278976,0.714817,0.307911,0.812861,-0.403497,-0.784382,-0.161823,-0.120835,0.323172,0.583739,0.732924,-0.220603,-0.594121,0.935093,-0.216736,0.659318,-0.750417,-0.284773,-0.271496,0.491731,-0.712174,-0.763681,0.0781023,0.951666,0.734031,0.826912,0.57488,-0.361951,-0.0739728,0.91438,-0.391653,0.0193887,0.412634,-0.169813,0.471794,0.660792,-0.350906,-0.612644,0.347876,0.112573,-0.501126,0.456761,-0.109004,0.289352,-0.566504,0.585042,0.584934,0.923676,0.895312,-0.161036,-0.995895,0.0853141,-0.583368,-0.157612,0.234119,0.875043,0.430805,0.706102,0.423887,0.296828,-0.265607,0.338806,-0.15829,0.642516,0.355126,0.174447,-0.975015,0.869905,-0.145285,-0.484002,-0.475966,-0.67704,0.996452,-0.0685748,-0.851985,0.416498,0.791047,-0.211323,-0.302819,0.640735,-0.317908,-0.116586,-0.896382,-0.817317,-0.948837,-0.597427,0.975863,-0.971371,-0.124115,0.4339,-0.254671,0.298068,-0.349803,-0.73185,0.488106,-0.0495073,0.253969,0.168116,0.148772,0.889593,-0.512213,-0.165437,0.666773,-0.976304,-0.170024,0.905794,0.473908,-0.855725,-0.0413591,-0.508661,0.443453,0.842925,-0.144503,0.936699,-0.443935,-0.182996,0.803564,0.960386,-0.0323329,0.638181,-0.895684,-0.360502,0.0646258,-0.202449,-0.717228,0.970489,0.404608,-0.0861868,-0.879417,-0.866462,-0.938336,-0.799685,0.213464,-0.932344,-0.668236,0.751366,-0.22712,-0.407783,0.657463,0.0970092,-0.579109,-0.868866,-0.504041,0.926483,0.169978,-0.00842563,-0.530324,0.282745,0.0255867,0.287686,0.410417,-0.766576,-0.536566,-0.628288,0.69665,0.820713,-0.506162,-0.404114,0.640099,-0.956123,-0.576586,0.435502,-0.470676,-0.367062,-0.831765,-0.294942,0.518991,0.922338,0.337886,-0.67474,-0.725667,0.916684,0.39175,0.759081,0.496979,-0.200691,0.0417966,-0.687391,0.438773,0.287357,0.316636,-0.262311,-0.0755541,-0.442313,0.621378,0.670105,0.060982,0.944162,0.643442,-0.750684,-0.639973,0.217424,0.592823,0.929094,-0.239135,-0.41628,0.570893,-0.0798988,-0.917135,-0.749545,-0.982047,0.0626998,-0.977963,0.660401,0.470569,-0.0528868,-0.00138645,0.931065,-0.748519,0.304188,-0.266153,0.672524,-0.105179,-0.874749,-0.154355,-0.774656,-0.69654,0.433098,0.615897,-0.387919,-0.429779,0.650202,0.122306,-0.237727,0.626817,-0.227929,0.405916,0.483328,0.282047,-0.262206,0.784123,0.83125,-0.662272,0.702768,0.875814,-0.701221,0.553793,0.471795,0.769147,0.059668,-0.841617,-0.191179,-0.972471,-0.825361,0.779826,-0.917201,0.43272,0.10301,0.358771,0.793448,-0.0379954,-0.870112,0.600442,-0.990603,0.549151,0.512146,-0.795843,0.490091,0.372046,-0.549437,0.0964285,0.753047,-0.86284,-0.589688,0.178612,-0.720358,};
  /*
  for(i = 0; i < m; i++){
    M[i] = (double *)malloc(sizeof(double)*m);
  }
  myFile = fopen("./data/L2ZDT3.txt", "r");

  for (i=0; i<m; i++)
  {
    for (j=0; j<m; j++)
    {
      if(!	fscanf(myFile, "%lf\n", &M[i][j]))
	     break;
    }
  }
  fclose(myFile);
*/
  for(i = 0; i < m; i++){
    if(x[i] >= 0 && x[i] <= 1)
      continue;
    else{
      printf("Invalid input range\n");
      return;
    }
  }

  y = (double*)malloc(sizeof(double)*m);
  for(i = 0; i < m; i++){
    y[i] = 0;
    for(j = 0; j < m; j++){
      y[i] = y[i] + ( M[i*m+j] * x[j] );
    }
  }

  ff1 = y[0]*y[0];
  gx = 0;

  for(i = 1; i < m; i++){
    gx = gx + y[i]*y[i];
  }
  gx = gx * (9.0/(m-1));
  gx = gx + 1;

  h = 1-sqrt(ff1/gx)-(ff1/gx)*sin(10*pi*ff1);

  f[0] = ff1;
  f[1] = gx*h;

  return;
}

/*
*
*   As described by Y. Jin, M. Olhofer and B. Sendhoff. "Dynamic weighted
*   aggregation for evolutionary multi-objective optimization: Why does it
*   work and how?", in Proceedings of Genetic and Evolutionary Computation
*   Conference, pp.1042-1049, San Francisco, USA, 2001.
*
*   Test function 2, F2.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
 */

void Jin2(double* f, double* x){
  double n = 2;
  double gx = 0;
  int i;
  for(i = 0; i < n; i++){
    if(x[i] >= 0 && x[i] <= 1)
      continue;
    else{
      printf("Invalid input range\n");
      return;
    }
  }

  for(i = 1; i < n; i++)
    gx = gx + x[i];

  gx=9.0*gx;
  gx = gx / (n-1);
  gx = gx + 1;

  f[1] = x[0];
  f[0] = gx * (1 - sqrt(x[0]/gx));

  return;
}

/*
*
*   As described by Huband et al. in "A Scalable Multi-objective Test Problem
*   Toolkit", in C. A. Coello Coello et al. (Eds.): EMO 2005, LNCS 3410,
*   pp. 280–295, 2005, Springer-Verlag Berlin Heidelberg 2005.
*
*   Example I5.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
 */

void I5(double* f, double* zz){
    int M = 3;
    double k = 4;
    double l = 4;
    int n = k+l;

    double pi2 = 2*atan(1);
    double AA = 0.98 / 49.98;
		double BB = 0.02;
		int CC = 50;

    double *S, *A;
    double *zmax, *x, *y, *t1, *t2, *t3, *h, *w, *r_sum, *z;
    int i, j, ii, jj;
    double sum1, sum2;

    z = (double*)malloc(sizeof(double)*(n+1));
    for(i=n; i>=1; --i)
      z[i] = zz[i-1];

    S = (double*)malloc(sizeof(double)*(M+1));
    A = (double*)malloc(sizeof(double)*(M));

    for(i=1; i<=M-1; ++i){
      S[i] = 1;
      A[i] = 1;
    }
    S[M] = 1;

    zmax = (double*)malloc(sizeof(double)*(n+1));
    y = (double*)malloc(sizeof(double)*(n+1));
    t1 = (double*)malloc(sizeof(double)*(n+1));
    t2 = (double*)malloc(sizeof(double)*(n+1));
    t3 = (double*)malloc(sizeof(double)*(M+1));
    x = (double*)malloc(sizeof(double)*(M+1));
    h = (double*)malloc(sizeof(double)*(M+1));

    zmax = (double*)malloc(sizeof(double)*(n+1));
		for(i=1 ; i<=n; i++)
			zmax[i] = 1;

		y = (double*)malloc(sizeof(double)*(n+1));
		for(i = 1; i <= n; i++)
			y[i] = z[i] / zmax[i];

		w = (double*)malloc(sizeof(double)*(n+1));
		for(i = 1; i <= n; i++)
			w[i] = 1;

		r_sum = (double*)malloc(sizeof(double)*(n+1));
		for(i = 2; i <= n; i++){
			sum1 = 0;
			sum2 = 0;
			for(j = 1; j <= i-1; j++){
				sum1 = sum1 + w[j]*y[j];
				sum2 = sum2 + w[j];
			}
			r_sum[i] = sum1 / sum2;
		}

		t1 = (double*)malloc(sizeof(double)*(n+1));
		for(i = 1; i <= n; i++){
			if(i < 2)
				t1[i] = y[i];
			else
				t1[i] = pow(y[i], (BB+(CC-BB)*(AA-(1-2*r_sum[i])*fabs(floor(0.5-r_sum[i])+AA))));
		}

    for(i=1; i<=n; ++i){
      if(i<=k)
        t2[i] = t1[i];
      else
        t2[i] = ((fabs(t1[i]-0.35)/fabs(floor(0.35-t1[i])+0.35)));
    }

    for(i=1; i<=n; ++i) {
        if( z[i]<0 || z[i]>zmax[i] ) {
            printf("%lf, Invalid input range\n", z[i]);
            return;
        }
    }

    for(i=1; i<=M; ++i){
      if(i<=(M-1)){
        sum1 = 0;

        for(ii=((i-1)*k/(M-1)+1); ii<=(i)*k/(M-1); ++ii){
          sum2 = 0;
          for(jj=0; jj<=(k/(M-1)-2); ++jj){
            sum2 += fabs(t2[ii] - t2[(int)(((i-1)*k/(M-1)+1)+((((int)(ii+jj-((i-1)*k/(M-1)+1)+1))) % (((int)(((i)*k/(M-1))-((i-1)*k/(M-1)+1)+1)))))]);
          }
          sum1 += (t2[ii] + sum2);
        }
        t3[i] = sum1/(((i*k/(M-1))-((i-1)*k/(M-1)+1)+1)/(k/(M-1))*ceil(k/(M-1)/2)*(1+2*k/(M-1)-2*ceil(k/(M-1)/2)));
      }
      else{
        sum1 = 0;
        sum2 = 0;
        for(ii=k+1; ii<=n; ++ii){
          sum2 = 0;
          for(jj=0; jj<=l-2; ++jj){
            sum2 += fabs(t2[ii] - t2[(int)(k+1+(((int)(ii+jj-(k+1)+1)) % ((int)(n-k))))]);
          }
          sum1 += t2[ii] + sum2;
        }
        t3[i] = sum1/(((n-k)/l)*ceil(l/2)*(1+2*l-2*ceil(l/2)));
      }
    }
    for(i=1; i<=M; ++i){
      if(i<=M-1){
        if(t3[M]>A[i])
          x[i] = t3[M]*(t3[i]-0.5)+0.5;
        else
          x[i] = A[i]*(t3[i]-0.5)+0.5;
      }
      else
        x[i] = t3[M];
    }

    for(i=0; i<=M; ++i){
      if(i==1){
        h[i] = 1;
        for(j=1; j<=M-1; ++j){
          h[i] *= sin(x[j]*pi2);
        }
      }
      else if(i<=M-1){
        h[i] = 1;
        for(j=1; j<=M-i; ++j)
          h[i] *= sin(x[j]*pi2);
        h[i] *= cos(x[M-i+1]*pi2);
      }
      else{
        h[i] = cos(x[1]*pi2);
      }
    }
    for(i=1; i<=M; ++i){
      f[i-1] = x[M]+S[i]*h[i];
    }

    return;
}

/*
*
*   As described by Huband et al. in "A Scalable Multi-objective Test Problem
*   Toolkit", in C. A. Coello Coello et al. (Eds.): EMO 2005, LNCS 3410,
*   pp. 280–295, 2005, Springer-Verlag Berlin Heidelberg 2005.
*
*   Example I1.
*
*   This file is part of a collection of problems developed for
*   the BMOBench platform and based on
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
 */

void I1(double* f, double* zz){
    int M = 3;
    double k = 4;
    double l = 4;
    int n = k+l;

    double pi2 = 2*atan(1);

    double *S, *A;
    double *zmax, *x, *y, *t1, *t2, *t3, *h, *w, *z;
    int i, j;
    double sum1, sum2;

    z = (double*)malloc(sizeof(double)*(n+1));
    for(i=n; i>=1; --i)
      z[i] = zz[i-1];

    S = (double*)malloc(sizeof(double)*(M+1));
    A = (double*)malloc(sizeof(double)*(M));

    for(i=1; i<=M-1; ++i){
      S[i] = 1;
      A[i] = 1;
    }
    S[M] = 1;

    zmax = (double*)malloc(sizeof(double)*(n+1));
    y = (double*)malloc(sizeof(double)*(n+1));
    t1 = (double*)malloc(sizeof(double)*(n+1));
    t2 = (double*)malloc(sizeof(double)*(n+1));
    t3 = (double*)malloc(sizeof(double)*(M+1));
    x = (double*)malloc(sizeof(double)*(M+1));
    h = (double*)malloc(sizeof(double)*(M+1));

    for(i=1; i<=n; ++i){
      zmax[i] = 1;
      y[i] = z[i]/zmax[i];
      t1[i] = y[i];
      if(i<=k)
        t2[i] = t1[i];
      else
        t2[i] = ((fabs(t1[i]-0.35)/fabs(floor(0.35-t1[i])+0.35)));
    }

    for(i=1; i<=n; ++i) {
        if( z[i]<0 || z[i]>zmax[i] ) {
            printf("%lf, Invalid input range\n", z[i]);
            return;
        }
    }

    w = (double*)malloc(sizeof(double)*(n+1));
    for(i = 0; i <= n; i++)
    	w[i] = 1;

    for(i = 1; i <= M; i++){
    	if(i <= M-1){
    		sum1 = 0;
    		for(j = ((i-1)*k/(M-1)+1); j <= (i*k/(M-1)); j++)
    			sum1 = sum1 + (w[j]*t2[j]);
    		sum2 = 0;
    		for(j = ((i-1)*k/(M-1)+1); j <= (i*k/(M-1)); j++)
    			sum2 = sum2 + w[j];
    		t3[i] = sum1 / sum2;
    	}else{
    		sum1 = 0;
    		for(j = k+1; j <= n; j++)
    			sum1 = sum1 + (w[j]*t2[j]);
    		sum2 = 0;
    		for(j = k+1; j <= n; j++)
    			sum2 = sum2 + w[j];
    		t3[i] = sum1 / sum2;
    	}
    }
    for(i=1; i<=M; ++i){
      if(i<=M-1){
        if(t3[M]>A[i])
          x[i] = t3[M]*(t3[i]-0.5)+0.5;
        else
          x[i] = A[i]*(t3[i]-0.5)+0.5;
      }
      else
        x[i] = t3[M];
    }

    for(i=0; i<=M; ++i){
      if(i==1){
        h[i] = 1;
        for(j=1; j<=M-1; ++j){
          h[i] *= sin(x[j]*pi2);
        }
      }
      else if(i<=M-1){
        h[i] = 1;
        for(j=1; j<=M-i; ++j)
          h[i] *= sin(x[j]*pi2);
        h[i] *= cos(x[M-i+1]*pi2);
      }
      else{
        h[i] = cos(x[1]*pi2);
      }
    }
    for(i=1; i<=M; ++i){
      f[i-1] = x[M]+S[i]*h[i];
    }

    return;
}

/*
 *
 *   As described by Huband et al. in "A review of multiobjective test problems
 *   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
 *   Computing 10(5): 477-506, 2006.
 *
 *   Example FES1, see the previous cited paper for the original reference.
 *
 *   In the above paper the number of variables was left undefined.
 *   We selected n=10 as default.
 *
 *   This file is part of a collection of problems developed for
 *   the BMOBench platform and based on
 *   derivative-free multiobjective optimization in
 *   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *   Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *   Written by the authors in June 1, 2010.
 *
 */

void FES1(double* f, double* x){
  double pi = 4*atan(1);
  double n = 10;
  int i;

  for(i = 0; i < n; i++){
    if(x[i] >= 0 && x[i] <= 1)
      continue;
    else{
      printf("Invalid input range\n");
      return;
    }
  }

  f[0] = 0;
  f[1] = 0;

  for(i = 0; i < n; i++){
    f[0] = f[0] + pow(fabs(x[i]-exp(((i+1)/n)*((i+1)/n))/3),0.5);
  }

  for(i = 0; i < n; i++){
    f[1] = f[1] + pow( (x[i]-0.5*cos(10*pi*(i+1)/n)-0.5) ,2);
  }

  return;
}

/*
 *
 *   As described by K. Deb, L. Thiele, M. Laumanns and E. Zitzler in "Scalable
 *   Multi-Objective Optimization Test Problems", Congress on Evolutionary
 *   Computation (CEC’2002): 825-830, 2002.
 *
 *   Example DTLZ6.
 *
 *   This file is part of a collection of problems developed for
 *   the BMOBench platform and based on
 *   derivative-free multiobjective optimization in
 *   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *   Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *   Written by the authors in June 1, 2010.
 *
 */

void DTLZ6(double* f, double* x){
  double m = 3;
  double n = 22;
  double k = n - m + 1;
  double pi = 4*atan(1);
  double gx = 0;
  double ffM = 0;
  int i;

  for(i = 0; i < n; i++){
    if(x[i] >= 0 && x[i] <= 1)
      continue;
    else{
      printf("Invalid input range\n");
      return;
    }
  }

  for(i = m-1; i < n; i++){
    gx = gx + x[i];
  }
  gx = gx*(9/k);
  gx += 1;

  for(i = 0; i < m-1; i++){
    ffM = ffM + (x[i]/(1+gx)*(1+sin(3*pi*x[i])));
  }
  ffM = m - ffM;
  ffM = ffM * (1 + gx);

  for(i = 1; i < m; i++)
    f[i] = x[i-1];
  f[0] = ffM;

  return;
}

/*
 *
 *   As described by K. Deb, L. Thiele, M. Laumanns and E. Zitzler in "Scalable
 *   Multi-Objective Optimization Test Problems", Congress on Evolutionary
 *   Computation (CEC’2002): 825-830, 2002.
 *
 *   Example DTLZ4.
 *
 *   This file is part of a collection of problems developed for
 *   the BMOBench platform and based on
 *   derivative-free multiobjective optimization in
 *   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *   Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *   Written by the authors in June 1, 2010.
 *
 */

void DTLZ4(double* f, double* x){
  int m = 3;
  int n = 12;
  double pi = 4*atan(1);
  int alpha = 100;
	double *y;
  int i, j;
  double gx = 0;

  for(i = 0; i < n; i++){
    if(x[i] >= 0 && x[i] <= 1)
      continue;
    else{
      printf("Invalid input range\n");
      return;
    }
  }

  y = (double*)malloc(sizeof(double)*n);
  for(i = 0; i < n; i++){
    y[i] = pow(x[i], alpha);
  }

  for(i = m-1; i < n; i++)
    gx = gx + ( (y[i] - 0.5)*(y[i] - 0.5) );

  f[0] = (1+gx);
  for(i = 0; i < m-1; i++){
    f[0] = f[0]*(cos(0.5*pi*y[i]));
  }

  for(i = 1; i < m; i++){
    f[i] = (1+gx)*(sin(0.5*pi*y[m-i-1]));
    for(j = 0; j < m-i-1; j++){
      f[i] = f[i] * (cos(0.5*pi*y[j]));
    }
  }

  return;
}

/*
 *
 *   As described by K. Deb, L. Thiele, M. Laumanns and E. Zitzler in "Scalable
 *   Multi-Objective Optimization Test Problems", Congress on Evolutionary
 *   Computation (CEC’2002): 825-830, 2002.
 *
 *   Example DTLZ2.
 *
 *   This file is part of a collection of problems developed for
 *   the BMOBench platform and based on
 *   derivative-free multiobjective optimization in
 *   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *   Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *   Written by the authors in June 1, 2010.
 *
 */

void DTLZ2(double* f, double* x){
  int m = 3;
  int n = 12;

  double pi = 4*atan(1);
  double gx = 0;
  int i, j;

  for(i = 0; i < n; i++){
    if(x[i] >= 0 && x[i] <= 1)
      continue;
    else{
      printf("Invalid input range\n");
      return;
    }
  }

  for(i = m-1; i < n; i++){
    gx = gx + ( (x[i] - 0.5) * (x[i] - 0.5) );
  }

  f[0] = (1+gx);
  for(i = 0; i < m - 1; i++ ){
    f[0] = f[0]*(cos(0.5*pi*x[i]));
  }

  for(i = 1; i < m ; i++){
    f[i] = (1+gx)*(sin(0.5*pi*x[m-i-1]));
    for(j = 0; j < m-i-1; j++){
      f[i] = f[i] * (cos(0.5*pi*x[j]));
    }
  }

  return;
}

/*
 *
 *   As described by Huband et al. in "A review of multiobjective test problems
 *   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
 *   Computing 10(5): 477-506, 2006.
 *
 *   Example DG01, see the previous cited paper for the original reference.
 *
 *   This file is part of a collection of problems developed for
 *   the BMOBench platform and based on
 *   derivative-free multiobjective optimization in
 *   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *   Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *   Written by the authors in June 1, 2010.
 *
 */

void DG01(double*f, double* x){
  if (x[0] >= -10 && x[0] <= 13);
  else{
    printf("Invalid input range\n");
    return;
  }

  f[0] = sin(x[0]);
  f[1] = sin(x[0] + 0.7);

  return;
}

/*
 *
 *   As described by K. Deb in "Multi-objective Genetic Algorithms: Problem
 *   Difficulties and Construction of Test Problems", Evolutionary Computation
 *   7(3): 205-230, 1999.
 *
 *   Example 5.1.2 (Non-convex local and convex global Pareto-optimal Front).
 *
 *   This file is part of a collection of problems developed for
 *   the BMOBench platform and based on
 *   derivative-free multiobjective optimization in
 *   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *   Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *   Written by the authors in June 1, 2010.
 *
 */


void Deb512c(double*f, double* x){
  double beta = 1;
  double ff1 = 4*x[0];
  double gx;
  double alpha;
  double h;

  if((x[0] >= 0 && x[0] <= 1) && (x[1] >= 0 && x[1] <= 1));
  else{
    printf("Invalid input range\n");
    return;
  }

  if(x[1] <= 0.4)
    gx = 4-3*exp(-pow(((x[1]-0.2)/0.02), 2));
  else
    gx = 4-2*exp(-pow(((x[1]-0.7)/0.2), 2));

  alpha = 0.25+3.75*(gx-1);
  if(ff1 <= beta*gx)
    h = (1-pow((ff1/(beta*gx)), alpha));
  else
    h = 0;

  f[0] = ff1;
  f[1] = gx*h;

  return;
}

/*
 *
 *   As described by K. Deb in "Multi-objective Genetic Algorithms: Problem
 *   Difficulties and Construction of Test Problems", Evolutionary Computation
 *   7(3): 205-230, 1999.
 *
 *   Example 4.1 (Multi-modal Multi-objective Problem).
 *
 *   This file is part of a collection of problems developed for
 *   the BMOBench platform and based on
 *   derivative-free multiobjective optimization in
 *   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *   Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *   Written by the authors in June 1, 2010.
 *
 */

void Deb41(double* f, double* x){
  double gx = 2-exp(-pow(((x[1]-0.2)/0.004), 2))-0.8*exp(-pow(((x[1]-0.6)/0.4), 2));

  if((x[0] >= 0.1 && x[0] <= 1) && (x[1] >= 0 && x[1] <= 1));
  else{
    printf("Invalid input range\n");
    return;
  }

  f[1] = x[0];
  f[0] = gx/x[0];

  return;
}

/*
 *   As described by Huband et al. in "A review of multiobjective test problems
 *   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
 *   Computing 10(5): 477-506, 2006.
 *
 *   Example BK1, see the previous cited paper for the original reference.
 *
 *   This file is part of a collection of problems developed for
 *   the BMOBench platform and based on
 *   derivative-free multiobjective optimization in
 *   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *   Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *   Written by the authors in June 1, 2010.
 *
 */

void BK1(double* f, double* x){
  if((x[0] >= -5 && x[0] <= 10) || (x[1] >= -5 && x[1] <= 10));
  else{
    printf("Invalid input range\n");
    return;
  }
  f[0] = pow(x[0], 2) + pow(x[1], 2);
  f[1] = pow((x[0] - 5), 2) + pow((x[1] - 5), 2);

  return;
}

/*
 *
 *  As described by K. Deb in "Multi-objective Genetic Algorithms: Problem
 *  Difficulties and Construction of Test Problems", Evolutionary Computation
 *  7(3): 205-230, 1999.
 *
 *  Example 5.1.3 (Discontinuous Pareto-optimal Front).
 *
 *  This file is part of a collection of problems developed for
 *  derivative-free multiobjective optimization in
 *  A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
  */
void Deb513(double *f, double *x){
	int n;
	double alpha;
	double q;
	double pi;
	double ff1;
	double gx;
	double h;
  int i;
	n = 2;

	alpha = 2.0;
	q = 4.0;
	pi = 4.0*atan(1.0);

	ff1 = x[0];
	gx = 1.0+10.0*x[1];
	h = 1.0-pow((ff1/gx), alpha)-(ff1/gx)*sin(2.0*pi*q*ff1);

  for(i=0; i<n; ++i){
  	if (!(x[i]<=1.0 && x[i]>=0.0))
  	{
  		printf("Invalid input range.\n");
  		return;
  	}
  }

  f[0] = ff1;
  f[1] = gx*h;

  return;
}

/*
 *
 *  As described by Huband et al. in "A review of multiobjective test problems
 *  and a scalable test problem toolkit", IEEE Transactions on Evolutionary
 *  Computing 10(5): 477-506, 2006.
 *
 *  Example DPAM1, see the previous cited paper for the original reference.
 *
 *  In the above paper the number of variables was left undefined.
 *  We selected n=10 as default.
 *
 *  This file is part of a collection of problems developed for
 *  derivative-free multiobjective optimization in
 *  A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
  */
void DPAM1(double *f, double *x){
	int n;
	double pi;
	int i, j;
	//double **A;
	double *y;
	double gx;
	//FILE *myFile;
double A[100] = {0.218418,-0.620254,0.843784,0.914311,-0.788548,0.428212,0.103064,-0.47373,-0.300792,-0.185507,0.330423,0.151614,0.884043,-0.272951,-0.993822,0.511197,-0.0997948,-0.659756,0.575496,0.675617,0.180332,-0.593814,-0.492722,0.0646786,-0.666503,-0.945716,-0.334582,0.611894,0.281032,0.508749,-0.0265389,-0.920133,0.308861,-0.0437502,-0.374203,0.207359,-0.219433,0.914104,0.184408,0.520599,-0.88565,-0.375906,-0.708948,-0.37902,0.576578,0.0194674,-0.470262,0.572576,0.351245,-0.480477,0.238261,-0.1596,-0.827302,0.669248,0.494475,0.691715,-0.198585,0.0492812,0.959669,0.884086,-0.218632,-0.865161,-0.715997,0.220772,0.692356,0.646453,-0.401724,0.615443,-0.0601957,-0.748176,-0.207987,-0.865931,0.613732,-0.525712,-0.995728,0.389633,-0.064173,0.662131,-0.707048,-0.340423,0.60624,0.0951648,-0.160446,-0.394585,-0.167581,0.0679849,0.449799,0.733505,-0.00918638,0.00446808,0.404396,0.449996,0.162711,0.294454,-0.563345,-0.114993,0.549589,-0.775141,0.677726,0.610715};

	n = 10;
	pi = 4.0*atan(1);
		/*
    A = (double**)malloc(sizeof(double*)*n);
    for (i=0; i<n; i++)
	{
		A[i] = (double*)malloc(sizeof(double)*n);
	}

    myFile = fopen("./data/DPAM1_Avalues.txt", "r");

	for (i=0; i<n; i++)
	{
		for (j=0; j<n; j++)
		{
			fscanf(myFile, "%lf\n", &A[i][j]);
		}
	}
	*/
	y = (double*)malloc(sizeof(double)*n);
	for (i=0; i<n; i++)
	{
		y[i] = 0.0;
		for (j=0; j<n; j++)
		{
			y[i]+=A[i*n+j]*x[j];
		}
	}

	gx = 1 + 10*(n-1);
	for (i=1; i<n; i++)
	{
		gx+=(pow(y[i], 2.0)-10.0*cos(4.0*pi*y[i]));
	}

	for (i=0; i<n; i++)
	{
		if (x[i]>0.3 || x[i]<-0.3)
		{
			printf("Invalid input range.\n");
			return;
		}
		else
		{
			f[0] = y[0];
			f[1] = gx*exp(-y[0]/gx);
		}
	}
	//fclose(myFile);
	//free(A);
	free(y);
	return;
}

/*
 *
 *  As described by K. Deb, L. Thiele, M. Laumanns and E. Zitzler in "Scalable
 *  Multi-Objective Optimization Test Problems", Congress on Evolutionary
 *  Computation (CEC-2002): 825-830, 2002.
 *
 *  Example DTLZ2 with M=2 and n=2.
 *
 *  This file is part of a collection of problems developed for
 *  derivative-free multiobjective optimization in
 *  A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
  */
void DTLZ2n2(double *f, double *x){
	int n;
	int M;
	double pi;
	int i, j;
	double gx;
	double ff1;
	double *ff;

	n = 2;
	M = 2;
	pi =  4.0*atan(1);

	gx = 0.0;
	for (i=M-1; i<n; i++)
	{
		gx+=pow((x[i]-0.5), 2.0);
	}

	ff1 = 1.0+gx;
	for (j=0; j<M-1; j++)
	{
		ff1*=cos(0.5*pi*x[j]);
	}

	ff = (double*)malloc(sizeof(double)*M);
	ff[0] = 0;
	for (i=1; i<M; i++)
	{
		ff[i] = (1.0+gx)*(sin(0.5*pi*x[M-i-1]));
		for (j=0; j<M-i-1; j++)
		{
			ff[i]*=(cos(0.5*pi*x[j]));
		}
	}

	for (i=0; i<n; i++)
	{
		if ((x[i]>1.0 || x[i]<0.0) || M<2 || n<M)
		{
			printf("Invalid input range.\n");
			return;
		}
		else
		{
			f[0] = ff1;
			for (i=1; i<M; i++)
			{
				f[i]=ff[i];
			}
		}
	}
	return;
}

/*
 *
 *  As described by K. Deb, L. Thiele, M. Laumanns and E. Zitzler in "Scalable
 *  Multi-Objective Optimization Test Problems", Congress on Evolutionary
 *  Computation (CEC-2002): 825-830, 2002.
 *
 *  Example DTLZ4 with M=2 and n=2.
 *
 *  This file is part of a collection of problems developed for
 *  derivative-free multiobjective optimization in
 *  A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
  */
void DTLZ4n2(double *f, double *x){
	int M;
	int n;
	double pi;
	double alpha;
	double *y;
	int i, j;
	double gx;
	double ff1;
	double *ff;

	M = 2;
	n = 2;
	alpha = 100.0;
	pi = 4.0*atan(1);

	y = (double*)malloc(sizeof(double)*n);
	for (i=0; i<n; i++)
	{
		y[i] = pow(x[i], alpha);
	}

	gx = 0.0;
	for (i=M-1; i<n; i++)
	{
		gx+=pow((y[i]-0.5), 2.0);
	}

	ff1 = 1.0+gx;
	for (j=0; j<M-1; j++)
	{
		ff1*=cos(0.5*pi*y[j]);
	}

	ff = (double*)malloc(sizeof(double)*M);
	ff[0] = 0.0;
	for (i=1; i<M; i++)
	{
		ff[i] = (1.0+gx)*(sin(0.5*pi*y[M-i-1]));
		for (j=0; j<M-i-1; j++)
		{
			ff[i]*=(cos(0.5*pi*y[j]));
		}
	}

	for (i=0; i<n; i++)
	{
		if ((x[i]>1.0 || x[i]<0.0) || M<2 || n<M)
		{
			printf("Invalid input range.\n");
			return;
		}
		else
		{
			f[0] = ff1;
			for (i=1; i<M; i++)
			{
				f[i]=ff[i];
			}
		}
	}
	return;
}

/*
 *
 *  As described by K. Deb, L. Thiele, M. Laumanns and E. Zitzler in "Scalable
 *  Multi-Objective Optimization Test Problems", Congress on Evolutionary
 *  Computation (CEC-2002): 825-830, 2002.
 *
 *  Example DTLZ6 with M=2 and n=2.
 *
 *  This file is part of a collection of problems developed for
 *  derivative-free multiobjective optimization in
 *  A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
  */
void DTLZ6n2(double *f, double *x){
	int n;
	int k;
	int M;
	double pi;
	double gx;
	int i;
	double ffM;
	double sum;

	n = 2;
	M = 2;
	k = n-M+1;
	pi = 4.0*atan(1);

	gx = 1.0;
	for (i=M-1; i<n; i++)
	{
		gx+=9.0/k*x[i];
	}

	sum = 0.0;
	for (i=0; i<M-1; i++)
	{
		sum+=(x[i]/(1.0+gx)*(1+sin(3*pi*x[i])));
	}
	ffM = (1.0+gx)*(M-sum);

	for (i=0; i<n; i++)
	{
		if ((x[i]>1.0 || x[i]<0.0) || M<2 || n<M)
		{
			printf("Invalid input range.\n");
			return;
		}
		else
		{
			for (i=0; i<M; i++)
			{
				f[i] = x[i-1];
			}
			f[0] = ffM;
		}
	}
	return;
}

/*
 *
 *  As described by Huband et al. in "A review of multiobjective test problems
 *  and a scalable test problem toolkit", IEEE Transactions on Evolutionary
 *  Computing 10(5): 477-506, 2006.
 *
 *  Example FES2, see the previous cited paper for the original reference.
 *
 *  In the above paper the number of variables was left undefined.
 *  We selected n=10 as default.
 *
 *  This file is part of a collection of problems developed for
 *  derivative-free multiobjective optimization in
 *  A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
  */
void FES2(double *f, double *x){
	double pi;
	int n;
	double f1, f2, f3;
	double sin_pow, cos_pow;
	int i;

	n = 10;
	pi = 4.0*atan(1);
	f1 = 0.0;
	f2 = 0.0;
	f3 = 0.0;

	for (i=0; i<n; i++)
	{
		f1+=pow((x[i]-0.5*cos(10.0*pi*(i+1)/n)-0.5), 2.0);
		sin_pow = pow(sin(i), 2.0);
		cos_pow = pow(cos(i), 2.0);
		f2+=pow(fabs(x[i]-sin_pow*cos_pow), 0.5);
		f3+=pow(fabs(x[i]-0.25*cos(i)*cos(2*i)-0.5), 0.5);
	}

	for (i=0; i<n; i++)
	{
		if ((x[i]>1.0 || x[i]<0.0))
		{
			printf("Invalid input range.\n");
			return;
		}
		else
		{
			f[0]=f1;
			f[1]=f2;
			f[2]=f3;
		}
	}
	return;
}

/*
 *
 *  As described by Huband et al. in "A review of multiobjective test problems
 *  and a scalable test problem toolkit", IEEE Transactions on Evolutionary
 *  Computing 10(5): 477-506, 2006.
 *
 *  Example LRS1, see the previous cited paper for the original reference.
 *
 *  This file is part of a collection of problems developed for
 *  derivative-free multiobjective optimization in
 *  A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
  */
void LRS1(double *f, double *x){
  int n, i;
  n = 2;

  for(i=0; i<n; ++i){
    if (!(x[i]<=50 && x[i]>=-50)){
    	printf("Invalid input range.\n");
    	return;
    }
  }

	f[0] = pow(x[0], 2)+pow(x[1], 2);
	f[1] = pow((x[0]+2), 2)+pow(x[1], 2);

	return;
}

/*
 *
 *  As described by Huband et al. in "A review of multiobjective test problems
 *  and a scalable test problem toolkit", IEEE Transactions on Evolutionary
 *  Computing 10(5): 477-506, 2006.
 *
 *  Example MHHM1, see the previous cited paper for the original reference.
 *
 *  This file is part of a collection of problems developed for
 *  derivative-free multiobjective optimization in
 *  A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
  */
void MHHM1(double *f, double *x){
    int n,i;
    n = 1;

    for(i=0; i<n; ++i){
      if(!(x[i]<=1 && x[i]>=0)){
	       printf("Invalid input range.\n");
         return;
       }
     }

     f[0] = pow((x[0]-0.8), 2);
     f[1] = pow((x[0]-0.85), 2);
     f[2] = pow((x[0]-0.9), 2);

     return;
}

/*
 *
 *  As described by Huband et al. in "A review of multiobjective test problems
 *  and a scalable test problem toolkit", IEEE Transactions on Evolutionary
 *  Computing 10(5): 477-506, 2006.
 *
 *  Example MHHM2, see the previous cited paper for the original reference.
 *
 *  This file is part of a collection of problems developed for
 *  derivative-free multiobjective optimization in
 *  A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
  */
void MHHM2(double *f, double *x){
    int n, i;

    n = 2;

    for(i=0; i<n; ++i){
      if(!(x[i]<=1 && x[i]>=0)){
		      printf("Invalid input range.\n");
          return;
      }
    }

    f[0] = pow((x[0]-0.8), 2)+pow((x[1]-0.6), 2);
    f[1] = pow((x[0]-0.85), 2)+pow((x[1]-0.7), 2);
    f[2] = pow((x[0]-0.9), 2)+pow((x[1]-0.6), 2);

  	return;
}

/*
 *
 *  As described by Huband et al. in "A review of multiobjective test problems
 *  and a scalable test problem toolkit", IEEE Transactions on Evolutionary
 *  Computing 10(5): 477-506, 2006.
 *
 *  Example WFG1.
 *
 *  This file is part of a collection of problems developed for
 *  derivative-free multiobjective optimization in
 *  A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
  */
void WFG1(double* f, double* zz){
    int M;
    int k;
    int l;
    int n;
    int i, j;
    double pi;
    double pi2;
    double *S;
    double *A;
    double *zmax;
    double *y;
    double *t1;
    double *t2;
    double AAA;
    double *t3;
    double AA, BB, CC;
    double *w;
    double *t4;
    double *x;
    double sum_wj1, sum_wj2, sum_wt31, sum_wt32;
    double *h;
    double alpha, AAAA;
    double *z;

    M = 3;
    k = 4;
    l = 4;
    n = k+l;
    pi = 4.0*atan(1);
    pi2 = 2.0*atan(1);
    AAA = 0.02;

    z = (double *)malloc(sizeof(double)*(n+1));
    for(i=1; i<=n; i++)
    {
        z[i] = zz[i-1];
    }

    S = (double*)malloc(sizeof(double)*(M+1));
    for(i = 1; i <= M; i++)
    {
        S[i] = 2.0*(double)i;
    }

    A = (double*)malloc(sizeof(double)*(M));
    for(i = 1; i < M; i++)
    {
        A[i] = 1.0;
    }

    zmax = (double*)malloc(sizeof(double)*(n+1));
    for(i = 1; i <= n; i++)
    {
        zmax[i] = 2.0*(double)i;
    }

    for(i = 1; i <= n; i++)
    {
        if(z[i] >= 0 && z[i] <= zmax[i])
        {
            continue;
        }
        else
        {
            printf("Invalid input range\n");
            return;
        }
    }

    y = (double*)malloc(sizeof(double)*(n+1));
    for(i = 1; i <= n; i++)
    {
        y[i] = z[i]/zmax[i];
    }

    t1 = (double*)malloc(sizeof(double)*(n+1));
    for(i = 1; i <= n; i++)
    {
        if(i <= k)
        {
            t1[i] = y[i];
        }
        else
        {
            t1[i] = (double)(fabs(y[i]-0.35))/(double)(fabs(floor(0.35-y[i])+0.35));
        }
    }

    AA = 0.8;
    BB = 0.75;
    CC = 0.85;
    t2 = (double*)malloc(sizeof(double)*(n+1));
    for (i=1; i<=n; i++)
    {
        if(i <= k)
        {
            t2[i] = t1[i];
        }
        else
        {
            if (0.0 < floor(t1[i]-BB))
            {
                if (0.0 < floor(CC-t1[i]))
                {
                    t2[i] =  AA;
                }
                else
                {
                    t2[i] =  AA - (double)(floor(CC-t1[i])*(1.0-AA)*(t1[i]-CC))/(double)(1.0-CC);
                }
            }
            else
            {
                if (0.0 < floor(CC-t1[i]))
                {
                    t2[i] =  AA + (double)(floor(t1[i]-BB)*(AA*(BB-t1[i])))/(double)BB;
                }
                else
                {
                    t2[i] =  AA + (double)(floor(t1[i]-BB)*(AA*(BB-t1[i])))/(double)BB-(double)(floor(CC-t1[i])*(1.0-AA)*(t1[i]-CC))/(double)(1.0-CC);
                }
            }
        }
    }

    t3 = (double*)malloc(sizeof(double)*(n+1));
    for (i=1; i<=n; i++)
    {
        t3[i] = pow(t2[i], AAA);
    }

    w = (double*)malloc(sizeof(double)*(n+1));
    for(i = 1; i <= n; i++)
    {
        w[i] = 2.0*(double)i;
    }


    t4 = (double*)malloc(sizeof(double)*(M+1));
    for (i=1; i<=M; i++)
    {
        sum_wt31 = 0.0;
        sum_wj1 = 0.0;
        sum_wt32 = 0.0;
        sum_wj2 = 0.0;
        if (i<=M-1)
        {
            for (j=((i-1)*k/(M-1)+1); j<=(i*k/(M-1)); j++)
            {
                sum_wj1+=w[j];
                sum_wt31+=w[j]*t3[j];
            }
            t4[i] = sum_wt31/sum_wj1;
        }
        else
        {
            for (j=k+1; j<=n; j++)
            {
                sum_wt32+=w[j]*t3[j];
                sum_wj2+=w[j];
            }
            t4[i] = sum_wt32/sum_wj2;
        }
    }

    x = (double*)malloc(sizeof(double)*(M+1));
    for(i = 1; i <= M; i++)
    {
        if(i <= M-1)
        {
            if (t4[M] > A[i])
                x[i] = t4[M]*(t4[i]-0.5)+0.5;
            else
                x[i] = A[i]*(t4[i]-0.5)+0.5;
        }
        else
        {
            x[i] = t4[M];
        }
    }

    alpha = 1.0;
    AAAA = 5.0;
    h = (double*)malloc(sizeof(double)*(M+1));
    for(i = 1; i <= M; i++)
    {
        if(i == 1)
        {
            h[i] = 1.0;
            for(j=1; j<=M-1; j++)
            {
                h[i] = h[i]*(1.0-cos(x[j]*pi2));
            }
        }
        else if(i<=(M-1))
        {
            h[i] = (1.0-sin(x[M-i+1]*pi2));
            for(j = 1; j <= M-i; j++)
            {
                h[i] = h[i]*(1.0-cos(x[j]*pi2));
            }
        }
        else
        {
            h[i] = pow((1.0-x[1]-(cos(2.0*AAAA*pi*x[1]+pi2))/(2.0*AAAA*pi)), alpha);
        }
    }

    for(i = 0; i < M; i++)
    {
        f[i] = x[M]+S[i+1]*h[i+1];
    }

    return;
}

/*
 *
 *  As described by Huband et al. in "A review of multiobjective test problems
 *  and a scalable test problem toolkit", IEEE Transactions on Evolutionary
 *  Computing 10(5): 477-506, 2006.
 *
 *  Example WFG2.
 *
 *  This file is part of a collection of problems developed for
 *  derivative-free multiobjective optimization in
 *  A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
  */
void WFG2(double* f, double* zz){
    int M;
    int k;
    int l;
    int n;
    int i, j;
    double pi;
    double pi2;
    double *S;
    double *A;
    double *zmax;
    double *y;
    double *t1;
    int ii, jj;
    double AA;
    double *t2;
    double sum1, sum2;
    double *w;
    double *t3;
    double sum_wt21, sum_wj1, sum_wt22, sum_wj2;
    double *x;
    double alpha;
    double beta;
    double AAAA;
    double *h;
	double *z;

    M = 3;
    k = 4;
    l = 4;
    n = k+l;
    AA = 2.0;
    pi = 4.0*atan(1);
    pi2 = 2.0*atan(1);
    sum_wj1 = 0.0;
    sum_wj2 = 0.0;
    sum_wt21 = 0.0;
    sum_wt22 = 0.0;
    alpha = 1.0;
    beta = 1.0;
    AAAA = 5.0;
	z = (double *)malloc(sizeof(double)*(n+1));

	for(i=1; i<=n; i++)
    {
        z[i] = zz[i-1];
    }

    S = (double*)malloc(sizeof(double)*(M+1));
    for(i = 1; i <= M; i++)
    {
        S[i] = 2.0*(double)i;
    }

    A = (double*)malloc(sizeof(double)*(M));
    for(i = 1; i <= M-1; i++)
    {
        A[i] = 1.0;
    }

    zmax = (double*)malloc(sizeof(double)*(n+1));
    for(i = 1; i <= n; i++)
    {
        zmax[i] = 2.0*(double)i;
    }

    for(i = 1; i <= n; i++)
    {
        if(z[i] >= 0 && z[i] <= zmax[i])
        {
            continue;
        }
        else
        {
            printf("Invalid input range\n");
            exit(1);
        }
    }

    y = (double*)malloc(sizeof(double)*(n+1));
    for(i = 1; i <= n; i++)
    {
        y[i] = z[i]/zmax[i];
    }

    t1 = (double*)malloc(sizeof(double)*(n+1));
    for(i=1; i<=n; i++)
    {
        if(i <= k)
        {
            t1[i] = y[i];
        }
        else
        {
            t1[i] = fabs(y[i]-0.35)/fabs(floor(0.35-y[i])+0.35);
        }
    }

    t2 = (double*)malloc(sizeof(double)*(k+l/2+1));
    for(i = 1; i <= k+l/2; i++)
    {
        if(i <= k)
        {
            t2[i] = t1[i];
        }
        else
        {
            sum1 = 0;
            for(ii = (k+2*(i-k)-1); ii <=(k+2*(i-k)); ii++)
            {
                sum2 = 0;
                for(jj = 0; jj <= AA-2; jj++)
                {
                    sum2+=fabs(t1[ii]-t1[(int)((k+2*(i-k)-1)+fmod((ii+jj-(k+2*(i-k)-1)+1), (k+2*(i-k)-(k+2*(i-k)-1)+1)))]);
                }
                sum1+=(t1[ii]+sum2)/((k+2*(i-k)-(k+2*(i-k)-1)+1)/AA*ceil(AA/2)*(1+2*AA-2*ceil(AA/2)));
            }
            t2[i] = sum1;
        }
    }

    w = (double*)malloc(sizeof(double)*(n+1));
    for(i = 1; i <= n; i++)
    {
        w[i] = 1.0;
    }

    t3 = (double*)malloc(sizeof(double)*(M+1));
    for (i=1; i<=M; i++)
    {
        if (i<=M-1)
        {
			sum_wj1 = 0.0;
			sum_wt21 = 0.0;
            for (j=((i-1)*k/(M-1)+1); j<=(i*k/(M-1)); j++)
            {
                sum_wj1+=w[j];
                sum_wt21+=w[j]*t2[j];
            }
            t3[i] = sum_wt21/sum_wj1;
        }
        else
        {
			sum_wt22 = 0.0;
			sum_wj2 = 0.0;
            for (j=k+1; j<=k+l/2; j++)
            {
                sum_wt22+=w[j]*t2[j];
                sum_wj2+=w[j];
            }
            t3[i] = sum_wt22/sum_wj2;
        }
    }

    x = (double*)malloc(sizeof(double)*(M+1));
    for(i = 1; i<=M; i++)
    {
        if(i<=M-1)
        {
            if (t3[M] > A[i])
                x[i] = t3[M]*(t3[i]-0.5)+0.5;
            else
                x[i] = A[i]*(t3[i]-0.5)+0.5;
        }
        else
        {
            x[i] = t3[M];
        }
    }

    h = (double*)malloc(sizeof(double)*(M+1));
    for(i = 1; i <= M; i++)
    {
        if(i == 1)
        {
            h[i] = 1.0;
            for(j=1; j<=M-1; j++)
            {
                h[i] = h[i]*(1.0-cos(x[j]*pi2));
            }
        }
        else if(i<=(M-1))
        {
            h[i] = (1.0-sin(x[M-i+1]*pi2));
            for(j = 1; j <= M-i; j++)
            {
                h[i] = h[i]*(1.0-cos(x[j]*pi2));
            }
        }
        else
        {
            h[i] = 1.0-pow((x[1]), alpha)*pow(cos(AAAA*pow((x[1]), beta)*pi), 2);
        }
    }

    for(i = 0; i < M; i++)
    {
        f[i] = x[M]+S[i+1]*h[i+1];
    }

    return;
}

/*
 *
 *  As described by Huband et al. in "A review of multiobjective test problems
 *  and a scalable test problem toolkit", IEEE Transactions on Evolutionary
 *  Computing 10(5): 477-506, 2006.
 *
 *  Example WFG4.
 *
 *  This file is part of a collection of problems developed for
 *  derivative-free multiobjective optimization in
 *  A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
  */
void WFG4(double* f, double* zz){
    int M;
    int k;
    int l;
    int n;
    int i, j;
    double pi;
    double pi2;
    double *S;
    double *A;
    double *zmax;
    double *y;
    double AA;
    double BB;
    double CC;
    double *t1;
    double *w;
    double *t2;
    double sum_wt21, sum_wj1, sum_wt22, sum_wj2;
    double *x;
    double *h;
	double *z;

    M = 3;
    k = 4;
    l = 4;
    n = k+l;
    pi = 4.0*atan(1);
    pi2 = 2.0*atan(1);
    AA = 30.0;
    BB = 10.0;
    CC = 0.35;
    sum_wt21 = 0.0;
    sum_wt22 = 0.0;
    sum_wj1 = 0.0;
    sum_wj2 = 0.0;
	z = (double *)malloc(sizeof(double)*(n+1));
    for(i=1; i<=n; i++)
    {
        z[i] = zz[i-1];
    }

	S = (double*)malloc(sizeof(double)*(M+1));
    for(i = 1; i <= M; i++)
    {
        S[i] = 2.0*(double)i;
    }

    A = (double*)malloc(sizeof(double)*(M));
    for(i = 1; i <= M-1; i++)
    {
        A[i] = 1.0;
    }

    zmax = (double*)malloc(sizeof(double)*(n+1));
    for(i = 1; i <= n; i++)
    {
        zmax[i] = 2.0*(double)i;
    }

    for(i = 1; i <= n; i++)
    {
        if(z[i] >= 0 && z[i] <= zmax[i])
        {
            continue;
        }
        else
        {
            printf("Invalid input range\n");
            return;
        }
    }

    y = (double*)malloc(sizeof(double)*(n+1));
    for(i = 1; i <= n; i++)
    {
        y[i] = z[i]/zmax[i];
    }

    t1 = (double*)malloc(sizeof(double)*(n+1));
    for(i=1; i<=n; i++)
    {
        t1[i] = (1+cos((4*AA+2)*pi*(0.5-fabs(y[i]-CC)/(2*(floor(CC-y[i])+CC))))+4*BB*pow((fabs(y[i]-CC)/(2*(floor(CC-y[i])+CC))), 2))/(BB+2);
    }

    w = (double*)malloc(sizeof(double)*(n+1));
    for(i = 1; i <= n; i++)
    {
        w[i] = 1.0;
    }

    t2 = (double*)malloc(sizeof(double)*(M+1));
    for (i=1; i<=M; i++)
    {
        if (i<=M-1)
        {
			sum_wj1 = 0.0;
			sum_wt21 = 0.0;
            for (j=((i-1)*k/(M-1)+1); j<=(i*k/(M-1)); j++)
            {
                sum_wj1+=w[j];
                sum_wt21+=w[j]*t1[j];
            }
            t2[i] = sum_wt21/sum_wj1;
        }
        else
        {
			sum_wt22 = 0.0;
			sum_wj2 = 0.0;
            for (j=k+1; j<=n; j++)
            {
                sum_wt22+=w[j]*t1[j];
                sum_wj2+=w[j];
            }
            t2[i] = sum_wt22/sum_wj2;
        }
    }

    x = (double*)malloc(sizeof(double)*(M+1));
    for(i = 1; i<=M; i++)
    {
        if(i<=M-1)
        {
            if (t2[M] > A[i])
                x[i] = t2[M]*(t2[i]-0.5)+0.5;
            else
                x[i] = A[i]*(t2[i]-0.5)+0.5;
        }
        else
        {
            x[i] = t2[M];
        }
    }

    h = (double*)malloc(sizeof(double)*(M+1));
    for(i = 1; i <= M; i++)
    {
        if(i == 1)
        {
			h[i] = 1.0;
			for (j=1; j<=M-1; j++)
			{
            	h[i] = h[i]*(sin(x[j]*pi2));
			}
        }
        else if(i<=(M-1))
        {
			h[i] = cos(x[M-i+1]*pi2);
            for(j = 1; j <= M-i; j++)
            {
                h[i] = h[i]*(sin(x[j]*pi2));
            }
        }
        else
        {
            h[i] = cos(x[1]*pi2);
        }
    }

    for(i = 0; i < M; i++)
    {
        f[i] = x[M]+S[i+1]*h[i+1];
    }
    return;
}

/*
 *
 *  As described by Huband et al. in "A review of multiobjective test problems
 *  and a scalable test problem toolkit", IEEE Transactions on Evolutionary
 *  Computing 10(5): 477-506, 2006.
 *
 *  Example WFG5.
 *
 *  This file is part of a collection of problems developed for
 *  derivative-free multiobjective optimization in
 *  A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
  */
void WFG5(double* f, double* zz){
    int M;
    int k;
    int l;
    int n;
    int i, j;
    double pi2;
    double *S;
    double *A;
    double *zmax;
    double *y;
    double AA;
    double BB;
    double CC;
    double *t1;
    double *w;
    double *t2;
    double sum_wt21, sum_wj1, sum_wt22, sum_wj2;
    double *x;
    double *h;
    double *z;

    M = 3;
    k = 4;
    l = 4;
    n = k+l;
    pi2 = 2.0*atan(1);
    AA = 0.35;
    BB = 0.001;
    CC = 0.05;
    sum_wt21 = 0.0;
    sum_wt22 = 0.0;
    sum_wj1 = 0.0;
    sum_wj2 = 0.0;

    z = (double *)malloc(sizeof(double)*(n+1));
    for(i=1; i<=n; i++)
    {
        z[i] = zz[i-1];
    }

    S = (double*)malloc(sizeof(double)*(M+1));
    for(i = 1; i <= M; i++)
    {
        S[i] = 2.0*(double)i;
    }

    A = (double*)malloc(sizeof(double)*(M));
    for(i = 1; i <= M-1; i++)
    {
        A[i] = 1.0;
    }

    zmax = (double*)malloc(sizeof(double)*(n+1));
    for(i = 1; i <= n; i++)
    {
        zmax[i] = 2.0*(double)i;
    }

    for(i = 1; i <= n; i++)
    {
        if(z[i] >= 0 && z[i] <= zmax[i])
        {
            continue;
        }
        else
        {
            printf("Invalid input range\n");
            return;
        }
    }

    y = (double*)malloc(sizeof(double)*(n+1));
    for(i = 1; i <= n; i++)
    {
        y[i] = z[i]/zmax[i];
    }

    t1 = (double*)malloc(sizeof(double)*(n+1));
    for(i=1; i<=n; i++)
    {
        t1[i] = 1.0+(fabs(y[i]-AA)-BB)*((double)(floor(y[i]-AA+BB)*(1.0-CC+(double)(AA-BB)/(double)BB))/(double)(AA-BB)+(floor(AA+BB-y[i])*(1.0-CC+(double)(1.0-AA-BB)/(double)BB))/(double)(1.0-AA-BB)+1.0/(double)BB);
    }

    w = (double*)malloc(sizeof(double)*(n+1));
    for(i = 1; i <= n; i++)
    {
        w[i] = 1.0;
    }

    t2 = (double*)malloc(sizeof(double)*(M+1));
    for (i=1; i<=M; i++)
    {
        if (i<=M-1)
        {
			sum_wj1 = 0.0;
			sum_wt21 = 0.0;
            for (j=((i-1)*k/(M-1)+1); j<=(i*k/(M-1)); j++)
            {
                sum_wj1+=w[j];
                sum_wt21+=w[j]*t1[j];
            }
            t2[i] = sum_wt21/sum_wj1;
        }
        else
        {
			sum_wt22 = 0.0;
			sum_wj2 = 0.0;
            for (j=k+1; j<=n; j++)
            {
                sum_wt22+=w[j]*t1[j];
                sum_wj2+=w[j];
            }
            t2[i] = sum_wt22/sum_wj2;
        }
    }

    x = (double*)malloc(sizeof(double)*(M+1));
    for(i = 1; i<=M; i++)
    {
        if(i<=M-1)
        {
            if (t2[M] > A[i])
                x[i] = t2[M]*(t2[i]-0.5)+0.5;
            else
                x[i] = A[i]*(t2[i]-0.5)+0.5;
        }
        else
        {
            x[i] = t2[M];
        }
    }

    h = (double*)malloc(sizeof(double)*(M+1));
    for(i = 1; i <= M; i++)
    {
        if(i == 1)
        {
			h[i] = 1.0;
			for (j=1; j<=M-1; j++)
			{
            	h[i] = h[i]*(sin(x[j]*pi2));
			}
        }
        else if(i<=(M-1))
        {
			h[i] = cos(x[M-i+1]*pi2);
            for(j = 1; j <= M-i; j++)
            {
                h[i] = h[i]*(sin(x[j]*pi2));
            }
        }
        else
        {
            h[i] = cos(x[1]*pi2);
        }
    }

    for(i = 0; i < M; i++)
    {
        f[i] = x[M]+S[i+1]*h[i+1];
    }
    return;
}

/*
 *
 *  As described by Huband et al. in "A review of multiobjective test problems
 *  and a scalable test problem toolkit", IEEE Transactions on Evolutionary
 *  Computing 10(5): 477-506, 2006.
 *
 *  Example IKK1, see the previous cited paper for the original reference.
 *
 *  This file is part of a collection of problems developed for
 *  derivative-free multiobjective optimization in
 *  A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
 *  Direct Multisearch for Multiobjective Optimization, 2010.
 *
 *  Written by the authors in June 1, 2010.
  */
void IKK1(double *f, double *x){
    int n, i;
    n = 2;

    for(i=0; i<n; ++i){
      if(!(x[i]>=-50 && x[i]<=50)){
            printf("Invalid input range\n");
            return;
      }
    }

    f[0] = pow(x[0], 2);
    f[1] = pow((x[0]-20), 2);
    f[2] = pow(x[1], 2);

    return;
}

/*
*
*   As described by Y. Jin, M. Olhofer and B. Sendhoff. "Dynamic weighted
*   aggregation for evolutionary multi-objective optimization: Why does it
*   work and how?", in Proceedings of Genetic and Evolutionary Computation
*   Conference, pp.1042-1049, San Francisco, USA, 2001.
*
*   Test function 3, F3.
*
*   This file is part of a collection of problems developed for
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
 */

void Jin3(double* f, double* x){
	int n = 2;
	int i;
	double gx;

	for(i = 0; i < n; i++){
		if(x[i] >= 0 && x[i] <= 1);
		else{
			printf("Invalid input range\n");
			return;
		}
	}

	gx = 0;
	for(i = 1; i < n; i++){
		gx = gx + x[i];
	}
	gx = 9*gx;
	gx = gx / (n-1);
	gx = gx + 1;

	f[1] = x[0];
	f[0] = gx*(1 - pow((x[0]/gx),2));

	return;
}


/*
*
*   As described by Huband et al. in "A review of multiobjective test problems
*   and a scalable test problem toolkit", IEEE Transactions on Evolutionary
*   Computing 10(5): 477-506, 2006.
*
*   Example WFG9.
*
*   This file is part of a collection of problems developed for
*   derivative-free multiobjective optimization in
*   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
*   Direct Multisearch for Multiobjective Optimization, 2010.
*
*   Written by the authors in June 1, 2010.
 */

void WFG9(double* f, double* zz){
    int M, l, n, i, j, ii, jj;
    double pi, pi2;
    double *S;
    double *A;
    double *zmax, *z;
    double *y;
    double k, AA, BB, CC, AAA, BBB, CCC, AAAA, BBBB, CCCC, sum_w, sum_wy;
    double  sum1, sum2;
    double *w;
    double *r_sum;
    double *t1;
    double *t2;
    double *t3;
    double *x;
    double *h;

    M = 3;
    k = 4;
    l = 4;
    n = k+l;
    pi = 4.0*atan(1);
    pi2 = 2.0*atan(1);
    AA = 0.98/49.98;
    BB = 0.02;
    CC = 50.0;
    sum_w = 0.0;
    sum_wy = 0.0;
    AAA = 0.35;
    BBB = 0.001;
    CCC = 0.05;
    AAAA = 30.0;
    BBBB = 95.0;
    CCCC = 0.35;

    S = (double*)malloc(sizeof(double)*(M+1));
    for(i = 1; i <= M; i++)
        S[i] = 2.0*i;

    A = (double*)malloc(sizeof(double)*(M));
    for(i = 1; i <= M-1; i++)
        A[i] = 1.0;

    zmax = (double*)malloc(sizeof(double)*(n+1));
    for(i = 1; i <= n; i++)
        zmax[i] = 2.0*i;

    z = (double*)malloc(sizeof(double)*(n+1));
    for(i=n; i>0; --i)
        z[i] = zz[i-1];

    for(i = 1; i <= n; i++){
        if(!(z[i] >= 0 && z[i] <= zmax[i])){
            printf("Invalid input range\n");
            return;
        }
    }
    y = (double*)malloc(sizeof(double)*(n+1));
    for(i = 1; i <= n; i++)
        y[i] = z[i]/zmax[i];

    w = (double*)malloc(sizeof(double)*(n+1));
    for(i = 1; i <= n; i++)
        w[i] = 1.0;

    r_sum = (double*)malloc(sizeof(double)*n);
    for (i=1; i<n; i++){
		    sum_w = 0.0;
		    sum_wy = 0.0;
        for (j=i+1; j<=n; j++){
            sum_w+=w[j];
            sum_wy+=w[j]*y[j];
        }
        r_sum[i] = sum_wy/sum_w;
    }

    t1 = (double*)malloc(sizeof(double)*(n+1));
    for(i = 1; i <= n; i++){
        if (i<=n-1)
            t1[i] = pow(y[i], (BB+(CC-BB)*(AA-(1-2*r_sum[i])*fabs(floor(0.5-r_sum[i])+AA))));
        else
            t1[i] = y[i];
    }

    t2 = (double*)malloc(sizeof(double)*(n+1));
    for(i=1; i<=n; i++){
        if(i <= k)
          t2[i] = 1+(fabs(t1[i]-AAA)-BBB)*((floor(t1[i]-AAA+BBB)*(1-CCC+(AAA-BBB)/BBB))/(AAA-BBB)+(floor(AAA+BBB-t1[i])*(1-CCC+(1-AAA-BBB)/BBB))/(1-AAA-BBB)+1/BBB);
        else
            t2[i] = (1+cos((4*AAAA+2)*pi*(0.5-fabs(t1[i]-CCCC)/(2*(floor(CCCC-t1[i])+CCCC))))+4*BBBB*pow((fabs(t1[i]-CCCC)/(2*(floor(CCCC-t1[i])+CCCC))), 2))/(BBBB+2);
    }

    t3 = (double*)malloc(sizeof(double)*(M+1));
    for(i = 1; i <= M; i++)    {
        if(i <= M-1){
            sum1 = 0;
            for(ii = ((i-1)*k/(M-1)+1); ii <=(i*k/(M-1)); ii++){
                sum2 = 0;
                for(jj = 0; jj <= (k/(M-1)-2); jj++)
                    sum2+=fabs(t2[ii]-t2[(int)(((i-1)*k/(M-1)+1)+(fmod((ii+jj-((i-1)*k/(M-1)+1)+1), (((i+1)*k/(M-1))-((i)*k/(M-1)+1)+1))))]);
                sum1+=(t2[ii]+sum2);
            }
            t3[i] = sum1/(((i*k/(M-1))-((i-1)*k/(M-1)+1)+1)/(k/(M-1))*ceil(k/(M-1)/2)*(1+2*k/(M-1)-2*ceil(k/(M-1)/2)));
        }
        else{
            sum1 = 0.0;
            for(ii = k+1; ii <= n; ii++){
                sum2 = 0.0;
                for(jj = 0; jj <= l-2; jj++)
                    sum2+=fabs(t2[ii]-t2[(int)(k+1+fmod((ii+jj-(k+1)+1), (n-k)))]);
                sum1+=(t2[ii]+sum2);
            }
            t3[i] = sum1/(((n-k)/l)*ceil(l/2)*(1+2*l-2*ceil(l/2)));
        }
    }

    x = (double*)malloc(sizeof(double)*(M+1));
    for(i = 1; i <= M; i++){
        if(i <= M-1)
            x[i] = ((t3[M]>A[i]) ? t3[M] : A[i])*(t3[i]-0.5)+0.5;
        else
            x[i] = t3[M];
    }

    h = (double*)malloc(sizeof(double)*(M+1));
    for(i = 1; i <= M; i++){
        if(i == 1){
            h[i] = 1.0;
            for(j=1; j<=M-1; j++)
                h[i] = h[i]*(sin(x[j]*pi2));
        }
        else if(i<=(M-1)){
            h[i] = (cos(x[M-i+1]*pi2));
            for(j = 1; j <= M-i; j++)
                h[i] = h[i]*(sin(x[j]*pi2));
        }
        else
            h[i] = cos(x[1]*pi2);
    }

    for(i = 0; i < M; i++){
        f[i] = x[M]+S[i+1]*h[i+1];
    }

    return;
}



int main(int argc, char const *argv[]) {
  /* code */
  return 0;
}
