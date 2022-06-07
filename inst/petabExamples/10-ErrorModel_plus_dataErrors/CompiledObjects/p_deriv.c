#include <R.h>
 #include <math.h>
 void p_deriv_pesnbnrl ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = exp(p[0]) ;
y[16+i**l] = exp(p[1]) ;
y[33+i**l] = exp(p[2]) ;
y[34+i**l] = exp(p[2])*exp(p[3])*exp(p[4])/exp(p[3]) ;
y[49+i**l] = exp(p[2])*exp(p[3])*exp(p[4])/exp(p[3])-exp(p[2])*exp(p[3])*exp(p[4])*exp(p[3])/pow(exp(p[3]),2.0) ;
y[54+i**l] = exp(p[3])*exp(p[4]) ;
y[55+i**l] = exp(p[3]) ;
y[64+i**l] = exp(p[2])*exp(p[3])*exp(p[4])/exp(p[3]) ;
y[69+i**l] = exp(p[3])*exp(p[4]) ;
y[81+i**l] = exp(p[5]) ;
y[97+i**l] = exp(p[6]) ;
y[113+i**l] = exp(p[7]) ;
y[131+i**l] = exp(p[8]) ;
y[147+i**l] = exp(p[9]) ;
y[164+i**l] = exp(p[10]) ; 
}
}