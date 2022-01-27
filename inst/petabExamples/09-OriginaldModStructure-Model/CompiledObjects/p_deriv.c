#include <R.h>
 #include <math.h>
 void p_deriv_pesnbnrl ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = exp(p[0]) ;
y[15+i**l] = exp(p[1]) ;
y[31+i**l] = exp(p[2]) ;
y[32+i**l] = exp(p[2])*exp(p[3])*exp(p[4])/exp(p[3]) ;
y[46+i**l] = exp(p[2])*exp(p[3])*exp(p[4])/exp(p[3])-exp(p[2])*exp(p[3])*exp(p[4])*exp(p[3])/pow(exp(p[3]),2.0) ;
y[51+i**l] = exp(p[3])*exp(p[4]) ;
y[52+i**l] = exp(p[3]) ;
y[60+i**l] = exp(p[2])*exp(p[3])*exp(p[4])/exp(p[3]) ;
y[65+i**l] = exp(p[3])*exp(p[4]) ;
y[76+i**l] = exp(p[5]) ;
y[91+i**l] = exp(p[6]) ;
y[106+i**l] = exp(p[7]) ;
y[123+i**l] = exp(p[8]) ;
y[138+i**l] = exp(p[9]) ;
y[153+i**l] = exp(p[10]) ; 
}
}