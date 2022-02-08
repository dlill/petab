#include <R.h>
 #include <math.h>
 void p_deriv_pesnbnrl ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = exp(p[0]) ;
y[17+i**l] = exp(p[1]) ;
y[35+i**l] = exp(p[2]) ;
y[36+i**l] = exp(p[2])*exp(p[3])*exp(p[4])/exp(p[3]) ;
y[52+i**l] = exp(p[2])*exp(p[3])*exp(p[4])/exp(p[3])-exp(p[2])*exp(p[3])*exp(p[4])*exp(p[3])/pow(exp(p[3]),2.0) ;
y[57+i**l] = exp(p[3])*exp(p[4]) ;
y[58+i**l] = exp(p[3]) ;
y[68+i**l] = exp(p[2])*exp(p[3])*exp(p[4])/exp(p[3]) ;
y[73+i**l] = exp(p[3])*exp(p[4]) ;
y[86+i**l] = exp(p[5]) ;
y[103+i**l] = exp(p[6]) ;
y[120+i**l] = exp(p[7]) ;
y[139+i**l] = exp(p[8]) ;
y[156+i**l] = exp(p[9]) ;
y[174+i**l] = exp(p[10]) ; 
}
}