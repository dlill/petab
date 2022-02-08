#include <R.h>
 #include <math.h>
 void P_petab_deriv_n9evusbd ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = exp(p[0]) ;
y[18+i**l] = exp(p[1]) ;
y[37+i**l] = exp(p[2]) ;
y[38+i**l] = exp(p[2])*exp(p[3])*exp(p[4])/exp(p[3]) ;
y[55+i**l] = exp(p[2])*exp(p[3])*exp(p[4])/exp(p[3])-exp(p[2])*exp(p[3])*exp(p[4])*exp(p[3])/pow(exp(p[3]),2.0) ;
y[60+i**l] = exp(p[3])*exp(p[4]) ;
y[61+i**l] = exp(p[3]) ;
y[72+i**l] = exp(p[2])*exp(p[3])*exp(p[4])/exp(p[3]) ;
y[77+i**l] = exp(p[3])*exp(p[4]) ;
y[91+i**l] = exp(p[5]) ;
y[109+i**l] = exp(p[6]) ;
y[127+i**l] = exp(p[7]) ;
y[149+i**l] = exp(p[8]) ;
y[167+i**l] = exp(p[9]) ;
y[186+i**l] = exp(p[10]) ; 
}
}