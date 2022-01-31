#include <R.h>
 #include <math.h>
 void P_petab_yppihbdm ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = exp(p[0]) ;
y[1+i**l] = exp(p[1]) ;
y[2+i**l] = p[2] ;
y[3+i**l] = exp(p[3]) ;
y[4+i**l] = exp(p[3])*exp(p[4])*exp(p[5])/exp(p[4]) ;
y[5+i**l] = p[6] ;
y[6+i**l] = exp(p[7]) ;
y[7+i**l] = exp(p[8]) ;
y[8+i**l] = exp(p[9]) ;
y[9+i**l] = exp(p[4])*exp(p[5]) ;
y[10+i**l] = exp(p[4]) ;
y[11+i**l] = p[10] ;
y[12+i**l] = exp(p[11]) ;
y[13+i**l] = exp(p[12]) ;
y[14+i**l] = exp(p[13]) ; 
}
}