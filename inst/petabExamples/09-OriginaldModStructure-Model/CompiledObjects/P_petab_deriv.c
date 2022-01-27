#include <R.h>
 #include <math.h>
 void P_petab_deriv_n9evusbd ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = exp(p[0]) ;
y[16+i**l] = exp(p[1]) ;
y[32+i**l] = 1.0 ;
y[48+i**l] = exp(p[3]) ;
y[49+i**l] = exp(p[3])*exp(p[4])*exp(p[5])/exp(p[4]) ;
y[64+i**l] = exp(p[3])*exp(p[4])*exp(p[5])/exp(p[4])-exp(p[3])*exp(p[4])*exp(p[5])*exp(p[4])/pow(exp(p[4]),2.0) ;
y[70+i**l] = exp(p[4]) ;
y[79+i**l] = exp(p[3])*exp(p[4])*exp(p[5])/exp(p[4]) ;
y[95+i**l] = 1.0 ;
y[111+i**l] = exp(p[7]) ;
y[127+i**l] = exp(p[8]) ;
y[143+i**l] = exp(p[9]) ;
y[159+i**l] = 1.0 ;
y[176+i**l] = 1.0 ;
y[192+i**l] = exp(p[12]) ;
y[208+i**l] = exp(p[13]) ;
y[224+i**l] = 1.0 ; 
}
}