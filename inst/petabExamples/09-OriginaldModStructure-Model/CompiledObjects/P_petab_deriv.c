#include <R.h>
 #include <math.h>
 void P_petab_deriv_84urzvmn ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = exp(p[0]) ;
y[16+i**l] = exp(p[1]) ;
y[32+i**l] = 1.0 ;
y[48+i**l] = exp(p[3]) ;
y[49+i**l] = exp(p[3])*exp(p[4])*exp(p[5])/exp(p[4]) ;
y[64+i**l] = exp(p[3])*exp(p[4])*exp(p[5])/exp(p[4])-exp(p[3])*exp(p[4])*exp(p[5])*exp(p[4])/pow(exp(p[4]),2.0) ;
y[69+i**l] = exp(p[4])*exp(p[5]) ;
y[70+i**l] = exp(p[4]) ;
y[79+i**l] = exp(p[3])*exp(p[4])*exp(p[5])/exp(p[4]) ;
y[84+i**l] = exp(p[4])*exp(p[5]) ;
y[95+i**l] = 1.0 ;
y[111+i**l] = exp(p[7]) ;
y[127+i**l] = exp(p[8]) ;
y[143+i**l] = exp(p[9]) ;
y[161+i**l] = 1.0 ;
y[177+i**l] = exp(p[11]) ;
y[193+i**l] = exp(p[12]) ;
y[209+i**l] = exp(p[13]) ; 
}
}