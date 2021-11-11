#include <R.h>
 #include <math.h>
 void P_petab_deriv_8miptvsx ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = pow(10.0,(p[0]))*log(10.0) ;
y[8+i**l] = pow(10.0,(p[1]))*log(10.0) ;
y[16+i**l] = pow(10.0,(p[2]))*log(10.0) ;
y[24+i**l] = pow(10.0,(p[3]))*log(10.0) ;
y[32+i**l] = pow(10.0,(p[4]))*log(10.0) ;
y[40+i**l] = 1.0 ;
y[48+i**l] = pow(10.0,(p[6]))*log(10.0) ; 
}
}