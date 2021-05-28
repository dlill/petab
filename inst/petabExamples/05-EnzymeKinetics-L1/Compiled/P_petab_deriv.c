#include <R.h>
 #include <math.h>
 void P_petab_deriv_uqmzfbim ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = pow(10.0,(p[0]))*log(10.0) ;
y[12+i**l] = pow(10.0,(p[1]))*log(10.0) ;
y[24+i**l] = pow(10.0,(p[2]))*log(10.0) ;
y[36+i**l] = pow(10.0,(p[3]))*log(10.0) ;
y[48+i**l] = 1.0 ;
y[60+i**l] = 1.0 ;
y[72+i**l] = pow(10.0,(p[6]))*log(10.0) ;
y[84+i**l] = pow(10.0,(p[7]))*log(10.0) ;
y[96+i**l] = pow(10.0,(p[8]))*log(10.0) ;
y[108+i**l] = pow(10.0,(p[9]))*log(10.0) ;
y[120+i**l] = pow(10.0,(p[10]))*log(10.0) ; 
}
}