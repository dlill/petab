#include <R.h>
 #include <math.h>
 void P_petab_deriv_02c5xuqq ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = pow(10.0,(p[0]))*log(10.0) ;
y[13+i**l] = pow(10.0,(p[1]))*log(10.0) ;
y[26+i**l] = pow(10.0,(p[2]))*log(10.0) ;
y[39+i**l] = pow(10.0,(p[3]))*log(10.0) ;
y[52+i**l] = pow(10.0,(p[4]))*log(10.0) ;
y[65+i**l] = pow(10.0,(p[5]))*log(10.0) ;
y[78+i**l] = pow(10.0,(p[6]))*log(10.0) ;
y[91+i**l] = 1.0 ;
y[104+i**l] = pow(10.0,(p[8]))*log(10.0) ;
y[117+i**l] = pow(10.0,(p[9]))*log(10.0) ;
y[130+i**l] = pow(10.0,(p[10]))*log(10.0) ;
y[143+i**l] = pow(10.0,(p[11]))*log(10.0) ; 
}
}