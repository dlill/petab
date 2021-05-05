#include <R.h>
 #include <math.h>
 void P_petab_deriv_104ovco3 ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = pow(10.0,(p[0]))*log(10.0) ;
y[14+i**l] = pow(10.0,(p[1]))*log(10.0) ;
y[28+i**l] = 1.0 ;
y[42+i**l] = 1.0 ;
y[56+i**l] = pow(10.0,(p[4]))*log(10.0) ;
y[70+i**l] = pow(10.0,(p[5]))*log(10.0) ;
y[84+i**l] = pow(10.0,(p[6]))*log(10.0) ;
y[98+i**l] = pow(10.0,(p[7]))*log(10.0) ;
y[112+i**l] = 1.0 ;
y[126+i**l] = pow(10.0,(p[9]))*log(10.0) ;
y[140+i**l] = pow(10.0,(p[10]))*log(10.0) ;
y[154+i**l] = pow(10.0,(p[11]))*log(10.0) ;
y[168+i**l] = pow(10.0,(p[12]))*log(10.0) ; 
}
}