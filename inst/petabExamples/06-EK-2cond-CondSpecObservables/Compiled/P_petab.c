#include <R.h>
 #include <math.h>
 void P_petab_ajhwpjtn ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = pow(10.0,(p[0])) ;
y[1+i**l] = pow(10.0,(p[1])) ;
y[2+i**l] = pow(10.0,(p[2])) ;
y[3+i**l] = pow(10.0,(p[3])) ;
y[4+i**l] = pow(10.0,(p[4])) ;
y[5+i**l] = p[5] ;
y[6+i**l] = pow(10.0,(p[6])) ;
y[7+i**l] = p[7] ;
y[8+i**l] = pow(10.0,(p[8])) ;
y[9+i**l] = pow(10.0,(p[9])) ;
y[10+i**l] = pow(10.0,(p[10])) ;
y[11+i**l] = pow(10.0,(p[11])) ;
y[12+i**l] = pow(10.0,(p[12])) ;
y[13+i**l] = pow(10.0,(p[13])) ;
y[14+i**l] = pow(10.0,(p[14])) ; 
}
}