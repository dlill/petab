#include <R.h>
 #include <math.h>
 void g_petab_deriv_5o55nvd1 ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = (1.0/x[2+i**k])*(x[6+i**k]) ;
y[1+i**l] = (1.0/x[2+i**k])*(x[10+i**k]) ;
y[2+i**l] = (1.0/x[2+i**k])*(x[14+i**k]) ;
y[3+i**l] = (1.0/x[2+i**k])*(x[18+i**k]) ;
y[4+i**l] = (1.0/x[2+i**k])*(x[22+i**k]) ;
y[5+i**l] = (1.0/x[2+i**k])*(x[26+i**k]) ; 
}
}