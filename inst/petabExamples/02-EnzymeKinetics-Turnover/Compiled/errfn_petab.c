#include <R.h>
 #include <math.h>
 void errfn_petab_llihsqzv ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = p[0]*x[4+i**k] ;
y[1+i**l] = p[1]*x[5+i**k] ;
y[2+i**l] = p[2]*x[6+i**k] ;
y[3+i**l] = p[3]*x[7+i**k] ; 
}
}