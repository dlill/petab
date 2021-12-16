#include <R.h>
 #include <math.h>
 void errfn_petab_h8x7ceux ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = p[0] ;
y[1+i**l] = p[1] ;
y[2+i**l] = p[2] ;
y[3+i**l] = p[3] ; 
}
}