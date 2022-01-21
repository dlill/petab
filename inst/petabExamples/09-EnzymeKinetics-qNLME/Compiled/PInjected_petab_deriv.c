#include <R.h>
 #include <math.h>
 void PInjected_petab_deriv_aejal9wc ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = p[1] ;
y[2+i**l] = p[0] ;
y[5+i**l] = p[3] ;
y[7+i**l] = p[2] ; 
}
}