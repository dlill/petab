#include <R.h>
 #include <math.h>
 void PInjected_petab_deriv_k0gcp17o ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = p[1] ;
y[1+i**l] = p[0] ; 
}
}