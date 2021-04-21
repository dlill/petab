#include <R.h>
 #include <math.h>
 void errfn_petab_deriv_niwmq22q ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[16+i**l] = 1.0 ;
y[21+i**l] = 1.0 ;
y[26+i**l] = 1.0 ;
y[31+i**l] = 1.0 ; 
}
}