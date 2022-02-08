#include <R.h>
 #include <math.h>
 void errfn_petab_deriv_cpx1igxa ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[4+i**l] = 1.0 ;
y[7+i**l] = 1.0 ; 
}
}