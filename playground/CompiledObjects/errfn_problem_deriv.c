#include <R.h>
 #include <math.h>
 void errfn_problem_deriv_t8u4byhb ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[1+i**l] = 1.0 ; 
}
}