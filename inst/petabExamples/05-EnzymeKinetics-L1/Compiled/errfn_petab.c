#include <R.h>
 #include <math.h>
 void errfn_petab_uia5ygrn ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = p[0]*x[4+i**k]+p[1] ;
y[1+i**l] = p[2] ;
y[2+i**l] = p[3] ;
y[3+i**l] = p[4] ; 
}
}