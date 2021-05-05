#include <R.h>
 #include <math.h>
 void preeq_petab_c0_dfdx_ijgvbzbm ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = -(p[5]+p[6]*p[1]) ;
y[2+i**l] = p[6]*p[1] ;
y[4+i**l] = -(p[6]*p[0]) ;
y[5+i**l] = 1.0 ;
y[6+i**l] = p[6]*p[0] ;
y[8+i**l] = p[7] ;
y[9+i**l] = 1.0 ;
y[10+i**l] = -(p[7]+p[8]) ;
y[11+i**l] = p[8] ; 
}
}