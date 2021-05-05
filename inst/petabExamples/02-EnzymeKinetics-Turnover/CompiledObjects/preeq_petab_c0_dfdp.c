#include <R.h>
 #include <math.h>
 void preeq_petab_c0_dfdp_qqznty1j ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = 1.0 ;
y[4+i**l] = -p[0] ;
y[8+i**l] = -(p[1]*p[0]) ;
y[10+i**l] = p[1]*p[0] ;
y[12+i**l] = p[2] ;
y[14+i**l] = -p[2] ;
y[18+i**l] = -p[2] ;
y[19+i**l] = p[2] ; 
}
}