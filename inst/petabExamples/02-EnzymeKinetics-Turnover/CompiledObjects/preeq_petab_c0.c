#include <R.h>
 #include <math.h>
 void preeq_petab_c0_1iwv8pmj ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = 1.0*(p[4]*1.0)-1.0*(p[5]*p[0]*1.0)-1.0*(p[6]*p[1]*p[0]*1.0)+1.0*(p[7]*p[2]*1.0) ;
y[1+i**l] = 1.0*p[1]+1.0*p[2]-1.0 ;
y[2+i**l] = 1.0*(p[6]*p[1]*p[0]*1.0)-1.0*(p[7]*p[2]*1.0)-1.0*(p[8]*p[2]*1.0) ;
y[3+i**l] = 1.0*(p[8]*p[2]*1.0) ; 
}
}