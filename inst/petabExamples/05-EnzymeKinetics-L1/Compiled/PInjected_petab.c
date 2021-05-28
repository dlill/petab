#include <R.h>
 #include <math.h>
 void PInjected_petab_2d982n99 ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = p[0]*p[1] ;
y[1+i**l] = p[2]*p[3] ; 
}
}