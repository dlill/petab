#include <R.h>
 #include <math.h>
 void e_coomjscl ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = p[0] ;
y[1+i**l] = p[1] ; 
}
}