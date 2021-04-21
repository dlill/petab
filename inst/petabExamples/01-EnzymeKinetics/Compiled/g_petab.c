#include <R.h>
 #include <math.h>
 void g_petab_8gtfb5mh ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = log(x[0+i**k]) ;
y[1+i**l] = log(x[1+i**k]) ;
y[2+i**l] = log(x[2+i**k]) ;
y[3+i**l] = log(x[3+i**k]) ; 
}
}