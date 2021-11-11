#include <R.h>
 #include <math.h>
 void g_petab_jy59jrjm ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = log(x[2+i**k]) ; 
}
}