#include <R.h>
 #include <math.h>
 void g_petab_s7fk3ang ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = log10(p[0]*x[4+i**k]+p[1]) ; 
}
}