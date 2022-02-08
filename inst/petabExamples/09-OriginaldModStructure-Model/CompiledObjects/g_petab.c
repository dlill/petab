#include <R.h>
 #include <math.h>
 void g_petab_wipy0sdl ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = log10((x[3+i**k]+x[4+i**k])+p[0]) ;
y[1+i**l] = log10(p[1]*x[4+i**k]+p[2]) ; 
}
}