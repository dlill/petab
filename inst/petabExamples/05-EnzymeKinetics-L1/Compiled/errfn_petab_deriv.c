#include <R.h>
 #include <math.h>
 void errfn_petab_deriv_71xrhf6f ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = (p[0])*(x[9+i**k]) ;
y[4+i**l] = (p[0])*(x[13+i**k]) ;
y[8+i**l] = (p[0])*(x[17+i**k]) ;
y[12+i**l] = (p[0])*(x[21+i**k]) ;
y[16+i**l] = (p[0])*(x[25+i**k])+x[4+i**k] ;
y[20+i**l] = (p[0])*(x[29+i**k])+1.0 ;
y[24+i**l] = (p[0])*(x[33+i**k]) ;
y[25+i**l] = 1.0 ;
y[28+i**l] = (p[0])*(x[37+i**k]) ;
y[30+i**l] = 1.0 ;
y[32+i**l] = (p[0])*(x[41+i**k]) ;
y[35+i**l] = 1.0 ;
y[36+i**l] = (p[0])*(x[45+i**k]) ;
y[40+i**l] = (p[0])*(x[49+i**k]) ;
y[44+i**l] = (p[0])*(x[53+i**k]) ; 
}
}