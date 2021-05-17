#include <R.h>
 #include <math.h>
 void errfn_Bachmann_MSB2011_32n0m4id ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = p[0] ;
y[1+i**l] = p[1] ;
y[2+i**l] = p[2] ;
y[3+i**l] = p[3] ;
y[4+i**l] = p[4] ;
y[5+i**l] = p[5] ;
y[6+i**l] = p[6] ;
y[7+i**l] = p[7] ;
y[8+i**l] = p[8] ;
y[9+i**l] = p[9] ;
y[10+i**l] = p[10] ;
y[11+i**l] = p[11] ;
y[12+i**l] = p[12] ;
y[13+i**l] = p[13] ;
y[14+i**l] = p[14] ;
y[15+i**l] = p[15] ;
y[16+i**l] = p[16] ;
y[17+i**l] = p[17]+p[18] ;
y[18+i**l] = p[19] ;
y[19+i**l] = p[20] ; 
}
}