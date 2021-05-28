#include <R.h>
 #include <math.h>
 void P_petab_deriv_dlfw84zv ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = pow(10.0,(p[0]))*log(10.0) ;
y[16+i**l] = pow(10.0,(p[1]))*log(10.0) ;
y[32+i**l] = pow(10.0,(p[2]))*log(10.0) ;
y[48+i**l] = pow(10.0,(p[3]))*log(10.0) ;
y[64+i**l] = 1.0 ;
y[80+i**l] = 1.0 ;
y[96+i**l] = pow(10.0,(p[6]))*log(10.0) ;
y[112+i**l] = pow(10.0,(p[7]))*log(10.0) ;
y[128+i**l] = pow(10.0,(p[8]))*log(10.0) ;
y[144+i**l] = pow(10.0,(p[9]))*log(10.0) ;
y[160+i**l] = pow(10.0,(p[10]))*log(10.0) ;
y[176+i**l] = pow(10.0,(p[11]))*log(10.0) ;
y[192+i**l] = pow(10.0,(p[12]))*log(10.0) ;
y[208+i**l] = pow(10.0,(p[13]))*log(10.0) ;
y[224+i**l] = pow(10.0,(p[14]))*log(10.0) ; 
}
}