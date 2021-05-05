#include <R.h>
 #include <math.h>
 void P_petab_deriv_2yoabz1q ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = pow(10.0,(p[0]))*log(10.0) ;
y[19+i**l] = 1.0 ;
y[38+i**l] = 1.0 ;
y[57+i**l] = 1.0 ;
y[76+i**l] = pow(10.0,(p[4]))*log(10.0) ;
y[95+i**l] = pow(10.0,(p[5]))*log(10.0) ;
y[114+i**l] = pow(10.0,(p[6]))*log(10.0) ;
y[133+i**l] = pow(10.0,(p[7]))*log(10.0) ;
y[152+i**l] = pow(10.0,(p[8]))*log(10.0) ;
y[171+i**l] = 1.0 ;
y[190+i**l] = pow(10.0,(p[10]))*log(10.0) ;
y[209+i**l] = 1.0 ;
y[228+i**l] = pow(10.0,(p[12]))*log(10.0) ;
y[247+i**l] = 1.0 ;
y[266+i**l] = pow(10.0,(p[14]))*log(10.0) ;
y[285+i**l] = 1.0 ;
y[304+i**l] = pow(10.0,(p[16]))*log(10.0) ;
y[323+i**l] = 1.0 ; 
}
}