#include <R.h>
 #include <math.h>
 void errfn_petab_deriv_ugo8kken ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = (p[0])*(x[9+i**k]) ;
y[1+i**l] = (p[1])*(x[10+i**k]) ;
y[2+i**l] = (p[2])*(x[11+i**k]) ;
y[3+i**l] = (p[3])*(x[12+i**k]) ;
y[4+i**l] = (p[0])*(x[13+i**k]) ;
y[5+i**l] = (p[1])*(x[14+i**k]) ;
y[6+i**l] = (p[2])*(x[15+i**k]) ;
y[7+i**l] = (p[3])*(x[16+i**k]) ;
y[8+i**l] = (p[0])*(x[17+i**k]) ;
y[9+i**l] = (p[1])*(x[18+i**k]) ;
y[10+i**l] = (p[2])*(x[19+i**k]) ;
y[11+i**l] = (p[3])*(x[20+i**k]) ;
y[12+i**l] = (p[0])*(x[21+i**k]) ;
y[13+i**l] = (p[1])*(x[22+i**k]) ;
y[14+i**l] = (p[2])*(x[23+i**k]) ;
y[15+i**l] = (p[3])*(x[24+i**k]) ;
y[16+i**l] = (p[0])*(x[25+i**k])+x[4+i**k] ;
y[17+i**l] = (p[1])*(x[26+i**k]) ;
y[18+i**l] = (p[2])*(x[27+i**k]) ;
y[19+i**l] = (p[3])*(x[28+i**k]) ;
y[20+i**l] = (p[0])*(x[29+i**k]) ;
y[21+i**l] = (p[1])*(x[30+i**k])+x[5+i**k] ;
y[22+i**l] = (p[2])*(x[31+i**k]) ;
y[23+i**l] = (p[3])*(x[32+i**k]) ;
y[24+i**l] = (p[0])*(x[33+i**k]) ;
y[25+i**l] = (p[1])*(x[34+i**k]) ;
y[26+i**l] = (p[2])*(x[35+i**k])+x[6+i**k] ;
y[27+i**l] = (p[3])*(x[36+i**k]) ;
y[28+i**l] = (p[0])*(x[37+i**k]) ;
y[29+i**l] = (p[1])*(x[38+i**k]) ;
y[30+i**l] = (p[2])*(x[39+i**k]) ;
y[31+i**l] = (p[3])*(x[40+i**k])+x[7+i**k] ;
y[32+i**l] = (p[0])*(x[41+i**k]) ;
y[33+i**l] = (p[1])*(x[42+i**k]) ;
y[34+i**l] = (p[2])*(x[43+i**k]) ;
y[35+i**l] = (p[3])*(x[44+i**k]) ;
y[36+i**l] = (p[0])*(x[45+i**k]) ;
y[37+i**l] = (p[1])*(x[46+i**k]) ;
y[38+i**l] = (p[2])*(x[47+i**k]) ;
y[39+i**l] = (p[3])*(x[48+i**k]) ;
y[40+i**l] = (p[0])*(x[49+i**k]) ;
y[41+i**l] = (p[1])*(x[50+i**k]) ;
y[42+i**l] = (p[2])*(x[51+i**k]) ;
y[43+i**l] = (p[3])*(x[52+i**k]) ;
y[44+i**l] = (p[0])*(x[53+i**k]) ;
y[45+i**l] = (p[1])*(x[54+i**k]) ;
y[46+i**l] = (p[2])*(x[55+i**k]) ;
y[47+i**l] = (p[3])*(x[56+i**k]) ;
y[48+i**l] = (p[0])*(x[57+i**k]) ;
y[49+i**l] = (p[1])*(x[58+i**k]) ;
y[50+i**l] = (p[2])*(x[59+i**k]) ;
y[51+i**l] = (p[3])*(x[60+i**k]) ; 
}
}