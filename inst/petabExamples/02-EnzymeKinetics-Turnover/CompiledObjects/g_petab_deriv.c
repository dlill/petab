#include <R.h>
 #include <math.h>
 void g_petab_deriv_drtulkcp ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = (1.0)*(x[6+i**k]) ;
y[1+i**l] = (1.0)*(x[7+i**k]) ;
y[2+i**l] = (1.0)*(x[8+i**k]) ;
y[3+i**l] = (1.0)*(x[5+i**k]) ;
y[4+i**l] = (1.0)*(x[11+i**k]) ;
y[5+i**l] = (1.0)*(x[12+i**k]) ;
y[6+i**l] = (1.0)*(x[13+i**k]) ;
y[7+i**l] = (1.0)*(x[10+i**k]) ;
y[8+i**l] = (1.0)*(x[16+i**k]) ;
y[9+i**l] = (1.0)*(x[17+i**k]) ;
y[10+i**l] = (1.0)*(x[18+i**k]) ;
y[11+i**l] = (1.0)*(x[15+i**k]) ;
y[12+i**l] = (1.0)*(x[21+i**k]) ;
y[13+i**l] = (1.0)*(x[22+i**k]) ;
y[14+i**l] = (1.0)*(x[23+i**k]) ;
y[15+i**l] = (1.0)*(x[20+i**k]) ;
y[16+i**l] = (1.0)*(x[26+i**k]) ;
y[17+i**l] = (1.0)*(x[27+i**k]) ;
y[18+i**l] = (1.0)*(x[28+i**k]) ;
y[19+i**l] = (1.0)*(x[25+i**k]) ;
y[20+i**l] = (1.0)*(x[31+i**k]) ;
y[21+i**l] = (1.0)*(x[32+i**k]) ;
y[22+i**l] = (1.0)*(x[33+i**k]) ;
y[23+i**l] = (1.0)*(x[30+i**k]) ;
y[24+i**l] = (1.0)*(x[36+i**k]) ;
y[25+i**l] = (1.0)*(x[37+i**k]) ;
y[26+i**l] = (1.0)*(x[38+i**k]) ;
y[27+i**l] = (1.0)*(x[35+i**k]) ;
y[28+i**l] = (1.0)*(x[41+i**k]) ;
y[29+i**l] = (1.0)*(x[42+i**k]) ;
y[30+i**l] = (1.0)*(x[43+i**k]) ;
y[31+i**l] = (1.0)*(x[40+i**k]) ;
y[32+i**l] = (1.0)*(x[46+i**k]) ;
y[33+i**l] = (1.0)*(x[47+i**k]) ;
y[34+i**l] = (1.0)*(x[48+i**k]) ;
y[35+i**l] = (1.0)*(x[45+i**k]) ;
y[36+i**l] = (1.0)*(x[51+i**k]) ;
y[37+i**l] = (1.0)*(x[52+i**k]) ;
y[38+i**l] = (1.0)*(x[53+i**k]) ;
y[39+i**l] = (1.0)*(x[50+i**k]) ; 
}
}