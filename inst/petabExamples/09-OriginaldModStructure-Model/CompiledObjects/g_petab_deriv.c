#include <R.h>
 #include <math.h>
 void g_petab_deriv_m7au6h93 ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = (p[0]/((p[0]*x[4+i**k]+p[1])*log(10.0)))*(x[11+i**k]) ;
y[1+i**l] = (p[0]/((p[0]*x[4+i**k]+p[1])*log(10.0)))*(x[18+i**k]) ;
y[2+i**l] = (p[0]/((p[0]*x[4+i**k]+p[1])*log(10.0)))*(x[25+i**k]) ;
y[3+i**l] = (p[0]/((p[0]*x[4+i**k]+p[1])*log(10.0)))*(x[32+i**k]) ;
y[4+i**l] = (p[0]/((p[0]*x[4+i**k]+p[1])*log(10.0)))*(x[39+i**k]) ;
y[5+i**l] = (p[0]/((p[0]*x[4+i**k]+p[1])*log(10.0)))*(x[46+i**k]) ;
y[6+i**l] = (p[0]/((p[0]*x[4+i**k]+p[1])*log(10.0)))*(x[53+i**k]) ;
y[7+i**l] = (p[0]/((p[0]*x[4+i**k]+p[1])*log(10.0)))*(x[60+i**k])+x[4+i**k]/((p[0]*x[4+i**k]+p[1])*log(10.0)) ;
y[8+i**l] = (p[0]/((p[0]*x[4+i**k]+p[1])*log(10.0)))*(x[67+i**k])+1.0/((p[0]*x[4+i**k]+p[1])*log(10.0)) ;
y[9+i**l] = (p[0]/((p[0]*x[4+i**k]+p[1])*log(10.0)))*(x[74+i**k]) ;
y[10+i**l] = (p[0]/((p[0]*x[4+i**k]+p[1])*log(10.0)))*(x[81+i**k]) ;
y[11+i**l] = (p[0]/((p[0]*x[4+i**k]+p[1])*log(10.0)))*(x[88+i**k]) ;
y[12+i**l] = (p[0]/((p[0]*x[4+i**k]+p[1])*log(10.0)))*(x[95+i**k]) ;
y[13+i**l] = (p[0]/((p[0]*x[4+i**k]+p[1])*log(10.0)))*(x[102+i**k]) ; 
}
}