#include <R.h>
 #include <math.h>
 void g_Bachmann_MSB2011_72ujbd7i ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = log10(x[16+i**k]*p[0]/p[1]+1.0) ;
y[1+i**l] = log10(x[16+i**k]*p[2]/p[1]+1.0) ;
y[2+i**l] = log10(x[16+i**k]*p[3]/p[1]+1.0) ;
y[3+i**l] = log10(x[17+i**k]) ;
y[4+i**l] = log10(p[4]+x[17+i**k]*p[5]/p[6]) ;
y[5+i**l] = log10(x[17+i**k]*p[7]/p[6]) ;
y[6+i**l] = log10(x[17+i**k]*p[8]/p[6]) ;
y[7+i**l] = log10(x[6+i**k]+x[7+i**k]) ;
y[8+i**l] = log10(x[23+i**k]*p[9]/p[10]+1.0) ;
y[9+i**l] = log10(x[23+i**k]*p[11]/p[10]+1.0) ;
y[10+i**l] = log10(x[23+i**k]*p[12]/p[10]+1.0) ;
y[11+i**l] = log10(x[24+i**k]) ;
y[12+i**l] = log10(p[13]+x[24+i**k]*p[14]/p[15]) ;
y[13+i**l] = log10(x[8+i**k]) ;
y[14+i**l] = log10(p[16]+p[17]*(16.0*x[4+i**k]+16.0*x[2+i**k]+16.0*x[3+i**k])/p[18]) ;
y[15+i**l] = log10(p[19]+p[20]*(2.0*x[1+i**k]+2.0*x[4+i**k]+2.0*x[2+i**k]+2.0*x[3+i**k])/p[18]) ;
y[16+i**l] = p[21]+100.0*x[9+i**k]/(x[8+i**k]+x[9+i**k]) ;
y[17+i**l] = log10(p[22]+x[9+i**k]*p[23]/p[24]) ;
y[18+i**l] = log10(p[25]*(x[6+i**k]+x[7+i**k])/p[26]) ;
y[19+i**l] = log10(p[27]*(x[8+i**k]+x[9+i**k])/p[24]) ; 
}
}