#include <R.h>
 #include <math.h>
 void P_Bachmann_MSB2011_sr4kzypn ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = pow(10.0,(p[0])) ;
y[1+i**l] = p[1] ;
y[2+i**l] = p[2] ;
y[3+i**l] = p[3] ;
y[4+i**l] = p[4] ;
y[5+i**l] = p[5] ;
y[6+i**l] = pow(10.0,(p[6]))*(p[7]*pow(10.0,(p[8]))+1.0) ;
y[7+i**l] = p[9] ;
y[8+i**l] = pow(10.0,(p[10])) ;
y[9+i**l] = p[11] ;
y[10+i**l] = p[12] ;
y[11+i**l] = p[13] ;
y[12+i**l] = p[14] ;
y[13+i**l] = p[15] ;
y[14+i**l] = p[16] ;
y[15+i**l] = p[17] ;
y[16+i**l] = p[18] ;
y[17+i**l] = p[19]*pow(10.0,(p[20]))*pow(10.0,(p[21])) ;
y[18+i**l] = p[22] ;
y[19+i**l] = p[23] ;
y[20+i**l] = p[24] ;
y[21+i**l] = p[25] ;
y[22+i**l] = p[26] ;
y[23+i**l] = p[27] ;
y[24+i**l] = p[28]*pow(10.0,(p[29]))*pow(10.0,(p[30])) ;
y[25+i**l] = p[31] ;
y[26+i**l] = pow(10.0,(p[32])) ;
y[27+i**l] = pow(10.0,(p[33])) ;
y[28+i**l] = pow(10.0,(p[30])) ;
y[29+i**l] = pow(10.0,(p[34])) ;
y[30+i**l] = pow(10.0,(p[6])) ;
y[31+i**l] = pow(10.0,(p[35])) ;
y[32+i**l] = pow(10.0,(p[36])) ;
y[33+i**l] = pow(10.0,(p[37])) ;
y[34+i**l] = pow(10.0,(p[0])) ;
y[35+i**l] = pow(10.0,(p[38])) ;
y[36+i**l] = pow(10.0,(p[39])) ;
y[37+i**l] = pow(10.0,(p[40])) ;
y[38+i**l] = pow(10.0,(p[41])) ;
y[39+i**l] = pow(10.0,(p[42])) ;
y[40+i**l] = pow(10.0,(p[20])) ;
y[41+i**l] = pow(10.0,(p[43])) ;
y[42+i**l] = pow(10.0,(p[44])) ;
y[43+i**l] = pow(10.0,(p[45])) ;
y[44+i**l] = pow(10.0,(p[46])) ;
y[45+i**l] = p[47] ;
y[46+i**l] = pow(10.0,(p[10])) ;
y[47+i**l] = pow(10.0,(p[48])) ;
y[48+i**l] = pow(10.0,(p[49])) ;
y[49+i**l] = pow(10.0,(p[21])) ;
y[50+i**l] = p[50] ;
y[51+i**l] = pow(10.0,(p[51])) ;
y[52+i**l] = pow(10.0,(p[52])) ;
y[53+i**l] = pow(10.0,(p[53])) ;
y[54+i**l] = pow(10.0,(p[54])) ;
y[55+i**l] = p[55] ;
y[56+i**l] = pow(10.0,(p[29])) ;
y[57+i**l] = p[56] ;
y[58+i**l] = p[57] ;
y[59+i**l] = pow(10.0,(p[58])) ;
y[60+i**l] = pow(10.0,(p[59])) ;
y[61+i**l] = pow(10.0,(p[60])) ;
y[62+i**l] = pow(10.0,(p[61])) ;
y[63+i**l] = pow(10.0,(p[62])) ;
y[64+i**l] = pow(10.0,(p[63])) ;
y[65+i**l] = pow(10.0,(p[64])) ;
y[66+i**l] = pow(10.0,(p[65])) ;
y[67+i**l] = pow(10.0,(p[66])) ;
y[68+i**l] = pow(10.0,(p[67])) ;
y[69+i**l] = pow(10.0,(p[68])) ;
y[70+i**l] = pow(10.0,(p[69])) ;
y[71+i**l] = pow(10.0,(p[70])) ;
y[72+i**l] = pow(10.0,(p[71])) ;
y[73+i**l] = pow(10.0,(p[72])) ;
y[74+i**l] = pow(10.0,(p[73])) ;
y[75+i**l] = pow(10.0,(p[74])) ;
y[76+i**l] = pow(10.0,(p[75])) ;
y[77+i**l] = pow(10.0,(p[76])) ;
y[78+i**l] = pow(10.0,(p[77])) ;
y[79+i**l] = pow(10.0,(p[78])) ;
y[80+i**l] = pow(10.0,(p[79])) ;
y[81+i**l] = pow(10.0,(p[80])) ;
y[82+i**l] = pow(10.0,(p[81])) ;
y[83+i**l] = pow(10.0,(p[82])) ;
y[84+i**l] = pow(10.0,(p[83])) ;
y[85+i**l] = pow(10.0,(p[84])) ;
y[86+i**l] = pow(10.0,(p[85])) ;
y[87+i**l] = pow(10.0,(p[86])) ;
y[88+i**l] = pow(10.0,(p[87])) ;
y[89+i**l] = pow(10.0,(p[88])) ;
y[90+i**l] = pow(10.0,(p[89])) ;
y[91+i**l] = pow(10.0,(p[90])) ;
y[92+i**l] = pow(10.0,(p[91])) ;
y[93+i**l] = pow(10.0,(p[92])) ;
y[94+i**l] = pow(10.0,(p[93])) ;
y[95+i**l] = pow(10.0,(p[94])) ;
y[96+i**l] = pow(10.0,(p[95])) ;
y[97+i**l] = pow(10.0,(p[96])) ;
y[98+i**l] = pow(10.0,(p[97])) ;
y[99+i**l] = pow(10.0,(p[98])) ;
y[100+i**l] = pow(10.0,(p[99])) ; 
}
}