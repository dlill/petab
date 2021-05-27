/** Code auto-generated by cOde 1.1.0 **/
#include <R.h> 
 #include <math.h> 

static double parms[32];
static double forc[0];
static double cons[0];
static double range[2];

#define nGridpoints 2 
#define nSplines 0 
#define precision 1e-05 

#define kon parms[0] 
 #define koff parms[1] 
 #define kcat parms[2] 
 #define y0_0 parms[3] 
 #define y1_0 parms[4] 
 #define y2_0 parms[5] 
 #define y3_0 parms[6] 
 #define y4_0 parms[7] 
 #define y5_0 parms[8] 
 #define y6_0 parms[9] 
 #define y7_0 parms[10] 
 #define y8_0 parms[11] 
 #define y9_0 parms[12] 
 #define y10_0 parms[13] 
 #define y11_0 parms[14] 
 #define y12_0 parms[15] 
 #define y13_0 parms[16] 
 #define y14_0 parms[17] 
 #define y15_0 parms[18] 
 #define y16_0 parms[19] 
 #define y17_0 parms[20] 
 #define y18_0 parms[21] 
 #define y19_0 parms[22] 
 #define y20_0 parms[23] 
 #define y21_0 parms[24] 
 #define y22_0 parms[25] 
 #define y23_0 parms[26] 
 #define y24_0 parms[27] 
 #define y25_0 parms[28] 
 #define y26_0 parms[29] 
 #define y27_0 parms[30] 
 #define y28_0 parms[31] 
#define tmin range[0]
#define tmax range[1]


void petab_s_initmod(void (* odeparms)(int *, double *)) {
	 int N=32;
	 odeparms(&N, parms);
}

void petab_s_initforc(void (* odeforcs)(int *, double *)) {
	 int N=0;
	 odeforcs(&N, forc);
}

/** Derivatives (ODE system) **/
void petab_s_derivs (int *n, double *t, double *y, double *ydot, double *RPAR, int *IPAR) {

	 double time = *t;

	 ydot[0] = -1.0*((kon)*y[0]*y[1])+1.0*(koff*y[2])+1.0*(kcat*y[2]);
 	 ydot[1] = -1.0*((kon)*y[0]*y[1])+1.0*(koff*y[2]);
 	 ydot[2] = 1.0*((kon)*y[0]*y[1])-1.0*(koff*y[2])-1.0*(kcat*y[2]);
 	 ydot[3] = 1.0*(kcat*y[2]);
 	 ydot[4] = (-((kon)*y[1]))*(y[4])+(-((kon)*y[0]))*(y[5])+(koff+kcat)*(y[6]);
 	 ydot[5] = (-((kon)*y[1]))*(y[4])+(-((kon)*y[0]))*(y[5])+(koff)*(y[6]);
 	 ydot[6] = ((kon)*y[1])*(y[4])+((kon)*y[0])*(y[5])+(-(koff+kcat))*(y[6]);
 	 ydot[7] = (kcat)*(y[6]);
 	 ydot[8] = (-((kon)*y[1]))*(y[8])+(-((kon)*y[0]))*(y[9])+(koff+kcat)*(y[10]);
 	 ydot[9] = (-((kon)*y[1]))*(y[8])+(-((kon)*y[0]))*(y[9])+(koff)*(y[10]);
 	 ydot[10] = ((kon)*y[1])*(y[8])+((kon)*y[0])*(y[9])+(-(koff+kcat))*(y[10]);
 	 ydot[11] = (kcat)*(y[10]);
 	 ydot[12] = (-((kon)*y[1]))*(y[12])+(-((kon)*y[0]))*(y[13])+(koff+kcat)*(y[14]);
 	 ydot[13] = (-((kon)*y[1]))*(y[12])+(-((kon)*y[0]))*(y[13])+(koff)*(y[14]);
 	 ydot[14] = ((kon)*y[1])*(y[12])+((kon)*y[0])*(y[13])+(-(koff+kcat))*(y[14]);
 	 ydot[15] = (kcat)*(y[14]);
 	 ydot[16] = (kcat)*(0.0);
 	 ydot[17] = (-((kon)*y[1]))*(y[17])+(-((kon)*y[0]))*(y[18])+(koff+kcat)*(y[19])-(y[0]*y[1]);
 	 ydot[18] = (-((kon)*y[1]))*(y[17])+(-((kon)*y[0]))*(y[18])+(koff)*(y[19])-(y[0]*y[1]);
 	 ydot[19] = ((kon)*y[1])*(y[17])+((kon)*y[0])*(y[18])+(-(koff+kcat))*(y[19])+y[0]*y[1];
 	 ydot[20] = (kcat)*(y[19]);
 	 ydot[21] = (-((kon)*y[1]))*(y[21])+(-((kon)*y[0]))*(y[22])+(koff+kcat)*(y[23])+y[2];
 	 ydot[22] = (-((kon)*y[1]))*(y[21])+(-((kon)*y[0]))*(y[22])+(koff)*(y[23])+y[2];
 	 ydot[23] = ((kon)*y[1])*(y[21])+((kon)*y[0])*(y[22])+(-(koff+kcat))*(y[23])-y[2];
 	 ydot[24] = (kcat)*(y[23]);
 	 ydot[25] = (-((kon)*y[1]))*(y[25])+(-((kon)*y[0]))*(y[26])+(koff+kcat)*(y[27])+y[2];
 	 ydot[26] = (-((kon)*y[1]))*(y[25])+(-((kon)*y[0]))*(y[26])+(koff)*(y[27]);
 	 ydot[27] = ((kon)*y[1])*(y[25])+((kon)*y[0])*(y[26])+(-(koff+kcat))*(y[27])-y[2];
 	 ydot[28] = (kcat)*(y[27])+y[2];

	 for(int i=  0 ; i <  3 ; ++i) RPAR[i] = 0;
}

