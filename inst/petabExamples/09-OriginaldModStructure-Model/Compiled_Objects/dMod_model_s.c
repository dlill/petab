/** Code auto-generated by cOde 1.1.0 **/
#include <R.h> 
 #include <math.h> 

static double parms[52];
static double forc[0];
static double cons[0];
static double range[2];

#define nGridpoints 2 
#define nSplines 0 
#define precision 1e-05 

#define k1 parms[0] 
 #define k2 parms[1] 
 #define k3 parms[2] 
 #define k4 parms[3] 
 #define k5 parms[4] 
 #define y0_0 parms[5] 
 #define y1_0 parms[6] 
 #define y2_0 parms[7] 
 #define y3_0 parms[8] 
 #define y4_0 parms[9] 
 #define y5_0 parms[10] 
 #define y6_0 parms[11] 
 #define y7_0 parms[12] 
 #define y8_0 parms[13] 
 #define y9_0 parms[14] 
 #define y10_0 parms[15] 
 #define y11_0 parms[16] 
 #define y12_0 parms[17] 
 #define y13_0 parms[18] 
 #define y14_0 parms[19] 
 #define y15_0 parms[20] 
 #define y16_0 parms[21] 
 #define y17_0 parms[22] 
 #define y18_0 parms[23] 
 #define y19_0 parms[24] 
 #define y20_0 parms[25] 
 #define y21_0 parms[26] 
 #define y22_0 parms[27] 
 #define y23_0 parms[28] 
 #define y24_0 parms[29] 
 #define y25_0 parms[30] 
 #define y26_0 parms[31] 
 #define y27_0 parms[32] 
 #define y28_0 parms[33] 
 #define y29_0 parms[34] 
 #define y30_0 parms[35] 
 #define y31_0 parms[36] 
 #define y32_0 parms[37] 
 #define y33_0 parms[38] 
 #define y34_0 parms[39] 
 #define y35_0 parms[40] 
 #define y36_0 parms[41] 
 #define y37_0 parms[42] 
 #define y38_0 parms[43] 
 #define y39_0 parms[44] 
 #define y40_0 parms[45] 
 #define y41_0 parms[46] 
 #define y42_0 parms[47] 
 #define y43_0 parms[48] 
 #define y44_0 parms[49] 
 #define y45_0 parms[50] 
 #define y46_0 parms[51] 
#define tmin range[0]
#define tmax range[1]


void dMod_model_s_initmod(void (* odeparms)(int *, double *)) {
	 int N=52;
	 odeparms(&N, parms);
}

void dMod_model_s_initforc(void (* odeforcs)(int *, double *)) {
	 int N=0;
	 odeforcs(&N, forc);
}

/** Derivatives (ODE system) **/
void dMod_model_s_derivs (int *n, double *t, double *y, double *ydot, double *RPAR, int *IPAR) {

	 double time = *t;

	 ydot[0] = -1.0*(k1*y[0]*y[1]*y[5])+1.0*(k2*y[2]);
 	 ydot[1] = -1.0*(k1*y[0]*y[1]*y[5])+1.0*(k2*y[2]);
 	 ydot[2] = 1.0*(k1*y[0]*y[1]*y[5])-1.0*(k2*y[2]);
 	 ydot[3] = -1.0*(k3*y[3]*y[2])-1.0*(k4*y[3])+1.0*(k5*y[4]);
 	 ydot[4] = 1.0*(k3*y[3]*y[2])+1.0*(k4*y[3])-1.0*(k5*y[4]);
 	 ydot[5] = 1.0*(0.0);
 	 ydot[6] = (-(k1*y[1]*y[5]))*(y[6])+(-(k1*y[0]*y[5]))*(y[7])+(k2)*(y[8])+(-(k1*y[0]*y[1]))*(0.0);
 	 ydot[7] = (-(k1*y[1]*y[5]))*(y[6])+(-(k1*y[0]*y[5]))*(y[7])+(k2)*(y[8])+(-(k1*y[0]*y[1]))*(0.0);
 	 ydot[8] = (k1*y[1]*y[5])*(y[6])+(k1*y[0]*y[5])*(y[7])+(-k2)*(y[8])+(k1*y[0]*y[1])*(0.0);
 	 ydot[9] = (-(k3*y[3]))*(y[8])+(-(k3*y[2]+k4))*(y[9])+(k5)*(y[10]);
 	 ydot[10] = (k3*y[3])*(y[8])+(k3*y[2]+k4)*(y[9])+(-k5)*(y[10]);
 	 ydot[11] = (-(k1*y[1]*y[5]))*(y[11])+(-(k1*y[0]*y[5]))*(y[12])+(k2)*(y[13])+(-(k1*y[0]*y[1]))*(0.0);
 	 ydot[12] = (-(k1*y[1]*y[5]))*(y[11])+(-(k1*y[0]*y[5]))*(y[12])+(k2)*(y[13])+(-(k1*y[0]*y[1]))*(0.0);
 	 ydot[13] = (k1*y[1]*y[5])*(y[11])+(k1*y[0]*y[5])*(y[12])+(-k2)*(y[13])+(k1*y[0]*y[1])*(0.0);
 	 ydot[14] = (-(k3*y[3]))*(y[13])+(-(k3*y[2]+k4))*(y[14])+(k5)*(y[15]);
 	 ydot[15] = (k3*y[3])*(y[13])+(k3*y[2]+k4)*(y[14])+(-k5)*(y[15]);
 	 ydot[16] = (-(k1*y[1]*y[5]))*(y[16])+(-(k1*y[0]*y[5]))*(y[17])+(k2)*(y[18])+(-(k1*y[0]*y[1]))*(0.0);
 	 ydot[17] = (-(k1*y[1]*y[5]))*(y[16])+(-(k1*y[0]*y[5]))*(y[17])+(k2)*(y[18])+(-(k1*y[0]*y[1]))*(0.0);
 	 ydot[18] = (k1*y[1]*y[5])*(y[16])+(k1*y[0]*y[5])*(y[17])+(-k2)*(y[18])+(k1*y[0]*y[1])*(0.0);
 	 ydot[19] = (-(k3*y[3]))*(y[18])+(-(k3*y[2]+k4))*(y[19])+(k5)*(y[20]);
 	 ydot[20] = (k3*y[3])*(y[18])+(k3*y[2]+k4)*(y[19])+(-k5)*(y[20]);
 	 ydot[21] = (-(k3*y[3]))*(0.0)+(-(k3*y[2]+k4))*(y[21])+(k5)*(y[22]);
 	 ydot[22] = (k3*y[3])*(0.0)+(k3*y[2]+k4)*(y[21])+(-k5)*(y[22]);
 	 ydot[23] = (-(k3*y[3]))*(0.0)+(-(k3*y[2]+k4))*(y[23])+(k5)*(y[24]);
 	 ydot[24] = (k3*y[3])*(0.0)+(k3*y[2]+k4)*(y[23])+(-k5)*(y[24]);
 	 ydot[25] = (-(k1*y[1]*y[5]))*(y[25])+(-(k1*y[0]*y[5]))*(y[26])+(k2)*(y[27])+(-(k1*y[0]*y[1]))*(y[30]);
 	 ydot[26] = (-(k1*y[1]*y[5]))*(y[25])+(-(k1*y[0]*y[5]))*(y[26])+(k2)*(y[27])+(-(k1*y[0]*y[1]))*(y[30]);
 	 ydot[27] = (k1*y[1]*y[5])*(y[25])+(k1*y[0]*y[5])*(y[26])+(-k2)*(y[27])+(k1*y[0]*y[1])*(y[30]);
 	 ydot[28] = (-(k3*y[3]))*(y[27])+(-(k3*y[2]+k4))*(y[28])+(k5)*(y[29]);
 	 ydot[29] = (k3*y[3])*(y[27])+(k3*y[2]+k4)*(y[28])+(-k5)*(y[29]);
 	 ydot[30] = 0;
 	 ydot[31] = (-(k1*y[1]*y[5]))*(y[31])+(-(k1*y[0]*y[5]))*(y[32])+(k2)*(y[33])+(-(k1*y[0]*y[1]))*(0.0)-(y[0]*y[1]*y[5]);
 	 ydot[32] = (-(k1*y[1]*y[5]))*(y[31])+(-(k1*y[0]*y[5]))*(y[32])+(k2)*(y[33])+(-(k1*y[0]*y[1]))*(0.0)-(y[0]*y[1]*y[5]);
 	 ydot[33] = (k1*y[1]*y[5])*(y[31])+(k1*y[0]*y[5])*(y[32])+(-k2)*(y[33])+(k1*y[0]*y[1])*(0.0)+y[0]*y[1]*y[5];
 	 ydot[34] = (-(k3*y[3]))*(y[33])+(-(k3*y[2]+k4))*(y[34])+(k5)*(y[35]);
 	 ydot[35] = (k3*y[3])*(y[33])+(k3*y[2]+k4)*(y[34])+(-k5)*(y[35]);
 	 ydot[36] = (-(k1*y[1]*y[5]))*(y[36])+(-(k1*y[0]*y[5]))*(y[37])+(k2)*(y[38])+(-(k1*y[0]*y[1]))*(0.0)+y[2];
 	 ydot[37] = (-(k1*y[1]*y[5]))*(y[36])+(-(k1*y[0]*y[5]))*(y[37])+(k2)*(y[38])+(-(k1*y[0]*y[1]))*(0.0)+y[2];
 	 ydot[38] = (k1*y[1]*y[5])*(y[36])+(k1*y[0]*y[5])*(y[37])+(-k2)*(y[38])+(k1*y[0]*y[1])*(0.0)-y[2];
 	 ydot[39] = (-(k3*y[3]))*(y[38])+(-(k3*y[2]+k4))*(y[39])+(k5)*(y[40]);
 	 ydot[40] = (k3*y[3])*(y[38])+(k3*y[2]+k4)*(y[39])+(-k5)*(y[40]);
 	 ydot[41] = (-(k3*y[3]))*(0.0)+(-(k3*y[2]+k4))*(y[41])+(k5)*(y[42])-(y[3]*y[2]);
 	 ydot[42] = (k3*y[3])*(0.0)+(k3*y[2]+k4)*(y[41])+(-k5)*(y[42])+y[3]*y[2];
 	 ydot[43] = (-(k3*y[3]))*(0.0)+(-(k3*y[2]+k4))*(y[43])+(k5)*(y[44])-y[3];
 	 ydot[44] = (k3*y[3])*(0.0)+(k3*y[2]+k4)*(y[43])+(-k5)*(y[44])+y[3];
 	 ydot[45] = (-(k3*y[3]))*(0.0)+(-(k3*y[2]+k4))*(y[45])+(k5)*(y[46])+y[4];
 	 ydot[46] = (k3*y[3])*(0.0)+(k3*y[2]+k4)*(y[45])+(-k5)*(y[46])-y[4];

	 for(int i=  0 ; i <  25 ; ++i) RPAR[i] = 0;
}
