/** Code auto-generated by cOde 1.1.0 **/
#include <R.h> 
 #include <math.h> 

static double parms[65];
static double forc[0];
static double cons[0];
static double eventcounter[1];
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
 #define y47_0 parms[52] 
 #define y48_0 parms[53] 
 #define y49_0 parms[54] 
 #define y50_0 parms[55] 
 #define y51_0 parms[56] 
 #define y52_0 parms[57] 
 #define y53_0 parms[58] 
 #define y54_0 parms[59] 
 #define y55_0 parms[60] 
 #define y56_0 parms[61] 
 #define y57_0 parms[62] 
 #define y58_0 parms[63] 
 #define y59_0 parms[64] 
#define tmin range[0]
#define tmax range[1]


void odemodel_petab_s_initmod(void (* odeparms)(int *, double *)) {
	 int N=65;
	 odeparms(&N, parms);
	 for(int i=0; i<1; ++i) eventcounter[i] = 0;
}

void odemodel_petab_s_initforc(void (* odeforcs)(int *, double *)) {
	 int N=0;
	 odeforcs(&N, forc);
}

/** Derivatives (ODE system) **/
void odemodel_petab_s_derivs (int *n, double *t, double *y, double *ydot, double *RPAR, int *IPAR) {

	 double time = *t;

	 ydot[0] = -1.0*(k1*y[0]*y[1]*y[5]*1.0)+1.0*(k2*y[2]*1.0);
 	 ydot[1] = -1.0*(k1*y[0]*y[1]*y[5]*1.0)+1.0*(k2*y[2]*1.0);
 	 ydot[2] = 1.0*(k1*y[0]*y[1]*y[5]*1.0)-1.0*(k2*y[2]*1.0);
 	 ydot[3] = -1.0*(k3*y[3]*y[2]*1.0)-1.0*(k4*y[3]*1.0)+1.0*(k5*y[4]*1.0);
 	 ydot[4] = 1.0*(k3*y[3]*y[2]*1.0)+1.0*(k4*y[3]*1.0)-1.0*(k5*y[4]*1.0);
 	 ydot[5] = 1.0*(0.0*1.0);
 	 ydot[6] = (-(k1*y[1]*y[5]))*(y[6])+(-(k1*y[0]*y[5]))*(y[7])+(k2)*(y[8])+(-(k1*y[0]*y[1]))*(y[11]);
 	 ydot[7] = (-(k1*y[1]*y[5]))*(y[6])+(-(k1*y[0]*y[5]))*(y[7])+(k2)*(y[8])+(-(k1*y[0]*y[1]))*(y[11]);
 	 ydot[8] = (k1*y[1]*y[5])*(y[6])+(k1*y[0]*y[5])*(y[7])+(-k2)*(y[8])+(k1*y[0]*y[1])*(y[11]);
 	 ydot[9] = (-(k3*y[3]))*(y[8])+(-(k3*y[2]+k4))*(y[9])+(k5)*(y[10]);
 	 ydot[10] = (k3*y[3])*(y[8])+(k3*y[2]+k4)*(y[9])+(-k5)*(y[10]);
 	 ydot[11] = 0;
 	 ydot[12] = (-(k1*y[1]*y[5]))*(y[12])+(-(k1*y[0]*y[5]))*(y[13])+(k2)*(y[14])+(-(k1*y[0]*y[1]))*(y[17]);
 	 ydot[13] = (-(k1*y[1]*y[5]))*(y[12])+(-(k1*y[0]*y[5]))*(y[13])+(k2)*(y[14])+(-(k1*y[0]*y[1]))*(y[17]);
 	 ydot[14] = (k1*y[1]*y[5])*(y[12])+(k1*y[0]*y[5])*(y[13])+(-k2)*(y[14])+(k1*y[0]*y[1])*(y[17]);
 	 ydot[15] = (-(k3*y[3]))*(y[14])+(-(k3*y[2]+k4))*(y[15])+(k5)*(y[16]);
 	 ydot[16] = (k3*y[3])*(y[14])+(k3*y[2]+k4)*(y[15])+(-k5)*(y[16]);
 	 ydot[17] = 0;
 	 ydot[18] = (-(k1*y[1]*y[5]))*(y[18])+(-(k1*y[0]*y[5]))*(y[19])+(k2)*(y[20])+(-(k1*y[0]*y[1]))*(y[23]);
 	 ydot[19] = (-(k1*y[1]*y[5]))*(y[18])+(-(k1*y[0]*y[5]))*(y[19])+(k2)*(y[20])+(-(k1*y[0]*y[1]))*(y[23]);
 	 ydot[20] = (k1*y[1]*y[5])*(y[18])+(k1*y[0]*y[5])*(y[19])+(-k2)*(y[20])+(k1*y[0]*y[1])*(y[23]);
 	 ydot[21] = (-(k3*y[3]))*(y[20])+(-(k3*y[2]+k4))*(y[21])+(k5)*(y[22]);
 	 ydot[22] = (k3*y[3])*(y[20])+(k3*y[2]+k4)*(y[21])+(-k5)*(y[22]);
 	 ydot[23] = 0;
 	 ydot[24] = (-(k1*y[1]*y[5]))*(y[24])+(-(k1*y[0]*y[5]))*(y[25])+(k2)*(y[26])+(-(k1*y[0]*y[1]))*(y[29]);
 	 ydot[25] = (-(k1*y[1]*y[5]))*(y[24])+(-(k1*y[0]*y[5]))*(y[25])+(k2)*(y[26])+(-(k1*y[0]*y[1]))*(y[29]);
 	 ydot[26] = (k1*y[1]*y[5])*(y[24])+(k1*y[0]*y[5])*(y[25])+(-k2)*(y[26])+(k1*y[0]*y[1])*(y[29]);
 	 ydot[27] = (-(k3*y[3]))*(y[26])+(-(k3*y[2]+k4))*(y[27])+(k5)*(y[28]);
 	 ydot[28] = (k3*y[3])*(y[26])+(k3*y[2]+k4)*(y[27])+(-k5)*(y[28]);
 	 ydot[29] = 0;
 	 ydot[30] = (-(k1*y[1]*y[5]))*(y[30])+(-(k1*y[0]*y[5]))*(y[31])+(k2)*(y[32])+(-(k1*y[0]*y[1]))*(y[35])-(y[0]*y[1]*y[5]);
 	 ydot[31] = (-(k1*y[1]*y[5]))*(y[30])+(-(k1*y[0]*y[5]))*(y[31])+(k2)*(y[32])+(-(k1*y[0]*y[1]))*(y[35])-(y[0]*y[1]*y[5]);
 	 ydot[32] = (k1*y[1]*y[5])*(y[30])+(k1*y[0]*y[5])*(y[31])+(-k2)*(y[32])+(k1*y[0]*y[1])*(y[35])+y[0]*y[1]*y[5];
 	 ydot[33] = (-(k3*y[3]))*(y[32])+(-(k3*y[2]+k4))*(y[33])+(k5)*(y[34]);
 	 ydot[34] = (k3*y[3])*(y[32])+(k3*y[2]+k4)*(y[33])+(-k5)*(y[34]);
 	 ydot[35] = 0;
 	 ydot[36] = (-(k1*y[1]*y[5]))*(y[36])+(-(k1*y[0]*y[5]))*(y[37])+(k2)*(y[38])+(-(k1*y[0]*y[1]))*(y[41])+y[2];
 	 ydot[37] = (-(k1*y[1]*y[5]))*(y[36])+(-(k1*y[0]*y[5]))*(y[37])+(k2)*(y[38])+(-(k1*y[0]*y[1]))*(y[41])+y[2];
 	 ydot[38] = (k1*y[1]*y[5])*(y[36])+(k1*y[0]*y[5])*(y[37])+(-k2)*(y[38])+(k1*y[0]*y[1])*(y[41])-y[2];
 	 ydot[39] = (-(k3*y[3]))*(y[38])+(-(k3*y[2]+k4))*(y[39])+(k5)*(y[40]);
 	 ydot[40] = (k3*y[3])*(y[38])+(k3*y[2]+k4)*(y[39])+(-k5)*(y[40]);
 	 ydot[41] = 0;
 	 ydot[42] = (-(k1*y[1]*y[5]))*(y[42])+(-(k1*y[0]*y[5]))*(y[43])+(k2)*(y[44])+(-(k1*y[0]*y[1]))*(y[47]);
 	 ydot[43] = (-(k1*y[1]*y[5]))*(y[42])+(-(k1*y[0]*y[5]))*(y[43])+(k2)*(y[44])+(-(k1*y[0]*y[1]))*(y[47]);
 	 ydot[44] = (k1*y[1]*y[5])*(y[42])+(k1*y[0]*y[5])*(y[43])+(-k2)*(y[44])+(k1*y[0]*y[1])*(y[47]);
 	 ydot[45] = (-(k3*y[3]))*(y[44])+(-(k3*y[2]+k4))*(y[45])+(k5)*(y[46])-(y[3]*y[2]);
 	 ydot[46] = (k3*y[3])*(y[44])+(k3*y[2]+k4)*(y[45])+(-k5)*(y[46])+y[3]*y[2];
 	 ydot[47] = 0;
 	 ydot[48] = (-(k1*y[1]*y[5]))*(y[48])+(-(k1*y[0]*y[5]))*(y[49])+(k2)*(y[50])+(-(k1*y[0]*y[1]))*(y[53]);
 	 ydot[49] = (-(k1*y[1]*y[5]))*(y[48])+(-(k1*y[0]*y[5]))*(y[49])+(k2)*(y[50])+(-(k1*y[0]*y[1]))*(y[53]);
 	 ydot[50] = (k1*y[1]*y[5])*(y[48])+(k1*y[0]*y[5])*(y[49])+(-k2)*(y[50])+(k1*y[0]*y[1])*(y[53]);
 	 ydot[51] = (-(k3*y[3]))*(y[50])+(-(k3*y[2]+k4))*(y[51])+(k5)*(y[52])-y[3];
 	 ydot[52] = (k3*y[3])*(y[50])+(k3*y[2]+k4)*(y[51])+(-k5)*(y[52])+y[3];
 	 ydot[53] = 0;
 	 ydot[54] = (-(k1*y[1]*y[5]))*(y[54])+(-(k1*y[0]*y[5]))*(y[55])+(k2)*(y[56])+(-(k1*y[0]*y[1]))*(y[59]);
 	 ydot[55] = (-(k1*y[1]*y[5]))*(y[54])+(-(k1*y[0]*y[5]))*(y[55])+(k2)*(y[56])+(-(k1*y[0]*y[1]))*(y[59]);
 	 ydot[56] = (k1*y[1]*y[5])*(y[54])+(k1*y[0]*y[5])*(y[55])+(-k2)*(y[56])+(k1*y[0]*y[1])*(y[59]);
 	 ydot[57] = (-(k3*y[3]))*(y[56])+(-(k3*y[2]+k4))*(y[57])+(k5)*(y[58])+y[4];
 	 ydot[58] = (k3*y[3])*(y[56])+(k3*y[2]+k4)*(y[57])+(-k5)*(y[58])-y[4];
 	 ydot[59] = 0;

}

/** Event function **/
void odemodel_petab_s_myevent(int *n, double *t, double *y) {

	 double time = *t;

	 if(*t == 0 & eventcounter[0] == 0) {
		y[11] = 0;
		y[17] = 0;
		y[23] = 0;
		y[29] = 0;
		y[35] = 0;
		y[41] = 0;
		y[47] = 0;
		y[53] = 0;
		y[59] = 0;
		y[5] = 1.0;
		eventcounter[0] = eventcounter[0] + 1.;
	 }


}
