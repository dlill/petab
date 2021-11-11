/** Code auto-generated by cOde 1.1.0 **/
#include <R.h> 
 #include <math.h> 

static double parms[16];
static double forc[0];
static double cons[0];
static double range[2];

#define nGridpoints 2 
#define nSplines 0 
#define precision 1e-05 

#define k_AB parms[0] 
 #define k_BC parms[1] 
 #define y0_0 parms[2] 
 #define y1_0 parms[3] 
 #define y2_0 parms[4] 
 #define y3_0 parms[5] 
 #define y4_0 parms[6] 
 #define y5_0 parms[7] 
 #define y6_0 parms[8] 
 #define y7_0 parms[9] 
 #define y8_0 parms[10] 
 #define y9_0 parms[11] 
 #define y10_0 parms[12] 
 #define y11_0 parms[13] 
 #define y12_0 parms[14] 
 #define y13_0 parms[15] 
#define tmin range[0]
#define tmax range[1]


void odemodel_petab_s_initmod(void (* odeparms)(int *, double *)) {
	 int N=16;
	 odeparms(&N, parms);
}

void odemodel_petab_s_initforc(void (* odeforcs)(int *, double *)) {
	 int N=0;
	 odeforcs(&N, forc);
}

/** Derivatives (ODE system) **/
void odemodel_petab_s_derivs (int *n, double *t, double *y, double *ydot, double *RPAR, int *IPAR) {

	 double time = *t;

	 ydot[0] = -1.0*(k_AB*y[0]*1.0);
 	 ydot[1] = 1.0*(k_AB*y[0]*1.0)-1.0*(k_BC*y[1]*1.0);
 	 ydot[2] = 1.0*(k_BC*y[1]*1.0);
 	 ydot[3] = (-k_AB)*(y[3]);
 	 ydot[4] = (k_AB)*(y[3])+(-k_BC)*(y[4]);
 	 ydot[5] = (k_BC)*(y[4]);
 	 ydot[6] = (k_AB)*(0.0)+(-k_BC)*(y[6]);
 	 ydot[7] = (k_BC)*(y[6]);
 	 ydot[8] = (k_BC)*(0.0);
 	 ydot[9] = (-k_AB)*(y[9])-y[0];
 	 ydot[10] = (k_AB)*(y[9])+(-k_BC)*(y[10])+y[0];
 	 ydot[11] = (k_BC)*(y[10]);
 	 ydot[12] = (k_AB)*(0.0)+(-k_BC)*(y[12])-y[1];
 	 ydot[13] = (k_BC)*(y[12])+y[1];

	 for(int i=  0 ; i <  4 ; ++i) RPAR[i] = 0;
}

