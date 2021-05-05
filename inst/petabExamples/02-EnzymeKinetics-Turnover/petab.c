/** Code auto-generated by cOde 1.1.0 **/
#include <R.h> 
 #include <math.h> 

static double parms[10];
static double forc[0];
static double cons[0];
static double eventcounter[1];
static double range[2];

#define nGridpoints 2 
#define nSplines 0 
#define precision 1e-05 

#define kproS parms[0] 
 #define kdegS parms[1] 
 #define kon parms[2] 
 #define koff parms[3] 
 #define kcat parms[4] 
 #define Eadd parms[5] 
 #define y0_0 parms[6] 
 #define y1_0 parms[7] 
 #define y2_0 parms[8] 
 #define y3_0 parms[9] 
#define tmin range[0]
#define tmax range[1]


void petab_initmod(void (* odeparms)(int *, double *)) {
	 int N=10;
	 odeparms(&N, parms);
	 for(int i=0; i<1; ++i) eventcounter[i] = 0;
}

void petab_initforc(void (* odeforcs)(int *, double *)) {
	 int N=0;
	 odeforcs(&N, forc);
}

/** Derivatives (ODE system) **/
void petab_derivs (int *n, double *t, double *y, double *ydot, double *RPAR, int *IPAR) {

	 double time = *t;

	 ydot[0] = 1.0*(kproS)-1.0*(kdegS*y[0])-1.0*((kon)*y[1]*y[0])+1.0*(koff*y[2]);
 	 ydot[1] = -1.0*((kon)*y[1]*y[0])+1.0*(koff*y[2])+1.0*(kcat*y[2]);
 	 ydot[2] = 1.0*((kon)*y[1]*y[0])-1.0*(koff*y[2])-1.0*(kcat*y[2]);
 	 ydot[3] = 1.0*(kcat*y[2]);

}

/** Event function **/
void petab_myevent(int *n, double *t, double *y) {

	 double time = *t;

	 if(*t == 0 & eventcounter[0] == 0) {
		y[1] = Eadd;
		eventcounter[0] = eventcounter[0] + 1.;
	 }


}
