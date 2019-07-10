#include <math.h>
/*
#define NX 91
#define NY 91
#define NZ 131

#define NX1 NX-1
#define NY1 NY-1
#define NZ1 NZ-1

#define mdim1 2
#define nstop 2500

#define ntest 2
*/
const int NX=111, NY=111, NZ=231;
const int NX1=NX-1, NY1=NY-1, NZ1=NZ-1;
const int mdim1=2;
const int nstop=1000;
const int ntest=2;

int IOBS[ntest], JOBS[ntest], KOBS[ntest];

double EXS[NX][NY][NZ], EYS[NX][NY][NZ], EZS[NX][NY][NZ];
double HXS[NX][NY][NZ], HYS[NX][NY][NZ], HZS[NX][NY][NZ];
int IDONE[NX][NY][NZ], IDTWO[NX][NY][NZ], IDTHREE[NX][NY][NZ];
double EXSY1[NX1][4][NZ1], EXSY2[NX1][4][NZ1];
double EXSZ1[NX1][NY1][4], EXSZ2[NX1][NY1][4];
double EYSX1[4][NY1][NZ1], EYSX2[4][NY1][NZ1];
double EYSZ1[NX1][NY1][4], EYSZ2[NX1][NY1][4];
double EZSX1[4][NY1][NZ1], EZSX2[4][NY1][NZ1];
double EZSY1[NX1][4][NZ1], EZSY2[NX1][4][NZ1];

double ESCTC[mdim1], EINCC[mdim1], EDEVCN[mdim1], ECRLX[mdim1], ECRLY[mdim1], ECRLZ[mdim1];

int nxc, nyc, nzc;



double delx=6.0e-3, dely=6.0e-3, delz=6.0e-3;
double radius2=3.0e-2;

double c=3.0e8;
double freq=300e6;
int iden=1145;
double dt, period, off, amp=1.0;
double ampx, ampy, ampz;
double ethinc=1.0, ephinc=0.0;
double xdisp, ydisp, zdisp;
double EPS[mdim1], SIGMA[mdim1], Ep[mdim1];
double dtedx, dtedy, dtedz, dtmdx, dtmdy, dtmdz;
double delay;
double cxd, cxu, cyd, cyu, czd, czu;
double cxx, cyy, czz, cxfyd, cxfzd, cyfxd, cyfzd, czfxd, czfyd;

double thinc=270.0, phinc=0.0;
double eps0=8.8542e-12, xmu0=1.2566306e-6;
double alpha=0.0, beta=110.0;
int nec;
int N;
double tau,t;

double sar[NX][NZ], erms[NX][NZ], s[NX], erms1d[NX];
double nrms,sar_Total,sar_wb,wholevol;
int cont,cont2;
const int wbsar=1,onegsar=1;

