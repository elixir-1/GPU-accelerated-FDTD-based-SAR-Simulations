const int NX=295, NY=295, NZ=1159;
const int NX1=NX-1, NY1=NY-1, NZ1=NZ-1;
const int mdim1=2;
const int nstop=2500;
const int ntest=2;

__device__ double ESCTC[mdim1], EINCC[mdim1], EDEVCN[mdim1], ECRLX[mdim1], ECRLY[mdim1], ECRLZ[mdim1];

__device__ int nxc, nyc, nzc;
__device__ double delx=0.9e-3, dely=0.9e-3, delz=0.9e-3;
double h_delx=0.9e-3, h_dely=0.9e-3, h_delz=0.9e-3;
__device__ double radius2=3.0e-2;

__device__ double c=3.0e8;
__device__ double freq=2400e6;
double frequency=2400e6;
int h_iden=1145;
__device__ int iden=1145;
__device__ double period, off, amp=1.0;
__device__ double ampx, ampy, ampz;
__device__ double ethinc=1.0, ephinc=0.0;
__device__ double xdisp, ydisp, zdisp;
__device__ double EPS[mdim1], SIGMA[mdim1], Ep[mdim1];
__device__ double dtedx, dtedy, dtedz, dtmdx, dtmdy, dtmdz;
__device__ double delay;
__device__ double cxd, cxu, cyd, cyu, czd, czu;
__device__ double cxx, cyy, czz, cxfyd, cxfzd, cyfxd, cyfzd, czfxd, czfyd;

__device__ double thinc=270.0, phinc=0.0;
__device__ double eps0=8.8542e-12, xmu0=1.2566306e-6;
__device__ double alpha=0.0, beta=110.0;
int nec;
int N;
double t;
double dt;
__device__ double tau;
__device__ int wbsar=1;
__device__ int onegsar=1;
