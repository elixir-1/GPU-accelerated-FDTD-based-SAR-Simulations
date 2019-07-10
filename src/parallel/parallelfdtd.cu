#include <stdio.h>
#include <math.h>
#include <cuda.h>
#include "parallel_var.h"
#include <iostream>
#include <thrust/device_vector.h>

using namespace std;

// #define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
// inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
// {
//    if (code != cudaSuccess) 
//    {
//       fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
//       if (abort) exit(code);
//    }
// }

#define pi 4.0*atan(1.0)

__global__ void zero_fields(double *EXS, double *EYS, double *EZS, double *HXS, double *HYS, double *HZS, int *IDONE, int *IDTWO, int *IDTHREE, int N){
  
  int tid=blockIdx.x*blockDim.x+threadIdx.x;

  if(tid<N){
    EXS[tid]=0.0;
    EYS[tid]=0.0;
    EZS[tid]=0.0;
    HXS[tid]=0.0;
    HYS[tid]=0.0;
    HZS[tid]=0.0;
    IDONE[tid]=0;
    IDTWO[tid]=0;
    IDTHREE[tid]=0;
  }

}

__global__ void zero_xplanes(double *EYSX1, double *EYSX2, double *EZSX1, double *EZSX2, int N){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;

  if(tid<N){
    EYSX1[tid]=0.0;
    EYSX2[tid]=0.0;
    EZSX1[tid]=0.0;
    EZSX2[tid]=0.0;
  }

}

__global__ void zero_yplanes(double *EXSY1, double *EXSY2, double *EZSY1, double *EZSY2, int N){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;

  if(tid<N){
    EXSY1[tid]=0.0;
    EXSY2[tid]=0.0;
    EZSY1[tid]=0.0;
    EZSY2[tid]=0.0;
  }

}

__global__ void zero_zplanes(double *EXSZ1, double *EXSZ2, double *EYSZ1, double *EYSZ2, int N){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;

  if(tid<N){
    EYSZ1[tid]=0.0;
    EYSZ2[tid]=0.0;
    EXSZ1[tid]=0.0;
    EXSZ2[tid]=0.0;
  }
  return;
}

__global__ void zeromdim(void){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;

  if(tid<mdim1){
    ESCTC[tid]=0.0;
    EINCC[tid]=0.0;
    EDEVCN[tid]=0.0;
    ECRLX[tid]=0.0;
    ECRLY[tid]=0.0;
    ECRLZ[tid]=0.0;
  }

  return;
}

__host__ __device__ long long int index(int i,int j, int k, int x, int y, int z){
  return (i*(y*z)+j*z+k);
}

__device__ void dcube(int **IDONE, int **IDTWO, int **IDTHREE, int istart, int jstart, int kstart, int nxwide, int nywide, int nzwide, int mtype){

  int imax, jmax, kmax;
  int i,j,k;
  
  imax = istart+nxwide-1;
  jmax = jstart+nywide-1;
  kmax = kstart+nzwide-1;

  if(nxwide==0){
    for(k=kstart;k<=kmax;k++){
      for(j=jstart;j<=jmax;j++){
	(*IDTWO)[index(istart,j,k,NX,NY,NZ)] = mtype;
	(*IDTWO)[index(istart,j,k+1,NX,NY,NZ)] = mtype;
	(*IDTHREE)[index(istart,j,k,NX,NY,NZ)] = mtype;
	(*IDTHREE)[index(istart,j,k+1,NX,NY,NZ)] = mtype;
      }
    }
  }

  else if(nywide==0){
    for(k=kstart;k<=kmax;k++){
      for(i=istart;i<=imax;i++){
	(*IDONE)[jstart+i*NY+k*(NX*NY)] = mtype;
	(*IDONE)[jstart+i*NY+(k+1)*(NX*NY)] = mtype;
	(*IDTHREE)[jstart+i*NY+k*(NX*NY)] = mtype;
	(*IDTHREE)[jstart+i*NY+(k+1)*(NX*NY)] = mtype;
      }
    }
  }

  else if(nzwide==0){
    for(j=jstart;j<=jmax;j++){
      for(i=istart;i<=imax;i++){
	(*IDONE)[index(i,j,kstart,NX,NY,NZ)] = mtype;
	(*IDONE)[index(i,j,kstart+1,NX,NY,NZ)] = mtype;
	(*IDTWO)[index(i,j,kstart,NX,NY,NZ)] = mtype;
	(*IDTWO)[index(i,j,kstart+1,NX,NY,NZ)] = mtype;
      }
    }
  }

  else{
    for(k=kstart;k<=kmax;k++){
      for(j=jstart;j<=jmax;j++){
	for(i=istart;i<=imax;i++){
	  (*IDONE)[index(i,j,k,NX,NY,NZ)] = mtype;
	  (*IDONE)[index(i,j,k+1,NX,NY,NZ)] = mtype;
	  (*IDONE)[index(i,j+1,k+1,NX,NY,NZ)] = mtype;
	  (*IDONE)[index(i,j+1,k,NX,NY,NZ)] = mtype;
	  (*IDTWO)[index(i,j,k,NX,NY,NZ)] = mtype;
	  (*IDTWO)[index(i+1,j,k,NX,NY,NZ)] = mtype;
	  (*IDTWO)[index(i+1,j,k+1,NX,NY,NZ)] = mtype;
	  (*IDTWO)[index(i,j,k+1,NX,NY,NZ)] = mtype;
	  (*IDTHREE)[index(i,j,k,NX,NY,NZ)] = mtype;
	  (*IDTHREE)[index(i+1,j,k,NX,NY,NZ)] = mtype;
	  (*IDTHREE)[index(i+1,j+1,k,NX,NY,NZ)] = mtype;
	  (*IDTHREE)[index(i,j+1,k,NX,NY,NZ)] = mtype;
	}
      }
    }
  }
  return;
  
}

__global__ void build(int NX, int NY, int NZ, int *IDONE, int *IDTWO, int *IDTHREE, long N){
  
  long tid=blockIdx.x*blockDim.x+threadIdx.x;
  
  if(tid<N){
    int mtype=2;
    double temp;
    double r2=0.5, r3=0.1;

    int i,j,k;

    k=(tid%(NY*NZ))%NZ;
    j=(tid%(NY*NZ))/NZ;
    i=tid/(NY*NZ);
    
    nxc=((NX+1)/2)-1;
    nyc=((NY+1)/2)-1;
    nzc=((NZ+1)/2)-1;

    temp=((pow((i-nxc)*delx,2)/pow(r3,2)) + (pow((j-nyc)*dely,2)/pow(r3,2)) + (pow((k-nzc)*delz,2)/pow(r2,2)));
    
    //temp=(pow((i-nxc)*delx,2) + pow((j-nyc)*dely,2) + pow((k-nzc)*delz,2));
    // r1=sqrt(temp);
    
    // if(radius2<delx || radius2<dely || radius2<delz)
    //   exit(-1);

    // if(r1<=radius2){
    //   dcube(&IDONE,&IDTWO,&IDTHREE,i,j,k,1,1,1,mtype);
    // }

    if(temp>0.0 && temp<=1.0){
      dcube(&IDONE,&IDTWO,&IDTHREE,i,j,k,1,1,1,mtype);
    }
  }
  return;
}

__global__ void setup(double *dt)
{

  // THIS SUBROUTINE INITIALIZES THE COMPUTATIONS
  double dtxi, dtyi, dtzi;
  int i;

  dtxi = c/delx;
  dtyi = c/dely;
  dtzi = c/delz;

  // CALCULATE DT--THE MAXIMUM TIME STEP ALLOWED BY THE COURANT STABILITY CONDITION
  *dt = (1.0/sqrt(dtxi*dtxi+dtyi*dtyi+dtzi*dtzi));
	
  // 	PARAMETER ALPHA IS THE DECAY RATE DETERMINED BY BETA.
  // TO CHANGE THE GAUSSIAN PULSE BY SINE WAVE WE HAVE TO MODIFY THE
  // CODE .WHERE EVER WE ARE USING THE ALPHA AND BETA WE HAVE
  // TO REMOVE IT .AND ALSO WE HAVE TO CHANGE THE TIME DURATION
  // ALSO TO EXIT THE SINE WAVE FOR THAT INTERVAL OF TIME.
	

  // in 3d fortran parameters.h
  period = 1e-6;
  
  // SET OFFSET FOR COMPUTING INCIDENT FIELDS

  off = 10;

  // THE FOLLOWING LINES ARE FOR SMOOTH COSINE INCIDENT FUNCTION
  // FIND DIRECTION COSINES FOR INCIDENT FIELD

  double costh = cos(pi*thinc/180.0);
  double sinth = sin(pi*thinc/180.0);
  double cosph = cos(pi*phinc/180.0);
  double sinph = sin(pi*phinc/180.0);
	
  // FIND AMPLITUDE OF INCIDENT FIELD COMPONENTS
	
  ampx = amp*(ethinc*cosph*costh - ephinc*sinph);
  ampy = amp*(ethinc*sinph*costh + ephinc*cosph);
  ampz = amp*(-ethinc*sinth);

  //printf("%lf %lf %lf\n",ampx,ampy,ampz);
  // FIND RELATIVE SPATIAL DELAY FOR X, Y, Z CELL DISPLACEMENT

  xdisp = -cosph*sinth;
  ydisp = -sinth*sinph;
  zdisp = -costh;
  //printf("%.15lf %.15lf %.15lf\n",xdisp,ydisp,zdisp);
  // for(i=1;i<mdim1;i++){
  //   EPS[i] = 0.0;
  //   SIGMA[i]=0.0;
  // }

  for(i=1;i<mdim1;i++){
    Ep[i]=45.5;
  }
  
  for(i = 1; i < mdim1; i++){
    EPS[i] = Ep[i]*eps0;
    SIGMA[i] = 1.75; 
    //printf("%.15lf %.15lf\n",EPS[i],SIGMA[i]);
  }
  
  // FREE SPACE:
  
  dtedx = (*dt)/(eps0*delx);
  dtedy = (*dt)/(eps0*dely);
  dtedz = (*dt)/(eps0*delz);
  dtmdx = (*dt)/(xmu0*delx);
  dtmdy = (*dt)/(xmu0*dely);
  dtmdz = (*dt)/(xmu0*delz);
  // lossy dielectrics
  
  for(i = 1; i < mdim1; i++)
    {
      ESCTC[i] = EPS[i]/(EPS[i]+SIGMA[i]*(*dt));
      EINCC[i] = SIGMA[i] *(*dt)/(EPS[i]+SIGMA[i]*(*dt));
      EDEVCN[i] = (*dt)*(EPS[i]-eps0)/(EPS[i]+SIGMA[i]*(*dt));
      ECRLX[i] = (*dt)/((EPS[i]+SIGMA[i]*(*dt))*delx);
      ECRLY[i] = (*dt)/((EPS[i]+SIGMA[i]*(*dt))*dely);
      ECRLZ[i] = (*dt)/((EPS[i]+SIGMA[i]*(*dt))*delz);
      //printf("%.15lf %.15lf %.15lf %.15lf %.15lf %.15lf\n",ESCTC[i],EINCC[i],EDEVCN[i],ECRLX[i],ECRLY[i],ECRLZ[i]);
    }

  // FIND MAXIMUM SPATIAL DELAY TO MAKE SURE PULSE PROPAGATES INTO SPACE PROPERLY.
  
  delay = 0.0;
  
  if(xdisp < 0.0)
    delay = delay-(xdisp*NX1*delx);
  if(ydisp < 0.0)
    delay = delay-(ydisp*NY1*dely);
  if(zdisp < 0.0)
    delay = delay-(zdisp*NZ1*delz);

  //printf("delay: %.25lf\n",delay);

  // COMPUTE OUTER RADIATION BOUNDARY CONDITION (ORBC) CONSTANTS:
  cxd =(c*(*dt)-delx)/(c*(*dt)+delx);
  cyd =(c*(*dt)-dely)/(c*(*dt)+dely);
  czd =(c*(*dt)-delz)/(c*(*dt)+delz);

  cxu = cxd;
  cyu = cyd;
  czu = czd;

  // COMPUTE 2ND ORDER ORBC CONSTANTS
  cxx = 2.0*delx/(c*(*dt)+delx);
  cyy = 2.0*dely/(c*(*dt)+dely);
  czz = 2.0*delz/(c*(*dt)+delz);

  cxfyd = delx*c*(*dt)*c*(*dt)/(2.0*dely*dely*(c*(*dt)+delx));
  cxfzd = delx*c*(*dt)*c*(*dt)/(2.0*delz*delz*(c*(*dt)+delx));
  cyfzd = dely*c*(*dt)*c*(*dt)/(2.0*delz*delz*(c*(*dt)+dely));
  cyfxd = dely*c*(*dt)*c*(*dt)/(2.0*delx*delx*(c*(*dt)+dely));
  czfxd = delz*c*(*dt)*c*(*dt)/(2.0*delx*delx*(c*(*dt)+delz));
  czfyd = delz*c*(*dt)*c*(*dt)/(2.0*dely*dely*(c*(*dt)+delz));	

}


__device__ double source(double dist, double t)
{
  double sourcev=0.0;
  tau=t-dist/c;
  if(tau<0.0)
    return sourcev;
  if(tau>period)
    return sourcev;
  double omega=2.0*pi*freq;
  sourcev=sin(omega*tau);
  return sourcev;
}

__device__ double dsrce(double dist, double t)
{
  double dsrcev;
  dsrcev=0.0;
  tau=t-dist/c;
  if(tau<0.0)
    return dsrcev;
  if(tau>period)
    return dsrcev;
  double omega=2.0*pi*freq;
  dsrcev=cos(omega*tau)*omega;
  return dsrcev;
}

__device__ double EXI(int i, int j, int k, double t){

  double dist;
  dist = ((i)*delx+0.5*delx*off)*xdisp+((j)*dely)*ydisp+((k)*delz)*zdisp + delay;
  double s = source(dist,t);
  return ampx*s;
}

__device__ double EYI(int i, int j, int k, double t)
{
  double dist;

  dist = ((i)*delx)*xdisp+((j)*dely+0.5*dely*off)*ydisp+((k)*delz)*zdisp + delay;

  return ampy*source(dist,t);
}

__device__ double EZI(int i, int j, int k, double t)
{
  double dist;

  dist = ((i)*delx)*xdisp+((j)*dely)*ydisp+((k)*delz+0.5*delz*off)*zdisp + delay;

  return ampz*source(dist,t);
}

__device__ double DEXI(int i, int j, int k, double t)
{
  double dist;

  dist = ((i)*delx+0.5*delx*off)*xdisp+((j)*dely)*ydisp+((k)*delz)*zdisp + delay;

  return ampx*dsrce(dist,t);
}

__device__ double DEYI(int i, int j, int k, double t)
{
  double dist;

  dist = ((i)*delx)*xdisp+((j)*dely+0.5*dely*off)*ydisp+((k)*delz)*zdisp+delay;

  return ampy*dsrce(dist,t);
}

__device__ double DEZI(int i, int j, int k, double t)
{
  double dist;

  dist = ((i)*delx)*xdisp+((j)*dely*off)*ydisp+((k)*delz+0.5*delz*off)*zdisp+delay;

  return ampz*dsrce(dist,t);
}

__global__ void eupdate(double *EXS, double *EYS, double *EZS, double *HXS, double *HYS, double *HZS, int *IDONE, int *IDTWO, int *IDTHREE, int NZ1, int NY1, int NX1, double t){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;
  int i,j,k;
  k=(tid%(NY*NZ))%NZ;
  j=(tid%(NY*NZ))/NZ;
  i=tid/(NY*NZ);

  if(k>0 && k<NZ1 && j>0 && j<NY1 && i<NX1){

    if(IDONE[tid]==0){
      EXS[tid]=EXS[tid]+(HZS[tid]-HZS[index(i,j-1,k,NX,NY,NZ)])*dtedy-(HYS[tid]-HYS[index(i,j,k-1,NX,NY,NZ)])*dtedz;
    }

    else if(IDONE[tid]==1){
      EXS[tid] = -EXI(i,j,k,t);
    }

    else if(IDONE[tid] == 2 || IDONE[tid] == 3 || IDONE[tid] == 4 || IDONE[tid] == 5){
      EXS[tid]=EXS[tid]*ESCTC[IDONE[tid]-1]-EINCC[IDONE[tid]-1]*EXI(i,j,k,t)-EDEVCN[IDONE[tid]-1]*DEXI(i,j,k,t)+(HZS[tid]-HZS[index(i,j-1,k,NX,NY,NZ)])*ECRLY[IDONE[tid]-1]-(HYS[tid]-HYS[index(i,j,k-1,NX,NY,NZ)])*ECRLZ[IDONE[tid]-1];
    }
    
  }

  if(k>0 && k<NZ1 && j<NY1 && i>0 && i<NX1){

    if(IDTWO[tid]==0){
      EYS[tid] = EYS[tid]+(HXS[tid]-HXS[index(i,j,k-1,NX,NY,NZ)])*dtedz-(HZS[tid]-HZS[index(i-1,j,k,NX,NY,NZ)])*dtedx;
    }

    else if(IDTWO[tid]==1){
      EYS[tid] = -EYI(i,j,k,t);
    }

    else if(IDTWO[tid]==2 || IDTWO[tid]==3 || IDTWO[tid]==4 || IDTWO[tid]==5){
      EYS[tid]=EYS[tid]*ESCTC[IDTWO[tid]-1]-EINCC[IDTWO[tid]-1]*EYI(i,j,k,t)-EDEVCN[IDTWO[tid]-1]*DEYI(i,j,k,t)+(HXS[tid]-HXS[index(i,j,k-1,NX,NY,NZ)])*ECRLZ[IDTWO[tid]-1]-(HZS[tid]-HZS[index(i-1,j,k,NX,NY,NZ)])*ECRLX[IDTWO[tid]-1];
    }
    
  }

  if(k<NZ1 && j>0 && j<NY1 && i>0 && i<NX1){

    if(IDTHREE[tid]==0){
      EZS[tid] = EZS[tid]+(HYS[tid]-HYS[index(i-1,j,k,NX,NY,NZ)])*dtedx-(HXS[tid]-HXS[index(i,j-1,k,NX,NY,NZ)])*dtedy;
    }

    else if(IDTHREE[tid]==1){
      EZS[tid] = -EZI(i,j,k,t);
    }

    else if(IDTHREE[tid]==2 || IDTHREE[tid]==3 || IDTHREE[tid]==4 || IDTHREE[tid]==5){
      EZS[tid]=EZS[tid]*ESCTC[IDTHREE[tid]-1]-EINCC[IDTHREE[tid]-1]*EZI(i,j,k,t)-EDEVCN[IDTHREE[tid]-1]*DEZI(i,j,k,t)+(HYS[tid]-HYS[index(i-1,j,k,NX,NY,NZ)])*ECRLX[IDTHREE[tid]-1]-(HXS[tid]-HXS[index(i,j-1,k,NX,NY,NZ)])*ECRLY[IDTHREE[tid]-1];
    }
    
  }
}

__global__ void radezx1(double *EZS, double *EZSX1, double *EZSX2, int NX, int NY, int NZ, int NX1, int NY1, int NZ1){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;
  int i,j,k;
  k=(tid%(NY*NZ))%NZ;
  j=(tid%(NY*NZ))/NZ;
  i=tid/(NY*NZ);

  if(i==0 && k<NZ1){
    j=1;
    EZS[index(i,j,k,NX,NY,NZ)] = EZSX1[index(1,j,k,4,NY1,NZ1)]+cxd*(EZS[index(1,j,k,NX,NY,NZ)]-EZSX1[index(i,j,k,4,NY1,NZ1)]);
    j = NY1-1;
    EZS[index(i,j,k,NX,NY,NZ)] = EZSX1[index(1,j,k,4,NY1,NZ1)]+ cxd*(EZS[index(1,j,k,NX,NY,NZ)] - EZSX1[index(i,j,k,4,NY1,NZ1)]);
  }

  if(i==NX-1 && k<NZ1){
    j=1;
    EZS[index(i,j,k,NX,NY,NZ)] = EZSX1[index(2,j,k,4,NY1,NZ1)] + cxu*(EZS[index(NX1-1,j,k,NX,NY,NZ)] - EZSX1[index(3,j,k,4,NY1,NZ1)]);
    j = NY1-1;
    EZS[index(i,j,k,NX,NY,NZ)] = EZSX1[index(2,j,k,4,NY1,NZ1)] + cxu*(EZS[index(NX1-1,j,k,NX,NY,NZ)] - EZSX1[index(3,j,k,4,NY1,NZ1)]);
  }

  if(i==0 && j>1 && j<NY1-1){
    k=0;
    EZS[index(i,j,k,NX,NY,NZ)] = EZSX1[index(1,j,k,4,NY1,NZ1)]+ cxd*(EZS[index(1,j,k,NX,NY,NZ)] - EZSX1[index(i,j,k,4,NY1,NZ1)]);
    k = NZ1-1;
    EZS[index(i,j,k,NX,NY,NZ)] = EZSX1[index(1,j,k,4,NY1,NZ1)]+ cxd*(EZS[index(1,j,k,NX,NY,NZ)] - EZSX1[index(i,j,k,4,NY1,NZ1)]);
  }

  if(i==NX-1 && j>1 && j<NY1-1){
    k=0;
    EZS[index(i,j,k,NX,NY,NZ)] = EZSX1[index(2,j,k,4,NY1,NZ1)] + cxu*(EZS[index(NX1-1,j,k,NX,NY,NZ)] - EZSX1[index(3,j,k,4,NY1,NZ1)]);
    k = NZ1-1;
    EZS[index(i,j,k,NX,NY,NZ)] = EZSX1[index(2,j,k,4,NY1,NZ1)] + cxu*(EZS[index(NX1-1,j,k,NX,NY,NZ)] - EZSX1[index(3,j,k,4,NY1,NZ1)]);
  }
}

__global__ void radezx2(double *EZS, double *EZSX1, double *EZSX2, int NX, int NY, int NZ, int NX1, int NY1, int NZ1){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;
  int i,j,k;
  k=(tid%(NY*NZ))%NZ;
  j=(tid%(NY*NZ))/NZ;
  i=tid/(NY*NZ);

  if(i==0 && k>0 && k<NZ1-1 && j>1 && j<NY1-1){
    EZS[index(i,j,k,NX,NY,NZ)] = -EZSX2[index(1,j,k,4,NY1,NZ1)]+cxd*(EZS[index(1,j,k,NX,NY,NZ)]+EZSX2[index(i,j,k,4,NY1,NZ1)])+cxx*(EZSX1[index(i,j,k,4,NY1,NZ1)]+EZSX1[index(1,j,k,4,NY1,NZ1)])+cxfyd*(EZSX1[index(i,j+1,k,4,NY1,NZ1)]-2.0*EZSX1[index(i,j,k,4,NY1,NZ1)]+EZSX1[index(i,j-1,k,4,NY1,NZ1)]+EZSX1[index(1,j+1,k,4,NY1,NZ1)]-2.0*EZSX1[index(1,j,k,4,NY1,NZ1)]+EZSX1[index(1,j-1,k,4,NY1,NZ1)])+cxfzd*(EZSX1[index(i,j,k+1,4,NY1,NZ1)]-2.0*EZSX1[index(i,j,k,4,NY1,NZ1)]+EZSX1[index(i,j,k-1,4,NY1,NZ1)]+EZSX1[index(1,j,k+1,4,NY1,NZ1)]-2.0*EZSX1[index(1,j,k,4,NY1,NZ1)]+EZSX1[index(1,j,k-1,4,NY1,NZ1)]);
  }

  if(i==NX-1 && k>0 && k<NZ1-1 && j>1 && j<NY1-1){
    EZS[index(i,j,k,NX,NY,NZ)] = -EZSX2[index(2,j,k,4,NY1,NZ1)]+cxd*(EZS[index(NX1-1,j,k,NX,NY,NZ)]+EZSX2[index(3,j,k,4,NY1,NZ1)])+cxx*(EZSX1[index(3,j,k,4,NY1,NZ1)]+EZSX1[index(2,j,k,4,NY1,NZ1)])+cxfyd*(EZSX1[index(3,j+1,k,4,NY1,NZ1)]-2.0*EZSX1[index(3,j,k,4,NY1,NZ1)]+EZSX1[index(3,j-1,k,4,NY1,NZ1)]+EZSX1[index(2,j+1,k,4,NY1,NZ1)]-2.0*EZSX1[index(2,j,k,4,NY1,NZ1)]+EZSX1[index(2,j-1,k,4,NY1,NZ1)])+cxfzd*(EZSX1[index(3,j,k+1,4,NY1,NZ1)]-2.0*EZSX1[index(3,j,k,4,NY1,NZ1)]+EZSX1[index(3,j,k-1,4,NY1,NZ1)]+EZSX1[index(2,j,k+1,4,NY1,NZ1)]-2.0*EZSX1[index(2,j,k,4,NY1,NZ1)]+EZSX1[index(2,j,k-1,4,NY1,NZ1)]);
  }

}

__global__ void radezx_save(double *EZS, double *EZSX1, double *EZSX2, int NX, int NY, int NZ, int NX1, int NY1, int NZ1){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;
  int i,j,k;
  k=(tid%(NY1*NZ1))%NZ1;
  j=(tid%(NY1*NZ1))/NZ1;
  i=tid/(NY1*NZ1);

  if(k<NZ1 && j>0 && j<NY1){
    EZSX2[index(i,j,k,4,NY1,NZ1)]=EZSX1[index(i,j,k,4,NY1,NZ1)];
  }

  if((i==0 || i==1) && k<NZ1 && j>0 && j<NY1){
    EZSX1[index(i,j,k,4,NY1,NZ1)]=EZS[index(i,j,k,NX,NY,NZ)];
  }

  if(i==2 && k<NZ1 && j>0 && j<NY1){
    EZSX1[index(i,j,k,4,NY1,NZ1)]=EZS[index(NX1-1,j,k,NX,NY,NZ)];
  }

  if(i==3 && k<NZ1 && j>0 && j<NY1){
    EZSX1[index(i,j,k,4,NY1,NZ1)]=EZS[index(NX-1,j,k,NX,NY,NZ)];
  }
}

__global__ void radeyx1(double *EYS, double *EYSX1, double *EYSX2, int NX, int NY, int NZ, int NX1, int NY1, int NZ1){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;
  int i,j,k;
  k=(tid%(NY*NZ))%NZ;
  j=(tid%(NY*NZ))/NZ;
  i=tid/(NY*NZ);

  if(i==0 && k>0 && k<NZ1){
    j = 0;
    EYS[index(i,j,k,NX,NY,NZ)] = EYSX1[index(1,j,k,4,NY1,NZ1)]+ cxd*(EYS[index(1,j,k,NX,NY,NZ)] - EYSX1[index(i,j,k,4,NY1,NZ1)]);
    j = NY1-1;
    EYS[index(i,j,k,NX,NY,NZ)] = EYSX1[index(1,j,k,4,NY1,NZ1)]+ cxd*(EYS[index(1,j,k,NX,NY,NZ)] - EYSX1[index(i,j,k,4,NY1,NZ1)]);
  }

  if(i==NX-1 && k>0 && k<NZ1){
    j = 0;
    EYS[index(i,j,k,NX,NY,NZ)] = EYSX1[index(2,j,k,4,NY1,NZ1)] + cxu*(EYS[index(NX1-1,j,k,NX,NY,NZ)] - EYSX1[index(3,j,k,4,NY1,NZ1)]);
    j = NY1-1;
    EYS[index(i,j,k,NX,NY,NZ)] = EYSX1[index(2,j,k,4,NY1,NZ1)] + cxu*(EYS[index(NX1-1,j,k,NX,NY,NZ)] - EYSX1[index(3,j,k,4,NY1,NZ1)]);
  }

  if(i==0 && j>0 && j<NY1-1){
    k = 1;
    EYS[index(i,j,k,NX,NY,NZ)] = EYSX1[index(1,j,k,4,NY1,NZ1)]+ cxd*(EYS[index(1,j,k,NX,NY,NZ)] - EYSX1[index(0,j,k,4,NY1,NZ1)]);
    k = NZ1-1;
    EYS[index(i,j,k,NX,NY,NZ)] = EYSX1[index(1,j,k,4,NY1,NZ1)]+ cxd*(EYS[index(1,j,k,NX,NY,NZ)] - EYSX1[index(0,j,k,4,NY1,NZ1)]);
  }

  if(i==NX-1 && j>0 && j<NY1-1){
    k = 1;
    EYS[index(i,j,k,NX,NY,NZ)] = EYSX1[index(2,j,k,4,NY1,NZ1)] + cxu*(EYS[index(NX1-1,j,k,NX,NY,NZ)] - EYSX1[index(3,j,k,4,NY1,NZ1)]);
    k = NZ1-1;
    EYS[index(i,j,k,NX,NY,NZ)] = EYSX1[index(2,j,k,4,NY1,NZ1)] + cxu*(EYS[index(NX1-1,j,k,NX,NY,NZ)] - EYSX1[index(3,j,k,4,NY1,NZ1)]);
  }
  
}

__global__ void radeyx2(double *EYS, double *EYSX1, double *EYSX2, int NX, int NY, int NZ, int NX1, int NY1, int NZ1){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;
  int i,j,k;
  k=(tid%(NY*NZ))%NZ;
  j=(tid%(NY*NZ))/NZ;
  i=tid/(NY*NZ);

  if(i==0 && k>1 && k<NZ1-1 && j>0 && j<NY1-1){
    EYS[index(i,j,k,NX,NY,NZ)] = -EYSX2[index(1,j,k,4,NY1,NZ1)]+cxd*(EYS[index(1,j,k,NX,NY,NZ)]+EYSX2[index(i,j,k,4,NY1,NZ1)])+cxx*(EYSX1[index(i,j,k,4,NY1,NZ1)]+EYSX1[index(1,j,k,4,NY1,NZ1)])+cxfyd*(EYSX1[index(i,j+1,k,4,NY1,NZ1)]-2.0*EYSX1[index(i,j,k,4,NY1,NZ1)]+EYSX1[index(i,j-1,k,4,NY1,NZ1)]+EYSX1[index(1,j+1,k,4,NY1,NZ1)]-2.0*EYSX1[index(1,j,k,4,NY1,NZ1)]+EYSX1[index(1,j-1,k,4,NY1,NZ1)])+cxfzd*(EYSX1[index(i,j,k+1,4,NY1,NZ1)]-2.0*EYSX1[index(i,j,k,4,NY1,NZ1)]+EYSX1[index(i,j,k-1,4,NY1,NZ1)]+EYSX1[index(1,j,k+1,4,NY1,NZ1)]-2.0*EYSX1[index(1,j,k,4,NY1,NZ1)]+EYSX1[index(1,j,k-1,4,NY1,NZ1)]);
  }

  if(i==NX-1 && k>1 && k<NZ1-1 && j>0 && j<NY1-1){
    EYS[index(i,j,k,NX,NY,NZ)] = -EYSX2[index(2,j,k,4,NY1,NZ1)]+cxd*(EYS[index(NX1-1,j,k,NX,NY,NZ)]+EYSX2[index(3,j,k,4,NY1,NZ1)])+cxx*(EYSX1[index(3,j,k,4,NY1,NZ1)]+EYSX1[index(2,j,k,4,NY1,NZ1)])+cxfyd*(EYSX1[index(3,j+1,k,4,NY1,NZ1)]-2.0*EYSX1[index(3,j,k,4,NY1,NZ1)]+EYSX1[index(3,j-1,k,4,NY1,NZ1)]+EYSX1[index(2,j+1,k,4,NY1,NZ1)]-2.0*EYSX1[index(2,j,k,4,NY1,NZ1)]+EYSX1[index(2,j-1,k,4,NY1,NZ1)])+cxfzd*(EYSX1[index(3,j,k+1,4,NY1,NZ1)]-2.0*EYSX1[index(3,j,k,4,NY1,NZ1)]+EYSX1[index(3,j,k-1,4,NY1,NZ1)]+EYSX1[index(2,j,k+1,4,NY1,NZ1)]-2.0*EYSX1[index(2,j,k,4,NY1,NZ1)]+EYSX1[index(2,j,k-1,4,NY1,NZ1)]);
  }
  
}

__global__ void radeyx_save(double *EYS, double *EYSX1, double *EYSX2, int NX, int NY, int NZ, int NX1, int NY1, int NZ1){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;
  int i,j,k;
  k=(tid%(NY1*NZ1))%NZ1;
  j=(tid%(NY1*NZ1))/NZ1;
  i=tid/(NY1*NZ1);

  if(k>0 && k<NZ1 && j<NY1){
    EYSX2[index(i,j,k,4,NY1,NZ1)]=EYSX1[index(i,j,k,4,NY1,NZ1)];
  }

  if((i==0 || i==1) && k>0 && k<NZ1 && j<NY1){
    EYSX1[index(i,j,k,4,NY1,NZ1)]=EYS[index(i,j,k,NX,NY,NZ)];
  }

  if(i==2 && k>0 && k<NZ1 && j<NY1){
    EYSX1[index(i,j,k,4,NY1,NZ1)]=EYS[index(NX1-1,j,k,NX,NY,NZ)];
  }

  if(i==3 && k>0 && k<NZ1 && j<NY1){
    EYSX1[index(i,j,k,4,NY1,NZ1)]=EYS[index(NX-1,j,k,NX,NY,NZ)];
  }

}

__global__ void radezy1(double *EZS, double *EZSY1, double *EZSY2, int NX, int NY, int NZ, int NX1, int NY1, int NZ1){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;
  int i,j,k;
  k=(tid%(NY*NZ))%NZ;
  j=(tid%(NY*NZ))/NZ;
  i=tid/(NY*NZ);

  if(j==0 && k<NZ1){
    i = 1;
    EZS[index(i,j,k,NX,NY,NZ)] = EZSY1[index(i,1,k,NX1,4,NZ1)]+ cyd*(EZS[index(i,1,k,NX,NY,NZ)] - EZSY1[index(i,j,k,NX1,4,NZ1)]);
    i = NX1-1;
    EZS[index(i,j,k,NX,NY,NZ)] = EZSY1[index(i,1,k,NX1,4,NZ1)]+ cyd*(EZS[index(i,1,k,NX,NY,NZ)] - EZSY1[index(i,j,k,NX1,4,NZ1)]);
  }

  if(j==(NY-1) && k<NZ1){
    i = 1;
    EZS[index(i,j,k,NX,NY,NZ)] = EZSY1[index(i,2,k,NX1,4,NZ1)] + cyd*(EZS[index(i,NY1-1,k,NX,NY,NZ)] - EZSY1[index(i,3,k,NX1,4,NZ1)]);
    i = NX1-1;
    EZS[index(i,j,k,NX,NY,NZ)] = EZSY1[index(i,2,k,NX1,4,NZ1)] + cyd*(EZS[index(i,NY1-1,k,NX,NY,NZ)] - EZSY1[index(i,3,k,NX1,4,NZ1)]);
  }

  
  if(j==0 && i>1 && i<NX1-1){
    k = 0;
    EZS[index(i,j,k,NX,NY,NZ)] = EZSY1[index(i,1,k,NX1,4,NZ1)]+ cyd*(EZS[index(i,1,k,NX,NY,NZ)] - EZSY1[index(i,j,k,NX1,4,NZ1)]);
    k = NZ1-1;
    EZS[index(i,j,k,NX,NY,NZ)] = EZSY1[index(i,1,k,NX1,4,NZ1)]+ cyd*(EZS[index(i,1,k,NX,NY,NZ)] - EZSY1[index(i,j,k,NX1,4,NZ1)]);
  }

  if(j==(NY-1) && i>1 && i<NX1-1){
    k = 0;
    EZS[index(i,j,k,NX,NY,NZ)] = EZSY1[index(i,2,k,NX1,4,NZ1)] + cyd*(EZS[index(i,NY1-1,k,NX,NY,NZ)] - EZSY1[index(i,3,k,NX1,4,NZ1)]);
    k = NZ1-1;
    EZS[index(i,j,k,NX,NY,NZ)] = EZSY1[index(i,2,k,NX1,4,NZ1)] + cyd*(EZS[index(i,NY1-1,k,NX,NY,NZ)] - EZSY1[index(i,3,k,NX1,4,NZ1)]);
  }
      
}

__global__ void radezy2(double *EZS, double *EZSY1, double *EZSY2, int NX, int NY, int NZ, int NX1, int NY1, int NZ1){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;
  int i,j,k;
  k=(tid%(NY*NZ))%NZ;
  j=(tid%(NY*NZ))/NZ;
  i=tid/(NY*NZ);

  if(j==0 && k>0 && k<NZ1-1 && i>1 && i<NX1-1){
    
    EZS[index(i,j,k,NX,NY,NZ)] = -EZSY2[index(i,1,k,NX1,4,NZ1)]+cyd*(EZS[index(i,1,k,NX,NY,NZ)]+EZSY2[index(i,j,k,NX1,4,NZ1)])+cyy*(EZSY1[index(i,j,k,NX1,4,NZ1)]+EZSY1[index(i,1,k,NX1,4,NZ1)])+cyfxd*(EZSY1[index(i+1,j,k,NX1,4,NZ1)]-2.0*EZSY1[index(i,j,k,NX1,4,NZ1)]+EZSY1[index(i-1,j,k,NX1,4,NZ1)]+EZSY1[index(i+1,1,k,NX1,4,NZ1)]-2.0*EZSY1[index(i,1,k,NX1,4,NZ1)]+EZSY1[index(i-1,1,k,NX1,4,NZ1)])+cyfzd*(EZSY1[index(i,j,k+1,NX1,4,NZ1)]-2.0*EZSY1[index(i,j,k,NX1,4,NZ1)]+EZSY1[index(i,j,k-1,NX1,4,NZ1)]+EZSY1[index(i,1,k+1,NX1,4,NZ1)]-2.0*EZSY1[index(i,1,k,NX1,4,NZ1)]+EZSY1[index(i,1,k-1,NX1,4,NZ1)]);
  }
  
  if(j==(NY-1) && k>0 && k<NZ1-1 && i>1 && i<NX1-1){
    EZS[index(i,j,k,NX,NY,NZ)] = -EZSY2[index(i,2,k,NX1,4,NZ1)]+cyd*(EZS[index(i,NY1-1,k,NX,NY,NZ)]+EZSY2[index(i,3,k,NX1,4,NZ1)])+cyy*(EZSY1[index(i,3,k,NX1,4,NZ1)]+EZSY1[index(i,2,k,NX1,4,NZ1)])+cyfxd*(EZSY1[index(i+1,3,k,NX1,4,NZ1)]-2.0*EZSY1[index(i,3,k,NX1,4,NZ1)]+EZSY1[index(i-1,3,k,NX1,4,NZ1)]+EZSY1[index(i+1,2,k,NX1,4,NZ1)]-2.0*EZSY1[index(i,2,k,NX1,4,NZ1)]+EZSY1[index(i-1,2,k,NX1,4,NZ1)])+ cyfzd*(EZSY1[index(i,3,k+1,NX1,4,NZ1)]-2.0*EZSY1[index(i,3,k,NX1,4,NZ1)]+EZSY1[index(i,3,k-1,NX1,4,NZ1)]+EZSY1[index(i,2,k+1,NX1,4,NZ1)]-2.0*EZSY1[index(i,2,k,NX1,4,NZ1)]+EZSY1[index(i,2,k-1,NX1,4,NZ1)]);
  }
}

__global__ void radezy_save(double *EZS, double *EZSY1, double *EZSY2, int NX, int NY, int NZ, int NX1, int NY1, int NZ1){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;
  int i,j,k;
  k=(tid%(4*NZ1))%NZ1;
  j=(tid%(4*NZ1))/NZ1;
  i=tid/(4*NZ1);

  if(k<NZ1 && i>0 && i<NX1){
    EZSY2[index(i,j,k,NX1,4,NZ1)]=EZSY1[index(i,j,k,NX1,4,NZ1)];
  }

  if((j==0 || j==1) && k<NZ1 && i>0 && i<NX1){
    EZSY1[index(i,j,k,NX1,4,NZ1)]=EZS[index(i,j,k,NX,NY,NZ)];
  }

  if(j==2 && k<NZ1 && i>0 && i<NX1){
    EZSY1[index(i,j,k,NX1,4,NZ1)]=EZS[index(i,NY1-1,k,NX,NY,NZ)];
  }

  if(j==3 && k<NZ1 && i>0 && i<NX1){
    EZSY1[index(i,j,k,NX1,4,NZ1)]=EZS[index(i,NY-1,k,NX,NY,NZ)];
  }

}

__global__ void radexy1(double *EXS, double *EXSY1, double *EXSY2, int NX, int NY, int NZ, int NX1, int NY1, int NZ1){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;
  int i,j,k;
  k=(tid%(NY*NZ))%NZ;
  j=(tid%(NY*NZ))/NZ;
  i=tid/(NY*NZ);

  if(j==0 && k>0 && k<NZ1){
    i = 0;
    EXS[index(i,j,k,NX,NY,NZ)] = EXSY1[index(i,1,k,NX1,4,NZ1)]+ cyd*(EXS[index(i,1,k,NX,NY,NZ)] - EXSY1[index(i,j,k,NX1,4,NZ1)]);
    i = NX1-1;
    EXS[index(i,j,k,NX,NY,NZ)] = EXSY1[index(i,1,k,NX1,4,NZ1)]+ cyd*(EXS[index(i,1,k,NX,NY,NZ)] - EXSY1[index(i,j,k,NX1,4,NZ1)]);
  }

  if(j==(NY-1) && k>0 && k<NZ1){
    i = 0;
    EXS[index(i,j,k,NX,NY,NZ)] = EXSY1[index(i,2,k,NX1,4,NZ1)] + cyd*(EXS[index(i,NY1-1,k,NX,NY,NZ)] - EXSY1[index(i,3,k,NX1,4,NZ1)]);
    i = NX1-1;
    EXS[index(i,j,k,NX,NY,NZ)] = EXSY1[index(i,2,k,NX1,4,NZ1)] + cyd*(EXS[index(i,NY1-1,k,NX,NY,NZ)] - EXSY1[index(i,3,k,NX1,4,NZ1)]);
  }

  
  if(j==0 && i>0 && i<NX1-1){
    
    k = 1;
    EXS[index(i,j,k,NX,NY,NZ)] = EXSY1[index(i,1,k,NX1,4,NZ1)]+ cyd*(EXS[index(i,1,k,NX,NY,NZ)] - EXSY1[index(i,j,k,NX1,4,NZ1)]);
    k = NZ1-1;
    EXS[index(i,j,k,NX,NY,NZ)] = EXSY1[index(i,1,k,NX1,4,NZ1)]+ cyd*(EXS[index(i,1,k,NX,NY,NZ)] - EXSY1[index(i,j,k,NX1,4,NZ1)]);
  }

  if(j==(NY-1) && i>0 && i<NX1-1){
    k = 1;
    EXS[index(i,j,k,NX,NY,NZ)] = EXSY1[index(i,2,k,NX1,4,NZ1)] + cyd*(EXS[index(i,NY1-1,k,NX,NY,NZ)] - EXSY1[index(i,3,k,NX1,4,NZ1)]);
    k = NZ1-1;
    EXS[index(i,j,k,NX,NY,NZ)] = EXSY1[index(i,2,k,NX1,4,NZ1)] + cyd*(EXS[index(i,NY1-1,k,NX,NY,NZ)] - EXSY1[index(i,3,k,NX1,3,NZ1)]);
  }
}

__global__ void radexy2(double *EXS, double *EXSY1, double *EXSY2, int NX, int NY, int NZ, int NX1, int NY1, int NZ1){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;
  int i,j,k;
  k=(tid%(NY*NZ))%NZ;
  j=(tid%(NY*NZ))/NZ;
  i=tid/(NY*NZ);

  if(j==0 && k>1 && k<NZ1-1 && i>0 && i<NX1-1){
    EXS[index(i,j,k,NX,NY,NZ)] = -EXSY2[index(i,1,k,NX1,4,NZ1)]+cyd*(EXS[index(i,1,k,NX,NY,NZ)]+EXSY2[index(i,j,k,NX1,4,NZ1)])+cyy*(EXSY1[index(i,j,k,NX1,4,NZ1)]+EXSY1[index(i,1,k,NX1,4,NZ1)])+cyfxd*(EXSY1[index(i+1,j,k,NX1,4,NZ1)]-2.0*EXSY1[index(i,j,k,NX1,4,NZ1)]+EXSY1[index(i-1,j,k,NX1,4,NZ1)]+EXSY1[index(i+1,1,k,NX1,4,NZ1)]-2.0*EXSY1[index(i,1,k,NX1,4,NZ1)]+EXSY1[index(i-1,1,k,NX1,4,NZ1)])+cyfzd*(EXSY1[index(i,j,k+1,NX1,4,NZ1)]-2.0*EXSY1[index(i,j,k,NX1,4,NZ1)]+EXSY1[index(i,j,k-1,NX1,4,NZ1)]+EXSY1[index(i,1,k+1,NX1,4,NZ1)]-2.0*EXSY1[index(i,1,k,NX1,4,NZ1)]+EXSY1[index(i,1,k-1,NX1,4,NZ1)]);
  }

  if(j==(NY-1) && k>1 && k<NZ1-1 && i>0 && i<NX1-1){
    EXS[index(i,j,k,NX,NY,NZ)] = -EXSY2[index(i,2,k,NX1,4,NZ1)]+cyd*(EXS[index(i,NY1-1,k,NX,NY,NZ)]+EXSY2[index(i,3,k,NX1,4,NZ1)])+cyy*(EXSY1[index(i,3,k,NX1,4,NZ1)]+EXSY1[index(i,2,k,NX1,4,NZ1)])+cyfxd*(EXSY1[index(i+1,3,k,NX1,4,NZ1)]-2.0*EXSY1[index(i,3,k,NX1,4,NZ1)]+EXSY1[index(i-1,3,k,NX1,4,NZ1)]+EXSY1[index(i+1,2,k,NX1,4,NZ1)]-2.0*EXSY1[index(i,2,k,NX1,4,NZ1)]+EXSY1[index(i-1,2,k,NX1,4,NZ1)])+cyfzd*(EXSY1[index(i,3,k+1,NX1,4,NZ1)]-2.0*EXSY1[index(i,3,k,NX1,4,NZ1)]+EXSY1[index(i,3,k-1,NX1,4,NZ1)]+EXSY1[index(i,2,k+1,NX1,4,NZ1)]-2.0*EXSY1[index(i,2,k,NX1,4,NZ1)]+EXSY1[index(i,2,k-1,NX1,4,NZ1)]);
  }
}

__global__ void radexy_save(double *EXS, double *EXSY1, double *EXSY2, int NX, int NY, int NZ, int NX1, int NY1, int NZ1){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;
  int i,j,k;
  k=(tid%(4*NZ1))%NZ1;
  j=(tid%(4*NZ1))/NZ1;
  i=tid/(4*NZ1);

  if(k>0 && k<NZ1 && i<NX1){
    EXSY2[index(i,j,k,NX1,4,NZ1)]=EXSY1[index(i,j,k,NX1,4,NZ1)];
  }

  if((j==0 || j==1) && k>0 && k<NZ1 && i<NX1){
    EXSY1[index(i,j,k,NX1,4,NZ1)]=EXS[index(i,j,k,NX,NY,NZ)];
  }

  if(j==2 && k>0 && k<NZ1 && i<NX1){
    EXSY1[index(i,j,k,NX1,4,NZ1)]=EXS[index(i,NY1-1,k,NX,NY,NZ)];
  }

  if(j==3 && k>0 && k<NZ1 && i<NX1){
    EXSY1[index(i,j,k,NX1,4,NZ1)]=EXS[index(i,NY-1,k,NX,NY,NZ)];
  }
}

__global__ void radexz1(double *EXS, double *EXSZ1, double *EXSZ2, int NX, int NY, int NZ, int NX1, int NY1, int NZ1){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;
  int i,j,k;
  k=(tid%(NY*NZ))%NZ;
  j=(tid%(NY*NZ))/NZ;
  i=tid/(NY*NZ);

  if(k==0 && j>0 && j<NY1){
    i = 0;
    EXS[index(i,j,k,NX,NY,NZ)] = EXSZ1[index(i,j,1,NX1,NY1,4)]+ czd*(EXS[index(i,j,1,NX,NY,NZ)] - EXSZ1[index(i,j,k,NX1,NY1,4)]);
    i = NX1-1;
    EXS[index(i,j,k,NX,NY,NZ)] = EXSZ1[index(i,j,1,NX1,NY1,4)]+ czd*(EXS[index(i,j,1,NX,NY,NZ)] - EXSZ1[index(i,j,k,NX1,NY1,4)]);
  }

  if(k==(NZ-1) && j>0 && j<NY1){
    i = 0;
    EXS[index(i,j,k,NX,NY,NZ)] = EXSZ1[index(i,j,2,NX1,NY1,4)] + czd*(EXS[index(i,j,NZ1-1,NX,NY,NZ)] - EXSZ1[index(i,j,3,NX1,NY1,4)]);
    i = NX1-1;
    EXS[index(i,j,k,NX,NY,NZ)] = EXSZ1[index(i,j,2,NX1,NY1,4)] + czd*(EXS[index(i,j,NZ1-1,NX,NY,NZ)] - EXSZ1[index(i,j,3,NX1,NY1,4)]);
  }

  if(k==0 && i>0 && i<NX1-1){
    j = 1;
    EXS[index(i,j,k,NX,NY,NZ)] = EXSZ1[index(i,j,1,NX1,NY1,4)]+ czd*(EXS[index(i,j,1,NX,NY,NZ)] - EXSZ1[index(i,j,k,NX1,NY1,4)]);
    j = NY1-1;
    EXS[index(i,j,k,NX,NY,NZ)] = EXSZ1[index(i,j,1,NX1,NY1,4)]+ czd*(EXS[index(i,j,1,NX,NY,NZ)] - EXSZ1[index(i,j,k,NX1,NY1,4)]);
  }

  if(k==(NZ-1) && i>0 && i<NX1-1){
    j = 1;
    EXS[index(i,j,k,NX,NY,NZ)] = EXSZ1[index(i,j,2,NX1,NY1,4)] + czd*(EXS[index(i,j,NZ1-1,NX,NY,NZ)] - EXSZ1[index(i,j,3,NX1,NY1,4)]);
    j = NY1-1;
    EXS[index(i,j,k,NX,NY,NZ)] = EXSZ1[index(i,j,2,NX1,NY1,4)] + czd*(EXS[index(i,j,NZ1-1,NX,NY,NZ)] - EXSZ1[index(i,j,3,NX1,NY1,4)]);
  }
}

__global__ void radexz2(double *EXS, double *EXSZ1, double *EXSZ2, int NX, int NY, int NZ, int NX1, int NY1, int NZ1){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;
  int i,j,k;
  k=(tid%(NY*NZ))%NZ;
  j=(tid%(NY*NZ))/NZ;
  i=tid/(NY*NZ);

  if(k==0 && j>1 && j<NY1-1 && i>0 && i<NX1-1){
    EXS[index(i,j,k,NX,NY,NZ)] = -EXSZ2[index(i,j,1,NX1,NY1,4)]+czd*(EXS[index(i,j,1,NX,NY,NZ)]+EXSZ2[index(i,j,k,NX1,NY1,4)])+czz*(EXSZ1[index(i,j,k,NX1,NY1,4)]+EXSZ1[index(i,j,1,NX1,NY1,4)])+czfxd*(EXSZ1[index(i+1,j,k,NX1,NY1,4)]-2.0*EXSZ1[index(i,j,k,NX1,NY1,4)]+EXSZ1[index(i-1,j,k,NX1,NY1,4)]+EXSZ1[index(i+1,j,1,NX1,NY1,4)]-2.0*EXSZ1[index(i,j,1,NX1,NY1,4)]+EXSZ1[index(i-1,j,1,NX1,NY1,4)])+czfyd*(EXSZ1[index(i,j+1,k,NX1,NY1,4)]-2.0*EXSZ1[index(i,j,k,NX1,NY1,4)]+EXSZ1[index(i,j-1,k,NX1,NY1,4)]+EXSZ1[index(i,j+1,1,NX1,NY1,4)]-2.0*EXSZ1[index(i,j,1,NX1,NY1,4)]+EXSZ1[index(i,j-1,1,NX1,NY1,4)]);
  }

  if(k==(NZ-1) && j>1 && j<NY1-1 && i>0 && i<NX1-1){
    EXS[index(i,j,k,NX,NY,NZ)] = -EXSZ2[index(i,j,2,NX1,NY1,4)]+czd*(EXS[index(i,j,NZ1-1,NX,NY,NZ)]+EXSZ2[index(i,j,3,NX1,NY1,4)])+czz*(EXSZ1[index(i,j,3,NX1,NY1,4)]+EXSZ1[index(i,j,2,NX1,NY1,4)])+czfxd*(EXSZ1[index(i+1,j,3,NX1,NY1,4)]-2.0*EXSZ1[index(i,j,3,NX1,NY1,4)]+EXSZ1[index(i-1,j,3,NX1,NY1,4)]+EXSZ1[index(i+1,j,2,NX1,NY1,4)]-2.0*EXSZ1[index(i,j,2,NX1,NY1,4)]+EXSZ1[index(i-1,j,2,NX1,NY1,4)])+czfyd*(EXSZ1[index(i,j+1,3,NX1,NY1,4)]-2.0*EXSZ1[index(i,j,3,NX1,NY1,4)]+EXSZ1[index(i,j-1,3,NX1,NY1,4)]+EXSZ1[index(i,j+1,2,NX1,NY1,4)]-2.0*EXSZ1[index(i,j,2,NX1,NY1,4)]+EXSZ1[index(i,j-1,2,NX1,NY1,4)]);
  }
}

__global__ void radexz_save(double *EXS, double *EXSZ1, double *EXSZ2, int NX, int NY, int NZ, int NX1, int NY1, int NZ1){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;
  int i,j,k;
  k=(tid%(NY1*4))%4;
  j=(tid%(NY1*4))/4;
  i=tid/(NY1*4);

  if(j>0 && j<NY1 && i<NX1){
    EXSZ2[index(i,j,k,NX1,NY1,4)]=EXSZ1[index(i,j,k,NX1,NY1,4)];
  }

  if((k==0 || k==1) && j>0 && j<NY1 && i<NX1){
    EXSZ1[index(i,j,k,NX1,NY1,4)]=EXS[index(i,j,k,NX,NY,NZ)];
  }

  if(k==2 && j>0 && j<NY1 && i<NX1){
    EXSZ1[index(i,j,k,NX1,NY1,4)]=EXS[index(i,j,NZ1-1,NX,NY,NZ)];
  }

  if(k==3 && j>0 && j<NY1 && i<NX1){
    EXSZ1[index(i,j,k,NX1,NY1,4)]=EXS[index(i,j,NZ-1,NX,NY,NZ)];
  }
}

__global__ void radeyz1(double *EYS, double *EYSZ1, double *EYSZ2, int NX, int NY, int NZ, int NX1, int NY1, int NZ1){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;
  int i,j,k;
  k=(tid%(NY*NZ))%NZ;
  j=(tid%(NY*NZ))/NZ;
  i=tid/(NY*NZ);

  if(k==0 && j<NY1){
    i = 1;
    EYS[index(i,j,k,NX,NY,NZ)] = EYSZ1[index(i,j,1,NX1,NY1,4)]+ czd*(EYS[index(i,j,1,NX,NY,NZ)] - EYSZ1[index(i,j,k,NX1,NY1,4)]);
    i = NX1-1;
    EYS[index(i,j,k,NX,NY,NZ)] = EYSZ1[index(i,j,1,NX1,NY1,4)]+ czd*(EYS[index(i,j,1,NX,NY,NZ)] - EYSZ1[index(i,j,k,NX1,NY1,4)]);
  }

  if(k==(NZ-1) && j<NY1){
    i = 1;
    EYS[index(i,j,k,NX,NY,NZ)] = EYSZ1[index(i,j,2,NX1,NY1,4)] + czd*(EYS[index(i,j,NZ1-1,NX,NY,NZ)] - EYSZ1[index(i,j,3,NX1,NY1,4)]);
    i = NX1-1;
    EYS[index(i,j,k,NX,NY,NZ)] = EYSZ1[index(i,j,2,NX1,NY1,4)] + czd*(EYS[index(i,j,NZ1-1,NX,NY,NZ)] - EYSZ1[index(i,j,3,NX1,NY1,4)]);
  }

  if(k==0 && i>1 && i<NX1-1){
    j = 0;
    EYS[index(i,j,k,NX,NY,NZ)] = EYSZ1[index(i,j,1,NX1,NY1,4)]+ czd*(EYS[index(i,j,1,NX,NY,NZ)] - EYSZ1[index(i,j,k,NX1,NY1,4)]);
    j = NY1-1;
    EYS[index(i,j,k,NX,NY,NZ)] = EYSZ1[index(i,j,1,NX1,NY1,4)]+ czd*(EYS[index(i,j,1,NX,NY,NZ)] - EYSZ1[index(i,j,k,NX1,NY1,4)]);
  }

  if(k==(NZ-1) && i>1 && i<NX1-1){
    j = 0;
    EYS[index(i,j,k,NX,NY,NZ)] = EYSZ1[index(i,j,2,NX1,NY1,4)] + czd*(EYS[index(i,j,NZ1-1,NX,NY,NZ)] - EYSZ1[index(i,j,3,NX1,NY1,4)]);
    j = NY1-1;
    EYS[index(i,j,k,NX,NY,NZ)] = EYSZ1[index(i,j,2,NX1,NY1,4)] + czd*(EYS[index(i,j,NZ1-1,NX,NY,NZ)] - EYSZ1[index(i,j,3,NX1,NY1,4)]);
  }
}

__global__ void radeyz2(double *EYS, double *EYSZ1, double *EYSZ2, int NX, int NY, int NZ, int NX1, int NY1, int NZ1){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;
  int i,j,k;
  k=(tid%(NY*NZ))%NZ;
  j=(tid%(NY*NZ))/NZ;
  i=tid/(NY*NZ);

  if(k==0 && j>0 && j<NY1-1 && i>1 && i<NX1-1){
    EYS[index(i,j,k,NX,NY,NZ)] = -EYSZ2[index(i,j,1,NX1,NY1,4)]+ czd*(EYS[index(i,j,1,NX,NY,NZ)]+EYSZ2[index(i,j,k,NX1,NY1,4)])+czz*(EYSZ1[index(i,j,k,NX1,NY1,4)]+EYSZ1[index(i,j,1,NX1,NY1,4)])+ czfxd*(EYSZ1[index(i+1,j,k,NX1,NY1,4)]-2.0*EYSZ1[index(i,j,k,NX1,NY1,4)]+EYSZ1[index(i-1,j,k,NX1,NY1,4)]+EYSZ1[index(i+1,j,1,NX1,NY1,4)]- 2.0*EYSZ1[index(i,j,1,NX1,NY1,4)]+EYSZ1[index(i-1,j,1,NX1,NY1,4)])+ czfyd*(EYSZ1[index(i,j+1,k,NX1,NY1,4)]-2.0*EYSZ1[index(i,j,k,NX1,NY1,4)]+EYSZ1[index(i,j-1,k,NX1,NY1,4)]+EYSZ1[index(i,j+1,1,NX1,NY1,4)]-2.0*EYSZ1[index(i,j,1,NX1,NY1,4)]+EYSZ1[index(i,j-1,1,NX1,NY1,4)]);
  }

  if(k==(NZ-1) && j>0 && j<NY1-1 && i>1 && i<NX1-1){
    EYS[index(i,j,k,NX,NY,NZ)] = -EYSZ2[index(i,j,2,NX1,NY1,4)]+czd*(EYS[index(i,j,NZ1-1,NX,NY,NZ)]+EYSZ2[index(i,j,3,NX1,NY1,4)])+czz*(EYSZ1[index(i,j,3,NX1,NY1,4)]+EYSZ1[index(i,j,2,NX1,NY1,4)])+czfxd*(EYSZ1[index(i+1,j,3,NX1,NY1,4)]-2.0*EYSZ1[index(i,j,3,NX1,NY1,4)]+EYSZ1[index(i-1,j,3,NX1,NY1,4)]+EYSZ1[index(i+1,j,2,NX1,NY1,4)]-2.0*EYSZ1[index(i,j,2,NX1,NY1,4)]+EYSZ1[index(i-1,j,2,NX1,NY1,4)])+czfyd*(EYSZ1[index(i,j+1,3,NX1,NY1,4)]-2.0*EYSZ1[index(i,j,3,NX1,NY1,4)]+EYSZ1[index(i,j-1,3,NX1,NY1,4)]+EYSZ1[index(i,j+1,2,NX1,NY1,4)]-2.0*EYSZ1[index(i,j,2,NX1,NY1,4)]+EYSZ1[index(i,j-1,2,NX1,NY1,4)]);
  }
}

__global__ void radeyz_save(double *EYS, double *EYSZ1, double *EYSZ2, int NX, int NY, int NZ, int NX1, int NY1, int NZ1){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;
  int i,j,k;
  k=(tid%(NY1*4))%4;
  j=(tid%(NY1*4))/4;
  i=tid/(NY1*4);

  if(j<NY1 && i>0 && i<NX1){
    EYSZ2[index(i,j,k,NX1,NY1,4)]=EYSZ1[index(i,j,k,NX1,NY1,4)];
  }

  if((k==0 || k==1) && j<NY1 && i>0 && i<NX1){
    EYSZ1[index(i,j,k,NX1,NY1,4)]=EYS[index(i,j,k,NX,NY,NZ)];
  }

  if(k==2 && j<NY1 && i>0 && i<NX1){
    EYSZ1[index(i,j,k,NX1,NY1,4)]=EYS[index(i,j,NZ1-1,NX,NY,NZ)];
  }

  if(k==3 && j<NY1 && i>0 && i<NX1){
    EYSZ1[index(i,j,k,NX1,NY1,4)]=EYS[index(i,j,NZ-1,NX,NY,NZ)];
  }
}

__global__ void hupdate(double *EXS, double *EYS, double *EZS, double *HXS, double *HYS, double *HZS, int NZ1, int NY1, int NX1){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;
  int i,j,k;
  k=(tid%(NY*NZ))%NZ;
  j=(tid%(NY*NZ))/NZ;
  i=tid/(NY*NZ);

  if(i>0 && i<NX1 && j<NY1 && k<NZ1){
    HXS[tid]=HXS[tid]-(EZS[index(i,j+1,k,NX,NY,NZ)]-EZS[tid])*dtmdy+(EYS[index(i,j,k+1,NX,NY,NZ)]-EYS[tid])*dtmdz;
  }

  if(i<NX1 && j>0 && j<NY1 && k<NZ1){
    HYS[tid]=HYS[tid]-(EXS[index(i,j,k+1,NX,NY,NZ)]-EXS[tid])*dtmdz+(EZS[index(i+1,j,k,NX,NY,NZ)]-EZS[tid])*dtmdx;
  }

  if(i<NX1 && j<NY1 && k>0 && k<NZ1){
    HZS[tid]=HZS[tid]-(EYS[index(i+1,j,k,NX,NY,NZ)]-EYS[tid])*dtmdx+(EXS[index(i,j+1,k,NX,NY,NZ)]-EXS[tid])*dtmdy;
  }
}

__global__ void datsav1(double *etimeavg, double *emax, double *EXS, double *EYS, double *EZS, int *count, int N, double t, double *dt, int NX, int NY, int NZ){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;
  double nrms,exc=0.0,eyc=0.0,ezc=0.0,exin,eyin,ezin,esq=0.0;
  nrms=(1.0/(freq*(*dt)));
  int i,j,k;
  k=(tid%(NY*NZ))%NZ;
  j=(tid%(NY*NZ))/NZ;
  i=tid/(NY*NZ);

  if(i<NX-2 && j<NY-2 && k<NZ-2){
    
    t=t-(*dt);
    exin=EXI(i,j,k,t)+EXI(i,j,k+1,t)+EXI(i,j+1,k,t)+EXI(i,j+1,k+1,t);
    eyin=EYI(i,j,k,t)+EYI(i+1,j,k,t)+EYI(i,j,k+1,t)+EYI(i+1,j,k+1,t);
    ezin=EZI(i,j,k,t)+EZI(i+1,j,k,t)+EZI(i+1,j+1,k,t)+EZI(i,j+1,k,t);
    t=t+(*dt);
    
    exc=(exin+EXS[index(i,j,k,NX,NY,NZ)]+EXS[index(i,j,k+1,NX,NY,NZ)]+EXS[index(i,j+1,k,NX,NY,NZ)]+EXS[index(i,j+1,k+1,NX,NY,NZ)])*(0.25);
    eyc=(eyin+EYS[index(i,j,k,NX,NY,NZ)]+EYS[index(i+1,j,k,NX,NY,NZ)]+EYS[index(i,j,k+1,NX,NY,NZ)]+EYS[index(i+1,j,k+1,NX,NY,NZ)])*(0.25);
    ezc=(ezin+EZS[index(i,j,k,NX,NY,NZ)]+EZS[index(i+1,j,k,NX,NY,NZ)]+EZS[index(i+1,j+1,k,NX,NY,NZ)]+EZS[index(i,j+1,k,NX,NY,NZ)])*(0.25);
    
    esq=((exc*exc)+(eyc*eyc)+(ezc*ezc));
    (etimeavg[index(i,j,k,NX,NY,NZ)])=(etimeavg[index(i,j,k,NX,NY,NZ)])+(esq/nrms);
    
    if(N==(nstop-nrms)) emax[index(i,j,k,NX,NY,NZ)]=esq;
    if(esq>(emax[index(i,j,k,NX,NY,NZ)])) emax[index(i,j,k,NX,NY,NZ)]=esq;
  }
  return;
}

__global__ void sar_cal(double *etimeavg, double *emax, int *IDONE, double *count2, double *sarx, int NX, int NY, int NZ){

  int tid=blockIdx.x*blockDim.x+threadIdx.x;
  int i,j,k;
  k=(tid%(NY*NZ))%NZ;
  j=(tid%(NY*NZ))/NZ;
  i=tid/(NY*NZ);
  double temp;
  double r2=0.5, r3=0.1;

  if(wbsar==1 && i<NX-2 && j<NY-2 && k<NZ-2){
    temp=((pow((i-nxc)*delx,2)/pow(r3,2)) + (pow((j-nyc)*dely,2)/pow(r3,2)) + (pow((k-nzc)*delz,2)/pow(r2,2)));

    // temp=((((i-nxc)*delx)*((i-nxc)*delx))+(((j-nyc)*dely)*((j-nyc)*dely))+(((k-nzc)*delz)*((k-nzc)*delz)));
    // r1=sqrt(temp);
    
    if(/*r1<=radius2*/ temp>0.0 && temp<=1.0){
      count2[tid]=1;
      sarx[tid]=((0.5)*SIGMA[IDONE[index(i,j,k,NX,NY,NZ)]-1]*(emax[index(i,j,k,NX,NY,NZ)])*(delx*dely*delz));
      // printf("%.25lf\n",sarx[tid]);
    }
  }
}


int main(){
  
  // double *EXS, *EYS, *EZS;
  // double *HXS, *HYS, *HZS;
  // double *EXSY1, *EXSY2, *EXSZ1, *EXSZ2, *EYSX1, *EYSX2, *EYSZ1, *EYSZ2, *EZSX1, *EZSX2, *EZSY1, *EZSY2;
  // int *IDONE, *IDTWO, *IDTHREE;
  // double *EZS;
  // int *IDONE;

  double *d_EXS, *d_EYS, *d_EZS;
  double *d_HXS, *d_HYS, *d_HZS;
  double *d_EXSY1, *d_EXSY2, *d_EXSZ1, *d_EXSZ2, *d_EYSX1, *d_EYSX2, *d_EYSZ1, *d_EYSZ2, *d_EZSX1, *d_EZSX2, *d_EZSY1, *d_EZSY2;
  int *d_IDONE, *d_IDTWO, *d_IDTHREE;

  long int domain_size=NX*NY*NZ;
  //double *emax;
  int *d_count;
  double *d_etimeavg, *d_emax;// , *d_s, *d_erms1d, *d_erms, *d_sar,

  // EXS=(double*)malloc(domain_size*sizeof(double));
  // EYS=(double*)malloc(domain_size*sizeof(double));
  // EZS=(double*)malloc(domain_size*sizeof(double));
  // HXS=(double*)malloc(domain_size*sizeof(double));
  // HYS=(double*)malloc(domain_size*sizeof(double));
  // HZS=(double*)malloc(domain_size*sizeof(double));
  // EXSY1=(double*)malloc(NX1*4*NZ1*sizeof(double));
  // EXSY2=(double*)malloc(NX1*4*NZ1*sizeof(double));
  // EXSZ1=(double*)malloc(NX1*NY1*4*sizeof(double));
  // EXSZ2=(double*)malloc(NX1*NY1*4*sizeof(double));
  // EYSX1=(double*)malloc(4*NY1*NZ1*sizeof(double));
  // EYSX2=(double*)malloc(4*NY1*NZ1*sizeof(double));
  // EYSZ1=(double*)malloc(NX1*NY1*4*sizeof(double));
  // EYSZ2=(double*)malloc(NX1*NY1*4*sizeof(double));
  // EZSX1=(double*)malloc(4*NY1*NZ1*sizeof(double));
  // EZSX2=(double*)malloc(4*NY1*NZ1*sizeof(double));
  // EZSY1=(double*)malloc(NX1*4*NZ1*sizeof(double));
  // EZSY2=(double*)malloc(NX1*4*NZ1*sizeof(double));
  // IDONE=(int*)malloc(domain_size*sizeof(int));
  // IDTWO=(int*)malloc(domain_size*sizeof(int));
  // IDTHREE=(int*)malloc(domain_size*sizeof(int));
  // emax=(double*)malloc(domain_size*sizeof(double));
  
  cudaMalloc((void **)&d_EXS, domain_size*sizeof(double));
  cudaMalloc((void **)&d_EYS, domain_size*sizeof(double));
  cudaMalloc((void **)&d_EZS, domain_size*sizeof(double));
  cudaMalloc((void **)&d_HXS, domain_size*sizeof(double));
  cudaMalloc((void **)&d_HYS, domain_size*sizeof(double));
  cudaMalloc((void **)&d_HZS, domain_size*sizeof(double));
  cudaMalloc((void **)&d_EXSY1, NX1*4*NZ1*sizeof(double));
  cudaMalloc((void **)&d_EXSY2, NX1*4*NZ1*sizeof(double));
  cudaMalloc((void **)&d_EXSZ1, NX1*NY1*4*sizeof(double));
  cudaMalloc((void **)&d_EXSZ2, NX1*NY1*4*sizeof(double));
  cudaMalloc((void **)&d_EYSX1, 4*NY1*NZ1*sizeof(double));
  cudaMalloc((void **)&d_EYSX2, 4*NY1*NZ1*sizeof(double));
  cudaMalloc((void **)&d_EYSZ1, NX1*NY1*4*sizeof(double));
  cudaMalloc((void **)&d_EYSZ2, NX1*NY1*4*sizeof(double));
  cudaMalloc((void **)&d_EZSX1, 4*NY1*NZ1*sizeof(double));
  cudaMalloc((void **)&d_EZSX2, 4*NY1*NZ1*sizeof(double));
  cudaMalloc((void **)&d_EZSY1, NX1*4*NZ1*sizeof(double));
  cudaMalloc((void **)&d_EZSY2, NX1*4*NZ1*sizeof(double));
  cudaMalloc((void **)&d_IDONE, domain_size*sizeof(int));
  cudaMalloc((void **)&d_IDTWO, domain_size*sizeof(int));
  cudaMalloc((void **)&d_IDTHREE, domain_size*sizeof(int));

  cudaMalloc((void **)&d_etimeavg,domain_size*sizeof(double));
  cudaMalloc((void **)&d_emax,domain_size*sizeof(double));

  double *d_dt;
  cudaMalloc((void **)&d_dt, sizeof(double));

  long int blocks_per_grid;
  int threads_per_block=512;
  blocks_per_grid=(domain_size+threads_per_block-1)/threads_per_block;

  cudaEvent_t start,stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start);

  zero_fields<<<blocks_per_grid, threads_per_block>>>(d_EXS, d_EYS, d_EZS, d_HXS, d_HYS, d_HZS, d_IDONE, d_IDTWO, d_IDTHREE, domain_size);

  zero_xplanes<<<(4*NY1*NZ1+threads_per_block-1)/threads_per_block,threads_per_block>>>(d_EYSX1, d_EYSX2, d_EZSX1, d_EZSX2, 4*NY1*NZ1);
  
  zero_yplanes<<<(NX1*4*NZ1+threads_per_block-1)/threads_per_block,threads_per_block>>>(d_EXSY1, d_EXSY2, d_EZSY1, d_EZSY2, NX1*4*NZ1);

  zero_zplanes<<<(NX1*NY1*4+threads_per_block-1)/threads_per_block,threads_per_block>>>(d_EXSZ1, d_EXSZ2, d_EYSZ1, d_EYSZ2, NX1*NY1*4);

  zeromdim<<<1,1>>>();

  build<<<(domain_size+threads_per_block-1)/threads_per_block, threads_per_block>>>(NX, NY, NZ, d_IDONE, d_IDTWO, d_IDTHREE, domain_size);

  cudaDeviceSynchronize();

  setup<<<1,1>>>(d_dt);

  cudaDeviceSynchronize();
  
  /*FILE *fpe;
  fpe=fopen("ezs.txt","w");*/
  cudaMemcpy(&dt, d_dt, sizeof(double), cudaMemcpyDeviceToHost);
  //printf("dt : %.15lf\n",dt);
  /*FILE *fsar;
  fsar=fopen("emax.txt", "w");*/
  double nrms=(1.0/(frequency*dt));
  
  for(N=1;N<=nstop;N++){
    //printf("%d\n", N);
    //printf("%d\n",N);
    //k<<<blocks_per_grid, threads_per_block>>>(d_EXS, d_EYS, d_EZS, t);
    //cudaDeviceSynchronize();
    //printf("kernel done\n");
    eupdate<<<blocks_per_grid, threads_per_block>>>(d_EXS, d_EYS, d_EZS, d_HXS, d_HYS, d_HZS, d_IDONE, d_IDTWO, d_IDTHREE, NZ1, NY1, NX1, t);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
    //cudaDeviceSynchronize();

    cudaDeviceSynchronize();

  // check for error
  /*cudaError_t error = cudaGetLastError();
  if(error != cudaSuccess)
  {
    // print the CUDA error message and exit
    printf("CUDA error: %s\n", cudaGetErrorString(error));
    exit(-1);
  }
    exit(0);*/
 
    radeyx1<<<blocks_per_grid,threads_per_block>>>(d_EYS, d_EYSX1, d_EYSX2, NX, NY, NZ, NX1, NY1, NZ1);
    //cudaDeviceSynchronize();
    radeyx2<<<blocks_per_grid,threads_per_block>>>(d_EYS, d_EYSX1, d_EYSX2, NX, NY, NZ, NX1, NY1, NZ1);
    cudaDeviceSynchronize();
    radeyx_save<<<((4*NY1*NZ1)+threads_per_block-1)/threads_per_block,threads_per_block>>>(d_EYS, d_EYSX1, d_EYSX2, NX, NY, NZ, NX1, NY1, NZ1);

    cudaDeviceSynchronize();
    
    radezx1<<<blocks_per_grid,threads_per_block>>>(d_EZS, d_EZSX1, d_EZSX2, NX, NY, NZ, NX1, NY1, NZ1);
    //cudaDeviceSynchronize();
    radezx2<<<blocks_per_grid,threads_per_block>>>(d_EZS, d_EZSX1, d_EZSX2, NX, NY, NZ, NX1, NY1, NZ1);
    cudaDeviceSynchronize();
    radezx_save<<<((4*NY1*NZ1)+threads_per_block-1)/threads_per_block,threads_per_block>>>(d_EZS, d_EZSX1, d_EZSX2, NX, NY, NZ, NX1, NY1, NZ1);

    cudaDeviceSynchronize();
    
    radezy1<<<blocks_per_grid,threads_per_block>>>(d_EZS, d_EZSY1, d_EZSY2, NX, NY, NZ, NX1, NY1, NZ1);
    radezy2<<<blocks_per_grid,threads_per_block>>>(d_EZS, d_EZSY1, d_EZSY2, NX, NY, NZ, NX1, NY1, NZ1);
    cudaDeviceSynchronize();
    radezy_save<<<((NX1*4*NZ1)+threads_per_block-1)/threads_per_block,threads_per_block>>>(d_EZS, d_EZSY1, d_EZSY2, NX, NY, NZ, NX1, NY1, NZ1);

    cudaDeviceSynchronize();
    
    radexy1<<<blocks_per_grid,threads_per_block>>>(d_EXS, d_EXSY1, d_EXSY2, NX, NY, NZ, NX1, NY1, NZ1);
    radexy2<<<blocks_per_grid,threads_per_block>>>(d_EXS, d_EXSY1, d_EXSY2, NX, NY, NZ, NX1, NY1, NZ1);
    cudaDeviceSynchronize();
    radexy_save<<<((NX1*4*NZ1)+threads_per_block-1)/threads_per_block,threads_per_block>>>(d_EXS, d_EXSY1, d_EXSY2, NX, NY, NZ, NX1, NY1, NZ1);

    cudaDeviceSynchronize();
    
    radexz1<<<blocks_per_grid,threads_per_block>>>(d_EXS, d_EXSZ1, d_EXSZ2, NX, NY, NZ, NX1, NY1, NZ1);
    radexz2<<<blocks_per_grid,threads_per_block>>>(d_EXS, d_EXSZ1, d_EXSZ2, NX, NY, NZ, NX1, NY1, NZ1);
    cudaDeviceSynchronize();
    radexz_save<<<((NX1*NY1*4)+threads_per_block-1)/threads_per_block,threads_per_block>>>(d_EXS, d_EXSZ1, d_EXSZ2, NX, NY, NZ, NX1, NY1, NZ1);
	
    cudaDeviceSynchronize();
    
    radeyz1<<<blocks_per_grid,threads_per_block>>>(d_EYS, d_EYSZ1, d_EYSZ2, NX, NY, NZ, NX1, NY1, NZ1);
    radeyz2<<<blocks_per_grid,threads_per_block>>>(d_EYS, d_EYSZ1, d_EYSZ2, NX, NY, NZ, NX1, NY1, NZ1);
    cudaDeviceSynchronize();
    radeyz_save<<<((NX1*NY1*4)+threads_per_block-1)/threads_per_block,threads_per_block>>>(d_EYS, d_EYSZ1, d_EYSZ2, NX, NY, NZ, NX1, NY1, NZ1);

    cudaDeviceSynchronize();
    
    t=t+((dt)/2.0);
    /*printf("t: %.15lf\n",t);
    cudaMemcpy(EZS,d_EZS, domain_size*sizeof(double), cudaMemcpyDeviceToHost);
    fprintf(fpe, "%.15lf %d %.25lf\n",t,N,EZS[index(50,50,50,NX,NY,NZ)]);
    */
    hupdate<<<blocks_per_grid, threads_per_block>>>(d_EXS, d_EYS, d_EZS, d_HXS, d_HYS, d_HZS, NZ1, NY1, NX1);

    cudaDeviceSynchronize();

    t=t+((dt)/2.0);

    if(N>=(nstop-nrms) && N<=(nstop-1)){
      datsav1<<<blocks_per_grid, threads_per_block>>>(d_etimeavg, d_emax, d_EXS, d_EYS, d_EZS, d_count, N, t, d_dt, NX, NY, NZ);
      cudaDeviceSynchronize();
      /*cudaMemcpy(emax, d_emax, domain_size*sizeof(double),cudaMemcpyDeviceToHost);
      fprintf(fsar, "%d %.15lf\n",N,emax[index(50,50,50,NX,NY,NZ)]);*/
    }
    cudaDeviceSynchronize();

  }

  //fclose(fpe);
  /* 
  cudaMemcpy(IDONE, d_IDONE, domain_size*sizeof(int), cudaMemcpyDeviceToHost);
  printf("copying scatterer data\n");
  FILE *fp;
  fp=fopen("log.txt","w");
  int i,j,k;
  for(long int it=0;it<domain_size;it++){
    k=(it%(NY*NZ))%NZ;
    j=(it%(NY*NZ))/NZ;
    i=it/(NY*NZ);
    if(IDONE[it]==2)
    fprintf(fp,"%d %d %d %d\n", i,j,k,IDONE[it]);
  }
  fclose(fp);
  */
  //free(IDONE);
  //free(EZS);
  
  cudaFree(d_EXS);
  cudaFree(d_EYS);
  cudaFree(d_EZS);
  cudaFree(d_HXS);
  cudaFree(d_HYS);
  cudaFree(d_HZS);
  cudaFree(d_EXSY1);
  cudaFree(d_EXSY2);
  cudaFree(d_EXSZ1);
  cudaFree(d_EXSZ2);
  cudaFree(d_EYSX1);
  cudaFree(d_EYSX2);
  cudaFree(d_EYSZ1);
  cudaFree(d_EYSZ2);
  cudaFree(d_EZSX1);
  cudaFree(d_EZSX2);
  cudaFree(d_EZSY1);
  cudaFree(d_EXSY2);
  // cudaFree(d_IDONE);
  cudaFree(d_IDTWO);
  cudaFree(d_IDTHREE);

  double *d_count2;
  double *d_sarx;

  cudaMalloc((void **)&d_count2, domain_size*sizeof(double));
  cudaMemset(d_count2,0,domain_size*sizeof(double));
  cudaMalloc((void **)&d_sarx, domain_size*sizeof(double));
  cudaMemset(d_sarx,0,domain_size*sizeof(double));
  
  sar_cal<<<blocks_per_grid, threads_per_block,threads_per_block*sizeof(double)>>>(d_etimeavg, d_emax, d_IDONE, d_count2, d_sarx, NX, NY, NZ);
  cudaDeviceSynchronize();

  thrust::device_ptr<double> cptr1=thrust::device_pointer_cast(d_count2);
  double count2=thrust::reduce(cptr1,cptr1+domain_size);
  thrust::device_ptr<double> cptr2=thrust::device_pointer_cast(d_sarx);
  double sarx=thrust::reduce(cptr2,cptr2+domain_size);

  printf("Total cell count: %lf\ncells used for SAR calculation: %lf\nTotal power: %.15lf\nWhole Body SAR: %.15lf\n",count2, count2, sarx, (sarx)/(count2*h_delx*h_dely*h_delz*h_iden));

  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  float runtime=0;
  cudaEventElapsedTime(&runtime, start, stop);
  printf("Run time: %lf\n",runtime);

  cudaFree(d_count2);
  cudaFree(d_sarx);  
  
  return 0;
}
