#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <time.h>
#include <malloc.h>
#include <math.h>
#include "var.h"

#define pi 4.0*atan(1.0)

/* Function Declarations */
void zero();
void build();
void dcube();
void setup();
void exsfld();
void eysfld();
void ezsfld();
void radeyx();
void radezx();
void radezy();
void radexy();
void radexz();
void radeyz();
void hxsfld();
void hysfld();
void hzsfld();
double EXI(int, int, int);
double EYI(int, int, int);
double EZI(int, int, int);
double source(double);
double DEXI(int, int, int);
double DEYI(int, int, int);
double DEZI(int, int, int);
double dsrce(double);

void datsav();
double *emax, *etimeavg;

long long offset(int i, int j, int k){
  return (i*NY*NZ+j*NZ+k);
}

//FILE *fp1;
//FILE *radfp;
/* Main */
int main(){
	float runtime=0;
	clock_t start=0;
	clock_t stop=0;
	start=clock();
	zero();
	build();
	setup();
	
	int i,j,k;
	/*fp2=fopen("build.txt","w");
	for(i=0;i<NX;i++){
	  for(j=0;j<NY;j++){
	    for(k=0;k<NZ;k++){
	      if(IDONE[i][j][k]==2)
	        fprintf(fp2,"%d %d %d %d\n",i,j,k,IDONE[i][j][k]);
	    }
	  }
	}*/
	
	FILE *fp;
	//fsar=fopen("eavg.txt", "w");
	fp=fopen("time.txt", "ab+");
	/*radfp=fopen("zsx2.txt","w");*/
	//fp1=fopen("etimeavg.txt", "w");
	emax=(double*)malloc(NX*NY*NZ*sizeof(double));
	etimeavg=(double*)malloc(NX*NY*NZ*sizeof(double));

	float e_runtime=0;
        clock_t e_start=0;
        clock_t e_stop=0;

        float h_runtime=0;
        clock_t h_start=0;
        clock_t h_stop=0;

        float rad_runtime=0;
        clock_t rad_start=0;
        clock_t rad_stop=0;

        float sar_runtime=0;
        clock_t sar_start=0;
        clock_t sar_stop=0;

	for(N=1;N<=nstop;N++){
	  printf("%d\n", N);
		
	  //advance scattered electric field
	  e_start=clock();
	  exsfld();
	  eysfld();
	  ezsfld();
	  e_stop=clock(); 
	  e_runtime+=((float)e_stop-e_start)/CLOCKS_PER_SEC;  
	//apply radiation boundary condition (second order)
	 rad_start=clock(); 
	  radeyx();
	  radezx();
	  radezy();
	  radexy();
	  radexz();
	  radeyz();
	 rad_stop=clock(); 
	  rad_runtime+=((float)rad_stop-rad_start)/CLOCKS_PER_SEC;
	//advance time by 1/2 time step
	  t=t+(dt/2.0);
	  //fprintf(fp, " %.15lf %d %.25lf\n",t, N,EZS[50][50][50]);
	  
	//advance scattered magentic field
	  
	h_start=clock();
	hxsfld();
	  hysfld();
	  hzsfld();
	  h_stop=clock();
	  h_runtime+=((float)h_stop-h_start)/CLOCKS_PER_SEC;
	
	t=t+(dt/2.0);
	sar_start=clock();  
	datsav();
	sar_stop=clock();	
	sar_runtime+=((float)sar_stop-sar_start)/CLOCKS_PER_SEC;
	}
	free(emax);
	free(etimeavg);

	stop=clock();
	runtime=((float)stop-start)/CLOCKS_PER_SEC;
//	fprintf(fp,,
	fprintf(fp,"size: %d time: %lf %lf %lf %lf %lf\n",NX,runtime,e_runtime,h_runtime,rad_runtime,sar_runtime);
	//fclose(fp1);
	fclose(fp);
	//fclose(radfp);
	return 0;
}

/* Function Definitions */
void zero(){
  
	int i, j, k;

	for(i = 0; i <NX; i++)
	{
		for(j = 0; j<NY; j++)
		{
			for(k = 0; k<NZ; k++)
			{
				EXS[i][j][k]=0.0;
				EYS[i][j][k]=0.0;
				EZS[i][j][k]=0.0;
				HXS[i][j][k]=0.0;
				HYS[i][j][k]=0.0;
				HZS[i][j][k]=0.0;
				IDONE[i][j][k]=0;
				IDTWO[i][j][k]=0;
				IDTHREE[i][j][k]=0;
			}
		}
	}

	for(i = 0; i <4; i++)
	{
		for(j = 0; j<NY1; j++)
		{
			for(k = 0; k<NZ1; k++)
			{
				EYSX1[i][j][k]=0.0;
				EYSX2[i][j][k]=0.0;
				EZSX1[i][j][k]=0.0;
				EZSX2[i][j][k]=0.0;
			}
		}
	}

	for(i = 0; i <NX1; i++)
	{
		for(j = 0; j<4; j++)
		{
			for(k = 0; k<NZ1; k++)
			{
				EXSY1[i][j][k]=0.0;
				EXSY2[i][j][k]=0.0;
				EZSY1[i][j][k]=0.0;
				EZSY2[i][j][k]=0.0;
			}
		}
	}

	for(i = 0; i <NX1; i++)
	{
		for(j = 0; j<NY1; j++)
		{
			for(k = 0; k<4; k++)
			{
				EYSZ1[i][j][k]=0.0;
				EYSZ2[i][j][k]=0.0;
				EXSZ1[i][j][k]=0.0;
				EXSZ2[i][j][k]=0.0;
			}
		}
	}

	for(i = 0; i<mdim1; i++)
	{
		ESCTC[i]=0.0;
        	EINCC[i]=0.0;
	        EDEVCN[i]=0.0;
	        ECRLX[i]=0.0;
        	ECRLY[i]=0.0;
	        ECRLZ[i]=0.0;
	}
}

void build(){
  int mtype = 0;
  int i,j,k;
  double temp;
  double r1,r2,r3;

	r2=50.0e-2,r3=10.0e-2;	 
 if(NX%2==0 && NY%2==0 && NZ%2==0){
    printf("Enter odd value for NX, NY, NZ");
  }

  else{
    nxc=((NX+1)/2)-1;
    nyc=((NY+1)/2)-1;
    nzc=((NZ+1)/2)-1;
    
    // spherical geometry

    //fp2=fopen("build.txt","w");
    mtype=2;
    for(i=0;i<NX;i++){
      for(j=0;j<NY;j++){
        for(k=0;k<NZ;k++){
          temp=((pow((i-nxc)*delx,2)/pow(r3,2)) + (pow((j-nyc)*dely,2)/pow(r3,2)) + (pow((k-nzc)*delz,2)/pow(r2,2)));
          //r1=sqrt(temp);

          if(radius2<delx || radius2<dely ||radius2<delz){
            printf("radius of sphere %f (%d, %d, %d) is less than cell size.\n",r1,i,j,k);
            exit(-1);
          }
          if(temp>0.0 && temp<=1.0){
	    dcube(i,j,k,1,1,1,mtype);
            //fprintf(fp2, "%d, %d, %d, %d\n", i, j, k, IDONE[i][j][k]);
          }
        }
      }
    }

    //GEOMETRY OF METAL PLATE
    /*mtype=1;
    for(i=0;i<NX;i++){
      for(j=0;j<NY;j++){
        for(k=0;k<NZ;k++){


      if((i==45) && (j>((nyc)-15) && j<((nyc)+15)) && (k>((nzc)-15) && k<((nzc)+15))){
          dcube(i,j,k,1,1,1,mtype);
	fprintf(fp2, "%f, %f, %f, %d\n", (i)*delx, (j)*dely, (k)*delz, IDONE[i][j][k]);
	     }   
	  }
       }
  }*/

  }

    for(k=0;k<NZ;k++){
      for(j=0;j<NY;j++){
        for(i=0;i<NX;i++){
          if(IDONE[i][j][k]>=10 || IDTWO[i][j][k]>=10 || IDTHREE[i][j][k]>=10){
            printf("error occured. illegal value for IDONE-IDTHREE at %d, %d, %d\n",i,j,k);
            exit(-1);
          }
        }
      }
    }
  
}

void dcube(int istart, int jstart, int kstart, int nxwide, int nywide, int nzwide, int mtype)
{
  	int imax, jmax, kmax;
	int i, j, k;

	imax = istart+nxwide-1;
	jmax = jstart+nywide-1;
	kmax = kstart+nzwide-1;

	if(nxwide == 0)
	{
		for(k = kstart; k<=kmax; k++)
		{
			for(j = jstart; j<=jmax; j++)
			{
				IDTWO[istart][j][k] = mtype;
				IDTWO[istart][j][k+1] = mtype;
				IDTHREE[istart][j][k] = mtype;
				IDTHREE[istart][j+1][k] = mtype;
			}
		}
	}

	else if(nywide == 0)
	{
		for(k = kstart; k<=kmax; k++)
		{
			for(i = istart; i<=imax; i++)
			{
				IDONE[i][jstart][k] = mtype;
				IDONE[i][jstart][k+1] = mtype;
				IDTHREE[i][jstart][k] = mtype;
				IDTHREE[i+1][jstart][k] = mtype;
			}
		}
	}

	else if(nzwide == 0)
	{
		for(j = jstart; j<=jmax; j++)
		{
			for(i = istart; i<=imax; i++)
			{
				IDONE[i][j][kstart] = mtype;
				IDONE[i][j+1][kstart] = mtype;
				IDTWO[i][j][kstart] = mtype;
				IDTWO[i+1][j][kstart] = mtype;
			}
		}
	}

	else
	{
		for(k = kstart; k<=kmax; k++)
		{
			for(j = jstart; j<=jmax; j++)
			{
				for(i = istart; i<=imax; i++)
				{
					IDONE[i][j][k] = mtype;
					IDONE[i][j][k+1] = mtype;
					IDONE[i][j+1][k+1] = mtype;
					IDONE[i][j+1][k] = mtype;
					IDTWO[i][j][k] = mtype;
					IDTWO[i+1][j][k] = mtype;
					IDTWO[i+1][j][k+1] = mtype;
					IDTWO[i][j][k+1] = mtype;
					IDTHREE[i][j][k] = mtype;
					IDTHREE[i+1][j][k] = mtype;
					IDTHREE[i+1][j+1][k] = mtype;
					IDTHREE[i][j+1][k] = mtype;
				}
			}
		}
	}

}

void setup()
{

	// THIS SUBROUTINE INITIALIZES THE COMPUTATIONS
	double dtxi, dtyi, dtzi;
	int i;

	dtxi = c/delx;
	dtyi = c/dely;
	dtzi = c/delz;

	// CALCULATE DT--THE MAXIMUM TIME STEP ALLOWED BY THE COURANT STABILITY CONDITION
	//printf("inside: %f\n",sqrt(pow(dtxi, 2.0)+pow(dtyi, 2.0)+pow(dtzi, 2.0)));
	dt = 1.0/sqrt(pow(dtxi, 2.0)+pow(dtyi, 2.0)+pow(dtzi, 2.0));
	printf("dt----- %0.15lf\n", dt);
	printf("%lf %lf %lf\n", delx, dely, delz);	
	/*
	PARAMETER ALPHA IS THE DECAY RATE DETERMINED BY BETA.
C  TO CHANGE THE GAUSSIAN PULSE BY SINE WAVE WE HAVE TO MODIFY THE
c  CODE .WHERE EVER WE ARE USING THE ALPHA AND BETA WE HAVE
c  TO REMOVE IT .AND ALSO WE HAVE TO CHANGE THE TIME DURATION
c  ALSO TO EXIT THE SINE WAVE FOR THAT INTERVAL OF TIME.
	*/

	// in 3d fortran parameters.h
	period = 1e-6;

	// SET OFFSET FOR COMPUTING INCIDENT FIELDS

	off = 10;

	// THE FOLLOWING LINES ARE FOR SMOOTH COSINE INCIDENT FUNCTION
	// FIND DIRECTION COSINES FOR INCIDENT FIELD
	//printf("arg: %lf\n", pi*thinc/180.0);
	double costh = cos(pi*thinc/180.0);
	double sinth = sin(pi*thinc/180.0);
	double cosph = cos(pi*phinc/180.0);
	double sinph = sin(pi*phinc/180.0);
	printf("costh, sinth, cosph, sinph\n%0.25lf, %.25lf, %.25lf, %.25lf\n", costh, sinth, cosph, sinph);
	
	// FIND AMPLITUDE OF INCIDENT FIELD COMPONENTS
	
	ampx = amp*(ethinc*cosph*costh - ephinc*sinph);
	ampy = amp*(ethinc*sinph*costh + ephinc*cosph);
	ampz = amp*(-ethinc*sinth);
	printf("ampx, ampy, ampz %.15lf, %.15lf, %.15lf", ampx, ampy, ampz);
	
	// FIND RELATIVE SPATIAL DELAY FOR X, Y, Z CELL DISPLACEMENT

	xdisp = -cosph*sinth;
	ydisp = -sinth*sinph;
	zdisp = -costh;
	printf("xdisp: %.15lf, ydisp: %.15lf, zdisp: %.25lf\n", xdisp, ydisp, zdisp);
	
	for(i=1;i<mdim1;i++){
		EPS[i] = 0.0;
		SIGMA[i]=0.0;
	}

	for(i=1;i<mdim1;i++){
		Ep[i]=57.0;
	}

	for(i = 1; i < mdim1; i++){
		EPS[i] = Ep[i]*eps0;
		SIGMA[i] =0.89;
		printf("%.15lf\n",EPS[i]);
	}

	// FREE SPACE:

	dtedx = dt/(eps0*delx);
	dtedy = dt/(eps0*dely);
	dtedz = dt/(eps0*delz);
	dtmdx = dt/(xmu0*delx);
	dtmdy = dt/(xmu0*dely);
	dtmdz = dt/(xmu0*delz);
	//printf("dtedx: %0.12lf, dtedy: %0.12lf, dtedz: %0.12lf, dtmdx: %0.12lf, dtmdy: %0.12lf, dtmdz: %0.12lf\n", dtedx, dtedy, dtedz, dtmdx, dtmdy, dtmdz);
	// lossy dielectrics

	for(i = 1; i < mdim1; i++)
	{
		ESCTC[i] = EPS[i]/(EPS[i]+SIGMA[i]*dt);
		EINCC[i] = SIGMA[i] *dt/(EPS[i]+SIGMA[i]*dt);
		EDEVCN[i] = dt*(EPS[i]-eps0)/(EPS[i]+SIGMA[i]*dt);
		ECRLX[i] = dt/((EPS[i]+SIGMA[i]*dt)*delx);
		ECRLY[i] = dt/((EPS[i]+SIGMA[i]*dt)*dely);
		ECRLZ[i] = dt/((EPS[i]+SIGMA[i]*dt)*delz);
		//printf("EPS %.15lf, SIGMA %.15lf, EP %.15lf\n", EPS[i],SIGMA[i], Ep[i]);
		printf("ESCTC, EINCC, EDEVCN, ECRLX, ECRLY, ECRLZ\n%.15lf, %.15lf, %.15lf, %.15lf, %.15lf, %.15lf", ESCTC[i], EINCC[i], EDEVCN[i], ECRLX[i], ECRLY[i], ECRLZ[i]);
	}

	// FIND MAXIMUM SPATIAL DELAY TO MAKE SURE PULSE PROPAGATES INTO SPACE PROPERLY.

	delay = 0.0;

	if(xdisp < 0.0)
		delay = delay-(xdisp*NX1*delx);
	if(ydisp < 0.0)
		delay = delay-(ydisp*NY1*dely);
	if(zdisp < 0.0)
		delay = delay-(zdisp*NZ1*delz);

	printf("delay: %.15lf\n", delay);	

	// COMPUTE OUTER RADIATION BOUNDARY CONDITION (ORBC) CONSTANTS:
	cxd =(c*dt-delx)/(c*dt+delx);
	cyd =(c*dt-dely)/(c*dt+dely);
	czd =(c*dt-delz)/(c*dt+delz);
	printf("\ncxd, cyd, czd\n%.14lf,%.14lf,%.14lf\n", cxd, cyd, czd);

	cxu = cxd;
	cyu = cyd;
	czu = czd;

	// COMPUTE 2ND ORDER ORBC CONSTANTS
	cxx = 2.0*delx/(c*dt+delx);
	cyy = 2.0*dely/(c*dt+dely);
	czz = 2.0*delz/(c*dt+delz);
	printf("\ncxx, cyy, czz\n%.14lf,%.14lf,%.14lf\n", cxx, cyy, czz);

	cxfyd = delx*c*dt*c*dt/(2.0*dely*dely*(c*dt+delx));
	cxfzd = delx*c*dt*c*dt/(2.0*delz*delz*(c*dt+delx));
	cyfzd = dely*c*dt*c*dt/(2.0*delz*delz*(c*dt+dely));
	cyfxd = dely*c*dt*c*dt/(2.0*delx*delx*(c*dt+dely));
	czfxd = delz*c*dt*c*dt/(2.0*delx*delx*(c*dt+delz));
	czfyd = delz*c*dt*c*dt/(2.0*dely*dely*(c*dt+delz));
	printf("\ncxfyd, cxfzd, cyfxd, cyfzd, czfxd, czfyd\n%.14lf,%.14lf,%.14lf,%.14lf,%.14lf,%.14lf\n", cxfyd, cxfzd, cyfxd, cyfzd, czfxd, czfyd);
	
}

void exsfld()
{
  int i, j, k, ii,jj,kk;
  for(k = 1; k<NZ1; k++)
    {
      kk=k;
      for(j = 1; j<NY1; j++)
	{
	  jj=j;
	  for(i = 0; i<NX1; i++)
	    {
	      if(IDONE[i][j][k] == 0)
		{
		  EXS[i][j][k] = EXS[i][j][k]+(HZS[i][j][k]-HZS[i][j-1][k])*dtedy-(HYS[i][j][k]-HYS[i][j][k-1])*dtedz;
		  //if(EXS[i][j][k]>100)exit(0);
		  //fprintf(fp,"%d %d %d %d %.15lf\n",N,5,nyc,nzc,EXS[i][j][k]);
		}
	      else if(IDONE[i][j][k] == 1)
		{
		  ii=i;
		  EXS[i][j][k] = -EXI(ii,jj,kk);
		  //if(EXS[i][j][k]>100)exit(0);
		  //fprintf(fp3,"%d %d %d %d %.15lf\n",N,i,j,k,EXS[i][j][k]);
		}
	      else if(IDONE[i][j][k] == 2 || IDONE[i][j][k] == 3 || IDONE[i][j][k] == 4 || IDONE[i][j][k] == 5)
		{
		  ii=i;
		  //double exii = EXI(i, j, k);

		  EXS[i][j][k]=EXS[i][j][k]*ESCTC[IDONE[i][j][k]-1]-EINCC[IDONE[i][j][k]-1]*EXI(ii,jj,kk)-EDEVCN[IDONE[i][j][k]-1]*DEXI(ii,jj,kk)+(HZS[i][j][k]-HZS[i][j-1][k])*ECRLY[IDONE[i][j][k]-1]-(HYS[i][j][k]-HYS[i][j][k-1])*ECRLZ[IDONE[i][j][k]-1];
		  //if(EXS[i][j][k]>100)exit(0);
		  //fprintf(fp3,"%d %d %d %d %.15lf\n",N,i,j,k,EXS[i][j][k]);
		}
	    }
	}
    }
}

void eysfld()
{
  int i, j, k,ii,jj,kk;
  for(k = 1; k<NZ1; k++)
    {
      kk=k;
      for(j = 0; j<NY1; j++)
	{
	  jj=j;
	  for(i = 1; i<NX1; i++)
	    {
	      if(IDTWO[i][j][k] == 0)
		{
		  EYS[i][j][k] = EYS[i][j][k]+(HXS[i][j][k]-HXS[i][j][k-1])*dtedz-(HZS[i][j][k]-HZS[i-1][j][k])*dtedx;
		}
	      else if(IDTWO[i][j][k] == 1)
		{
		  ii=i;
		  EYS[i][j][k] = -EYI(ii,jj,kk);
		}
	      else if(IDTWO[i][j][k] == 2 || IDTWO[i][j][k] == 3 || IDTWO[i][j][k] == 4 || IDTWO[i][j][k] == 5)
		{
		  ii=i;
		  //double eyii = EYI(i, j, k);

		  EYS[i][j][k]=EYS[i][j][k]*ESCTC[IDTWO[i][j][k]-1]-EINCC[IDTWO[i][j][k]-1]*EYI(ii,jj,kk)-EDEVCN[IDTWO[i][j][k]-1]*DEYI(ii,jj,kk)+(HXS[i][j][k]-HXS[i][j][k-1])*ECRLZ[IDTWO[i][j][k]-1]-(HZS[i][j][k]-HZS[i-1][j][k])*ECRLX[IDTWO[i][j][k]-1];
		}
	    }
	}
    }
}

void ezsfld()
{
  int i, j, k,ii,jj,kk;
  for(k = 0; k<NZ1; k++)
    {
      kk=k;
      for(j = 1; j<NY1; j++)
	{
	  jj=j;
	  for(i = 1; i<NX1; i++)
	    {
	      if(IDTHREE[i][j][k] == 0)
		{
		  EZS[i][j][k] = EZS[i][j][k]+(HYS[i][j][k]-HYS[i-1][j][k])*dtedx-(HXS[i][j][k]-HXS[i][j-1][k])*dtedy;
		  //fprintf(fp,"%d %d %d %d %.25lf %.25lf %.25lf\n",N,5,nyc,nzc,EZS[5][nyc][nzc],(HYS[i][j][k]-HYS[i-1][j][k])*dtedx,-(HXS[i][j][k]-HXS[i][j-1][k])*dtedy);
		}
	      else if(IDTHREE[i][j][k] == 1)
		{
		  ii=i;
		  EZS[i][j][k] = -EZI(ii,jj,kk);
		}
	      else if(IDTHREE[i][j][k] == 2 || IDTHREE[i][j][k] == 3 || IDTHREE[i][j][k] == 4 || IDTHREE[i][j][k] == 5)
		{
		  ii=i;
		  //double ezii = EZI(i, j, k);

		  EZS[i][j][k]=EZS[i][j][k]*ESCTC[IDTHREE[i][j][k]-1]-EINCC[IDTHREE[i][j][k]-1]*EZI(ii,jj,kk)-EDEVCN[IDTHREE[i][j][k]-1]*DEZI(ii,jj,kk)+(HYS[i][j][k]-HYS[i-1][j][k])*ECRLX[IDTHREE[i][j][k]-1]-(HXS[i][j][k]-HXS[i][j-1][k])*ECRLY[IDTHREE[i][j][k]-1];
		}
	    }
	}
    }
}


void radezx()
{
  int k, j;
  // DO EDGES WITH FIRST ORDER ORBC
  for(k = 0; k<NZ1; k++)
    {
      j = 1;
      EZS[0][j][k] = EZSX1[1][j][k]+ cxd*(EZS[1][j][k] - EZSX1[0][j][k]);
      EZS[NX-1][j][k] = EZSX1[2][j][k] + cxu*(EZS[NX1-1][j][k] - EZSX1[3][j][k]);
      j = NY1-1;
      EZS[0][j][k] = EZSX1[1][j][k]+ cxd*(EZS[1][j][k] - EZSX1[0][j][k]);
      EZS[NX-1][j][k] = EZSX1[2][j][k] + cxu*(EZS[NX1-1][j][k] - EZSX1[3][j][k]);
    }
  for(j = 2; j<NY1-1; j++)
    {
      k = 0;
      EZS[0][j][k] = EZSX1[1][j][k]+ cxd*(EZS[1][j][k] - EZSX1[0][j][k]);
      EZS[NX-1][j][k] = EZSX1[2][j][k] + cxu*(EZS[NX1-1][j][k] - EZSX1[3][j][k]);
      k = NZ1-1;
      EZS[0][j][k] = EZSX1[1][j][k]+ cxd*(EZS[1][j][k] - EZSX1[0][j][k]);
      EZS[NX-1][j][k] = EZSX1[2][j][k] + cxu*(EZS[NX1-1][j][k] - EZSX1[3][j][k]);
    }
  // NOW DO 2ND ORDER ORBC ON REMAINING PORTIONS OF FACES
  for(k = 1; k<NZ1-1; k++)
    {
      for(j = 2; j<NY1-1; j++)
	{
	  EZS[0][j][k] = -EZSX2[1][j][k]+cxd*(EZS[1][j][k]+EZSX2[0][j][k])+cxx*(EZSX1[0][j][k]+EZSX1[1][j][k])+cxfyd*(EZSX1[0][j+1][k]-2.0*EZSX1[0][j][k]+EZSX1[0][j-1][k]+EZSX1[1][j+1][k]-2.0*EZSX1[1][j][k]+EZSX1[1][j-1][k])+cxfzd*(EZSX1[0][j][k+1]-2.0*EZSX1[0][j][k]+EZSX1[0][j][k-1]+EZSX1[1][j][k+1]-2.0*EZSX1[1][j][k]+EZSX1[1][j][k-1]);

	  EZS[NX-1][j][k] = -EZSX2[2][j][k]+cxd*(EZS[NX1-1][j][k]+EZSX2[3][j][k])+cxx*(EZSX1[3][j][k]+EZSX1[2][j][k])+cxfyd*(EZSX1[3][j+1][k]-2.0*EZSX1[3][j][k]+EZSX1[3][j-1][k]+EZSX1[2][j+1][k]-2.0*EZSX1[2][j][k]+EZSX1[2][j-1][k])+cxfzd*(EZSX1[3][j][k+1]-2.0*EZSX1[3][j][k]+EZSX1[3][j][k-1]+EZSX1[2][j][k+1]-2.0*EZSX1[2][j][k]+EZSX1[2][j][k-1]);
	}
    }
  
  // NOW SAVE PAST VALUES

  for(k = 0; k<NZ1; k++)
    {
      for(j = 1; j<NY1; j++)
	{
	  EZSX2[0][j][k]=EZSX1[0][j][k];
	  EZSX2[1][j][k]=EZSX1[1][j][k];
	  EZSX2[2][j][k]=EZSX1[2][j][k];
	  EZSX2[3][j][k]=EZSX1[3][j][k];
	  EZSX1[0][j][k]=EZS[0][j][k];
	  EZSX1[1][j][k]=EZS[1][j][k];
	  EZSX1[2][j][k]=EZS[NX1-1][j][k];
	  EZSX1[3][j][k]=EZS[NX-1][j][k];
	}
    }
    //fprintf(radfp,"%d %.25lf\n",N,EZSX2[0][20][20]);
}

void radeyx()
{
  int k, j;
  // DO EDGES WITH FIRST ORDER ORBC
  for(k = 1; k<NZ1; k++)
    {
      j = 0;
      EYS[0][j][k] = EYSX1[1][j][k]+ cxd*(EYS[1][j][k] - EYSX1[0][j][k]);
      EYS[NX-1][j][k] = EYSX1[2][j][k] + cxu*(EYS[NX1-1][j][k] - EYSX1[3][j][k]);

      j = NY1-1;
      EYS[0][j][k] = EYSX1[1][j][k]+ cxd*(EYS[1][j][k] - EYSX1[0][j][k]);
      EYS[NX-1][j][k] = EYSX1[2][j][k] + cxu*(EYS[NX1-1][j][k] - EYSX1[3][j][k]);
    }

  for(j = 1; j<NY1-1; j++)
    {
      k = 1;
      EYS[0][j][k] = EYSX1[1][j][k]+ cxd*(EYS[1][j][k] - EYSX1[0][j][k]);
      EYS[NX-1][j][k] = EYSX1[2][j][k] + cxu*(EYS[NX1-1][j][k] - EYSX1[3][j][k]);

      k = NZ1-1;
      EYS[0][j][k] = EYSX1[1][j][k]+ cxd*(EYS[1][j][k] - EYSX1[0][j][k]);
      EYS[NX-1][j][k] = EYSX1[2][j][k] + cxu*(EYS[NX1-1][j][k] - EYSX1[3][j][k]);
    }

  // NOW DO 2ND ORDER ORBC ON REMAINING PORTIONS OF FACES

  for(k = 2; k<NZ1-1; k++)
    {
      for(j = 1; j<NY1-1; j++)
	{
	  EYS[0][j][k] = -EYSX2[1][j][k]+cxd*(EYS[1][j][k]+EYSX2[0][j][k])+cxx*(EYSX1[0][j][k]+EYSX1[1][j][k])+cxfyd*(EYSX1[0][j+1][k]-2.0*EYSX1[0][j][k]+EYSX1[0][j-1][k]+EYSX1[1][j+1][k]-2.0*EYSX1[1][j][k]+EYSX1[1][j-1][k])+cxfzd*(EYSX1[0][j][k+1]-2.0*EYSX1[0][j][k]+EYSX1[0][j][k-1]+EYSX1[1][j][k+1]-2.0*EYSX1[1][j][k]+EYSX1[1][j][k-1]);

	  EYS[NX-1][j][k] = -EYSX2[2][j][k]+cxd*(EYS[NX1-1][j][k]+EYSX2[3][j][k])+cxx*(EYSX1[3][j][k]+EYSX1[2][j][k])+cxfyd*(EYSX1[3][j+1][k]-2.0*EYSX1[3][j][k]+EYSX1[3][j-1][k]+EYSX1[2][j+1][k]-2.0*EYSX1[2][j][k]+EYSX1[2][j-1][k])+cxfzd*(EYSX1[3][j][k+1]-2.0*EYSX1[3][j][k]+EYSX1[3][j][k-1]+EYSX1[2][j][k+1]-2.0*EYSX1[2][j][k]+EYSX1[2][j][k-1]);
	}
    }
  // NOW SAVE PAST VALUES
  for(k = 1; k<NZ1; k++)
    {
      for(j = 0; j<NY1; j++)
	{
	  EYSX2[0][j][k]=EYSX1[0][j][k];
	  EYSX2[1][j][k]=EYSX1[1][j][k];
	  EYSX2[2][j][k]=EYSX1[2][j][k];
	  EYSX2[3][j][k]=EYSX1[3][j][k];
	  EYSX1[0][j][k]=EYS[0][j][k];
	  EYSX1[1][j][k]=EYS[1][j][k];
	  EYSX1[2][j][k]=EYS[NX1-1][j][k];
	  EYSX1[3][j][k]=EYS[NX-1][j][k];
	}
    }
}

void radezy()
{
  int k, i;

  // DO EDGES WITH FIRST ORDER ORBC

  for(k = 0; k<NZ1; k++)
    {
      i = 1;

      EZS[i][0][k] = EZSY1[i][1][k]+ cyd*(EZS[i][1][k] - EZSY1[i][0][k]);
      EZS[i][NY-1][k] = EZSY1[i][2][k] + cyd*(EZS[i][NY1-1][k] - EZSY1[i][3][k]);

      i = NX1-1;

      EZS[i][0][k] = EZSY1[i][1][k]+ cyd*(EZS[i][1][k] - EZSY1[i][0][k]);
      EZS[i][NY-1][k] = EZSY1[i][2][k] + cyd*(EZS[i][NY1-1][k] - EZSY1[i][3][k]);
    }

  for(i = 2; i<NX1-1; i++)
    {
      k = 0;

      EZS[i][0][k] = EZSY1[i][1][k]+ cyd*(EZS[i][1][k] - EZSY1[i][0][k]);
      EZS[i][NY-1][k] = EZSY1[i][2][k] + cyd*(EZS[i][NY1-1][k] - EZSY1[i][3][k]);

      k = NZ1-1;

      EZS[i][0][k] = EZSY1[i][1][k]+ cyd*(EZS[i][1][k] - EZSY1[i][0][k]);
      EZS[i][NY-1][k] = EZSY1[i][2][k] + cyd*(EZS[i][NY1-1][k] - EZSY1[i][3][k]);
    }

  // NOW DO 2ND ORDER ORBC ON REMAINING PORTIONS OF FACES

  for(k = 1; k<NZ1-1; k++)
    {
      for(i = 2; i<NX1-1; i++)
	{
	  EZS[i][0][k] = -EZSY2[i][1][k]+cyd*(EZS[i][1][k]+EZSY2[i][0][k])+cyy*(EZSY1[i][0][k]+EZSY1[i][1][k])+cyfxd*(EZSY1[i+1][0][k]-2.0*EZSY1[i][0][k]+EZSY1[i-1][0][k]+EZSY1[i+1][1][k]-2.0*EZSY1[i][1][k]+EZSY1[i-1][1][k])+cyfzd*(EZSY1[i][0][k+1]-2.0*EZSY1[i][0][k]+EZSY1[i][0][k-1]+EZSY1[i][1][k+1]-2.0*EZSY1[i][1][k]+EZSY1[i][1][k-1]);

	  EZS[i][NY-1][k] = -EZSY2[i][2][k]+cyd*(EZS[i][NY1-1][k]+EZSY2[i][3][k])+cyy*(EZSY1[i][3][k]+EZSY1[i][2][k])+cyfxd*(EZSY1[i+1][3][k]-2.0*EZSY1[i][3][k]+EZSY1[i-1][3][k]+EZSY1[i+1][2][k]-2.0*EZSY1[i][2][k]+EZSY1[i-1][2][k])+ cyfzd*(EZSY1[i][3][k+1]-2.0*EZSY1[i][3][k]+EZSY1[i][3][k-1]+EZSY1[i][2][k+1]-2.0*EZSY1[i][2][k]+EZSY1[i][2][k-1]);
	}
    }
  // NOW SAVE PAST VALUES

  for(k = 0; k<NZ1; k++)
    {
      for(i = 1; i<NX1; i++)
	{
	  EZSY2[i][0][k]=EZSY1[i][0][k];
	  EZSY2[i][1][k]=EZSY1[i][1][k];
	  EZSY2[i][2][k]=EZSY1[i][2][k];
	  EZSY2[i][3][k]=EZSY1[i][3][k];
	  EZSY1[i][0][k]=EZS[i][0][k];
	  EZSY1[i][1][k]=EZS[i][1][k];
	  EZSY1[i][2][k]=EZS[i][NY1-1][k];
	  EZSY1[i][3][k]=EZS[i][NY-1][k];
	}
    }
}

void radexy()
{
  int k, i;
  // DO EDGES WITH FIRST ORDER ORBC

  for(k = 1; k<NZ1; k++)
    {
      i = 0;
      EXS[i][0][k] = EXSY1[i][1][k]+ cyd*(EXS[i][1][k] - EXSY1[i][0][k]);

      EXS[i][NY-1][k] = EXSY1[i][2][k] + cyd*(EXS[i][NY1-1][k] - EXSY1[i][3][k]);

      i = NX1-1;

      EXS[i][0][k] = EXSY1[i][1][k]+ cyd*(EXS[i][1][k] - EXSY1[i][0][k]);
      EXS[i][NY-1][k] = EXSY1[i][2][k] + cyd*(EXS[i][NY1-1][k] - EXSY1[i][3][k]);
    }

  for(i = 1; i<NX1-1; i++)
    {
      k = 1;
      EXS[i][0][k] = EXSY1[i][1][k]+ cyd*(EXS[i][1][k] - EXSY1[i][0][k]);
      EXS[i][NY-1][k] = EXSY1[i][2][k] + cyd*(EXS[i][NY1-1][k] - EXSY1[i][3][k]);

      k = NZ1-1;

      EXS[i][0][k] = EXSY1[i][1][k]+ cyd*(EXS[i][1][k] - EXSY1[i][0][k]);
      EXS[i][NY-1][k] = EXSY1[i][2][k] + cyd*(EXS[i][NY1-1][k] - EXSY1[i][3][k]);

    }

  // NOW DO 2ND ORDER ORBC ON REMAINING PORTIONS OF FACES

  for(k = 2; k<NZ1-1; k++)
    {
      for(i = 1; i<NX1-1; i++)
	{
	  EXS[i][0][k] = -EXSY2[i][1][k]+cyd*(EXS[i][1][k]+EXSY2[i][0][k])+cyy*(EXSY1[i][0][k]+EXSY1[i][1][k])+cyfxd*(EXSY1[i+1][0][k]-2.0*EXSY1[i][0][k]+EXSY1[i-1][0][k]+EXSY1[i+1][1][k]-2.0*EXSY1[i][1][k]+EXSY1[i-1][1][k])+cyfzd*(EXSY1[i][0][k+1]-2.0*EXSY1[i][0][k]+EXSY1[i][0][k-1]+EXSY1[i][1][k+1]-2.0*EXSY1[i][1][k]+EXSY1[i][1][k-1]);

	  EXS[i][NY-1][k] = -EXSY2[i][2][k]+cyd*(EXS[i][NY1-1][k]+EXSY2[i][3][k])+cyy*(EXSY1[i][3][k]+EXSY1[i][2][k])+cyfxd*(EXSY1[i+1][3][k]-2.0*EXSY1[i][3][k]+EXSY1[i-1][3][k]+EXSY1[i+1][2][k]-2.0*EXSY1[i][2][k]+EXSY1[i-1][2][k])+cyfzd*(EXSY1[i][3][k+1]-2.0*EXSY1[i][3][k]+EXSY1[i][3][k-1]+EXSY1[i][2][k+1]-2.0*EXSY1[i][2][k]+EXSY1[i][2][k-1]);

	}
    }

  // NOW SAVE PAST VALUES

  for(k = 1; k<NZ1; k++)
    {
      for(i = 0; i<NX1; i++)
	{
	  EXSY2[i][0][k]=EXSY1[i][0][k];
	  EXSY2[i][1][k]=EXSY1[i][1][k];
	  EXSY2[i][2][k]=EXSY1[i][2][k];
	  EXSY2[i][3][k]=EXSY1[i][3][k];
	  EXSY1[i][0][k]=EXS[i][0][k];
	  EXSY1[i][1][k]=EXS[i][1][k];
	  EXSY1[i][2][k]=EXS[i][NY1-1][k];
	  EXSY1[i][3][k]=EXS[i][NY-1][k];
	}
    }
}

void radexz()
{
  int j, i;
  // DO EDGES WITH FIRST ORDER ORBC

  for(j = 1; j<NY1; j++)
    {
      i = 0;
      EXS[i][j][0] = EXSZ1[i][j][1]+ czd*(EXS[i][j][1] - EXSZ1[i][j][0]);
      EXS[i][j][NZ-1] = EXSZ1[i][j][2] + czd*(EXS[i][j][NZ1-1] - EXSZ1[i][j][3]);

      i = NX1-1;

      EXS[i][j][0] = EXSZ1[i][j][1]+ czd*(EXS[i][j][1] - EXSZ1[i][j][0]);
      EXS[i][j][NZ-1] = EXSZ1[i][j][2] + czd*(EXS[i][j][NZ1-1] - EXSZ1[i][j][3]);
    }

  for(i = 1; i<NX1-1; i++)
    {
      j = 1;
      EXS[i][j][0] = EXSZ1[i][j][1]+ czd*(EXS[i][j][1] - EXSZ1[i][j][0]);
      EXS[i][j][NZ-1] = EXSZ1[i][j][2] + czd*(EXS[i][j][NZ1-1] - EXSZ1[i][j][3]);

      j = NY1-1;

      EXS[i][j][0] = EXSZ1[i][j][1]+ czd*(EXS[i][j][1] - EXSZ1[i][j][0]);
      EXS[i][j][NZ-1] = EXSZ1[i][j][2] + czd*(EXS[i][j][NZ1-1] - EXSZ1[i][j][3]);

    }

  // NOW DO 2ND ORDER ORBC ON REMAINING PORTIONS OF FACES

  for(j = 2; j<NY1-1; j++)
    {
      for(i = 1; i<NX1-1; i++)
	{
	  EXS[i][j][0] = -EXSZ2[i][j][1]+czd*(EXS[i][j][1]+EXSZ2[i][j][0])+czz*(EXSZ1[i][j][0]+EXSZ1[i][j][1])+czfxd*(EXSZ1[i+1][j][0]-2.0*EXSZ1[i][j][0]+EXSZ1[i-1][j][0]+EXSZ1[i+1][j][1]-2.0*EXSZ1[i][j][1]+EXSZ1[i-1][j][1])+czfyd*(EXSZ1[i][j+1][0]-2.0*EXSZ1[i][j][0]+EXSZ1[i][j-1][0]+EXSZ1[i][j+1][1]-2.0*EXSZ1[i][j][1]+EXSZ1[i][j-1][1]);

	  EXS[i][j][NZ-1] = -EXSZ2[i][j][2]+czd*(EXS[i][j][NZ1-1]+EXSZ2[i][j][3])+czz*(EXSZ1[i][j][3]+EXSZ1[i][j][2])+czfxd*(EXSZ1[i+1][j][3]-2.0*EXSZ1[i][j][3]+EXSZ1[i-1][j][3]+EXSZ1[i+1][j][2]-2.0*EXSZ1[i][j][2]+EXSZ1[i-1][j][2])+czfyd*(EXSZ1[i][j+1][3]-2.0*EXSZ1[i][j][3]+EXSZ1[i][j-1][3]+EXSZ1[i][j+1][2]-2.0*EXSZ1[i][j][2]+EXSZ1[i][j-1][2]);
	}
    }
  // NOW SAVE PAST VALUES

  for(j = 1; j<NY1; j++)
    {
      for(i = 0; i<NX1; i++)
	{
	  EXSZ2[i][j][0]=EXSZ1[i][j][0];
	  EXSZ2[i][j][1]=EXSZ1[i][j][1];
	  EXSZ2[i][j][2]=EXSZ1[i][j][2];
	  EXSZ2[i][j][3]=EXSZ1[i][j][3];
	  EXSZ1[i][j][0]=EXS[i][j][0];
	  EXSZ1[i][j][1]=EXS[i][j][1];
	  EXSZ1[i][j][2]=EXS[i][j][NZ1-1];
	  EXSZ1[i][j][3]=EXS[i][j][NZ-1];
	}
    }
}

void radeyz()
{
  int i, j;
  // DO EDGES WITH FIRST ORDER ORBC

  for(j = 0; j<NY1; j++)
    {
      i = 1;
      EYS[i][j][0] = EYSZ1[i][j][1]+ czd*(EYS[i][j][1] - EYSZ1[i][j][0]);
      EYS[i][j][NZ-1] = EYSZ1[i][j][2] + czd*(EYS[i][j][NZ1-1] - EYSZ1[i][j][3]);

      i = NX1-1;
      EYS[i][j][0] = EYSZ1[i][j][1]+ czd*(EYS[i][j][1] - EYSZ1[i][j][0]);
      EYS[i][j][NZ-1] = EYSZ1[i][j][2] + czd*(EYS[i][j][NZ1-1] - EYSZ1[i][j][3]);
    }

  for(i = 2; i<NX1-1; i++)
    {
      j = 0;
      EYS[i][j][0] = EYSZ1[i][j][1]+ czd*(EYS[i][j][1] - EYSZ1[i][j][0]);
      EYS[i][j][NZ-1] = EYSZ1[i][j][2] + czd*(EYS[i][j][NZ1-1] - EYSZ1[i][j][3]);

      j = NY1-1;
      EYS[i][j][0] = EYSZ1[i][j][1]+ czd*(EYS[i][j][1] - EYSZ1[i][j][0]);
      EYS[i][j][NZ-1] = EYSZ1[i][j][2] + czd*(EYS[i][j][NZ1-1] - EYSZ1[i][j][3]);
    }

  // NOW DO 2ND ORDER ORBC ON REMAINING PORTIONS OF FACES

  for(j = 1; j<NY1-1; j++)
    {
      for(i = 2; i<NX1-1; i++)
	{
	  EYS[i][j][0] = -EYSZ2[i][j][1]+ czd*(EYS[i][j][1]+EYSZ2[i][j][0])+ czz*(EYSZ1[i][j][0]+EYSZ1[i][j][1])+ czfxd*(EYSZ1[i+1][j][0]- 2*EYSZ1[i][j][0]+EYSZ1[i-1][j][0]+EYSZ1[i+1][j][1]- 2*EYSZ1[i][j][1]+EYSZ1[i-1][j][1])+ czfyd*(EYSZ1[i][j+1][0]- 2*EYSZ1[i][j][0]+EYSZ1[i][j-1][0]+EYSZ1[i][j+1][1]- 2*EYSZ1[i][j][1]+EYSZ1[i][j-1][1]);

	  EYS[i][j][NZ-1] = -EYSZ2[i][j][2]+czd*(EYS[i][j][NZ1-1]+EYSZ2[i][j][3])+czz*(EYSZ1[i][j][3]+EYSZ1[i][j][2])+czfxd*(EYSZ1[i+1][j][3]-2.0*EYSZ1[i][j][3]+EYSZ1[i-1][j][3]+EYSZ1[i+1][j][2]-2.0*EYSZ1[i][j][2]+EYSZ1[i-1][j][2])+czfyd*(EYSZ1[i][j+1][3]-2.0*EYSZ1[i][j][3]+EYSZ1[i][j-1][3]+EYSZ1[i][j+1][2]-2.0*EYSZ1[i][j][2]+EYSZ1[i][j-1][2]);

	}
    }

  // NOW SAVE PAST VALUES

  for(j = 0; j<NY1; j++)
    {
      for(i = 1; i<NX1; i++)
	{
	  EYSZ2[i][j][0]=EYSZ1[i][j][0];
	  EYSZ2[i][j][1]=EYSZ1[i][j][1];
	  EYSZ2[i][j][2]=EYSZ1[i][j][2];
	  EYSZ2[i][j][3]=EYSZ1[i][j][3];
	  EYSZ1[i][j][0]=EYS[i][j][0];
	  EYSZ1[i][j][1]=EYS[i][j][1];
	  EYSZ1[i][j][2]=EYS[i][j][NZ1-1];
	  EYSZ1[i][j][3]=EYS[i][j][NZ-1];
	}
    }
}

void hxsfld()
{
  int i, j, k;

  for(k = 0; k<NZ1; k++)
    {
      for(j = 0; j<NY1; j++)
	{
	  for(i = 1; i<NX1; i++)
	    {
	      HXS[i][j][k]=HXS[i][j][k]-(EZS[i][j+1][k]-EZS[i][j][k])*dtmdy+(EYS[i][j][k+1]-EYS[i][j][k])*dtmdz;
	      //printf("%.11f\n", EZS[i][j+1][k]);
	    }
	}
    }
}

void hysfld()
{
  int i, j, k;
  for(k = 0; k<NZ1; k++)
    {
      for(j = 1; j<NY1; j++)
	{
	  for(i = 0; i<NX1; i++)
	    {
	      HYS[i][j][k]=HYS[i][j][k]-(EXS[i][j][k+1]-EXS[i][j][k])*dtmdz+(EZS[i+1][j][k]-EZS[i][j][k])*dtmdx;
	    }
	}
    }
}

void hzsfld()
{
  int i, j, k;

  for(k = 1; k<NZ1; k++)
    {
      for(j = 0; j<NY1; j++)
	{
	  for(i = 0; i<NX1; i++)
	    {
	      HZS[i][j][k]=HZS[i][j][k]-(EYS[i+1][j][k]-EYS[i][j][k])*dtmdx+(EXS[i][j+1][k]-EXS[i][j][k])*dtmdy;
	    }
	}
    }
}

double EXI(int i, int j, int k)
{
	// THIS FUNCTION COMPUTES THE X COMPONENT OF THE INCIDENT ELECTRIC FIELD
	//if(i == 0 && j == 0 && k == 0)
	//printf("\ndist = %.11f\n", t);
	double dist;

	// First calc. the distance to the specified cell (i,j,k):

	dist = ((i)*delx+0.5*delx*off)*xdisp+((j)*dely)*ydisp+((k)*delz)*zdisp + delay;
	//printf("\ndist = %.11f\n", t);

	double s = source(dist);
	//printf("\nsource = %.11f\n", s);
    return ampx*s;
}

double EYI(int i, int j, int k)
{
	// THIS FUNCTION COMPUTES THE Y COMPONENT OF THE INCIDENT ELECTRIC FIELD

	double dist;
	// First calc. the distance to the specified cell (i,j,k):

	dist = ((i)*delx)*xdisp+((j)*dely+0.5*dely*off)*ydisp+((k)*delz)*zdisp + delay;
    return ampy*source(dist);
}

double EZI(int i, int j, int k)
{

	double dist;

	dist = ((i)*delx)*xdisp+((j)*dely)*ydisp+((k)*delz+0.5*delz*off)*zdisp + delay;
    return ampz*source(dist);
}

double DEXI(int i, int j, int k)
{
	double dist;

	dist = ((i)*delx+0.5*delx*off)*xdisp+((j)*dely)*ydisp+((k)*delz)*zdisp + delay;
    return ampx*dsrce(dist);
}

double DEYI(int i, int j, int k)
{
	// THIS FUNCTION COMPUTES THE Y COMPONENT OF THE INCIDENT ELECTRIC FIELD

	double dist;

	// First calc. the distance to the specified cell (i,j,k):

	dist = ((i)*delx)*xdisp+((j)*dely+0.5*dely*off)*ydisp+((k)*delz)*zdisp+delay;
    return ampy*dsrce(dist);
}

double DEZI(int i, int j, int k)
{
	// THIS FUNCTION COMPUTES THE Z COMPONENT OF THE INCIDENT ELECTRIC FIELD

	double dist;

	// First calc. the distance to the specified cell (i,j,k):

	dist = ((i)*delx)*xdisp+((j)*dely*off)*ydisp+((k)*delz+0.5*delz*off)*zdisp+delay;
    return ampz*dsrce(dist);
}


double source(double dist)
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

double dsrce(double dist)
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

void datsav(){
  
  //printf("in datsav\n");
  int i,j,k;
  double exc,eyc,ezc,exin,eyin,ezin,temp,r1,esq,sarx,r2,r3;
  nrms=(1.0/(freq*dt));
  //printf("nrms: %lf\n",nrms);
  r2=50.0e-2,r3=10.0e-2;
  if(N==1){

    for(i=0;i<NX;i++){
      s[i]=0.0;
      erms1d[i]=0.0;
    }
  }

  for(i=0;i<NX;i++){
    for(j=0;j<NZ;j++){
      sar[i][j]=0.0;
      erms[i][j]=0.0;
    }
  }

  if(N==1){
    memset(etimeavg,0,sizeof(double)*NX*NY*NZ);
    memset(emax,0,sizeof(double)*NX*NY*NZ);
    sar_Total=0.0;
    sar_wb=0.0;
  }

  cont=0;

  if(N>=(nstop-nrms) && N<=(nstop-1)){
    for(i=0;i<NX-2;i++){
      for(j=0;j<NY-2;j++){
	for(k=0;k<NZ-2;k++){
	  esq=0.0;
	  t=t-dt;
	  exc=0.0;
	  eyc=0.0;
	  ezc=0.0;

	  exin=EXI(i,j,k)+EXI(i,j,k+1)+EXI(i,j+1,k)+EXI(i,j+1,k+1);
	  eyin=EYI(i,j,k)+EYI(i+1,j,k)+EYI(i,j,k+1)+EYI(i+1,j,k+1);
	  ezin=EZI(i,j,k)+EZI(i+1,j,k)+EZI(i+1,j+1,k)+EZI(i,j+1,k);
	  t=t+dt;
	  exc=(exin+EXS[i][j][k]+EXS[i][j][k+1]+EXS[i][j+1][k]+EXS[i][j+1][k+1])*(0.25);
	  eyc=(eyin+EYS[i][j][k]+EYS[i+1][j][k]+EYS[i][j][k+1]+EYS[i+1][j][k+1])*(0.25);
	  ezc=(ezin+EZS[i][j][k]+EZS[i+1][j][k]+EZS[i+1][j+1][k]+EZS[i][j+1][k])*(0.25);

	  esq=((exc*exc)+(eyc*eyc)+(ezc*ezc));
	  (etimeavg[offset(i,j,k)])=(etimeavg[offset(i,j,k)])+(esq/nrms);
	  if(N==(nstop-nrms)) (emax[offset(i,j,k)])=esq;
	  if(esq>(emax[offset(i,j,k)])) (emax[offset(i,j,k)])=esq;
	  cont=cont+1;
	}
      }
    }
    //printf("done cal\n");
    //fprintf(fsar,"%d %.15lf\n",N,emax[offset(50,50,50)]);
  }
  //fprintf(fsar,"%d %.15lf\n",N,etimeavg[offset(100,100,100)]);

  if(N==(nstop-1)){
    printf("Number of cells: %d\n",cont);
    j=nyc;
    k=nzc;

    for(i=0;i<NX-2;i++){
      s[i]=(0.5)*SIGMA[IDONE[i][j][k]-1]*(emax[offset(i,j,k)])*(1.0/iden);

      erms1d[i]=sqrt(etimeavg[offset(i,j,k)])*(0.707);
    }
  }
//FILE *fp1;
//fp1=fopen("etimeavg.txt", "w");
  if(N==(nstop-1)){
    j=nyc;
    for(k=0;k<NZ-2;k++){
      for(i=0;i<NX-2;i++){
	sar[i][k]=(0.5)*SIGMA[IDONE[i][j][k]-1]*(emax[offset(i,j,k)])*(1.0/iden);

	erms[i][k]=sqrt(etimeavg[offset(i,j,k)])*(0.707);
   //fprintf(fp1,"%d %d %.15lf\n",i,k,erms[i][k]);
      }
    }
  }
//fclose(fp1);
  if(N==(nstop-1)){
    cont2=0;
    temp=0.0;
    sarx=0.0;
    wholevol=0.0;

    if(wbsar==1){
      for(i=0;i<NX-2;i++){
	for(j=0;j<NY-2;j++){
	  for(k=0;k<NZ-2;k++){
	    //temp=((((i-nxc)*delx)*((i-nxc)*delx))+(((j-nyc)*dely)*((j-nyc)*dely))+(((k-nzc)*delz)*((k-nzc)*delz)));
           temp=((pow((i-nxc)*delx,2)/pow(r3,2)) + (pow((j-nyc)*dely,2)/pow(r3,2)) + (pow((k-nzc)*delz,2)/pow(r2,2)));
	   // r1=sqrt(temp);

	    //if(r1<=radius2){
          if(temp>0.0 && temp<=1.0){
	      cont2++;
	      sarx=sarx+((0.5)*SIGMA[IDONE[i][j][k]-1]*(emax[offset(i,j,k)])*(delx*dely*delz));
	      wholevol+=(delx*dely*delz);
	    }
	  }
	}
      }
      printf("number of cell in whole body sar calculation: %d\n",cont2);
      printf("sarx: %.15lf\n wholevol: %.15lf\n",sarx,wholevol);
      sar_wb=sarx/(wholevol*iden);
      printf("total power: %lf\nsar_wb: %lf\n",sarx,sar_wb);
    }
  }

  /* t=t-dt; */
  /* dum1=EXI(IOBS[0],JOBS[0],KOBS[0]); */
  /* t=t+dt; */
  return;
}

