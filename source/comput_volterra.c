#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "params.h"
#include "comput_volterra.h"
#include "comput_fourierb.h"

void p1init( int nx, int nn, double dx, double ome, double czr, double zzer, double *p1zr, double *p1zim){

     
       int i=0, n=0, mm=0;

       double ak = ome/czr;
       double al = (nn-1)*dx;
       double zz = zzer; 
       double akn=0.0, akn1=0.0, delta=0.0;

        int ka = nn*nn;

         // intialization
              for (i=0;i<ka;i++){
                    p1zr[i]=0.0;
                    p1zim[i]=0.0;
              }


       // intialization is diagonal !!! cf eq a-25 
       // only propagating wave
       for (i=0;i<nn;i++){

            n = -nx + i;
            akn = 2.0*M_PI*n/al;
            akn *= -akn;
            akn += ak*ak;  

            mm = i*nn +i; 
            akn1=0.0;  

	        if(akn>=0.0){

                akn1 = sqrt(akn); 
		       p1zr[mm]  = cos(akn1*zz);
		       p1zim[mm] = sin(akn1*zz);

		}
                /*else{

                akn1 = -sqrt(fabs(akn)); 
		       p1zr[mm]  = exp(akn1*zz);
		       p1zim[mm] = 0.0;

                }*/

       }

}


void aknmat( int nx, double ome, double czr, double dx, double *aknr, double *aknim, double *akninvr, double *akninvim){

      int nn = 2*nx + 1;
      int i=0, i1=0;

      double ak = ome/czr;
      double al= (nn-1)*dx;
      double akn, akn1;


      FILE* file2;
      char fname2[100];
        sprintf(fname2,"akn.dat");
	file2 = fopen(fname2,"w");

	for (i=0;i<nn;i++){
             i1 = -nx + i;
             akn = 2.0*M_PI*i1/al;
             akn *= akn;
             akn = ak*ak - akn;
             akn1 = akn;
             if(akn>=0.0){
               akn1 = sqrt(akn1);
	       aknr[i]  = akn1;
	       aknim[i] = 0.0;
	       akninvr[i] = 1.0/akn1;
	       akninvim[i] = 0.0;
             }else{
               akn1 = sqrt(fabs(akn1));
	       aknr[i]    = 0.0;
	       aknim[i]   = 0.0;
//	       aknim[i]   = akn1;
	       akninvr[i] = 0.0;
//	       akninvim[i] =-1.0/akn1;
	       akninvim[i] =0.0;
             }
          fprintf(file2," %d %lf %lf \n", i, aknr[i], akninvr[i]);
        }
       fclose(file2);

}

void compute_volterrap1( double *p1r, double *p1im,int nx, int jz, double dz, double ak2, int iw, int *itrv,
                  double *vnnmr, double *vnnmim, double *aknr, double *akninvr,
                  double *rst1pr, double *rst1pim, double *rst1mr, double *rst1mim){

      int i=0, j=0, k=0;
      int ii=0, jj=0;
      int mm=0, mm1=0, mm2=0;
      int nn = 2*nx+1; 
      int ka = nn*nn;   
      double intp1zr=0.0, intp1zim=0.0;
      double z=jz*dz;
      double sum1, sum2;
      double delta=0.0;

      double* expozpr = malloc(nn*sizeof(double));
      double* expozpim = malloc(nn*sizeof(double));
      double* expozmr = malloc(nn*sizeof(double));
      double* expozmim = malloc(nn*sizeof(double));

      double* prodvp1r = malloc(ka*sizeof(double));
      double* prodvp1im = malloc(ka*sizeof(double));

      double* rst1ppr  = malloc(ka*sizeof(double));
      double* rst1ppim = malloc(ka*sizeof(double));
      double* rst1mpr  = malloc(ka*sizeof(double));
      double* rst1mpim = malloc(ka*sizeof(double));


      FILE* file2;
      char fname2[100];

      comput_expoz( nx, jz-1, dz, iw, itrv, aknr, expozpr, expozpim, expozmr, expozmim);


     for(i=0;i<nn;i++){

	      for(j=0;j<nn;j++){

                  mm = i*nn + j;
                  sum1 = 0.0;
                  sum2 = 0.0;

	      	for(k=0;k<nn;k++){

                   mm1 = i*nn+k;
                   mm2 = k*nn+j;

                   sum1 += vnnmr[mm1]*p1r[mm2] - vnnmim[mm1]*p1im[mm2];
                   sum2 += vnnmr[mm1]*p1im[mm2] + vnnmim[mm1]*p1r[mm2];
                     
              	} 

		prodvp1r[mm]  = sum1;        
		prodvp1im[mm] = sum2;   
     
                  // int e+ v p1
                  intp1zr  = expozpr[i]*prodvp1r[mm]  - expozpim[i]*prodvp1im[mm];
                  intp1zim = expozpr[i]*prodvp1im[mm] + expozpim[i]*prodvp1r[mm];
		  rst1ppr[mm]  = rst1pr[mm]  + dz*intp1zr; 
                  rst1ppim[mm] = rst1pim[mm] + dz*intp1zim; 
     
                  // int e- v p1
                  intp1zr  = expozmr[i]*prodvp1r[mm]  - expozmim[i]*prodvp1im[mm];
                  intp1zim = expozmr[i]*prodvp1im[mm] + expozmim[i]*prodvp1r[mm];
		  rst1mpr[mm]  = rst1mr[mm]  + dz*intp1zr; 
                  rst1mpim[mm] = rst1mim[mm] + dz*intp1zim; 

	      }
         
      }

    

      comput_expoz( nx, jz, dz, iw, itrv, aknr, expozpr, expozpim, expozmr, expozmim);
       
      for(i=0;i<nn;i++){

	      for(j=0;j<nn;j++){

                  mm = i*nn + j;
                  delta = func_delta( i, j);

               
			
		// -i e+ int e- v p1 + i e- int e+ v p1
                  intp1zr = expozpr[i]*rst1mpim[mm]+expozpim[i]*rst1mpr[mm];
                  intp1zr -= expozmr[i]*rst1ppim[mm]+expozmim[i]*rst1ppr[mm];

                  intp1zim = -expozpr[i]*rst1mpr[mm]+expozpim[i]*rst1mpim[mm];
                  intp1zim += expozmr[i]*rst1ppr[mm]-expozmim[i]*rst1ppim[mm];

                  intp1zr  *= ak2*0.5*akninvr[i];
                  intp1zim *= ak2*0.5*akninvr[i];

                  p1r[mm]  = delta*expozpr[i];
                  p1im[mm] = delta*expozpim[i];

                  p1r[mm]  += intp1zr;
                  p1im[mm] += intp1zim; 

              }
      } 

	for(i=0;i<nn;i++){

	      for(j=0;j<nn;j++){

                  mm = i*nn + j;
                  rst1pr[mm]  = rst1ppr[mm]  ; 
                  rst1pim[mm] = rst1ppim[mm] ; 
     
                  rst1mr[mm]  = rst1mpr[mm] ; 
                  rst1mim[mm] = rst1mpim[mm] ; 

	      }
         
      }
      free(expozpr);
      free(expozpim);
      free(expozmr);
      free(expozmim);
      free(prodvp1r);
      free(prodvp1im);
      free(rst1ppr);
      free(rst1ppim);
      free(rst1mpr);
      free(rst1mpim);


}


void compute_volterrap2( double *p2r, double *p2im,int nx, int jz, double dz, double ak2, int iw, int *itrv,
                  double *vnnmr, double *vnnmim,  double *aknr, double *akninvr,
                  double *rst1pr, double *rst1pim, double *rst1mr, double *rst1mim){

      int i=0, j=0, k=0, ii=0, jj=0;
      int mm=0, mm1=0, mm2=0;
      int nn = 2*nx+1; 
      int ka = nn*nn;   
      double intp2zr=0.0, intp2zim=0.0;
      double z=jz*dz;
      double sum1, sum2;
      double delta=0.0;

      double* expozpr = malloc(nn*sizeof(double));
      double* expozpim = malloc(nn*sizeof(double));
      double* expozmr = malloc(nn*sizeof(double));
      double* expozmim = malloc(nn*sizeof(double));

      double* prodvp2r = malloc(ka*sizeof(double));
      double* prodvp2im = malloc(ka*sizeof(double));

      double* rst1ppr  = malloc(ka*sizeof(double));
      double* rst1ppim = malloc(ka*sizeof(double));
      double* rst1mpr  = malloc(ka*sizeof(double));
      double* rst1mpim = malloc(ka*sizeof(double));


      FILE* file2;
      char fname2[100];

      comput_expoz( nx, jz-1, dz, iw, itrv, aknr, expozpr, expozpim, expozmr, expozmim);

     for(i=0;i<nn;i++){

	      for(j=0;j<nn;j++){

                  mm = i*nn + j;
                  sum1 = 0.0;
                  sum2 = 0.0;

	      	for(k=0;k<nn;k++){

                   mm1 = i*nn+k;
                   mm2 = k*nn+j;

                   sum1 += vnnmr[mm1]*p2r[mm2] - vnnmim[mm1]*p2im[mm2];
                   sum2 += vnnmr[mm1]*p2im[mm2] + vnnmim[mm1]*p2r[mm2];
                     
              	} 

		prodvp2r[mm]  = sum1;        
		prodvp2im[mm] = sum2;   
     
                  // int e+ v p2
                  intp2zr  = expozpr[i]*prodvp2r[mm]  - expozpim[i]*prodvp2im[mm];
                  intp2zim = expozpr[i]*prodvp2im[mm] + expozpim[i]*prodvp2r[mm];
                  rst1ppr[mm]  = rst1pr[mm]  + dz*intp2zr; 
                  rst1ppim[mm] = rst1pim[mm] + dz*intp2zim; 
     
                  // int e- v p2
                  intp2zr  = expozmr[i]*prodvp2r[mm]  - expozmim[i]*prodvp2im[mm];
                  intp2zim = expozmr[i]*prodvp2im[mm] + expozmim[i]*prodvp2r[mm];
                  rst1mpr[mm]  = rst1mr[mm]  + dz*intp2zr; 
                  rst1mpim[mm] = rst1mim[mm] + dz*intp2zim; 

	      }
         
      }


    

      comput_expoz( nx, jz, dz, iw, itrv, aknr, expozpr, expozpim, expozmr, expozmim);

      for(i=0;i<nn;i++){

	      for(j=0;j<nn;j++){

                  mm = i*nn + j;
                  delta = func_delta( i, j);

                  // -i e+ int e- v p2 + i e- int e+ v p2
                  intp2zr = expozpr[i]*rst1mpim[mm]+expozpim[i]*rst1mpr[mm];
                  intp2zr -= expozmr[i]*rst1ppim[mm]+expozmim[i]*rst1ppr[mm];

                  intp2zim = -expozpr[i]*rst1mpr[mm]+expozpim[i]*rst1mpim[mm];
                  intp2zim += expozmr[i]*rst1ppr[mm]-expozmim[i]*rst1ppim[mm];

                  intp2zr  *= ak2*0.5*akninvr[i];
                  intp2zim *= ak2*0.5*akninvr[i];

                  p2r[mm]  = delta*expozmr[i];
                  p2im[mm] = delta*expozmim[i];

                  p2r[mm]  += intp2zr;
                  p2im[mm] += intp2zim; 

              }
      } 

	
	for(i=0;i<nn;i++){

	      for(j=0;j<nn;j++){

                  mm = i*nn + j;
                  rst1pr[mm]  = rst1ppr[mm]  ; 
                  rst1pim[mm] = rst1ppim[mm] ; 
     
                  rst1mr[mm]  = rst1mpr[mm] ; 
                  rst1mim[mm] = rst1mpim[mm] ; 

	      }
         
      }
     
      free(expozpr);
      free(expozpim);
      free(expozmr);
      free(expozmim);

      free(prodvp2r);
      free(prodvp2im);
      free(rst1ppr);
      free(rst1ppim);
      free(rst1mpr);
      free(rst1mpim);

}



void comput_expoz( int nx, int jz, double dz, int iw, int *itrv, double *aknr, double *expozpr, double *expozpim, double *expozmr, double *expozmim){

      int i=0, i1=0;
      int nn = 2*nx+1;
      double z=jz*dz; 

//      fprintf(stdout," z %lf \n", z); 

      for(i=0;i<nn;i++){
        expozpr[i]  = 0.0;
        expozpim[i] = 0.0;
        expozmr[i]  = 0.0;
        expozmim[i] = 0.0;
      } 
     
      for(i=0;i<itrv[iw];i++){
         i1 = nx - itrv[iw]/2 + i;
        expozpr[i1]  = cos(aknr[i1]*z);
        expozpim[i1] = sin(aknr[i1]*z);
        expozmr[i1]  = cos(aknr[i1]*z);
        expozmim[i1] = -1.0*sin(aknr[i1]*z);
        //fprintf(stdout," i %d i1 %d akn %lf \n", i, i1, aknr[i1]);
      }
      //exit(1);
}
