#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "params.h"
#include "comput_fourierb.h"


double func_delta( int i, int j){

      double delta;

      if(i==j){
         delta=1.0;
      }else{
         delta=0.0;
      } 

      return delta;

}


//   Compute the Fourier components for V matrix (at z=iz*dz)
//   similar to B-16   
void computevnn( int jz, int nx, int nz, double czr, double ome, double dx, double *vnnr, double *vnnim, double *vel){

       int j, k, l;
       int j1, k1, l1;
       int mm, mm1;
       int nn = 2*nx+1; 
       double sum1=0.0, sum2=0.0,sum3=0.0,sum4=0.0;
       double c = (nn-1)*dx;
       double delta = 0.0;
       double* cmin= malloc(2*sizeof(double));;
       double* xnr   = malloc(nn*sizeof(double));
       double* xnim   = malloc(nn*sizeof(double));     
       c = 1.0/sqrt(c);
       vsn( jz, nx, nn, dx, xnr, xnim, vel, czr,cmin);

   /*   FILE* file2;
      char fname2[100];
      sprintf(fname2,"vnnc.dat");
      file2 = fopen(fname2,"w");*/


       for (j=0;j<nn;j++){

            j1 = -nx + j;

	       for (k=0;k<nn;k++){

		    k1 = -nx + k;
                    mm = j*nn + k;
                    sum1 = 0.0;
                    sum2 = 0.0;
                  

	            for (l=0;l<nn;l++){

			    l1  = -nx + l;
		            mm1 = l;

                            delta = func_delta( j1, k1+l1);

                            sum1 = sum1 +  xnr[mm1] * delta;
                            sum2 = sum2 + xnim[mm1] * delta;
                          

                            //fprintf(file2," %d %lf %lf %lf %lf %lf \n", l, delta, sum1, sum2, sum3, sum4); 

                    }

                    vnnr[mm]  = sum1*c;
                    vnnim[mm] = sum2*c;
                            //fprintf(file2," %.12lf %.12lf\n", vnnabcr[mm], vnnabcim[mm]); 

                    //fclose(file2);

                    

	       }
       }
	czr=cmin[0];
	
	//fprintf(stdout,"%d  %lf \n",jz, czr);
	
	double* projppv   = malloc(nn*sizeof(double));
	projp( ome, nx, dx, czr, projppv);
	applyprojv( nx, vnnr, vnnim, projppv);
       free(xnr);
       free(xnim); 
     

}


void vsn( int jz, int nx, int nn, double dx, double *xnr, double *xnim, double *vseg, double czr,double* cmin){

       int j, nq, nq1;
       int mm, mm1;
       double sum1=0.0, sum2=0.0;
       double vel=0.0, xx=0.0, c=0.0;
       double al = (nn-1)*dx;

       c = 1.0/sqrt(al);
       cmin[0]=czr;
	 for (j=0;j<nn;j++){
                 mm  = jz*nn+j;
	 	if (vseg[mm]>cmin[0]){
			cmin[0]=vseg[mm];
		 }
	}
       for (nq=0;nq<nn;nq++){

            mm1  = nq;
//            nq1  = -nx + nq - 1;
            nq1  = -nx + nq;
            sum1 = 0.0;
            sum2 = 0.0;

            	 j=0;
                 mm  = jz*nn+j;
                
                 vel = 1.0/vseg[mm];
                 vel = 1.0 - vel*vel*czr*czr;

                 xx = j*dx; 

                 sum1 = sum1 + 0.5*dx*cos(2.0*M_PI*nq1*xx/al)*vel;                  
                 sum2 = sum2 - 0.5*dx*sin(2.0*M_PI*nq1*xx/al)*vel;                  

            for (j=1;j<nn-1;j++){

                 mm  = jz*nn+j;
//                 vel = czr/vseg[mm];
//                 vel = 1.0 - vel*vel;
                 vel = 1.0/vseg[mm];
                 vel = 1.0 - vel*vel*czr*czr;

                 xx = j*dx; 

                 sum1 = sum1 + dx*cos(2.0*M_PI*nq1*xx/al)*vel;                  
                 sum2 = sum2 - dx*sin(2.0*M_PI*nq1*xx/al)*vel;                  

            }

            	 j=nn-1;
                 mm  = jz*nn+j;
//                 vel = czr/vseg[mm];
//                 vel = 1.0 - vel*vel;
                 vel = 1.0/vseg[mm];
                 vel = 1.0 - vel*vel*czr*czr;

                 xx = j*dx; 

                 sum1 = sum1 + 0.5*dx*cos(2.0*M_PI*nq1*xx/al)*vel;                  
                 sum2 = sum2 - 0.5*dx*sin(2.0*M_PI*nq1*xx/al)*vel;                  

            xnr[mm1]  = sum1*c;
            xnim[mm1] = sum2*c;

       }

}




void projp( double ome, int nx, double dx, double czr, double *projpp){

      int nn = 2*nx+1;
      int kd = nn*nn;
      int i=0, j=0, k=0, j1=0;
      int mm=0;

      double am0j = 0.0;
      double al = (nn-1)*dx;

      // diagonal matrices
      double* amr   = malloc(nn*sizeof(double));     // lambda_j
    //  fprintf(stdout," nn %d al %lf czr %lf ome %lf pi %lf \n", nn, al, czr, ome, M_PI);      

     for (j=0;j<nn;j++){

          j1 = -nx + j;

//      amr = S = lambda_j (+ integration Fourier basis)
//      eq B-15
//      if lambda_j >0, eigenvalue of M is real
        am0j = (2.0*M_PI*j1)*(2.0*M_PI*j1)/(al*al) - (ome*ome)/(czr*czr);
        amr[j] = am0j;

        if(am0j<0.0){
            projpp[j] = 1.0;
            //projpq[j] = 0.0; 
            //fprintf(stdout," j1 %d am0j %lf sqrt %lf \n", j1, am0j, sqrt(fabs(am0j)));

        }else{
            projpp[j] = 0.0;
//            projpp[j] = 0.5;
            //projpq[j] = 0.5/sqrt(abs(am0j));
            //fprintf(stdout," j1 %d am0j %lf sqrt %lf \n", j1, am0j, sqrt(fabs(am0j)));

        }
    }


}

void applyproj( int nx, double *p1zr, double *p1zim, double *p2zr, double *p2zim, double *projpp){

    int i=0, j=0, mm=0;
    int nn = 2*nx+1;

     for (i=0;i<nn;i++){
     	for (j=0;j<nn;j++){
            mm = i*nn+j;
            p1zr[mm] = projpp[i]*p1zr[mm];
            p1zim[mm] = projpp[i]*p1zim[mm];
            p2zr[mm] = projpp[i]*p2zr[mm];
            p2zim[mm] = projpp[i]*p2zim[mm];
     	}
     }

}
void applyprojv( int nx, double *vnnr, double *vnnim, double *projpp){

    int i=0, j=0, mm=0;
    int nn = 2*nx+1;

     for (i=0;i<nn;i++){
     	for (j=0;j<nn;j++){
            mm = i*nn+j;
            vnnr[mm] = projpp[i]*vnnr[mm];
            vnnim[mm] = projpp[i]*vnnim[mm];
           
     	}
     }

}
