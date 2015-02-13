#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "params.h"
#include "comput_reflection.h"
#include "matrice.h"

void mult( int nx, int iom, double ak2, double *rst1pr, double *rst1pim, double *rst2pr, double *rst2pim, int *itrv, double *akninvr){

      int i, j, ii, jj, mm;
      int nn = 2*nx+1;
     double st1=0, st2=0;

	for (i=0;i<itrv[iom];i++){
               ii = nx - (itrv[iom]-1)/2+ i;
		for (j=0;j<itrv[iom];j++){
                	jj = nx - (itrv[iom]-1)/2+ j;
                        mm = ii*nn+jj;

                        st1 = akninvr[ii]*rst1pr[mm]*0.5*ak2;
                        st2 = akninvr[ii]*rst1pim[mm]*0.5*ak2;
                        rst1pr[mm] = -1.0*st2;
                        rst1pim[mm] = st1;

                        st1 = akninvr[ii]*rst2pr[mm]*0.5*ak2;
                        st2 = akninvr[ii]*rst2pim[mm]*0.5*ak2;
                        rst2pr[mm] = -1.0*st2;
                        rst2pim[mm] = st1;

                }
        }

}

void reflec( FILE *file1, int nx, double *rst1pr, double *rst1pim, double *rst2pr, double *rst2pim, double *rfampr, double *rfampim, int iom, int *itrv){

      int nn = 2*nx+1;
      int ka = nn*nn;
      int mm=0, mm1=0, mm2=0, mm3=0;
      int jc=0, jc1=0, i=0, j=0, k=0;
      int ii=0, jj=0, kk=0;
      int ret=0; 

      double sum1=0.0, sum2=0.0;

    //  double* rs1 = malloc(ka*sizeof(double));
    //  double* rs2 = malloc(ka*sizeof(double));
      double* rs3 = malloc(ka*sizeof(double));
      double* rs4 = malloc(ka*sizeof(double));
      double* rs5 = malloc(ka*sizeof(double));
      double* rs6 = malloc(ka*sizeof(double));


  //FILE* file2;
  //char fname2[100];

        //sprintf(fname2,"rstp2rf.dat");
	//file2 = fopen(fname2,"w");

	for (j=0;j<itrv[iom];j++){
                     ii = nx - (itrv[iom]-1)/2+ j;
		for (k=0;k<itrv[iom];k++){
                     mm = j*itrv[iom] + k;
                     jj = nx - (itrv[iom]-1)/2+ k;
                     mm1 = ii*nn+jj; 
		       rs3[mm] = rst2pr[mm1];
		       rs4[mm] = rst2pim[mm1];
		}
	}

        for (i=0;i<itrv[iom];i++){
             mm = i*itrv[iom] + i;
	     rs3[mm] = 1.0 + rs3[mm];
        }

        /*for (i=0;i<itrv[iom];i++){
	        for (j=0;j<itrv[iom];j++){
                   mm = i*itrv[iom] + j;
                   fprintf(file2," %d %d %.12lf \n", i, j, rs3[mm]);
                }
                fprintf(file2," \n");
        }
        fclose(file2);

        sprintf(fname2,"rstp2if.dat");
	file2 = fopen(fname2,"w");

        for (i=0;i<itrv[iom];i++){
	        for (j=0;j<itrv[iom];j++){
                   mm = i*itrv[iom] + j;
                   fprintf(file2," %d %d %.12lf \n", i, j, rs4[mm]);
                }
                fprintf(file2," \n");
        }
        fclose(file2);*/

       jc  = itrv[iom];
       jc1 = jc*jc;


       // code matrix inversion 
       complexmatrixinv( jc, rs3, rs4, rs5, rs6);


       /* sprintf(fname2,"rstp2invr.dat");
	file2 = fopen(fname2,"w");

        for (i=0;i<itrv[iom];i++){
	        for (j=0;j<itrv[iom];j++){
                   mm = i*itrv[iom] + j;
                   fprintf(file2," %d %d %.12lf \n", i, j, rs5[mm]);
                }
                fprintf(file2," \n");
        }
        fclose(file2);

        sprintf(fname2,"rstp2invi.dat");
	file2 = fopen(fname2,"w");

        for (i=0;i<itrv[iom];i++){
	        for (j=0;j<itrv[iom];j++){
                   mm = i*itrv[iom] + j;
                   fprintf(file2," %d %d %.12lf \n", i, j, rs6[mm]);
                }
                fprintf(file2," \n");
        }
        fclose(file2);*/

	for (i=0;i<itrv[iom];i++){
		for (j=0;j<itrv[iom];j++){

                     mm = i*itrv[iom] + j;
                     jj = nx - (itrv[iom]-1)/2+ j;
                     sum1 = 0.0;
                     sum2 = 0.0; 

			for (k=0;k<itrv[iom];k++){

	                     kk = nx - (itrv[iom]-1)/2+ k;
			       mm2 = i*itrv[iom]+k;
			       mm3 = kk*nn+jj;
                    /*      fprintf(stdout," i %d j %d k %d \n", i, j, k);
                          fprintf(stdout," jj %d kk %d mm2 %d mm3 %d \n", jj, kk, mm2, mm3);
                         fprintf(stdout," rst1pr %lf \n", rst1pr[mm3]);
                         fprintf(stdout," rst1pim %lf \n", rst1pim[mm3]);
                         fprintf(stdout," rs5 %lf \n", rs5[mm2]);
                         fprintf(stdout," rs6 %lf \n", rs6[mm2]);*/

			       sum1 = sum1+(rst1pr[mm3]*rs5[mm2]-rst1pim[mm3]*rs6[mm2]);
			       sum2 = sum2+(rst1pr[mm3]*rs6[mm2]+rst1pim[mm3]*rs5[mm2]);
			}
		       rfampr[mm]  = -1.0*sum1;
		       rfampim[mm] = -1.0*sum2;
                }
       } 



}


// code matrix inversion like Kaushik
void complexmatrixinv2( int n, double *matr, double *mati, double *matinvr, double *matinvi){

   int i, j, k;
   int mm1, mm2, mm3, mm4;

   int nmat = n*n;
   double *a1   =  malloc(nmat*sizeof(double));
   double *a2   =  malloc(nmat*sizeof(double));
   double sum=1.0;


  FILE* file2;
  char fname2[100];

   for (i=0;i<nmat;i++) a1[i] = matr[i];
   sum = inverse_matrice( a1, n, a2);

        sprintf(fname2,"matrinv.dat");
	file2 = fopen(fname2,"w");

    for (i=0;i<n;i++){                
	for (j=0;j<n;j++){
             mm1 = i*n+j;
             fprintf(file2," %d %d %.12lf \n", i, j, a2[mm1]);
        }
        fprintf(file2,"\n");
   }
   fclose(file2);

}

// code matrix inversion 
void complexmatrixinv( int n, double *matr, double *mati, double *matinvr, double *matinvi){

   int i, j, k;
   int mm1, mm2, mm3, mm4;

   double sum1, sum2, mataij, matcij;
   int nmat = n*n;

   double *verifr   =  malloc(nmat*sizeof(double));
   double *verifi   =  malloc(nmat*sizeof(double));

   double *bigm      =  malloc(4*nmat*sizeof(double));
   double *bigminv   =  malloc(4*nmat*sizeof(double));

  FILE* file2;
  char fname2[100];

   //     sprintf(fname2,"bigm.dat");
	//file2 = fopen(fname2,"w");

   //fprintf(stdout," invert size n %d \n", n);

    for (i=0;i<n;i++){                
	for (j=0;j<n;j++){

            mm1 = i*2*n+j; 
            mm2 = i*2*n+j+n; 
            mm3 = (i+n)*2*n+j; 
            mm4 = (i+n)*2*n+j+n; 

            mataij = matr[i*n+j];
            matcij = mati[i*n+j];

            bigm[mm1] =  mataij;
            bigm[mm2] =  matcij;
            bigm[mm3] = -matcij;
            bigm[mm4] =  mataij;
	}
     }

    /*for (i=0;i<2*n;i++){                
	for (j=0;j<2*n;j++){
             mm1 = i*2*n+j;
             fprintf(file2," %d %d %.12lf \n", i, j, bigm[mm1]);
        }
        fprintf(file2,"\n");
   }
   fclose(file2);*/

//   matrix_inv( bigm, bigminv, 2*n);
   sum1 = inverse_matrice( bigm, 2*n, bigminv);


       /* sprintf(fname2,"bigminv.dat");
	file2 = fopen(fname2,"w");

    for (i=0;i<2*n;i++){                
	for (j=0;j<2*n;j++){
             mm1 = i*2*n+j;
             fprintf(file2," %d %d %.12lf \n", i, j, bigminv[mm1]);
        }
        fprintf(file2,"\n");
   }
   fclose(file2);*/

    for (i=0;i<n;i++){                
	for (j=0;j<n;j++){
            matinvr[i*n+j] = bigminv[i*2*n+j];
            matinvi[i*n+j] = bigminv[i*2*n+j+n];
        }
    }

   free(bigm);
   free(bigminv);
//   free(matinvr);
//   free(matinvi);
   free(verifr);
   free(verifi);


}


// subroutine to calculate the inverse of a matrix using Gauss-Jordan elimination
// the inverse of matrix a(n,n) is calculated and stored in the matrix b(n,n)
void matrix_inv( double *a, double *b, int n){

    int mmi, mmji, mm, mm1, mm2;
    int i, j, k, irow;
    double big=0.0, dum=0.0;

	for (i=0;i<n;i++){

		for (j=0;j<n;j++){
                     mm = i*n + j;
                     b[mm] = 0.0;
		}

                mm = i*n + i;
                b[mm] = 1.0;

	}

	for (i=0;i<n;i++){

// this is the big loop over all the columns of a(n,n)
// in case the entry a(i,i) is zero, we need to find a good pivot; this pivot
// is chosen as the largest value on the column i from a(j,i) with j = 1,n

                mmi = i*n + i;
                big = a[mmi]; 

		for (j=i;j<n;j++){
                    mmji = j*n + i;
                    if(a[mmji]>big){
                       big = a[mmji];
                       irow = j;
                    }
		}
// interchange lines i with irow for both a() and b() matrices
                if(big>a[mmi]){

                    for (j=0;j<n;j++){

                         mm1 = i*n + j;
                         mm2 = irow*n + j;

                         dum = a[mm1];
                         a[mm1] = a[mm2];
                         a[mm2] = dum;    

                         dum = b[mm1];
                         b[mm1] = b[mm2];
                         b[mm2] = dum;    
                    
                    }

                }

// divide all entries in line i from a(i,j) by the value a(i,i);
// same operation for the identity matrix

                 dum = a[mmi];
                 for (j=0;j<n;j++){
                         mm1 = i*n + j;
                         a[mm1] = a[mm1]/dum;
                         b[mm1] = b[mm1]/dum;
                 }

// make zero all entries in the column a(j,i); same operation for indent()
                 for (j=i+1;j<n;j++){
	                 mmji = j*n + i;
                         dum = a[mmji];
	                 for (k=0;k<n;k++){
        	                 mm1  = i*n + k;
        	                 mm2  = j*n + k;
                             a[mm2] -= dum*a[mm1];
                             b[mm2] -= dum*b[mm1];
                         }
                 }

        }

// substract appropiate multiple of row j from row j-1
	for (i=0;i<n-1;i++){                
		for (j=i+1;j<n;j++){
                       dum = a[i*n+j];
                       for (k=0;k<n;k++){
                            mm1 = i*n+k; 
                            mm2 = j*n+k; 
                            a[mm1] -= dum*a[mm2];
                            b[mm1] -= dum*b[mm2];
                       }
                }
         } 

}



