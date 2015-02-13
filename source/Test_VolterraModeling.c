//!########################################################################
//!                                                                       #
//! Copyright (C) University of Houston
//!                                                                       #
//! Anne-Cecile Lesage, Nelka Wijesinghe Bandara                          # 
//!                                                                       #
//!
//! Research code for 2D acoustic modeling
//! Full Wave Phase Shift method
//!
//!########################################################################
#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define  MASTER		0

#include "params.h"
#include "comput_volterra.h"
#include "comput_fourierb.h"
#include "comput_reflection.h"

int main( int argc, char *argv[] )
{

   // parameters
   // z-depth discretization 
   int nz;                // z grid point number 
   double dz;              // z spatial step
   double lz;              // z domain length

   // x discretization 
   int nx;                   
   int nn;              // number of grid point in x
   double dx;            // x spatial step
   double lx;            // x domain length

   // time discretization
   int nit;             // number of time iterations
   double at;            // time variable
   double delt;            // time step  
   int nom;             // number of discrete frequencies
   int nfr=11;           // selected sole frequency/ starting frequency
   double delw;          // frequency step
   double aw;            // frequency variable
   double sum1, sum2;

   // source variables
   double fmax=10.0, sig=1.5, gam=8.0, tau=1.0, tt=0.0, arg=0.0;
   //sig *= sig;
   
   // abc variables
   int ne;
    
   
  FILE* file1;
  char fname1[100];
  int ret=0;

  const char * NameWeekDay[] = {"Sunday", "Monday", "Tuesday", "Wenesday", "Thursday", "Friday", "Saturday"};

  time_t timestamp;
  struct tm * tc;
    
  timestamp = time(NULL);
  tc = localtime(&timestamp);

  sprintf(fname1,"test_fwps_modeling.out");
  file1 = fopen(fname1,"w");

  
  fprintf(file1,"--!\t                                     \n");
  fprintf(file1,"--!\tStart Time : %s, ", NameWeekDay[tc->tm_wday]);
  fprintf(file1,"%02u/%02u/%04u, ", tc->tm_mon+1, tc->tm_mday, 1900 + tc->tm_year);
  fprintf(file1,"Clock : %02uh %02umn %02usec.\n", tc->tm_hour, tc->tm_min, tc->tm_sec);

  ret = fflush(file1);

  /*fprintf(stdout, "--!\t                                \n");
  fprintf(stdout, "--!\t Full Wave Phase Shift modeling \n");
  fprintf(stdout, "--!\t                                \n");*/

  fprintf(file1, "--!\t                                 \n");
  fprintf(file1, "--!\t Full Wave Phase Shift modeling  \n");
  fprintf(file1, "--!\t                                 \n");

  ret = fflush(file1);

  // velocity data 
  char file_vp[MAX_LINE_LEN];     // vpz file name
  file_vp[0]   = '\0';
  int homogeneous;  
  double c0;                       // homogeneous velocity

  int dim_mod[2], dimt_mod[2];
  double length_mod[2], dh_mod[2], dt_mod[2], coefabc[3];
  //int ret=0;
  
  // Read main parameter file 
  params_readfwps_mod( argv[1], file1, length_mod, dim_mod, dh_mod, dimt_mod, dt_mod, &nfr, &ne, coefabc, &homogeneous, &c0, file_vp); 

  lz = length_mod[0];
  lx = length_mod[1];

  nz = dim_mod[0];
  nx = dim_mod[1];
  nn = 2*nx + 1; 

  dz = dh_mod[0];
  dx = dh_mod[1];

  nit = dimt_mod[0];
  nom = dimt_mod[1];

  delt = dt_mod[0];
  delw = dt_mod[1];
//  delw = dt_mod[1]/(2.0*M_PI);

  int i, j, mm;
  int nze, nxe;
  nze = nz + 2*ne;
  nxe = nn + 2*ne;
  int ndims = nz*nn, km=ndims;  
  int ka = nn*nn;
  int kc = 4*nn*nn;
  int kd = 2*ka;
  int kb = 2*nn; 

//  int kdisc = ka*500;
  int kdisc = ka*nz;

    int iw, iwi, iwf;

  double* velr    = malloc(nz*nn*sizeof(double));   // velocity from file
  double* vele    = malloc(nze*nxe*sizeof(double));   // extended velocity
  double* pwave   = malloc(nz*nn*sizeof(double));
  double* vpot   = malloc(nz*nn*sizeof(double));

  double* ag      =  malloc(nit*sizeof(double));
  double* afr     =  malloc(nom*sizeof(double));
  double* agomer  =  malloc(nom*sizeof(double));
  double* agomeim =  malloc(nom*sizeof(double));
  int* itrv      =  malloc(nom*sizeof(int));   // how many propagating waves by frequency !!!

  fprintf(file1, "--!\tCompute S(w)...\n");
  sig *= fmax;

  FILE* file2;
  char fname2[100];

  // todo 
  // calculations for the source function
  // equation 34 in paper
  // calculates the source functions S(t) for different times
 sprintf(fname2,"signal_t.dat");
  file2 = fopen(fname2,"w");
  for (i=0;i<nit;i++){
           at = i*delt;
           ag[i]=-sqrt(2.0/M_PI)*sig*gam;
           ag[i]=ag[i]*(sig-2.0*sig*gam*(sig*at-tau)*(sig*at-tau));
           ag[i]=ag[i]*exp(-gam*(sig*at-tau)*(sig*at-tau));
      	   fprintf(file2," %lf %lf  \n", at,  ag[i]);

  }
fclose(file2);

  fprintf(file1, "--!\tCompute Itrv...\n");

//   if lambda_k imaginary
//   counting the number of propagating waves for each frequency
//   eq 23
  int m1=0, j1=0;
  double a1; 
  double al = (nn-1)*dx;
  double czr = 2000.0;
  double xx=0.0, zz=0.0;
//  double czr = 2654.39882;
  czr = c0;
  double v0=0.5, vexs=0.0;
  double temp, temp1, temp2, temp3, temp4;



    iwi = nfr;
    iwf = nfr+nom; 


  sprintf(fname2,"signal_iwi_%d.dat", iwi);
  file2 = fopen(fname2,"w");
  //   loop over the discrete frequencies
  //   aw: frequency value
  for (iw=iwi;iw<iwf;iw++){

            aw   = iw*delw+0.001;
            sum1 = 0.0;
            sum2 = 0.0;

                  sum1 = sum1+0.50*delt*cos(arg)*ag[0];
                  sum2 = sum2-0.50*delt*sin(arg)*ag[0];
          
            for (j=1;j<nit-1;j++){
//            for (j=0;j<nit;j++){

               at  = j*delt;
               arg = aw*at;

                  sum1 = sum1+delt*cos(arg)*ag[j];
                  sum2 = sum2-delt*sin(arg)*ag[j];

  //   approximated values of the fourier transform S^(w)
  //   S^(w) = 0.5*S(t_init)*dt+ sum(i=2,i=nit-1)*S(i*dt)*dt + 0.5*S(t_final)*dt

            }

                  sum1 = sum1+0.5*delt*cos(arg)*ag[nit-1];
                  sum2 = sum2-0.5*delt*sin(arg)*ag[nit-1];

  //  real and imaginary S(w)
            agomer[iw-iwi]  = sum1;
            agomeim[iw-iwi] = sum2;
            afr[iw-iwi]     = aw;

       fprintf(file2," %lf %lf %lf \n", aw, agomer[iw-iwi], agomeim[iw-iwi]);

  }
  fclose(file2);

  sprintf(fname2,"itrv.dat");
  file2 = fopen(fname2,"w");


  for (i=0;i<nom;i++){

        m1 = 0;
  	for (j=0;j<nn;j++){
             j1 = -nx + j;
             a1 = afr[i]/czr;
             a1 = a1*a1;
             a1 = a1 - (2.0*M_PI*j1/al)*(2.0*M_PI*j1/al);
             if(a1>=0.0){
                m1++;
             }
             itrv[i] = m1;
 	} 

        fprintf(file2," %d %d %lf \n", i, itrv[i], afr[i]);  
  }
  fclose(file2);


  //exit(1);

  fprintf(file1, " Initialize velocity \n");

  // read ascii file for velocity
  int sizevel=0;
  FILE *fic; 
  int i1=0, i2=0;
  double val=0.0;
  
  	if(homogeneous==0){
		fic = fopen(file_vp, "r"); // open to read
	  	if(fic == NULL){
	  	    fprintf(stdout,"ERROR : Impossible to open file %s\n", file_vp);
	  	    exit(1);
	  	}
	  	fprintf(file1,".............. Open file %s\n", file_vp);

	  	// count and get time sequence values
	  	i=0;
	  	
	  	while(fscanf(fic, "%f", &val) != EOF){
	  	  velr[i] = val;
		  i++; 
	  	}

		  // fermeture du fichier
	  	if(fclose(fic) == EOF) {
	      		fprintf(file1,"Problem to close file %s", file_vp);
	      		exit(1);
	  	}
	  	fprintf(file1,".............. Close file %s\n", file_vp);

        }else if(homogeneous==1){

               fprintf(file1, "--!\t...HOMOGENEOUS TEST CASE\n");

               for (i=0;i<ndims;i++) velr[i] = c0;

        }else if(homogeneous==2){

               for (i=0;i<nz;i++){
		       for (j=0;j<nn;j++){
		          mm = i*nn+j;
                          val = 0.0;
                          if((i>29)&&(i<80)){
				if((j>100)&&(j<nn-100)){ 
                            	val=0.5;
				}
                          } 
                          vpot[mm] = val; 
                          velr[mm] = c0/(sqrt(1.0-val));
                       }
               }


        // write steep model
        sprintf(fname2,"layermod.dat");
	file2 = fopen(fname2,"w");

	for(j=0;j<nn;j++){

 	  		for(i=0;i<nz;i++){
		          mm = (nz-1-i)*nn+j;
                          xx = j*dx;
                          zz = (nz-1-i)*dz;
                          fprintf(file2," %lf %lf %lf %lf \n", xx, -zz, velr[mm],vpot[mm]);
			}
                          fprintf(file2," \n");

	}

        fclose(file2);


        }else if(homogeneous==3){

             // steep model
          
                for (i=0;i<nz;i++){
		       for (j=0;j<nn;j++){
		          mm = i*nn+j;
                          val = 0.0;
                          i1=40+(j-nx)*15/nx;
                          i2=60+(j-nx)*15/nx;
                            if((i<i2)&&(i>i1)) val = 0.5;
                            vpot[mm] = val;
                            velr[mm] = c0/(sqrt(1.0-val));
		       }
                }


        // write steep model
        sprintf(fname2,"steepmod.dat");
	file2 = fopen(fname2,"w");

	for(j=0;j<nn;j++){

 	  		for(i=0;i<nz;i++){
		          mm = (nz-1-i)*nn+j;
                          xx = j*dx;
                          zz = (nz-1-i)*dz;
                          fprintf(file2," %lf %lf %lf %lf \n", xx, -zz, velr[mm], vpot[mm]);
			}
                          fprintf(file2," \n");

	}

        fclose(file2);


        }else if(homogeneous==4){


               for (i=0;i<nz;i++){

                       zz = i*7.0/(nz-1);
          temp = 2.0*zz-6.0;
//          temp = 2.0*z-M_PI;
          temp1 = temp*temp;
          temp2 = temp1*temp1;
          temp3 = temp2*temp1;
          temp4 = temp2*temp2; 

          vexs = 1;
          vexs += temp1;
          vexs += temp2/2.0;
          vexs += temp3/6.0;
          vexs += temp4/24.0;
          vexs = v0*exp(-1.0*temp1)*vexs;

		       for (j=0;j<nn;j++){
		          mm = i*nn+j;
                          velr[mm] = c0/sqrt(1.0-vexs);
                          vpot[mm] = vexs;
                       }
               } 


        // write steep model
        sprintf(fname2,"smoothlayer.dat");
	file2 = fopen(fname2,"w");

	for(j=0;j<nn;j++){

 	  		for(i=0;i<nz;i++){
		          mm = (nz-1-i)*nn+j;
                          xx = j*dx;
                          zz = (nz-1-i)*dz;
                          fprintf(file2," %lf %lf %lf %lf \n", xx, -zz, velr[mm], vpot[mm]);
			}
                          fprintf(file2," \n");

	}

        fclose(file2);

   }else if(homogeneous==5){

                for (i=0;i<nz;i++){
		       for (j=0;j<nn;j++){
		          mm = i*nn+j;
                          val = 0.0;
                          i1=29;
                          i2=80; 

                       if((i<i2)&&(i>i1)){
				val =0;
		               if((j<nn-100)&&(j>100)){
				 val =0.9;
		               }
                       }
                            vpot[mm] = val;
                            velr[mm] = c0/(sqrt(1.0-val));
		       }
                }

        // write steep model
        sprintf(fname2,"rectangle.dat");
	file2 = fopen(fname2,"w");

	for(j=0;j<nn;j++){

 	  		for(i=0;i<nz;i++){
		          mm = (nz-1-i)*nn+j;
                          xx = j*dx;
                          zz = (nz-1-i)*dz;
                          fprintf(file2," %lf %lf %lf %lf \n", xx, -zz, velr[mm], vpot[mm]);
			}
                          fprintf(file2," \n");

	}

        fclose(file2);


   }

 // exit(1);


  ret=fflush(file1);

    double ome;
    double zzer = 0.0; 
    double sum3=0.0, xj=0;//(nn-1)*dx/2;
    int iom1, ic, ic2;
    int iz, k, mm3, mm1, jj1, ic1; 
        iw=4;

    double* p1zr   = malloc(ka*sizeof(double));
    double* p1zim  = malloc(ka*sizeof(double));
    double* p2zr   = malloc(ka*sizeof(double));
    double* p2zim  = malloc(ka*sizeof(double));

    double* plsr   = malloc(ka*sizeof(double));
    double* plsim  = malloc(ka*sizeof(double));
    double* plsr1   = malloc(ka*sizeof(double));
    double* plsim1  = malloc(ka*sizeof(double));

    double* plsr2   = malloc(km*sizeof(double));   // Re P(x,z,w)
    double* plsim2  = malloc(km*sizeof(double));   // Im P(x,z,w)

    // reflection amplitude r(n',n,w)
    double* rwgr   = malloc(ka*sizeof(double));
    double* rwgim  = malloc(ka*sizeof(double));
    double* rfampr   = malloc(ka*sizeof(double));
    double* rfampim  = malloc(ka*sizeof(double));
    // reflection amplitude for inversion r(n',0,w)
    double* rinvr  = malloc(nom*nn*sizeof(double));
    double* rinvim = malloc(nom*nn*sizeof(double));
 

    // fourier basis transform of ome^2(1/c0^2-a/c^2) 
    double* vnnmr   = malloc(ka*sizeof(double));
    double* vnnmim   = malloc(ka*sizeof(double));
  
    //sum dz*exp(+)*v(z)*p1(z)
    double* rst1pr   = malloc(ka*sizeof(double));
    double* rst1pim   = malloc(ka*sizeof(double));
    double* rst1mr   = malloc(ka*sizeof(double));
    double* rst1mim   = malloc(ka*sizeof(double));
    //sum dz*exp(+)*v(z)*p2(z)
    double* rst2pr   = malloc(ka*sizeof(double));
    double* rst2pim   = malloc(ka*sizeof(double));
    double* rst2mr   = malloc(ka*sizeof(double));
    double* rst2mim   = malloc(ka*sizeof(double));

    double* projpp   = malloc(nn*sizeof(double));
    double* akninvr   = malloc(nn*sizeof(double));
    double* akninvim  = malloc(nn*sizeof(double));
    double* aknr   = malloc(nn*sizeof(double));
    double* aknim  = malloc(nn*sizeof(double));

    // saving P1(n',z,n,w) to file
    double* p1zrdisc   = malloc(kdisc*sizeof(double)); // Re P1(n',z,n,w)
    double* p1zimdisc  = malloc(kdisc*sizeof(double)); // Im P1(n',z,n,w)
    double* p2zrdisc   = malloc(kdisc*sizeof(double)); // Re P2(n',z,n,w)
    double* p2zimdisc  = malloc(kdisc*sizeof(double)); // Im P2(n',z,n,w)

   //saving P(x,0,w) to file
    double* pxz0r   = malloc(nn*nom*sizeof(double));   // Re P(x,z,w)
    double* pxz0i  = malloc(nn*nom*sizeof(double));   // Im P(x,z,w)

    double ak2=0.0;
    double c1=4000.0;

//
//   pq1init, pq2init - calculates the Volterra initial conditions
//   Eq A-25 & 26 in paper
    FILE* file3;  
    char fname3[100];

    int fresize, frei, fref, nbfreid;
    int numtasks, taskid;
    MPI_Init(&argc, &argv);
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);  /*total processors used */
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);     /* the number of current processor */
     sprintf(fname3,"proc_%d.dat", taskid);
     file3 = fopen(fname3,"w");
      fresize = (nom/numtasks)+1; /* total number of source for each processor */
     frei = taskid*fresize;  /*starting source number for each processor */
     fref = (taskid+1)*fresize;     /*ending source number for each processor */
     if(frei>nom) frei = nom;
     if(fref>nom) fref = nom;
     nbfreid = fref-frei;
     fprintf (file3, "MPI task %d has started... frei %d fref %d \n", taskid, frei, fref);
        ret = fflush(file3);
   
  for(iw=frei;iw<fref;iw++){

     	    ome = afr[iw];  
            iom1 = iw-frei;
            ak2 = ome*ome/(czr*czr);
            fprintf(file3," ome %lf ak2 %lf ", ome, ak2);


	    for(i=0;i<nn;i++){
			for(j=0;j<nn;j++){
		          mm = i*nn+j;
		          rwgr[mm] = 0.0;
		          rwgim[mm] = 0.0;

			}
	   }

            for (i=0;i<ka;i++){

		          rst1pr[i]=0.0;
		          rst1pim[i]=0.0;
		          rst1mr[i]=0.0;
		          rst1mim[i]=0.0;
		          rst2pr[i]=0.0;
		          rst2pim[i]=0.0;
		          rst2mr[i]=0.0;
		          rst2mim[i]=0.0;

            }


     	    fprintf(file3, " Compute pq1init aw %lf \n", ome);
     	    p1init( nx, nn, dx, ome, czr, zzer, p1zr, p1zim);


     	    iz = 0;
     	    projp( ome, nx, dx, czr, projpp);

     	    aknmat( nx, ome, czr, dx, aknr, aknim, akninvr, akninvim);
     
     	for(i=0;i<ka;i++){

        	p2zr[i]  = p1zr[i];
         	p2zim[i] = -1.0*p1zim[i];

         	p1zrdisc[i+iz*ka]  = p1zr[i];
         	p1zimdisc[i+iz*ka] = p1zim[i];
         	p2zrdisc[i+iz*ka]  = p2zr[i];
         	p2zimdisc[i+iz*ka] = p2zim[i];

     	}

       
     	for(iz=1;iz<nz;iz++){
		 computevnn( iz-1, nx, nz, czr, ome, dx, vnnmr, vnnmim, velr);
     		compute_volterrap1( p1zr, p1zim, nx, iz, dz, ak2, iw, itrv, vnnmr, vnnmim,
                  aknr, akninvr, rst1pr, rst1pim, rst1mr, rst1mim);
     		compute_volterrap2( p2zr, p2zim, nx, iz, dz, ak2, iw, itrv, vnnmr, vnnmim, 
                  aknr, akninvr, rst2pr, rst2pim, rst2mr, rst2mim);


       	        applyproj( nx, p1zr, p1zim, p2zr, p2zim, projpp);

	   	if((iz%50)==0) {
	    	fprintf(file3, " Compute P1 P2 iz %d \n", iz);
	    	ret=fflush(file3);
		if ((iw%100)==0){
		sprintf(fname2,"vnn_%d_fre_%d.dat", iz,iw);
		file2 = fopen(fname2,"w");
		for(j=0;j<nn;j++){
 	  		for(i=0;i<nn;i++){
		           mm = i*nn+j;
                          
                          fprintf(file2," %d %d %lf %lf\n", i,j, vnnmr[mm],vnnmim[mm]);
			}
                          fprintf(file2," \n");
		}

        	fclose(file2);
		}

		}
     		 for(i=0;i<ka;i++){

             

		 p1zrdisc[i+iz*ka]  = p1zr[i];
		 p1zimdisc[i+iz*ka] = p1zim[i];
		 p2zrdisc[i+iz*ka]  = p2zr[i];
		 p2zimdisc[i+iz*ka] = p2zim[i];

	     	}

    	 }
	
	    fprintf(file3, " Compute reflection \n");
	    ret=fflush(file3);

         mult( nx, iw, ak2, rst1pr, rst1pim, rst2pr, rst2pim, itrv, akninvr);

     	 reflec( file1, nx, rst1pr, rst1pim, rst2pr, rst2pim, rfampr, rfampim, iw, itrv);

	 for(i=0;i<itrv[iw];i++){
                  ic=nx-(itrv[iw]-1)/2+i;
	          for(j=0;j<itrv[iw];j++){
                          ic2=nx-(itrv[iw]-1)/2+j;
		          mm  = i*itrv[iw]+j;
		          mm1 = ic*nn+ic2;
		          rwgr[mm1]  = rfampr[mm];
		          rwgim[mm1] = rfampim[mm];

		}
	}

     	 for(j=0;j<nn;j++){
             rinvr[j+iw*nn] = rwgr[j*nn+nx];
             rinvim[j+iw*nn] = rwgim[j*nn+nx];
        } 
      	for(iz=0;iz<nz;iz++){

    		// Read P1(n',z,n,w) P2(n',z,n,w)
		for(j=0;j<ka;j++){
             		 p1zr[j]  =  p1zrdisc[j+iz*ka];
              		 p1zim[j] = p1zimdisc[j+iz*ka];
              		 p2zr[j]  =  p2zrdisc[j+iz*ka];
              		 p2zim[j] = p2zimdisc[j+iz*ka];
       		}

	  	 if((iz%50)==0) {
	    		fprintf(file3, " Compute P+(n',z,n,w) iz %d \n", iz);
	   		 ret=fflush(file3);
	   	}


     		 // P+(n',z,n,w) = P1(n',z,n,w) + r(n'',n,w) P2(n',z,n'',w)
      		// P2 = P*1
		for(i=0;i<nn;i++){
		     for(j=0;j<nn;j++){
		          mm = i*nn + j;
                          sum1 = 0.0;
                          sum2 = 0.0;
		   	  for(k=0;k<nn;k++){
		                  mm3 = i*nn + k;                                
		                  mm1 = k*nn + j;                  
		                  sum1=sum1+(p2zr[mm3]*rwgr[mm1] - p2zim[mm3]*rwgim[mm1]);
		                  sum2=sum2+(p2zim[mm3]*rwgr[mm1] + p2zr[mm3]*rwgim[mm1]);
		          }
			  plsr[mm] = p1zr[mm] + sum1;
                          plsim[mm] = p1zim[mm] + sum2;
			// plsr[mm] = p1zr[mm];
                        //  plsim[mm] =p1zim[mm];
                     }
                } 
		
	

  //		fprintf(file1, " Compute P+(x,z,n,w) \n");


		// A-10 P+(n',z,n,w) -> P+(x,z,n,w)
		for(k=0;k<nn;k++){
               	 	xx = (k)*dx; 
			for(i=0;i<nn;i++){
                     		mm = k*nn + i;
                     		sum1 = 0.0;
                     		sum2 = 0.0; 
		     		for(j=0;j<nn;j++){
                         		jj1 = -nx + j;
                         		mm3 = j*nn + i;
                         		arg = 2.0*M_PI*jj1*xx/al;
                         		sum1 += cos(arg)*plsr[mm3]  - sin(arg)*plsim[mm3];   
                         		sum2 += cos(arg)*plsim[mm3] + sin(arg)*plsr[mm3];   
                    		 }
                     		plsr1[mm]   = sum1/sqrt(al);
                     		plsim1[mm] = sum2/sqrt(al);
                	} 
        	}
		for(k=0;k<nn;k++){
                	xx = (-nx+k)*dx;
                	mm = iz*nn+k; 
                 	sum1=0.0;
                 	sum2=0.0;
                 	sum3=0.0;
                 	for(j=0;j<itrv[iw];j++){
                         	ic1 = nx - (itrv[iw]-1)/2 + j;
                         	jj1 = j - (itrv[iw]-1)/2;
                         	mm1 = k*nn + ic1;
                         	arg = 2.0*M_PI*jj1*(xj+nx*dx)/al;
                         	a1 = afr[iw]/czr;
                         	a1 *= a1;
                         	a1 = a1 - (2.0*M_PI*jj1*2.0*M_PI*jj1)/(al*al);
                         	a1 = sqrt(fabs(a1));  
                         	a1 = 1.0/(sqrt(al)*a1*2.0);
        
                         	sum1 += (cos(arg)*plsim1[mm1]-sin(arg)*plsr1[mm1])*a1;
                         	sum2 -= (cos(arg)*plsr1[mm1]+sin(arg)*plsim1[mm1])*a1;
                	 }


              		plsr2[mm]=sum1* agomer[iw]-sum2*agomeim[iw];
              		plsim2[mm]=sum1* agomeim[iw]+sum2*agomer[iw];
			

       	 	}
     		

      } 
	 for(j=0;j<nn;j++){
		mm=j;
		pxz0r[iw*nn+j]=plsr2[2*nn+mm];
		pxz0i[iw*nn+j]=plsim2[2*nn+mm];
	}
	  if((iw%100)==0){
		sprintf(fname2,"reflecr_iw_%d.dat", iw);
		file2 = fopen(fname2,"w");

		for(i=0;i<nn;i++){
			for(j=0;j<nn;j++){
		          mm = i*nn+j;
                          fprintf(file2," %d %d %.12lf \n", i, j, rwgr[mm]);
			}
                          fprintf(file2," \n");
		}

       		 fclose(file2);

        	sprintf(fname2,"refleci_iw_%d.dat", iw);
		file2 = fopen(fname2,"w");

		for(i=0;i<nn;i++){
			for(j=0;j<nn;j++){
		          mm = i*nn+j;
                          fprintf(file2," %d %d %.12lf \n", i, j, rwgim[mm]);
			}
                          fprintf(file2," \n");
		}

        	fclose(file2);
        	fprintf(file1, " Print pressure result iw %d aw %lf \n", iw, afr[iw]);
        	sprintf(fname2,"pxzr_iw_%d.dat", iw);
		file2 = fopen(fname2,"w");

		for(j=0;j<nn;j++){
 	  		for(i=0;i<nz;i++){
		          mm = (nz-1-i)*nn+j;
                          xx = j*dx;
                          zz = (nz-1-i)*dz;
                          fprintf(file2," %lf %lf %lf \n", xx, -zz, plsr2[mm]);
			}
                          fprintf(file2," \n");
		}

        	fclose(file2);


        	sprintf(fname2,"pxzim_iw_%d.dat", iw);
		file2 = fopen(fname2,"w");

		for(i=0;i<nz;i++){
			for(j=0;j<nn;j++){
		          mm = i*nn+j;
                          xx = j*dx;
                          zz = i*dz;
                          fprintf(file2," %lf %lf %lf \n", zz, xx, plsim2[mm]);
			}
                          fprintf(file2," \n");
		}

        	fclose(file2);

        	sprintf(fname2,"pxzr_iw_%d.bin", iw);
		file2 = fopen(fname2,"a");
    		 if( fwrite( plsr2, sizeof(double), ndims, file2 ) != (size_t) ndims)
     		{
			printf(" Can't output data to file \n");
			exit(1);
     		}

     		fclose(file2);


        	sprintf(fname2,"pxzim_iw_%d.bin", iw);
		file2 = fopen(fname2,"a");
     		if( fwrite( plsr2, sizeof(double), ndims, file2 ) != (size_t) ndims)
     		{
			printf(" Can't output data to file \n");
			exit(1);
     		}

    		 fclose(file2);

      	}

	 

    }

     int offset, lsize;
     int tag1=2, tag2=1, tag3=3, tag4=4;
     int dest;

     	offset = frei*nn;
     	lsize =(fref-frei)*nn;
	 if (taskid > MASTER) {
	 
	  dest = MASTER;
      	  ret=fflush(file3);
	  MPI_Send(&pxz0r[offset], lsize, MPI_DOUBLE, MASTER, tag1, MPI_COMM_WORLD);
	  MPI_Send(&pxz0i[offset], lsize, MPI_DOUBLE, MASTER, tag2, MPI_COMM_WORLD);
	  MPI_Send(&rinvr[offset], lsize, MPI_DOUBLE, MASTER, tag3, MPI_COMM_WORLD);
	  MPI_Send(&rinvim[offset], lsize, MPI_DOUBLE, MASTER, tag4, MPI_COMM_WORLD);
     	}
    fprintf(file3," after sending \n");
     ret = fflush(file3);

  	
     
     if (taskid == MASTER){

	  
	  for (i=1; i<numtasks; i++) {

	     frei = i*fresize;
	     fref = (i+1)*fresize;     
	     if(frei>nom) frei = nom;
	     if(fref>nom) fref = nom;
    		 offset = frei*nn;
     		lsize =(fref-frei)*nn;


        fprintf(file3," master receive slave %d \n", i);
        ret = fflush(file3);
	MPI_Recv(&pxz0r[offset], lsize, MPI_DOUBLE, i, tag1,MPI_COMM_WORLD, &status);
        MPI_Recv(&pxz0i[offset], lsize, MPI_DOUBLE, i, tag2,MPI_COMM_WORLD, &status);
	MPI_Recv(&rinvr[offset], lsize, MPI_DOUBLE, i, tag3,MPI_COMM_WORLD, &status);
        MPI_Recv(&rinvim[offset], lsize, MPI_DOUBLE, i, tag4,MPI_COMM_WORLD, &status);		

         }
     }
    if(taskid==MASTER){
    double w;
    sprintf(fname2,"reflection_normal.dat");
    file2 = fopen(fname2,"w");
    for(i=0;i<nom;i++){
		for(j=0;j<nn;j++){
		          mm = i*nn+j;
                          w = i*delw;
                          fprintf(file2," %d %lf %lf %lf\n", i, w, rinvr[i*nn+j],rinvim[i*nn+j]);
		}
                          fprintf(file2," \n");
	}

     fclose(file2);
    sprintf(fname2,"pxz0r.dat");
    file2 = fopen(fname2,"w");
    for(i=0;i<nom;i++){
		for(j=0;j<nn;j++){
		          mm = i*nn+j;
                          xx = j*dx;
                          w = i*delw;
                          fprintf(file2," %lf %lf %lf \n", xx, w, pxz0r[i*nn+j]);
		}
                          fprintf(file2," \n");
	}

     fclose(file2);
     sprintf(fname2,"pxz0i.dat");
     file2 = fopen(fname2,"w");
     for(i=0;i<nom;i++){
		for(j=0;j<nn;j++){
		          mm = i*nn+j;
                          xx = j*dx;
                          w = i*delw;
                          fprintf(file2," %lf %lf %lf \n", xx, w, pxz0i[i*nn+j]);
		}
                          fprintf(file2," \n");
	}

     fclose(file2);
  double* pxz0t   = malloc(nn*nit*sizeof(double));  
  double pi=3.1415926;
  sprintf(fname2,"pxz0t.dat");
    
 file2 = fopen(fname2,"w");
     for (j=0;j<nn;j++){
	xx=j*dx;
	for (i=0;i<nit;i++){
		at=i*delt-0.2;
		pxz0t[j*nit+i]=0;
		for (k=0;k<nom;k++){
			arg=-(afr[k])*at;
			pxz0t[j*nit+i]+=cos(arg)*pxz0r[k*nn+j]-sin(arg)*pxz0i[k*nn+j];
		}
		pxz0t[j*nit+i]*=delw*2/2/pi;
		fprintf(file2," %f %f %f \n", xx, at,pxz0t[j*nit+i]);
	}
	fprintf(file2,"  \n");
    }
     fclose(file2);
		
  } 
   MPI_Finalize();


	
}
