#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "params.h"

/* Read main parameter file */
void params_readfwps_mod( char *file, FILE *file1, double *length_mod, int *dim_mod, double *dh_mod,  int *dimt_mod, double *dt_mod, int *nfr, int *ne, double *coefabc, int *homogeneous, double *c0, char *file_vp){

  fprintf(stdout, "--!\tRead parameter file for 2D FWPS Modeling...\n");

  FILE *file_id;
  char *line, *string, *head;
  int  done = 0, i;

  // open file
  if ( (file_id = fopen(file, "r")) == NULL )
  {                                         
    fprintf(stderr, "\n Problem opening file %s", file);
    exit(EXIT_FAILURE);                                
  }

  string = (char *) malloc(MAX_LINE_LEN);
  line = (char *) malloc(MAX_LINE_LEN);
  head = (char *) malloc(MAX_LINE_LEN);


  fgets(line, MAX_LINE_LEN, file_id);

  
  strcpy(string, "read_fwps");
  while ( fgets(line, MAX_LINE_LEN, file_id) && !done )
  {
    if ( strstr(line, string) != NULL )
    {
      while ( fgets(line, MAX_LINE_LEN, file_id) && !done )
      {
        if ( strstr(line, "physdims") != NULL )
        {
          strcpy(head, strtok(line, "'"));
          strcpy(head, strtok(NULL, "'"));
          length_mod[0] = (double) atof(strtok(head, " "));
          length_mod[1] = (double) atof(strtok(NULL, " "));
          fgets(line, MAX_LINE_LEN, file_id);
          strcpy(head, strtok(line, "'"));
          strcpy(head, strtok(NULL, "'"));
          dim_mod[0] = atoi(strtok(head, " "));
          dim_mod[1] = atoi(strtok(NULL, " "));
          fgets(line, MAX_LINE_LEN, file_id);
          strcpy(head, strtok(line, "'"));
          strcpy(head, strtok(NULL, "'"));
          dh_mod[0] = (double) atof(strtok(head, " "));
          dh_mod[1] = (double) atof(strtok(NULL, " "));
          for (i=0;i<3;i++) fgets(line, MAX_LINE_LEN, file_id);
          strcpy(head, strtok(line, "'"));
          strcpy(head, strtok(NULL, "'"));
          dimt_mod[0] = atoi(strtok(head, " "));
          fgets(line, MAX_LINE_LEN, file_id);
          strcpy(head, strtok(line, "'"));
          strcpy(head, strtok(NULL, "'"));
          dimt_mod[1] = atoi(strtok(head, " "));
          fgets(line, MAX_LINE_LEN, file_id);
          strcpy(head, strtok(line, "'"));
          strcpy(head, strtok(NULL, "'"));
          dt_mod[0] = (double) atof(strtok(head, " "));
          fgets(line, MAX_LINE_LEN, file_id);
          strcpy(head, strtok(line, "'"));
          strcpy(head, strtok(NULL, "'"));
          dt_mod[1] = (double) atof(strtok(head, " "));
          fgets(line, MAX_LINE_LEN, file_id);
          strcpy(head, strtok(line, "'"));
          strcpy(head, strtok(NULL, "'"));
          *nfr = atoi(strtok(head, " "));
          for (i=0;i<3;i++) fgets(line, MAX_LINE_LEN, file_id);
          strcpy(head, strtok(line, "'"));
          strcpy(head, strtok(NULL, "'"));
          *ne = atoi(strtok(head, " "));
          fgets(line, MAX_LINE_LEN, file_id);
          strcpy(head, strtok(line, "'"));
          strcpy(head, strtok(NULL, "'"));
          coefabc[0] = (double) atof(strtok(head, " "));
          fgets(line, MAX_LINE_LEN, file_id);
          strcpy(head, strtok(line, "'"));
          strcpy(head, strtok(NULL, "'"));
          coefabc[1] = (double) atof(strtok(head, " "));
          fgets(line, MAX_LINE_LEN, file_id);
          strcpy(head, strtok(line, "'"));
          strcpy(head, strtok(NULL, "'"));
          coefabc[2] = (double) atof(strtok(head, " "));
          for (i=0;i<3;i++) fgets(line, MAX_LINE_LEN, file_id);
          strcpy(head, strtok(line, "'"));
          strcpy(head, strtok(NULL, "'"));
          *homogeneous = atoi(strtok(head, " "));
          fgets(line, MAX_LINE_LEN, file_id);
          strcpy(head, strtok(line, "'"));
          strcpy(head, strtok(NULL, "'"));
          *c0 = (double) atof(strtok(head, " "));
	  fgets(line, MAX_LINE_LEN, file_id);
          strcpy(head, strtok(line, "'"));
          strcpy(file_vp, strtok(NULL, "'"));
	  done = 1;
	}
      }
      done = 0;
    }
  }
  

  free(line);
  free(head);
  free(string);

  fprintf(file1, "--!\tRead parameter file for 2D FWPS Modeling...\n");
  fprintf(file1, "--!\t\tFull Model length [m]: lz %lf lx %lf \n", length_mod[0], length_mod[1]);
  fprintf(file1, "--!\t\tFull Model dimension : nz %d nx %d \n", dim_mod[0], dim_mod[1]);
  fprintf(file1, "--!\t\tFull Model discretization [m]: dz %lf dx %lf \n", dh_mod[0], dh_mod[1]);
  fprintf(file1, "--!\t\tTime step number %d\n", dimt_mod[0]);
  fprintf(file1, "--!\t\tFrequency number %d. Starting Frequency %d \n", dimt_mod[1], *nfr);
  fprintf(file1, "--!\t\tTime step %lf\n", dt_mod[0]);
  fprintf(file1, "--!\t\tFrequency step %lf\n", dt_mod[1]);
  fprintf(file1, "--!\t\tTaper dimension: ne %d\n", *ne);
  fprintf(file1, "--!\t\tCoefs ABC: aaz %lf  aax %lf alpha %lf \n", coefabc[0], coefabc[1], coefabc[2]);
  if(*homogeneous==1) fprintf(file1, "--!\t\tP-wave velocity %lf\n", *c0);
  if(*homogeneous==0) fprintf(file1, "--!\t\tP-wave velocity file %s\n", file_vp);

}
