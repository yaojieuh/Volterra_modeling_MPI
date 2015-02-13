//   Compute the Fourier components for V matrix (at z=iz*dz)
//   similar to B-16   
void computevnn( int jz, int nx, int nz, double czr, double ome, double dx, double *vnnr, double *vnnim, double *vel);
void vsn( int jz, int nx, int nn, double dx, double *xnr, double *xnim, double *vseg, double czr,double* cmin);
void projp( double, int, double, double, double *);
double func_delta( int, int);
void applyproj( int, double *, double *, double *, double *, double *);
void applyprojv( int nx, double *vnnr, double *vnnim, double *projpp);
