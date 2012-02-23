
char **hlines; //header lines for vac config files 
hlines=(char **)calloc(5, sizeof(char*));
// Define time-domain
real dt;


real *d_w;
real *d_wnew;

real *d_wmod,  *d_dwn1,  *d_dwn2,  *d_dwn3,  *d_dwn4,  *d_wd;

real *w,*wnew,*wd, *temp2,*wmod;
real *d_wtemp,*d_wtemp1,*d_wtemp2;

#ifdef USE_MPI
  real *gmpivisc, *gmpiw, *gmpiwmod;
  real *d_gmpivisc, *d_gmpiw, *d_gmpiwmod;
#endif


#ifdef USE_MPI
//buffers to use on GPU
  real *d_gmpisendbuffer;
  real *d_gmpirecvbuffer;

   
  real *d_gmpisrcbufferl;
  real *d_gmpisrcbufferr;
  real *d_gmpitgtbufferl;
  real *d_gmpitgtbufferr;
#endif

