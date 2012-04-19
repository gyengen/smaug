
#include "iotypes.h"
#include "iobparams.h"




int cuinit(struct params **p, struct bparams **bp,real **w, real **wnew, real **wd, struct state **state, struct params **d_p, struct bparams **d_bp,real **d_w, real **d_wnew, real **d_wmod, real **d_dwn1, real **d_wd, struct state **d_state, real **d_wtemp, real **d_wtemp1, real **d_wtemp2);
int initgrid(struct params **p, real **w, real **wnew,   struct state **state, real **wd, struct params **d_p, real **d_w, real **d_wnew, real **d_wmod, real **d_dwn1, real **d_wd, struct state **d_state, real **d_wtemp, real **d_wtemp1, real **d_wtemp2);

int cufinish(struct params **p, real **w, real **wnew,   struct state **state, struct params **d_p, struct bparams **d_bp,real **d_w, real **d_wnew, real **d_wmod, real **d_dwn1, real **d_wd, struct state **d_state, real **d_wtemp, real **d_wtemp1, real **d_wtemp2);

int cuboundary(struct params **p, struct bparams **bp,struct params **d_p, struct bparams **d_bp,struct state **d_s,  real **d_wmod,  int order,int idir,int field);
int cuupdate(struct params **p, real **w, real **wmod,real **wtemp2,  struct state **state,struct params **d_p, real **d_w,  real **d_wmod,  real **d_wtemp2,  struct state **d_state,int step);

int cucentdiff1(struct params **p, struct params **d_p, struct state **d_s,real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero, real dt, int field, int dir);
int cucentdiff2(struct params **p, struct params **d_p, struct state **d_s, real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd, int order,int ordero, real dt, int field,int dir);

int cugrav(struct params **p, struct params **d_p, struct state **d_s, real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd, int order,int ordero, real dt);
int cusource(struct params **p, struct params **d_p, struct state **d_s, real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd, int order,int ordero, real dt);

int cuadvance(struct params **p, struct params **d_p,    real **d_wmod, real **d_w,int order);
int cucomputedervfields(struct params **p, struct params **d_p,  real **d_wmod, real **d_wd, int order);
int cucomputevels(struct params **p, struct params **d_p,  real **d_wmod, real **d_wd, int order,int dir);
int cucomputemaxc(struct params **p, struct params **d_p,  real **d_wmod, real **d_wd, int order,int dir, real **wd, real **d_wtemp);
int cucomputemaxcourant(struct params **p, struct params **d_p,  real **d_wmod, real **d_wd, int order,int dir, real **wd, real **d_wtemp);
int cucomputec(struct params **p, struct params **d_p,  real **d_wmod, real **d_wd, int order,int dir);
int cucomputept(struct params **p, struct params **d_p,  real **d_wmod, real **d_wd, int order,int dir);
int cucomputepk(struct params **p, struct params **d_p,  real **d_wmod, real **d_wd, int order,int dir);
int cucomputepbg(struct params **p, struct params **d_p,  real **d_wmod, real **d_wd, int order,int dir);
int cudivb(struct params **p,struct params **d_p, real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd, int order,int ordero, real dt);

int cuhyperdifmomsource1(struct params **p,  struct params **d_p,   real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero, real **d_wtemp, int field, int dim, int ii, int ii0, real dt);

int cuhyperdifmomsourcene1(struct params **p,  struct params **d_p,   real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero, real **d_wtemp, int field, int dim, int ii, int ii0, real dt);

int cuhyperdifesource1(struct params **p,  struct params **d_p,   real **d_wmod, real **d_dwn1, real **d_wd, int order,int ordero, real **d_wtemp, int field, int dim, real dt);

int cuhyperdifbsource1(struct params **p,  struct params **d_p,   real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero, real **d_wtemp, int field, int dim, int jj, int ii0,int mm,real sb,real dt);

int cuhyperdifbsourcene1(struct params **p, struct params **d_p,   real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero, real **d_wtemp, int field, int dim, int jj, int ii0,int mm,real sb,real dt);


int cuhyperdifrhosource1(struct params **p,  struct params **d_p,   real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero, real **d_wtemp, int field, int dim, real dt);



int cunushk1(struct params **p,  struct params **d_p,   real **d_wmod,  real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2);
int cugetdtvisc1(struct params **p,  struct params **d_p,   real **d_wmod, real **wd, real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2);

int cuhyperdifvisc1(struct params **p,  struct params **d_p,   real **d_wmod,  real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2, int field, int dim,int hand);
int cuhyperdifvisc1r(struct params **p,  struct params **d_p,   real **d_wmod,real **wd,  real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2, int field, int dim);
int cuhyperdifvisc1l(struct params **p,  struct params **d_p,   real **d_wmod, real **wd, real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2, int field, int dim);
int cuhyperdifvisc1ir(struct params **p,  struct params **d_p,   real **d_wmod,  real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2, int field, int dim);
int cuhyperdifvisc1il(struct params **p,  struct params **d_p,   real **d_wmod,  real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2, int field, int dim);


int cuhyperdifviscmax(struct params **p, real **w, real **wnew, struct params **d_p, real **d_w, real **d_wnew,  real **d_wmod, real **d_dwn1, real **d_wd, int order, real **d_wtemp, int field, int dim);

#ifdef USE_MPI
 int cuinitmpibuffers(struct params **p,real **w, real **wmod, real **temp2, real **gmpivisc,   real **gmpiw, real **gmpiwmod, struct params **d_p,   real **d_w, real **d_wmod,real **d_wtemp2,    real **d_gmpivisc,   real **d_gmpiw, real **d_gmpiwmod);
 int cufinishmpi(struct params **p,real **w, real **wmod, real **temp2, real **gmpivisc,   real **gmpiw, real **gmpiwmod, struct params **d_p,   real **d_w, real **d_wmod,real **d_wtemp2,    real **d_gmpivisc,   real **d_gmpiw, real **d_gmpiwmod);
int cucopywtompiw(struct params **p,real **w, real **wmod,    real **gmpiw, real **gmpiwmod, struct params **d_p  ,real **d_w, real **d_wmod,   real **d_gmpiw, real **d_gmpiwmod, int order);
int cucopywfrommpiw(struct params **p,real **w, real **wmod,    real **gmpiw, real **gmpiwmod, struct params **d_p  ,real **d_w, real **d_wmod,   real **d_gmpiw, real **d_gmpiwmod, int order);
int cucopytompivisc(struct params **p,real **temp2, real **gmpivisc,  struct params **d_p,real **d_wtemp2,    real **d_gmpivisc);
int cucopyfrommpivisc(struct params **p,real **temp2, real **gmpivisc,  struct params **d_p,real **d_wtemp2,    real **d_gmpivisc);


#endif
