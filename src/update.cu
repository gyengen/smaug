#include "../include/cudapars.h"
#include "../include/paramssteeringtest1.h"

/////////////////////////////////////
// standard imports
/////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include "../include/step.h"

/////////////////////////////////////
// kernel function (CUDA device)
/////////////////////////////////////
#include "../include/gradops_u.cuh"


__device__ __host__
int updatestate (struct params *p, struct state *s, real *w ,int *ii, int field) {

  int status=0;
                      // atomicExch(&(p->cmax),(wd[fencode3_pre(p,ii,soundspeed)]));
                    switch(field)
                    {
                      case rho:
                    	s->rho=s->rho+(w[fencode3_u(p,ii,field)]);
		      break;
                      case mom1:
                    	s->m1=s->m1+(w[fencode3_u(p,ii,field)]);
		      break;
                      case mom2:
                    	s->m2=s->m2+(w[fencode3_u(p,ii,field)]);
		      break;
                      /*case mom3:
                    	s->m3=s->m3+(w[fencode3_u(p,ii,field)]);
		      break;*/
                      case energy:
                    	s->e=s->e+(w[fencode3_u(p,ii,field)]);
		      break;
                      case b1:
                    	s->b1=s->b1+(w[fencode3_u(p,ii,field)]);
		      break;
                      case b2:
                    	s->b2=s->b2+(w[fencode3_u(p,ii,field)]);
		      break;
                      /*case b3:
                    	s->b3=s->b3+(w[fencode3_u(p,ii,field)]);
		      break;*/
                    };
  return status;
}



__global__ void update_parallel(struct params *p, struct state *s, real *w, real *wmod)
{

   int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j,f;
  int index,k;
  __shared__ int ntot;

  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  //real g=p->g;
  real *u,  *v,  *h;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

   int ip,jp,ipg,jpg;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   //int ip,jp,ipg,jpg;
  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     



//int shift=order*NVAR*dimp;

  h=w+dimp*rho;
  u=w+dimp*mom1;
  v=w+dimp*mom2;


     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];
           for( f=rho; f<=b3; f++)
     #else
           for( f=rho; f<=b2; f++)
     #endif
             {  
         #ifdef USE_SAC_3D
      if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
     #else
       if(i<((p->n[0])) && j<((p->n[1])))
     #endif           
	{
            
                  w[fencode3_u(p,iia,f)]=wmod[fencode3_u(p,iia,f)];

	}


}

__syncthreads(); 







  
}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_u(char *label)
{
  // we need to synchronise first to catch errors due to
  // asynchroneous operations that would otherwise
  // potentially go unnoticed

  cudaError_t err;

  err = cudaThreadSynchronize();
  if (err != cudaSuccess)
  {
    char *e = (char*) cudaGetErrorString(err);
    fprintf(stderr, "CUDA Error: %s (at %s)", e, label);
  }

  err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    char *e = (char*) cudaGetErrorString(err);
    fprintf(stderr, "CUDA Error: %s (at %s)", e, label);
  }
}
int cuupdate(struct params **p, real **w, real **wmod,real **wtemp2, struct state **state,struct params **d_p, real **d_w, real **d_wmod, real ** d_wtemp2, struct state **d_state, int step)
//int cuupdate(struct params **p, real **w, real **wmod, real **wd, real **temp2, struct state **state,
//             struct params **d_p, real **d_w, real **d_wmod, real **d_wtemp2, struct state **d_state, int step)
{
  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 
    dim3 dimBlock(dimblock, 1);
 
    dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;
  // cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyHostToDevice);
cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);
     update_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_state,*d_w,*d_wmod);
	    //printf("called update\n"); 
    cudaThreadSynchronize();
//following comments removed from if def pragmas  if
//using MPI and copying all cell data to host (how slow!?)
//#ifdef USE_MPI

//#else
    if((step%((*p)->cfgsavefrequency))==0)
//#endif
    {

//following commentes removed from section if
//using MPI and copying all cell data to host (how slow!?)
/*#ifdef USE_MPI
    cudaMemcpy(*wmod, *d_w, NVAR*dimp*sizeof(real), cudaMemcpyDeviceToHost);
    #ifdef USE_SAC_3D  
           cudaMemcpy(*wtemp2, *d_wtemp2,NTEMP2*(((*p)->n[0])+2)* (((*p)->n[1])+2)* (((*p)->n[2])+2)*sizeof(real), cudaMemcpyDeviceToHost);
    #else
       cudaMemcpy(*wtemp2, *d_wtemp2,NTEMP2*(((*p)->n[0])+2)* (((*p)->n[1])+2)*sizeof(real), cudaMemcpyDeviceToHost);
    #endif

#endif */ 
    cudaMemcpy(*w, *d_w, NVAR*dimp*sizeof(real), cudaMemcpyDeviceToHost);

    //cudaMemcpy(*wnew, *d_wd, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);

   cudaMemcpy(*state, *d_state, sizeof(struct state), cudaMemcpyDeviceToHost);
    }

//cudaMemcpy(*wnew, *d_wnew, 8*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_u, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

  //checkErrors("copy data from device");


 


}


int cufinish(struct params **p, real **w, real **wnew, struct state **state, struct params **d_p,struct bparams **d_bp, real **d_w, real **d_wnew, real **d_wmod, real **d_dwn1, real **d_wd, struct state **d_state, real **d_wtemp, real **d_wtemp1, real **d_wtemp2)
{
  

 //cudaMemcpy(*w, *d_w, 8*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*wnew, *d_wnew, 8*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_u, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

  checkErrors_u("copy data from device");


  cudaFree(*d_p);
  cudaFree(*d_bp);
//  cudaFree(*d_state);

  cudaFree(*d_w);
  cudaFree(*d_wnew);
 // cudaFree(*d_u);

  cudaFree(*d_wmod);
  cudaFree(*d_dwn1);
  cudaFree(*d_wd);
  cudaFree(*d_wtemp);
  cudaFree(*d_wtemp1);
  cudaFree(*d_wtemp2);
  




}

  #ifdef USE_MPI

int cufinishmpi(struct params **p,real **w, real **wmod, real **temp2, real **gmpivisc,   real **gmpiw, real **gmpiwmod, struct params **d_p,   real **d_w, real **d_wmod,real **d_wtemp2,    real **d_gmpivisc,   real **d_gmpiw, real **d_gmpiwmod)
{
  

 //cudaMemcpy(*w, *d_w, 8*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*wnew, *d_wnew, 8*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_u, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

  checkErrors_u("copy data from device");


  cudaFree(*d_gmpiw);
  cudaFree(*d_gmpiwmod);
  cudaFree(*d_gmpivisc);

  free(*gmpiw);
  free(*gmpiwmod);
  free(*gmpivisc);
  free(*temp2);
}
#endif
