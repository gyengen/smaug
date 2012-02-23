#include "../include/cudapars.h"
#include "../include/paramssteeringtest1.h"
#include "../include/iobparams.h"
/////////////////////////////////////
// standard imports
/////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include "../include/step.h"

/////////////////////////////////////
// kernel function (CUDA device)
/////////////////////////////////////
#include "../include/gradops_cd1.cuh"
#include "../include/dervfields_cd1.cuh"
#include "../include/usersource_cd1.cuh"

__device__ __host__
int divflux1(real *dw, real *wd, real *w, struct params *p,int *ii,int field,int dir) {

  int direction;
  int status=0;
  real divflux=0;

dw[fencode3_cd1(p,ii,field)]+= grad3dn_cd1(wd,wd,p,ii,flux,dir); 
//dw[fencode3_cd1(p,ii,field)]+= grad3d_cd1(wd,p,ii,flux,dir); 
//dw[fencode3_cd1(p,ii,field)]=0.0;
 /*switch(field)
  {

     case mom1:
       dw[fencode3_cd1(p,ii,field)]-= (p->g[dir])*w[fencode3_cd1(p,ii,rho)];
      break;
    case mom2:
      dw[fencode3_cd1(p,ii,field)]-= (p->g[dir])*w[fencode3_cd1(p,ii,rho)];
      break;
#if defined USE_SAC_3D
    case mom3:
      dw[fencode3_cd1(p,ii,field)]-= (p->g[dir])*w[fencode3_cd1(p,ii,rho)];
      break;
#endif
    case rho:
     ;// dw[fencode3_cd1(p,ii,field)]+= ix/800;
      break;

  } */   
 // dw[fencode3_cd1(p,ii,field)]= gradd0_cd1(wd,p,ii,f1,0)+gradd1_cd1(wd,p,ii,f2,1);    
  return ( status);
}






__device__ __host__
real transportflux (real *dw, real *wd, real *w, struct params *p,int *ii,int field, int direction) {

  

   //real fluxt=0;

   //transport flux
   //this will work without the switch as follows
        #if defined USE_SAC || defined USE_SAC_3D
     return(w[fencode3_cd1(p,ii,mom1+direction)]*w[fencode3_cd1(p,ii,field)]/(w[fencode3_cd1(p,ii,rho)]+w[fencode3_cd1(p,ii,rhob)]));

        #else
     return(w[fencode3_cd1(p,ii,mom1+direction)]*w[fencode3_cd1(p,ii,field)]/w[fencode3_cd1(p,ii,rho)]);

        #endif


  
 
}






__device__ __host__
real fluxmom1 (real *dw, real *wd, real *w, struct params *p,int *ii,int field, int direction) {


         #if defined USE_SAC || defined USE_SAC_3D
     		return( -(w[fencode3_cd1(p,ii,field+(NDIM+1))]*w[fencode3_cd1(p,ii,b1b+direction)]+w[fencode3_cd1(p,ii,field+(2*NDIM+3))]*w[fencode3_cd1(p,ii,b1+direction)])-w[fencode3_cd1(p,ii,field+(NDIM+1))]*w[fencode3_cd1(p,ii,b1+direction)]);
        #endif


}




__device__ __host__
real fluxmom10 (real *dw, real *wd, real *w, struct params *p,int *ii,int field, int direction) {


         #if defined USE_SAC || defined USE_SAC_3D
         return(direction==0?wd[fencode3_cd1(p,ii,pressuret)]-(w[fencode3_cd1(p,ii,field+(NDIM+1))]*w[fencode3_cd1(p,ii,b1b+direction)]+w[fencode3_cd1(p,ii,field+(2*NDIM+3))]*w[fencode3_cd1(p,ii,b1+direction)])-w[fencode3_cd1(p,ii,field+(NDIM+1))]*w[fencode3_cd1(p,ii,b1+direction)]:-(w[fencode3_cd1(p,ii,field+(NDIM+1))]*w[fencode3_cd1(p,ii,b1b+direction)]+w[fencode3_cd1(p,ii,field+(2*NDIM+3))]*w[fencode3_cd1(p,ii,b1+direction)])-w[fencode3_cd1(p,ii,field+(NDIM+1))]*w[fencode3_cd1(p,ii,b1+direction)]);
        #endif


}

__device__ __host__
real fluxmom11 (real *dw, real *wd, real *w, struct params *p,int *ii,int field, int direction) {


         #if defined USE_SAC || defined USE_SAC_3D
         return(direction==1?wd[fencode3_cd1(p,ii,pressuret)]-(w[fencode3_cd1(p,ii,field+(NDIM+1))]*w[fencode3_cd1(p,ii,b1b+direction)]+w[fencode3_cd1(p,ii,field+(2*NDIM+3))]*w[fencode3_cd1(p,ii,b1+direction)])-w[fencode3_cd1(p,ii,field+(NDIM+1))]*w[fencode3_cd1(p,ii,b1+direction)]:-(w[fencode3_cd1(p,ii,field+(NDIM+1))]*w[fencode3_cd1(p,ii,b1b+direction)]+w[fencode3_cd1(p,ii,field+(2*NDIM+3))]*w[fencode3_cd1(p,ii,b1+direction)])-w[fencode3_cd1(p,ii,field+(NDIM+1))]*w[fencode3_cd1(p,ii,b1+direction)]);
        #endif


}


__device__ __host__
real fluxmom12 (real *dw, real *wd, real *w, struct params *p,int *ii,int field, int direction) {


         #if defined USE_SAC || defined USE_SAC_3D
         return(direction==2?wd[fencode3_cd1(p,ii,pressuret)]-(w[fencode3_cd1(p,ii,field+(NDIM+1))]*w[fencode3_cd1(p,ii,b1b+direction)]+w[fencode3_cd1(p,ii,field+(2*NDIM+3))]*w[fencode3_cd1(p,ii,b1+direction)])-w[fencode3_cd1(p,ii,field+(NDIM+1))]*w[fencode3_cd1(p,ii,b1+direction)]:-(w[fencode3_cd1(p,ii,field+(NDIM+1))]*w[fencode3_cd1(p,ii,b1b+direction)]+w[fencode3_cd1(p,ii,field+(2*NDIM+3))]*w[fencode3_cd1(p,ii,b1+direction)])-w[fencode3_cd1(p,ii,field+(NDIM+1))]*w[fencode3_cd1(p,ii,b1+direction)]);
        #endif


}







__device__ __host__
int computefluxrho (real *dw, real *wd, real *w, struct params *p,int *ii,int direction) {

  int field;
  int status=0;
      wd[fencode3_cd1(p,ii,flux)]=0.0;
  // for(direction=0;direction<2;direction++)
         #if defined USE_SAC || defined USE_SAC_3D
	      wd[fencode3_cd1(p,ii,flux)]= transportflux(dw,wd,w,p,ii,rho,direction)+(w[fencode3_cd1(p,ii,rhob)]*w[fencode3_cd1(p,ii,mom1+direction)])/(w[fencode3_cd1(p,ii,rhob)]+w[fencode3_cd1(p,ii,rho)]);
         #else
             wd[fencode3_cd1(p,ii,flux)]= transportflux(dw,wd,w,p,ii,rho,direction);
         #endif
  
  return ( status);
}


__device__ __host__
int computefluxmom3 (real *dw, real *wd, real *w, struct params *p,int *ii, int field,int direction) {

 
  int status=0;

#ifdef USE_SAC_3D
               wd[fencode3_cd1(p,ii,flux)]=0.0;


    		wd[fencode3_cd1(p,ii,flux)]+=transportflux(dw,wd,w,p,ii,field,direction)+fluxmom12(dw,wd,w,p,ii,field,direction);
               

#endif

  return ( status);
}




__device__ __host__
int computefluxmom2 (real *dw, real *wd, real *w, struct params *p,int *ii, int field,int direction) {

 
  int status=0;

               wd[fencode3_cd1(p,ii,flux)]=0.0;
 
        #ifdef USE_SAC
    		wd[fencode3_cd1(p,ii,flux)]+= transportflux(dw,wd,w,p,ii,field,direction)+fluxmom11(dw,wd,w,p,ii,field,direction);

 
        #endif
        #ifdef USE_SAC_3D
    		wd[fencode3_cd1(p,ii,flux)]+= transportflux(dw,wd,w,p,ii,field,direction)+fluxmom11(dw,wd,w,p,ii,field,direction);
 
 
        #endif

  return ( status);
}







__device__ __host__
int computefluxmom1 (real *dw, real *wd, real *w, struct params *p,int *ii, int field,int direction) {

 
  int status=0;

               wd[fencode3_cd1(p,ii,flux)]=0.0;



        #ifdef ADIABHYDRO
     		wd[fencode3_cd1(p,ii,flux)]+= transportflux(dw,wd,w,p,ii,field,direction);
        #endif
        #ifdef USE_SAC
    		wd[fencode3_cd1(p,ii,flux)]+= transportflux(dw,wd,w,p,ii,field,direction)+fluxmom10(dw,wd,w,p,ii,field,direction);
 
        #endif
        #ifdef USE_SAC_3D
    		wd[fencode3_cd1(p,ii,flux)]= transportflux(dw,wd,w,p,ii,field,direction)+fluxmom10(dw,wd,w,p,ii,field,direction);
 
        #endif
        
  return ( status);
}







//rho, mom1, mom2, mom3, energy, b1, b2, b3
__device__ __host__
void computeflux (real *dw, real *wd, real *w, struct params *p,int *ii, int field,int dir) {

  //int status=0;
  switch(field)
  {
     case rho:
      //computevel3_cd1(w,wd,p,ii);
      //computept3_cd1(w,wd,p,ii);
      computefluxrho(dw,wd,w,p,ii,dir);
     break;
     case mom1:
      computefluxmom1(dw,wd,w,p,ii,field,dir);
      //wd[fencode3_cd1(p,ii,f1)]+=wd[fencode3_cd1(p,ii,pressuret)];
     break;
     case mom2:
       computefluxmom2(dw,wd,w,p,ii,field,dir);
       //wd[fencode3_cd1(p,ii,f2)]+=wd[fencode3_cd1(p,ii,pressuret)];
     break;
     #ifdef USE_SAC_3D
       case mom3:
        computefluxmom3(dw,wd,w,p,ii,field,dir);
        //wd[fencode3_cd1(p,ii,f3)]+=wd[fencode3_cd1(p,ii,pressuret)];
       break;
     #endif
  }
  //return ( status);
}




__global__ void centdiff1init_parallel(struct params *p, struct state *s, real *w, real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real dt, int f, int dir)
{
  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int fid;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dy=p->dx[1];
  real dx=p->dx[0];

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     
 

   fid=0;
   

     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
                        {
                            dwn1[fencode3_cd1(p,ii,f)]=0.0;
                  	    //for(fid=0;fid<2;fid++)
                               wd[fencode3_cd1(p,ii,flux)]=0.0;
                               //wmod[fencode_cd1(p,i,j,flux)+order*NVAR*(p->n[0])*(p->n[1])]=0.0;
                               /*wmod[fencode3_cd1(p,ii,b1b)+(ordero*NVAR*dimp)]=0.0;
                               wmod[fencode3_cd1(p,ii,b1b)+(order*NVAR*dimp)]=0.0;
                               wmod[fencode3_cd1(p,ii,b2b)+(ordero*NVAR*dimp)]=0.0;
                               wmod[fencode3_cd1(p,ii,b2b)+(order*NVAR*dimp)]=0.0;
                               wmod[fencode3_cd1(p,ii,b3b)+(ordero*NVAR*dimp)]=0.0;
                               wmod[fencode3_cd1(p,ii,b3b)+(order*NVAR*dimp)]=0.0;
                               wmod[fencode3_cd1(p,ii,energyb)+(ordero*NVAR*dimp)]=0.0;
                               wmod[fencode3_cd1(p,ii,energyb)+(order*NVAR*dimp)]=0.0;*/
                        }

   
 __syncthreads();                       




}



__global__ void centdiff1_parallel(struct params *p, struct state *s, real *w, real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real dt, int f, int dir)
{
  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int fid;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dy=p->dx[1];
  real dx=p->dx[0];

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     


   fid=0;




     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif


//1. 11/1/11 could swap cases below
                        switch(dir)
                        {
                         case 0:
                          #ifdef USE_SAC_3D
       				if(ii[0]<p->n[0] && ii[1]>1 && ii[1]<(p->n[1]-2) && ii[2]>1 && ii[2]<(p->n[2]-2))
     			  #else
       				if(ii[0]<p->n[0] && ii[1]>1 && ii[1]<(p->n[1]-2))
     			  #endif
                         //if(i<(ni)  && j >1 &&  j<(nj-1))
                            computeflux(dwn1,wd,wmod+order*NVAR*dimp,p,ii,f,0); 
                         break;
                         case 1:
                          #ifdef USE_SAC_3D
       				if(ii[1]<p->n[1] && ii[0]>1 && ii[0]<(p->n[0]-2) && ii[2]>1 && ii[2]<(p->n[2]-2))
     			  #else
       				if(ii[1]<p->n[1] && ii[0]>1 && ii[0]<(p->n[0]-2))
     			  #endif
                         //if(i>1 &&  i<(ni-1) && j<(nj))
                            computeflux(dwn1,wd,wmod+order*NVAR*dimp,p,ii,f,1); 
                         break;
                          #ifdef USE_SAC_3D
                         case 2:

       				if(ii[2]<p->n[2] && ii[0]>1 && ii[0]<(p->n[0]-2) && ii[1]>1 && ii[1]<(p->n[1]-2))

                         //if(i>1 &&  i<(ni-1) && j<(nj))
                            computeflux(dwn1,wd,wmod+order*NVAR*dimp,p,ii,f,2); 
                         break;
                         #endif
                        }
              //  }
                        //might need to set boundaries correctly
 

__syncthreads();                        



}










__global__ void centdiff1a_parallel(struct params *p, struct state *s, real *w, real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real dt, int f, int dir)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  //int i = blockIdx.x * blockDim.x + threadIdx.x;
  //int j = blockIdx.y * blockDim.y + threadIdx.y;

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int fid;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];

  //real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
  //real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
    int nk=p->n[2];
    real dz=p->dx[2];
#endif
 #ifdef USE_SAC_3D
   int kp;
   
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     

   fid=0;





     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif

			// if(i>1 && j >1 && i<(ni-2) && j<(nj-2))
			     #ifdef USE_SAC
				   if(ii[0]>1 && ii[1] >1 && ii[0]<(ni-2) && ii[1]<(nj-2))
			     #endif
			     #ifdef USE_SAC_3D
				  if(ii[0]>1 && ii[1] >1 && ii[2] >1 && ii[0]<(ni-2) && ii[1]<(nj-2) && ii[2]<(nk-2))
			     #endif                        
                               divflux1(dwn1,wd,wmod+order*NVAR*dimp,p,ii,f,dir);  


 __syncthreads();


}

__global__ void centdiff1af_parallel(struct params *p, struct state *s, real *w, real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real dt, int f, int dir)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  //int i = blockIdx.x * blockDim.x + threadIdx.x;
  //int j = blockIdx.y * blockDim.y + threadIdx.y;

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int fid;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];

  //real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
  //real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
    int nk=p->n[2];
    real dz=p->dx[2];
#endif
 #ifdef USE_SAC_3D
   int kp;
   
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     


   fid=0;

             // for(int f=rho; f<=mom3; f++)
              // {



     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif


                        switch(dir)
                        {
                         case 0:

                         //if(i<(ni)  && j >1 &&  j<(nj-2))
			     #ifdef USE_SAC
				   if(ii[1]>1 && ii[1] <(nj-2) && ii[0]<(ni) )
			     #endif
			     #ifdef USE_SAC_3D
				   if(ii[1]>1 && ii[1] <(nj-2) && ii[0]<(ni) &&  ii[2]>1 && ii[2] <(nk-2) )
			     #endif                          
                              wmod[fencode3_cd1(p,ii,f)+(ordero*NVAR*dimp)]=wmod[fencode3_cd1(p,ii,f)+(ordero*NVAR*dimp)]-dt*dwn1[fencode3_cd1(p,ii,f)]; 
                         break;
                         case 1:
			     #ifdef USE_SAC
				   if(ii[0]>1 && ii[1] <(nj) && ii[0]<(ni-2) )
			     #endif
			     #ifdef USE_SAC_3D
				   if(ii[0]>1 && ii[1] <(nj) && ii[0]<(ni-2) &&  ii[2]>1 && ii[2] <(nk-2) )
			     #endif 
                         
                         //if(i>1 &&  i<(ni-2) && j<(nj))
                              wmod[fencode3_cd1(p,ii,f)+(ordero*NVAR*dimp)]=wmod[fencode3_cd1(p,ii,f)+(ordero*NVAR*dimp)]-dt*dwn1[fencode3_cd1(p,ii,f)];
                         break;
                         case 2:

 
			     #ifdef USE_SAC_3D
				   if(ii[0]>1 &&  ii[0]<(ni-2)  && ii[1]>1 &&  ii[1]<(nj-2) && ii[2] <(nk) )
                               wmod[fencode3_cd1(p,ii,f)+(ordero*NVAR*dimp)]=wmod[fencode3_cd1(p,ii,f)+(ordero*NVAR*dimp)]-dt*dwn1[fencode3_cd1(p,ii,f)];
			     #endif                         
                         //if(i>1 &&  i<(ni-2) && j<(nj))
                             
                         break;
                        }


              //  }
	

  __syncthreads();


}


__global__ void centdiff1binit_parallel(struct params *p, struct state *s, real *w, real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real dt, int f, int dir)
{
  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int fid;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dy=p->dx[1];
  real dx=p->dx[0];

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     
 

   fid=0;
   

     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
                        {
                            dwn1[fencode3_cd1(p,ii,f)]=0.0;

                        }

   
 __syncthreads();                       




}


__global__ void centdiff1b_parallel(struct params *p, struct state *s, real *w, real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real dt, int f, int dir)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  //int i = blockIdx.x * blockDim.x + threadIdx.x;
  //int j = blockIdx.y * blockDim.y + threadIdx.y;

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int fid;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];

  //real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
  //real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
    int nk=p->n[2];
    real dz=p->dx[2];
#endif
 #ifdef USE_SAC_3D
   int kp;
   
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     

   fid=0;





     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif


#if(defined(USE_USERSOURCE))
   {

     ii[0]=ip;
     ii[1]=jp;
#endif
     #if(defined(USE_SAC_3D) && defined(USE_USERSOURCE))
	   ii[2]=kp;
     #endif


     #if(defined(USE_SAC_3D) && defined(USE_USERSOURCE))
       //if(ii[0]<((p->n[0])-2) && ii[1]<((p->n[1])-2) && ii[2]<((p->n[2])-2)     && ii[0]>1    &&  ii[1]>1   && ii[2]>1   )
       if(ii[0]<((p->n[0])) && ii[1]<((p->n[1])) && ii[2]<((p->n[2]))    )
     #endif
     #if(defined(USE_SAC) && defined(USE_USERSOURCE))
       //if(ii[0]<(p->n[0])-2 && ii[1]<(p->n[1])-2)
      if(ii[0]<(p->n[0]) && ii[1]<(p->n[1]))
     #endif

                     #ifdef USE_USERSOURCE
                            addsourceterms1_cd1(dwn1,wd,wmod+ordero*NVAR*dimp,p,s,ii,f,dir); 


                      }
                    __syncthreads();
                     #endif



               // }
    


}



__global__ void centdiff1bf_parallel(struct params *p, struct state *s, real *w, real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real dt, int f, int dir)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  //int i = blockIdx.x * blockDim.x + threadIdx.x;
  //int j = blockIdx.y * blockDim.y + threadIdx.y;

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int fid;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];

  //real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
  //real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
    int nk=p->n[2];
    real dz=p->dx[2];
#endif
 #ifdef USE_SAC_3D
   int kp;
   
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     


   fid=0;

             // for(int f=rho; f<=mom3; f++)
              // {



     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif


			     #ifdef USE_SAC
				   if( ii[1] <(nj) && ii[0]<(ni) )
			     #endif
			     #ifdef USE_SAC_3D
				   if(ii[1] <(nj) && ii[0]<(ni) &&   ii[2] <(nk) )
			     #endif                          
                              wmod[fencode3_cd1(p,ii,f)+(ordero*NVAR*dimp)]=wmod[fencode3_cd1(p,ii,f)+(ordero*NVAR*dimp)]-dt*dwn1[fencode3_cd1(p,ii,f)];

            

  __syncthreads();


}



/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_cd1(char *label)
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

int cucentdiff1(struct params **p, struct params **d_p,struct state **d_s, real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero, real dt, int field, int dir)
{
 int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 
 //dim3 dimBlock(dimblock, 1);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   // dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;
 //  cudaMemcpy(*w, *d_w, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
 // if(order==0)
    cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);
    //printf("gamma %g\n", (*p)->gamma);
     centdiff1init_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_s,*d_w,*d_wmod, *d_dwn1,  *d_wd, order, ordero,dt,field,dir);
     //prop_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_b,*d_u,*d_v,*d_h);
	    //printf("called prop\n"); 
     cudaThreadSynchronize();
     centdiff1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_s,*d_w,*d_wmod, *d_dwn1,  *d_wd, order, ordero,dt,field,dir);
     //prop_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_b,*d_u,*d_v,*d_h);
	    //printf("called prop\n"); 
     cudaThreadSynchronize();
     centdiff1a_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_s,*d_w,*d_wmod, *d_dwn1,  *d_wd, order, ordero,dt,field,dir);
     cudaThreadSynchronize();


     centdiff1af_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_s,*d_w,*d_wmod, *d_dwn1,  *d_wd, order, ordero,dt,field,dir);
     cudaThreadSynchronize();
     
    // centdiff1binit_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_s,*d_w,*d_wmod, *d_dwn1,  *d_wd, order, ordero,dt,field,dir);
    // cudaThreadSynchronize();
     
   //  centdiff1b_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_s,*d_w,*d_wmod, *d_dwn1,  *d_wd, order, ordero,dt,field,dir);
   //  cudaThreadSynchronize();
     
   //  centdiff1bf_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_s,*d_w,*d_wmod, *d_dwn1,  *d_wd, order, ordero,dt,field,dir);
   //  cudaThreadSynchronize();
     
}


