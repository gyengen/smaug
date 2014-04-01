#include "../include/cudapars.h"
#include "../include/paramssteeringtest1.h"

/////////////////////////////////////
// standard imports
/////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include "../include/smaugcukernels.h"

/////////////////////////////////////
// kernel function (CUDA device)
/////////////////////////////////////
#include "../include/gradops_b.cuh"

__global__ void boundary_parallel_lower0(struct params *p, struct bparams *bp, struct state *s,  real *wmod, int order, int dir, int field)
{

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;

  int ni=p->n[0];
  int nj=p->n[1];
  
  int ip,jp,ipg,jpg;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));


 #if defined USE_SAC
    j=iindex/(ni);
   i=iindex-(jp*ni);
#endif  


int shift=order*NVAR*dimp;



     iia[0]=i;
     iia[1]=j;
     k=0;

      for( f=rho; f<=b2; f++)
      {  
	       if(i<((p->n[0])) && j<((p->n[1]))) 
		{
		  	//wmod[order*NVAR*dimp+encode3_b(p,i,j,k,f)]=((( ) )?((i==0 ||
			//i==1) && dir==0)*wmod[encode3_b(p,(p->n[0])-4+i,j,k,f)]+(((i==((p->n[0])-1)) || (i==((p->n[0])-2))) &&
			//dir==0)*wmod[order*NVAR*dimp+encode3_b(p,4-(p->n[0])+i,j,k,f)]:wmod[order*NVAR*dimp+encode3_b(p,i,j,k,f)]);

			wmod[order*NVAR*dimp+encode3_b(p,i,j,k,f)]=( i==0 || i==1  ?  wmod[encode3_b(p,(p->n[0])-4+i,j,k,f)] :wmod[order*NVAR*dimp+encode3_b(p,i,j,k,f)]);
		}
       }

 __syncthreads();





  
}


__global__ void boundary_parallel_lower1(struct params *p, struct bparams *bp, struct state *s,  real *wmod, int order, int dir, int field)
{

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;

  int ni=p->n[0];
  int nj=p->n[1];
  
  int ip,jp,ipg,jpg;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));


 #if defined USE_SAC
    j=iindex/(ni);
   i=iindex-(jp*ni);
#endif  


int shift=order*NVAR*dimp;



     iia[0]=i;
     iia[1]=j;
     k=0;

      for( f=rho; f<=b2; f++)
      {  
	       if(i<((p->n[0])) && j<((p->n[1]))) 
		{
		  	//wmod[order*NVAR*dimp+encode3_b(p,i,j,k,f)]=((( i==((p->n[0])-1) || i==((p->n[0])-2)) && dir==0)?((i==0 ||
			//i==1) && dir==0)*wmod[encode3_b(p,(p->n[0])-4+i,j,k,f)]+(((i==((p->n[0])-1)) || (i==((p->n[0])-2))) &&
			//dir==0)*wmod[order*NVAR*dimp+encode3_b(p,4-(p->n[0])+i,j,k,f)]:wmod[order*NVAR*dimp+encode3_b(p,i,j,k,f)]);

		      // bc3_periodic1_dir_b(wmod+order*NVAR*dimp,p,iia,f,dir);  //for OZT
		       //  bc3_cont_cd4_b(wmod+order*NVAR*dimp,p,iia,f);  //for BW
		       //  bc3_fixed_b(wmod+order*NVAR*dimp,p,iia,f,0.0);
		       
		       	wmod[order*NVAR*dimp+encode3_b(p,i,j,k,f)]=( j==0 || j==1  ?  wmod[encode3_b(p,i,(p->n[1])-4+j,k,f)] :wmod[order*NVAR*dimp+encode3_b(p,i,j,k,f)]);

			//wmod[order*NVAR*dimp+encode3_b(p,i,j,k,f)]=( i==0 || i==1  ?  wmod[encode3_b(p,(p->n[0])-4+i,j,k,f)] :wmod[order*NVAR*dimp+encode3_b(p,i,j,k,f)]);

		       
		       
		}
       }

 __syncthreads();





  
}




__global__ void boundary_parallel_upper0(struct params *p, struct bparams *bp, struct state *s,  real *wmod, int order, int dir, int field)
{

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;

  int ni=p->n[0];
  int nj=p->n[1];
  
  int ip,jp,ipg,jpg;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));


 #if defined USE_SAC
    j=iindex/(ni);
   i=iindex-(jp*ni);
#endif  


int shift=order*NVAR*dimp;



     iia[0]=i;
     iia[1]=j;
     k=0;

      for( f=rho; f<=b2; f++)
      {  
	       if(i<((p->n[0])) && j<((p->n[1]))) 
		{
		       //bc3_periodic1_dir_b(wmod+order*NVAR*dimp,p,iia,f,dir);  //for OZT
		       //  bc3_cont_cd4_b(wmod+order*NVAR*dimp,p,iia,f);  //for BW
		       //  bc3_fixed_b(wmod+order*NVAR*dimp,p,iia,f,0.0);
		       
		       	wmod[order*NVAR*dimp+encode3_b(p,i,j,k,f)]=( ((i==((p->n[0])-1)) || (i==((p->n[0])-2)))?wmod[order*NVAR*dimp+encode3_b(p, 4-(p->n[0])+i  ,j,k,f)] :wmod[order*NVAR*dimp+encode3_b(p,i,j,k,f)]);

		}
       }

 __syncthreads();





  
}



__global__ void boundary_parallel_upper1(struct params *p, struct bparams *bp, struct state *s,  real *wmod, int order, int dir, int field)
{

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;

  int ni=p->n[0];
  int nj=p->n[1];
  
  int ip,jp,ipg,jpg;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));


 #if defined USE_SAC
    j=iindex/(ni);
   i=iindex-(jp*ni);
#endif  


int shift=order*NVAR*dimp;



     iia[0]=i;
     iia[1]=j;
     k=0;

      for( f=rho; f<=b2; f++)
      {  
	       if(i<((p->n[0])) && j<((p->n[1]))) 
		{
		       //bc3_periodic1_dir_b(wmod+order*NVAR*dimp,p,iia,f,dir);  //for OZT
		       //  bc3_cont_cd4_b(wmod+order*NVAR*dimp,p,iia,f);  //for BW
		       //  bc3_fixed_b(wmod+order*NVAR*dimp,p,iia,f,0.0);
		       wmod[order*NVAR*dimp+encode3_b(p,i,j,k,f)]=( ((j==((p->n[1])-1)) || (j==((p->n[1])-2)))?wmod[order*NVAR*dimp+encode3_b(p, i,4-(p->n[1])+j  ,k,f)] :wmod[order*NVAR*dimp+encode3_b(p,i,j,k,f)]);

		}
       }

 __syncthreads();





  
}








int cuboundary(struct params **p, struct bparams **bp,struct params **d_p, struct bparams **d_bp, struct state **d_s,  real **d_wmod,  int order,int idir,int field)
{


 dim3 dimBlock(dimblock, 1);


int numBlocks = ((dimproduct_b(*p)+numThreadsPerBlock-1)) / numThreadsPerBlock;



if(idir==0 &&  ((*p)->boundtype[field][0][0])==0)
    boundary_parallel_lower0<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);

    cudaThreadSynchronize();
if(idir==1 && ((*p)->boundtype[field][1][0])==0)
    boundary_parallel_lower1<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
    cudaThreadSynchronize();
    
if(idir ==0 && ((*p)->boundtype[field][0][1])==0)
    boundary_parallel_upper0<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);

    cudaThreadSynchronize();
if(idir == 1 && ((*p)->boundtype[field][1][1])==0)
    boundary_parallel_upper1<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
    cudaThreadSynchronize();
    
    
    
    
}
