#include "../include/cudapars.h"
#include "../include/iotypes.h"
#include "../include/iobparams.h"
/////////////////////////////////////
// standard imports
/////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include "../include/smaugcukernels.h"

/////////////////////////////////////
// kernel function (CUDA device)
/////////////////////////////////////
#include "../include/gradops_mpiu.cuh"
//#include "../include/init_user_mpiu.cuh"





__device__ __host__
int encodempiw (struct params *p,int ix, int iy, int iz, int field,int bound,int dim) {
  #ifdef USE_SAC_3D
    return (dim*(    4*NVAR*(         ((p->n[0])*(p->n[1]))+((p->n[1])*(p->n[2]))+((p->n[0])*(p->n[2]))   )           )+4*field*(         ((p->n[0])*(p->n[1]))+((p->n[1])*(p->n[2]))+((p->n[0])*(p->n[2]))   )+
bound*(         (dim==2)*((p->n[0])*(p->n[1]))   +  (dim==0)*((p->n[1])*(p->n[2]))  +   (dim==1)*((p->n[0])*(p->n[2]))    )+   (  (ix+iz*(p->n[0]))*(dim==1)+(iy+iz*(p->n[1]))*(dim==0)+(iz+ix*(p->n[2]))*(dim==2)    ));
  #else
    return (dim*(4*NVAR*((p->n[0])+(p->n[1])))+4*field*((p->n[0])+(p->n[1]))+bound*((dim==1)*(p->n[0])+(dim==0)*(p->n[1]))  +   (ix*(dim==1)+iy*(dim==0)));
  #endif
}

__device__ __host__
int encodempiw0 (struct params *p,int ix, int iy, int iz, int field,int bound) {
  #ifdef USE_SAC_3D
    return (4*field*(         ((p->n[1])*(p->n[2]))   )+
bound*(            +  ((p->n[1])*(p->n[2]))      )+   (  (iy+iz*(p->n[1]))    ));
  #else
    return (   4*field*(p->n[1]) +bound*((p->n[1]))  +   (iy)   );
  #endif
}


__device__ __host__
int encodempiw1 (struct params *p,int ix, int iy, int iz, int field,int bound) {
  #ifdef USE_SAC_3D
    return (4*field*(         ((p->n[0])*(p->n[2]))   )+
bound*(            +  ((p->n[0])*(p->n[2]))      )+   (  (ix+iz*(p->n[0]))    ));
  #else
    return (4*field*(p->n[0]) +bound*((p->n[0]))  +   (ix));
  #endif
}

__device__ __host__
int encodempiw2 (struct params *p,int ix, int iy, int iz, int field,int bound) {
  #ifdef USE_SAC_3D
    return (4*field*(         ((p->n[0])*(p->n[1]))   )+
bound*(            +  ((p->n[0])*(p->n[1]))      )+   (  (ix+iy*(p->n[0]))    ));
  #endif
}


__device__ __host__
int encodempivisc (struct params *p,int ix, int iy, int iz, int bound,int dim) {
  #ifdef USE_SAC_3D
    return (dim*(    2*(         (((p->n[0])+2)*((p->n[1])+2))+(((p->n[1])+2)*((p->n[2])+2))+(((p->n[0])+2)*((p->n[2])+2))   )           )+
bound*(         (dim==2)*(((p->n[0])+2)*((p->n[1])+2))   +  (dim==0)*(((p->n[1])+2)*((p->n[2])+2))  +   (dim==1)*(((p->n[0])+2)*((p->n[2])+2))    )+   (  (ix+iz*((p->n[0])+2))*(dim==1)+(iy+iz*((p->n[1])+2))*(dim==0)+(iz+ix*((p->n[2])+2))*(dim==2)    ));
  #else
    return (   dim*(2*(  ((p->n[0])+2)+((p->n[1])+2)   ))      +bound*(    (dim==1)*((p->n[0])+2)+(dim==0)*((p->n[1])+2)  )  +   (ix*(dim==1)+iy*(dim==0))     );
  #endif
}


__device__ __host__
int encodempivisc0 (struct params *p,int ix, int iy, int iz, int bound,int dim) {
  #ifdef USE_SAC_3D
    return (
bound*(           (((p->n[1])+2)*((p->n[2])+2))      )+   (  (iy+iz*((p->n[1])+2))    ));
  #else
    return (   bound*(    ((p->n[1])+2)  )  +   iy     );
  #endif
}


__device__ __host__
int encodempivisc1 (struct params *p,int ix, int iy, int iz, int bound,int dim) {
  #ifdef USE_SAC_3D
    return (
bound*(           (((p->n[0])+2)*((p->n[2])+2))      )+   (  (ix+iz*((p->n[0])+2))    ));
  #else
    return (   bound*(    ((p->n[0])+2)  )  +   ix     );
  #endif
}

__device__ __host__
int encodempivisc2 (struct params *p,int ix, int iy, int iz, int bound,int dim) {
  #ifdef USE_SAC_3D
    return (
bound*(           (((p->n[0])+2)*((p->n[1])+2))      )+   (  (ix+iy*((p->n[0])+2))    ));
  #endif
}


#ifdef USE_MPI

     __device__ __host__ void mpiwtogpu(struct params *p,real *d_w,real *d_wmod,real *d_mpiw0,real *d_mpiwmod0,real *d_mpiw1,real *d_mpiwmod1,real *d_mpiw2,real *d_mpiwmod2,int *ii, int var, int dim)
    {

             int i,j,k,bound;
i=ii[0];
j=ii[1];
k=0;




                if((i==0 || i==1) && dim==0)
                {              
                    bound=i;
                    d_w[encode3_mpiu(p,i,j,k,var)]=d_mpiw0[encodempiw0(p,i,j,k,var,bound)];
                   // d_wmod[encode3_mpiu(p,i,j,k,var)]=d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)];
                    // if(var==4 && ((p)->ipe)==0)                        
                    //    printf(" %d %d %d %d actual %d  mpi data%d %g\n",i,j,bound,dim,var,encodempiw0(p,i,j,k,var,bound),d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]);   


 

     
      
                }
                else if((( i>=((p->n[0])-2)   ))  && dim==0)               
                {
                    bound=2*(i==((p->n[0])-1))+(p->n[0])-i;
                    d_w[encode3_mpiu(p,i,j,k,var)]=d_mpiw0[encodempiw0(p,i,j,k,var,bound)];
                  //  d_wmod[encode3_mpiu(p,i,j,k,var)]=d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)];    


                          // if(p->ipe==0    && var==rho && dim==0 )
                          //   {
                             // for(int bound=0;bound<=1;bound++)
                             //   printf("mpiw0 %d %d %d %d %lg \n",dim,bound,i,j,d_mpiw0[encodempiw0(p,i,j,k,var,bound)]);
			     	//printf("mpiwmod0 %d %d %d %d %lg %lg\n\n",dim,bound,i,j,d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)],d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)]);
                           //  }
          
                }

              

                if((j==0 || j==1) && dim==1)              
                {              
                    bound=j;
                    d_w[encode3_mpiu(p,i,j,k,var)]=d_mpiw1[encodempiw1(p,i,j,k,var,bound)];
                   // d_wmod[encode3_mpiu(p,i,j,k,var)]=d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)];              
                }            
                 else if((( j>=((p->n[1])-2)   ))  && dim==1)               
                {
                   bound=2*(j==((p->n[1])-1))+(p->n[1])-j;
                    d_w[encode3_mpiu(p,i,j,k,var)]=d_mpiw1[encodempiw1(p,i,j,k,var,bound)];
                  //  d_wmod[encode3_mpiu(p,i,j,k,var)]=d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)];              
                }

       #ifdef USE_SAC_3D
               k=ii[2];
                if((k==0 || k==1) && dim==2)              
                {              
                    bound=k;
                    d_w[encode3_mpiu(p,i,j,k,var)]=d_mpiw2[encodempiw2(p,i,j,k,var,bound)];
                  //  d_wmod[encode3_mpiu(p,i,j,k,var)]=d_mpiwmod2[encodempiw2(p,i,j,k,var,bound)];              
                }        
                 else if((( k>=((p->n[2])-2)   ))  && dim==2)               
                {
                    bound=2*(k==((p->n[2])-1))+(p->n[2])-k;
                    d_w[encode3_mpiu(p,i,j,k,var)]=d_mpiw2[encodempiw2(p,i,j,k,var,bound)];
                  //  d_wmod[encode3_mpiu(p,i,j,k,var)]=d_mpiwmod2[encodempiw2(p,i,j,k,var,bound)];              
                }

     #endif
 
// if( (p->ipe)==0  && ((p)->it)==0 && (isnan(d_wmod[fencode3_mpiu(p,ii,rho)]) || d_wmod[fencode3_mpiu(p,ii,rho)]==0      ))
//        { 
//    	printf("nant %d %d %d %lg\n",ii[0],ii[1], dim, d_wmod[fencode3_mpiu(p,ii,rho)] );
//}

    }

     __device__ __host__ void mpiwmodtogpu(struct params *p,real *d_w,real *d_wmod,real *d_mpiw0,real *d_mpiwmod0,real *d_mpiw1,real *d_mpiwmod1,real *d_mpiw2,real *d_mpiwmod2,int *ii, int var, int dim, int order)
    {

             int i,j,k,bound;
i=ii[0];
j=ii[1];
k=0;

 #ifdef USE_SAC_3D
	int ntot=((p->n[0]))*((p->n[1]))*((p->n[2]))*NVAR;
 #else
	int ntot=((p->n[0]))*((p->n[1]))*NVAR;
 #endif
                //remember only update the boundaries if they are mpiupper boundaries 
                //or an mpi period 


                if((i==0 || i==1) && dim==0 /* && ( ((p->mpilowerb[dim])==1) || ((p->boundtype[0][dim][0])==2))  */  )
                {              
                    bound=i;
                    d_wmod[order*ntot+encode3_mpiu(p,i,j,k,var)]=d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)];
                     //if(var==0 && ((p)->ipe)==0)                        
                     //   printf(" %d %d %d %d actual mpi data %d %g\n",i,j,bound,dim,encodempiw0(p,i,j,k,var,bound),d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]); 

       
      
                }
                else if((( i>=((p->n[0])-2)   ))  && dim==0 /* && ( ((p->mpiupperb[dim])==1) || ((p->boundtype[0][dim][0])==2)) */)               
                {
                    bound=2*(i==((p->n[0])-1))+(p->n[0])-i;
                    d_wmod[order*ntot+encode3_mpiu(p,i,j,k,var)]=d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)];  

                             //if(p->ipe==0    && var==rho && dim==0 )
                             //{
                             // for(int bound=0;bound<=1;bound++)
                             //   printf("mpiw0 %d %d %d %d %lg \n",dim,bound,i,j,d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]);
			     //	printf("mpiwmod0 %d %d %d %d %lg %lg\n\n",dim,bound,i,j,d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)],d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)]);
                             //}
          
                }

              

                if((j==0 || j==1) && dim==1  /*&& ( ((p->mpilowerb[dim])==1) || ((p->boundtype[0][0][0])==2))*/   )              
                {              
                    bound=j;
                    d_wmod[order*ntot+encode3_mpiu(p,i,j,k,var)]=d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)]; 
//if(i>=0 && i<10  &&  var==rhob)
//printf("nani %d %d %d %d  %lg %lg %d  \n",p->ipe,order,i,j, d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)],d_wmod[order*ntot+encode3_mpiu(p,i,j,k,var)],bound );

             
                }            
                 else if((( j>=((p->n[1])-2)   ))  && dim==1  /* && ( ((p->mpiupperb[dim])==1) || ((p->boundtype[0][0][0])==2)) */   )               
                {
                   bound=2*(j==((p->n[1])-1))+(p->n[1])-j;
                    d_wmod[order*ntot+encode3_mpiu(p,i,j,k,var)]=d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)]; 
//if(i>=0 && i<10  &&  var==rhob)
//printf("nani %d %d %d %d  %lg %lg %d  \n",p->ipe,order,i,j, d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)],d_wmod[order*ntot+encode3_mpiu(p,i,j,k,var)],bound );

  

             
                }

       #ifdef USE_SAC_3D
               k=ii[2];
                if((k==0 || k==1) && dim==2  /* && ( ((p->mpilowerb[dim])==1) || ((p->boundtype[0][dim][0])==2))  */ )              
                {              
                    bound=k;
                    d_wmod[order*ntot+encode3_mpiu(p,i,j,k,var)]=d_mpiwmod2[encodempiw2(p,i,j,k,var,bound)];              
                }        
                 else if((( k>=((p->n[2])-2)   ))  && dim==2  /* && ( ((p->mpiupperb[dim])==1) || ((p->boundtype[0][dim][0])==2))  */ )               
                {
                    bound=2*(k==((p->n[2])-1))+(p->n[2])-k;
                    d_wmod[order*ntot+encode3_mpiu(p,i,j,k,var)]=d_mpiwmod2[encodempiw2(p,i,j,k,var,bound)];              
                }

     #endif
 
// if( (p->ipe)==0  && ((p)->it)==0 && (isnan(d_wmod[fencode3_mpiu(p,ii,rho)]) || d_wmod[fencode3_mpiu(p,ii,rho)]==0      ))
//        { 
//    	printf("nant %d %d %d %lg\n",ii[0],ii[1], dim, d_wmod[fencode3_mpiu(p,ii,rho)] );
//}

    }



     __device__ __host__ void mpiwdtogpu(struct params *p,real *d_wd,real *d_mpiw0,real *d_mpiw1,real *d_mpiw2,int *ii, int var, int dim)
    {

             int i,j,k,bound;
i=ii[0];
j=ii[1];
k=0;

                //remember only update the boundaries if they are mpiupper boundaries 
                //or an mpi period 



 
                if((i==0 || i==1) && dim==0  /* &&  ((p->mpilowerb[dim])==1)*/)
                {              
                    bound=i;
                    d_wd[encode3_mpiu(p,i,j,k,var)]=d_mpiw0[encodempiw0(p,i,j,k,var,bound)];
                    // if(var==4 && ((p)->ipe)==0)                        
                    //    printf(" %d %d %d %d actual %d  mpi data%d %g\n",i,j,bound,dim,var,encodempiw0(p,i,j,k,var,bound),d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]);        
      
                }
                else if((( i>=((p->n[0])-2)   ))  && dim==0 /* &&  ((p->mpiupperb[dim])==1)*/)               
                {
                    bound=2*(i==((p->n[0])-1))+(p->n[0])-i;
                    d_wd[encode3_mpiu(p,i,j,k,var)]=d_mpiw0[encodempiw0(p,i,j,k,var,bound)];
                }

              

                if((j==0 || j==1) && dim==1  /* &&  ((p->mpilowerb[dim])==1)*/)              
                {              
                    bound=j;
                    d_wd[encode3_mpiu(p,i,j,k,var)]=d_mpiw1[encodempiw1(p,i,j,k,var,bound)];
                }            
                 else if((( j>=((p->n[1])-2)   ))  && dim==1 /* &&  ((p->mpiupperb[dim])==1)*/)               
                {
                   bound=2*(j==((p->n[1])-1))+(p->n[1])-j;
                    d_wd[encode3_mpiu(p,i,j,k,var)]=d_mpiw1[encodempiw1(p,i,j,k,var,bound)];
                }


               /* if((i==0 || i==1) && (j==0 || j==1))
                {              
                    bound=i;
                    d_wd[encode3_mpiu(p,i,j,k,var)]=d_mpiw0[encodempiw0(p,i,j,k,var,bound)];
                    // if(var==4 && ((p)->ipe)==0)                        
                    //    printf(" %d %d %d %d actual %d  mpi data%d %g\n",i,j,bound,dim,var,encodempiw0(p,i,j,k,var,bound),d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]);        
      
                }
                if((( j>=((p->n[1])-2)   ))  && (( i>=((p->n[0])-2)   )))               
                {
                   bound=2*(j==((p->n[1])-1))+(p->n[1])-j;
                    d_wd[encode3_mpiu(p,i,j,k,var)]=d_mpiw1[encodempiw0(p,i,j,k,var,bound)];
                }*/
                



       #ifdef USE_SAC_3D
               k=ii[2];
                if((k==0 || k==1) && dim==2)              
                {              
                    bound=k;
                    d_wd[encode3_mpiu(p,i,j,k,var)]=d_mpiw2[encodempiw2(p,i,j,k,var,bound)];
                }        
                 else if((( k>=((p->n[2])-2)   ))  && dim==2)               
                {
                    bound=2*(k==((p->n[2])-1))+(p->n[2])-k;
                    d_wd[encode3_mpiu(p,i,j,k,var)]=d_mpiw2[encodempiw2(p,i,j,k,var,bound)];
                }

     #endif
 


    }




__device__ __host__ void   mpivisctogpu(struct params *p,real *d_wtemp2,real *d_gmpivisc0,real *d_gmpivisc1,real *d_gmpivisc2,int *ii,  int dim)
{
                                
               int i,j,k,bound,var;
              var=0;
i=ii[0];
j=ii[1];
k=0;
                //remember only update the boundaries if they are mpiupper boundaries 
                //or an mpi period 
 
                if((i==0 ) && dim==0 /* && ( ((p->mpilowerb[dim])==1) || ((p->boundtype[0][dim][0])==2))*/)
                {              
                    bound=i;
                    d_wtemp2[encode3p2_mpiu(p,i,j,k,var)]=d_gmpivisc0[encodempivisc0(p,i,j,k,bound,dim)];
                    
                }
                else if((( i==((p->n[0])+1)   ))  && dim==0  /* && ( ((p->mpiupperb[dim])==1) || ((p->boundtype[0][dim][0])==2)) */ )               
                {
                    bound=1;
                    d_wtemp2[encode3p2_mpiu(p,i,j,k,var)]=d_gmpivisc0[encodempivisc0(p,i,j,k,bound,dim)];
                }

              

                if((j==0) && dim==1 /* && ( ((p->mpilowerb[dim])==1) || ((p->boundtype[0][dim][0])==2))*/)              
                {              
                    bound=j;
                    d_wtemp2[encode3p2_mpiu(p,i,j,k,var)]=d_gmpivisc1[encodempivisc1(p,i,j,k,bound,dim)];
                }            
                 else if((( j==((p->n[1])+1)   ))  && dim==1   /*&& ( ((p->mpiupperb[dim])==1) || ((p->boundtype[0][dim][0])==2))*/)               
                {
                    bound=1;
                    d_wtemp2[encode3p2_mpiu(p,i,j,k,var)]=d_gmpivisc1[encodempivisc1(p,i,j,k,bound,dim)];
             
                }

       #ifdef USE_SAC_3D
               k=ii[2];
                if((k==0 ) && dim==2  /*&& ( ((p->mpilowerb[dim])==1) || ((p->boundtype[0][dim][0])==2))*/)              
                {              
                    bound=k;
                    d_wtemp2[encode3p2_mpiu(p,i,j,k,var)]=d_gmpivisc2[encodempivisc2(p,i,j,k,bound,dim)];
                }        
                 else if(((k==((p->n[2])+1)   ))  && dim==2   /* && ( ((p->mpiupperb[dim])==1) || ((p->boundtype[0][dim][0])==2))*/)               
                {
                    bound=1;
                    d_wtemp2[encode3p2_mpiu(p,i,j,k,var)]=d_gmpivisc2[encodempivisc2(p,i,j,k,bound,dim)];
                }

     #endif
                               
                                
}

__device__ __host__ void   gputompivisc(struct params *p,real *d_wtemp2,real *d_gmpivisc0,real *d_gmpivisc1,real *d_gmpivisc2,int *ii,  int dim)
{
                                
              int i,j,k,bound,var;
              var=0;
i=ii[0];
j=ii[1];
k=0;
 
 
                if((i==0 ) && dim==0)
                {              
                    bound=i;
                    d_gmpivisc0[encodempivisc0(p,i,j,k,bound,dim)]=d_wtemp2[encode3p2_mpiu(p,i,j,k,var)];
                    
                }
                else if((( i==((p->n[0])+1)   ))  && dim==0)               
                {
                    bound=1;
                    d_gmpivisc0[encodempivisc0(p,i,j,k,bound,dim)]=d_wtemp2[encode3p2_mpiu(p,i,j,k,var)];
                }

              

                if((j==0) && dim==1)              
                {              
                    bound=j;
                    d_gmpivisc1[encodempivisc1(p,i,j,k,bound,dim)]=d_wtemp2[encode3p2_mpiu(p,i,j,k,var)];
                }            
                 else if((( j==((p->n[1])+1)   ))  && dim==1)               
                {
                    bound=1;
                    d_gmpivisc1[encodempivisc1(p,i,j,k,bound,dim)]=d_wtemp2[encode3p2_mpiu(p,i,j,k,var)];
             
                }

       #ifdef USE_SAC_3D
               k=ii[2];
                if((k==0 ) && dim==2)              
                {              
                    bound=k;
                    d_gmpivisc2[encodempivisc2(p,i,j,k,bound,dim)]=d_wtemp2[encode3p2_mpiu(p,i,j,k,var)];
                }        
                 else if(((k==((p->n[2])+1)   ))  && dim==2)               
                {
                    bound=1;
                    d_gmpivisc2[encodempivisc2(p,i,j,k,bound,dim)]=d_wtemp2[encode3p2_mpiu(p,i,j,k,var)];
                }

     #endif
                               
                                
}

     __device__ __host__ void gputompiw(struct params *p,real *d_w,real *d_wmod,real *d_mpiw0,real *d_mpiwmod0,real *d_mpiw1,real *d_mpiwmod1,real *d_mpiw2,real *d_mpiwmod2,int *ii, int var, int dim)
    {
             int i,j,k,bound;
i=ii[0];
j=ii[1];
k=0;

/*  for(int field=rho;field<=rho ; field++)
if(  (p->ipe)==0  && ((p)->it)==1 && ( isnan(d_wmod[fencode3_mpiu(p,ii,field)])|| d_wmod[fencode3_mpiu(p,ii,field)]==0 ))
        { 
    				printf("nant %d %d %d %d %lg %lg \n",ii[0],ii[1],field,dim,d_wmod[fencode3_mpiu(p,ii,rho)],d_wmod[fencode3_mpiu(p,ii,field)] );
}*/
 
                if((i==0 || i==1) && dim==0)
                {              
                    bound=i;
                    d_mpiw0[encodempiw0(p,i,j,k,var,bound)]=d_w[encode3_mpiu(p,i+2,j,k,var)];
                   // d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]=d_wmod[encode3_mpiu(p,i+2,j,k,var)];



              
                }
                else if((( i>=((p->n[0])-2)   ))  && dim==0)               
                {
                    bound=2*(i==((p->n[0])-1))+(p->n[0])-i;
                    d_mpiw0[encodempiw0(p,i,j,k,var,bound)]=d_w[encode3_mpiu(p,i-2,j,k,var)];
                   // d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]=d_wmod[encode3_mpiu(p,i-2,j,k,var)];  




                }

              

                if((j==0 || j==1) && dim==1)              
                {              
                    bound=j;
                    d_mpiw1[encodempiw1(p,i,j,k,var,bound)]=d_w[encode3_mpiu(p,i,j+2,k,var)];
                   // d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)]=d_wmod[encode3_mpiu(p,i,j+2,k,var)];   








			/*  if( var==0 && (p)->ipe==3 && ((p)->it)==2 && bound==3)
			    {
				 //printf("ipe3 mpiwmod \n");
				 //for(int iii=0; iii<4*((p)->n[0]);iii++)
				     printf(" %lg %d \n",d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)],encodempiw1(p,i,j,k,var,bound));
				 //printf("\n");
			     }*/


           
                }            
                 else if((( j>=((p->n[1])-2)   ))  && dim==1)               
                {
                    bound=2*(j==((p->n[1])-1))+(p->n[1])-j;
                    d_mpiw1[encodempiw1(p,i,j,k,var,bound)]=d_w[encode3_mpiu(p,i,j-2,k,var)];
                   // d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)]=d_wmod[encode3_mpiu(p,i,j-2,k,var)];  

                          /*   if(p->ipe==0    && var==rho && dim==1 )
                             {
                             // for(int bound=0;bound<=1;bound++)
                                printf("mpiw0 %d %d %d %d %lg \n",dim,bound,i,j,d_mpiw0[encodempiw0(p,i,j,k,var,bound)]);
			     	//printf("mpiwmod0 %d %d %d %d %lg %lg\n\n",dim,bound,i,j,d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)],d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)]);
                             }*/



             
                }

       #ifdef USE_SAC_3D
               k=ii[2];
                if((k==0 || k==1) && dim==2)              
                {              
                    bound=k;
                    d_mpiw2[encodempiw2(p,i,j,k,var,bound)]=d_w[encode3_mpiu(p,i,j,k+2,var)];
                   // d_mpiwmod2[encodempiw2(p,i,j,k,var,bound)]=d_wmod[encode3_mpiu(p,i,j,k+2,var)];              
                }        
                 else if((( k>=((p->n[2])-2)   ))  && dim==2)               
                {
                    bound=2*(k==((p->n[2])-1))+(p->n[2])-k;
                    d_mpiw2[encodempiw2(p,i,j,k,var,bound)]=d_w[encode3_mpiu(p,i,j,k-2,var)];
                   // d_mpiwmod2[encodempiw2(p,i,j,k,var,bound)]=d_wmod[encode3_mpiu(p,i,j,k-2,var)];               
                }

     #endif



    /*if((p)->ipe==3 && ((p)->it)==2  && i==((p->n[0])-1) && j==((p->n[1])-1))
    {
         printf("ipe3 mpiwmod \n");
         for(int iii=0; iii<4*((p)->n[0]);iii++)
             printf(" %lg ",d_mpiwmod1[iii]);
         printf("\n");
     }*/

 
 
 }


     __device__ __host__ void gputompiwmod(struct params *p,real *d_w,real *d_wmod,real *d_mpiw0,real *d_mpiwmod0,real *d_mpiw1,real *d_mpiwmod1,real *d_mpiw2,real *d_mpiwmod2,int *ii, int var, int dim, int order)
    {
             int i,j,k,bound;
i=ii[0];
j=ii[1];
k=0;

 #ifdef USE_SAC_3D
	int ntot=((p->n[0]))*((p->n[1]))*((p->n[2]))*NVAR;
 #else
	int ntot=((p->n[0]))*((p->n[1]))*NVAR;
 #endif

/*  for(int field=rho;field<=rho ; field++)
if(  (p->ipe)==0  && ((p)->it)==1 && ( isnan(d_wmod[fencode3_mpiu(p,ii,field)])|| d_wmod[fencode3_mpiu(p,ii,field)]==0 ))
        { 
    				printf("nant %d %d %d %d %lg %lg \n",ii[0],ii[1],field,dim,d_wmod[fencode3_mpiu(p,ii,rho)],d_wmod[fencode3_mpiu(p,ii,field)] );
}*/
 
                if((i==0 || i==1) && dim==0)
                {              
                    bound=i;
                    d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]=d_wmod[(order*ntot)+encode3_mpiu(p,i+2,j,k,var)];
                   // d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]=d_wmod[encode3_mpiu(p,i+2,j,k,var)];

                            // if(p->ipe==0    && var==rho && dim==0 )
                             //{
                             // for(int bound=0;bound<=1;bound++)
                              //  printf("mpiw0 %d %d %d %d %lg \n",dim,bound,i,j,d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]);
			     	//printf("mpiwmod0 %d %d %d %d %lg %lg\n\n",dim,bound,i,j,d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)],d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)]);
                             //}
              
                }
                else if((( i>=((p->n[0])-2)   ))  && dim==0)               
                {
                    bound=2*(i==((p->n[0])-1))+(p->n[0])-i;
                    //d_mpiw0[encodempiw0(p,i,j,k,var,bound)]=d_w[encode3_mpiu(p,i-2,j,k,var)];
                    d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]=d_wmod[(order*ntot)+encode3_mpiu(p,i-2,j,k,var)];  

                             //if(p->ipe==0    && var==rho && dim==0 )
                             //{
                             // for(int bound=0;bound<=1;bound++)
                             //   printf("mpiw0 %d %d %d %d %lg \n",dim,bound,i,j,d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]);
			     	//printf("mpiwmod0 %d %d %d %d %lg %lg\n\n",dim,bound,i,j,d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)],d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)]);
                             //}





                }

              

                if((j==0 || j==1) && dim==1)              
                {              
                    bound=j;
                    //d_mpiw1[encodempiw1(p,i,j,k,var,bound)]=d_w[encode3_mpiu(p,i,j+2,k,var)];
                    d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)]=d_wmod[(order*ntot)+encode3_mpiu(p,i,j+2,k,var)];   






			/*  if( var==0 && (p)->ipe==3 && ((p)->it)==2 && bound==3)
			    {
				 //printf("ipe3 mpiwmod \n");
				 //for(int iii=0; iii<4*((p)->n[0]);iii++)
				     printf(" %lg %d \n",d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)],encodempiw1(p,i,j,k,var,bound));
				 //printf("\n");
			     }*/


           
                }            
                 else if((( j>=((p->n[1])-2)   ))  && dim==1)               
                {
                    bound=2*(j==((p->n[1])-1))+(p->n[1])-j;
                   // d_mpiw1[encodempiw1(p,i,j,k,var,bound)]=d_w[encode3_mpiu(p,i,j-2,k,var)];
                    d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)]=d_wmod[(order*ntot)+encode3_mpiu(p,i,j-2,k,var)];  





             
                }

       #ifdef USE_SAC_3D
               k=ii[2];
                if((k==0 || k==1) && dim==2)              
                {              
                    bound=k;
                   // d_mpiw2[encodempiw2(p,i,j,k,var,bound)]=d_w[encode3_mpiu(p,i,j,k+2,var)];
                    d_mpiwmod2[encodempiw2(p,i,j,k,var,bound)]=d_wmod[(order*ntot)+encode3_mpiu(p,i,j,k+2,var)];              
                }        
                 else if((( k>=((p->n[2])-2)   ))  && dim==2)               
                {
                    bound=2*(k==((p->n[2])-1))+(p->n[2])-k;
                    //d_mpiw2[encodempiw2(p,i,j,k,var,bound)]=d_w[encode3_mpiu(p,i,j,k-2,var)];
                    d_mpiwmod2[encodempiw2(p,i,j,k,var,bound)]=d_wmod[(order*ntot)+encode3_mpiu(p,i,j,k-2,var)];               
                }

     #endif



    /*if((p)->ipe==3 && ((p)->it)==2  && i==((p->n[0])-1) && j==((p->n[1])-1))
    {
         printf("ipe3 mpiwmod \n");
         for(int iii=0; iii<4*((p)->n[0]);iii++)
             printf(" %lg ",d_mpiwmod1[iii]);
         printf("\n");
     }*/

 
 
 }





     __device__ __host__ void gputompiwd(struct params *p,real *d_wd,real *d_mpiw0,real *d_mpiw1,real *d_mpiw2,int *ii, int var, int dim)
    {
             int i,j,k,bound;
i=ii[0];
j=ii[1];
k=0;
 
 
                if((i==0 || i==1) && dim==0)
                {              
                    bound=i;
                    d_mpiw0[encodempiw0(p,i,j,k,var,bound)]=d_wd[encode3_mpiu(p,i+2,j,k,var)];
                 }
                else if((( i>=((p->n[0])-2)   ))  && dim==0)               
                {
                    bound=2*(i==((p->n[0])-1))+(p->n[0])-i;
                    d_mpiw0[encodempiw0(p,i,j,k,var,bound)]=d_wd[encode3_mpiu(p,i-2,j,k,var)];
                  }

              

                if((j==0 || j==1) && dim==1)              
                {              
                    bound=j;
                    d_mpiw1[encodempiw1(p,i,j,k,var,bound)]=d_wd[encode3_mpiu(p,i,j+2,k,var)];
 

			/*  if( var==0 && (p)->ipe==3 && ((p)->it)==2 && bound==3)
			    {
				 //printf("ipe3 mpiwmod \n");
				 //for(int iii=0; iii<4*((p)->n[0]);iii++)
				     printf(" %lg %d \n",d_mpiwmod1[encodempiw1(p,i,j,k,var,bound)],encodempiw1(p,i,j,k,var,bound));
				 //printf("\n");
			     }*/


           
                }            
                 else if((( j>=((p->n[1])-2)   ))  && dim==1)               
                {
                    bound=2*(j==((p->n[1])-1))+(p->n[1])-j;
                    d_mpiw1[encodempiw1(p,i,j,k,var,bound)]=d_wd[encode3_mpiu(p,i,j-2,k,var)];
                 }

       #ifdef USE_SAC_3D
               k=ii[2];
                if((k==0 || k==1) && dim==2)              
                {              
                    bound=k;
                    d_mpiw2[encodempiw2(p,i,j,k,var,bound)]=d_wd[encode3_mpiu(p,i,j,k+2,var)];
                 }        
                 else if((( k>=((p->n[2])-2)   ))  && dim==2)               
                {
                    bound=2*(k==((p->n[2])-1))+(p->n[2])-k;
                    d_mpiw2[encodempiw2(p,i,j,k,var,bound)]=d_wd[encode3_mpiu(p,i,j,k-2,var)];
                 }

     #endif



    /*if((p)->ipe==3 && ((p)->it)==2  && i==((p->n[0])-1) && j==((p->n[1])-1))
    {
         printf("ipe3 mpiwmod \n");
         for(int iii=0; iii<4*((p)->n[0]);iii++)
             printf(" %lg ",d_mpiwmod1[iii]);
         printf("\n");
     }*/

 
 
 }






__global__ void  mpiwtogpu_parallel(struct params *p,real *d_w, real *d_wmod, real *d_mpiw0, real *d_mpiwmod0, real *d_mpiw1, real *d_mpiwmod1, real *d_mpiw2, real *d_mpiwmod2, int idir)
{

int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;

  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[0];
  real dx=p->dx[1];
                real val=0;
  
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

//if(iindex==0)
//         printf("in mpiwtogpu\n");

     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];
      //for(int dim=0; dim<NDIM;dim++)
           //for( f=rho; f<=b3; f++)
            for( f=rho; f<NVAR; f++)
     #else
     //for(int dim=0; dim<NDIM;dim++)
           //for( f=rho; f<=b2; f++)
	  for( f=rho; f<NVAR; f++)
     #endif     
         #ifdef USE_SAC_3D
           if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
         #else
           if(i<((p->n[0])) && j<((p->n[1])))
         #endif     
{      
    // if(iindex==0)
    //     printf("calling  mpiwtogpu %d %d\n",dim,f);

                 mpiwtogpu(p,d_w,d_wmod,d_mpiw0,d_mpiwmod0,d_mpiw1,d_mpiwmod1,d_mpiw2,d_mpiwmod2,iia,f,idir);

}


 __syncthreads();

           
               
}



__global__ void  mpiwmodtogpu_parallel(struct params *p,real *d_w, real *d_wmod, real *d_mpiw0, real *d_mpiwmod0, real *d_mpiw1, real *d_mpiwmod1, real *d_mpiw2, real *d_mpiwmod2, int idir, int order)
{

int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;

//int var,bound;

  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[0];
  real dx=p->dx[1];
                real val=0;
  
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

//if(iindex==0)
//         printf("in mpiwtogpu\n");

     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];
      //for(int dim=0; dim<NDIM;dim++)
           //for( f=rho; f<=b3; f++)
            for( f=rho; f<NVAR; f++)
     #else
     //for(int dim=0; dim<NDIM;dim++)
           //for( f=rho; f<=b2; f++)
	  for( f=rho; f<NVAR; f++)
     #endif     
         #ifdef USE_SAC_3D
           if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
         #else
           if(i<((p->n[0])) && j<((p->n[1])))
         #endif     
{      
    // if(iindex==0)
    //     printf("calling  mpiwtogpu %d %d\n",dim,f);

                    //bound=2*(i==((p->n[0])-1))+(p->n[0])-i;
                    /*bound=i;
                    var=f;
                    d_mpiwmod0[encodempiw0(p,i,j,k,var,bound)]=d_wmod[encode3_mpiu(p,i+2,j,k,var)];*/


                 //if( f==rho && (p->ipe)==3 && (iia[0]==0  || iia[0]==1) )
                      // if(idir==0)
		      // {
		       //d_mpiwmod0[encodempiw0(p,iia[0],iia[1],iia[2],f,iia[0])]=j+(p->ipe)*(1000);
    		 		//printf("nani0 %d %d %d  %lg  \n",p->ipe, iia[0],iia[1], d_mpiwmod0[encodempiw0(p,iia[0],iia[1],iia[2],f,iia[0])]); 
				//}


 mpiwmodtogpu(p,d_w,d_wmod,d_mpiw0,d_mpiwmod0,d_mpiw1,d_mpiwmod1,d_mpiw2,d_mpiwmod2,iia,f,idir,order);

                
		/*int bound;
                 if( f==rhob /*&& (p->ipe)==1  && (  d_wmod[fencode3_mpiu(p,iia,f)]==0 )*/ /*  && (j==0 || j==513)   && i>=0 && i<20)
                       if(idir==1)
                       {

			if(j==0 || j==1)
				bound=j;

                           if(j>=((p->n[1])-2))
 				bound=2*(j==((p->n[1])-1))+(p->n[1])-j;


    				printf("nani %d %d %d  %lg %lg %d  \n",p->ipe,iia[0],iia[1], d_mpiwmod1[encodempiw1(p,i,j,k,f,bound)],d_wmod[encode3_mpiu(p,i+2,j,k,f)],bound );
				
                       }*/


                 //if( f==rho && (p->ipe)==3 && (iia[0]==0  || iia[0]==1) )
                       //if(idir==0)
    		 		//printf("nani0 %d %d  %lg %lg \n",iia[0],iia[1], d_wmod[fencode3_mpiu(p,iia,rho)],d_wmod[fencode3_mpiu(p,iia,f)+dimp*NVAR] );

}


 __syncthreads();

           
               
}




__global__ void  mpiwdtogpu_parallel(struct params *p,  int dir, int var, real *d_wd, real *d_mpiw0,  real *d_mpiw1,  real *d_mpiw2, int idir)
{

int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;

  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[0];
  real dx=p->dx[1];
                real val=0;
  
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

//if(iindex==0)
//         printf("in mpiwtogpu\n");

     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];
      //for(int dim=0; dim<NDIM;dim++)
           //for( f=rho; f<=b3; f++)
            for( f=pos1; f<=delx3; f++)
     #else
     //for(int dim=0; dim<NDIM;dim++)
           //for( f=rho; f<=b2; f++)
	  for( f=pos1; f<=delx2; f++)
     #endif     
         #ifdef USE_SAC_3D
           if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
         #else
           if(i<((p->n[0])) && j<((p->n[1])))
         #endif     
{      
    // if(iindex==0)
    //     printf("calling  mpiwtogpu %d %d\n",dim,f);

                 mpiwdtogpu(p,d_wd,d_mpiw0,d_mpiw1,d_mpiw2,iia,f,idir);

}


 __syncthreads();

           
               
}



     __global__ void gputompiwmod_parallel(struct params *p,real *d_w,real *d_wmod,real *d_mpiw0,real *d_mpiwmod0,real *d_mpiw1,real *d_mpiwmod1,real *d_mpiw2,real *d_mpiwmod2,int order, int idir)
    {

 int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;
int dim;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[0];
  real dx=p->dx[1];
                real val=0;
  
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


     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];
      //for(dim=0; dim<NDIM;dim++)
           //for( f=rho; f<=b3; f++)
           for( f=rho; f<NVAR; f++)
     #else
           //for(dim=0; dim<NDIM;dim++)
           //for( f=rho; f<=b2; f++)
	   for( f=rho; f<NVAR; f++)
     #endif
             {
            
         #ifdef USE_SAC_3D
      if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
     #else
       if(i<((p->n[0])) && j<((p->n[1])))
     #endif           
	{

 

                  gputompiwmod(p,d_w,d_wmod,d_mpiw0,d_mpiwmod0,d_mpiw1,d_mpiwmod1,d_mpiw2,d_mpiwmod2,iia,f,idir,order);

                             /*if( f==rho && idir==0 )
                                          if((i==0 || i==1) )
						{
                                 		
							d_mpiwmod0[encodempiw0(p,i,j,k,f,i)]=1000*(p->ipe)+j;
							//printf("mpiw0 %d %d %d %d %d %lg %lg\n",p->ipe,idir,i,iia[0],iia[1],d_mpiwmod0[encodempiw0(p,i,j,k,f,i)],d_mpiwmod1[encodempiw1(p,i,j,k,f,i)]);
                                                  }*/


                            // if(p->ipe==3    && f==rho && idir==0 )
                                       //   if((i==0 || i==1) )
						//{
                                 		
							//d_mpiwmod0[encodempiw0(p,i,j,k,f,i)]=1000*(p->ipe)+j;
							//printf("mpiw0 %d %d %d %d %d %lg %lg\n",p->ipe,idir,i,iia[0],iia[1],d_mpiwmod0[encodempiw0(p,i,j,k,f,i)],d_mpiwmod1[encodempiw1(p,i,j,k,f,i)]);
                                                 // }

	}

               }





 __syncthreads();

}

     __global__ void gputompiw_parallel(struct params *p,real *d_w,real *d_wmod,real *d_mpiw0,real *d_mpiwmod0,real *d_mpiw1,real *d_mpiwmod1,real *d_mpiw2,real *d_mpiwmod2,int order, int idir)
    {

 int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;
int dim;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[0];
  real dx=p->dx[1];
                real val=0;
  
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


     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];
      //for(dim=0; dim<NDIM;dim++)
           //for( f=rho; f<=b3; f++)
           for( f=rho; f<NVAR; f++)
     #else
           //for(dim=0; dim<NDIM;dim++)
           //for( f=rho; f<=b2; f++)
	   for( f=rho; f<NVAR; f++)
     #endif
             {
            
         #ifdef USE_SAC_3D
      if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
     #else
       if(i<((p->n[0])) && j<((p->n[1])))
     #endif           
	{

 

                  gputompiw(p,d_w,d_wmod,d_mpiw0,d_mpiwmod0,d_mpiw1,d_mpiwmod1,d_mpiw2,d_mpiwmod2,iia,f,idir);

                             //if(p->ipe==0    && f==rho && idir==0 )
                             // for(int bound=0;bound<=1;bound++)
                             //   printf("mpiw0 %d %d %d %d %lg %lg\n",idir,bound,iia[0],iia[1],d_mpiw0[encodempiw0(p,i,j,k,f,bound)],d_mpiw1[encodempiw1(p,i,j,k,f,bound)]);

	}

               }





 __syncthreads();

}


     __global__ void gputompiwd_parallel(struct params *p,real *d_wd,real *d_mpiw0,real *d_mpiw1,real *d_mpiw2,int order, int idir)
    {

 int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;
int dim;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[0];
  real dx=p->dx[1];
                real val=0;
  
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


     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];
      //for(dim=0; dim<NDIM;dim++)
           //for( f=rho; f<=b3; f++)
           for( f=pos1; f<=delx3; f++)
     #else
        //   for(dim=0; dim<NDIM;dim++)
           //for( f=rho; f<=b2; f++)
	   for( f=pos1; f<=delx2; f++)
     #endif
             {
            
         #ifdef USE_SAC_3D
      if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
     #else
       if(i<((p->n[0])) && j<((p->n[1])))
     #endif           
	{

 

                  gputompiwd(p,d_wd,d_mpiw0,d_mpiw1,d_mpiw2,iia,f,idir);


	}

               }





 __syncthreads();

}


     __global__ void gputompivisc_parallel(struct params *p,real *d_wtemp2,real *d_gmpivisc0,real *d_gmpivisc1,real *d_gmpivisc2)
     {
               
  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;
int dim;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[0];
  real dx=p->dx[1];
                real val=0;
  
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
   kp=iindex/((nj+2)*(ni+2));
   jp=(iindex-(kp*((nj+2)*(ni+2))))/(ni+2);
   ip=iindex-(kp*(nj+2)*(ni+2))-(jp*(ni+2));
#else
    jp=iindex/(ni+2);
   ip=iindex-(jp*(ni+2));
#endif     


//int shift=order*NVAR*dimp;


     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];

     #else

     #endif
           for(dim=0; dim<NDIM;dim++)
             {
            
         #ifdef USE_SAC_3D
      if(i<(((p->n[0])+2)) && j<(((p->n[1])+2))  && k<(((p->n[2])+2)))
     #else
       if(i<(((p->n[0])+2)) && j<(((p->n[1])+2)))
     #endif           
	{

 

                  gputompivisc(p,d_wtemp2,d_gmpivisc0,d_gmpivisc1,d_gmpivisc2,iia,dim);

	}

               }

 __syncthreads();
              
               }    
     
     
    __global__ void  mpivisctogpu_parallel(struct params *p,real *d_wtemp2,real *d_gmpivisc0,real *d_gmpivisc1,real *d_gmpivisc2)
    {
               
  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;
int dim;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[0];
  real dx=p->dx[1];
                real val=0;
  
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
   kp=iindex/((nj+2)*(ni+2));
   jp=(iindex-(kp*((nj+2)*(ni+2))))/(ni+2);
   ip=iindex-(kp*(nj+2)*(ni+2))-(jp*(ni+2));
#else
    jp=iindex/(ni+2);
   ip=iindex-(jp*(ni+2));
#endif     


//int shift=order*NVAR*dimp;


     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];

     #else

     #endif
           for(dim=0; dim<NDIM;dim++)
             {
            
         #ifdef USE_SAC_3D
      if(i<(((p->n[0])+2)) && j<(((p->n[1])+2))  && k<(((p->n[2])+2)))
     #else
       if(i<(((p->n[0])+2)) && j<(((p->n[1])+2)))
     #endif           
	{

 

                  mpivisctogpu(p,d_wtemp2,d_gmpivisc0,d_gmpivisc1,d_gmpivisc2,iia,dim);

	}

               }

 __syncthreads();
               
               
}

#endif



/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_mpiu(char *label)
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



#ifdef USE_MULTIGPU




//prepare data buffers used to copy data between gpu and cpu
//this will update only the ghost cells transferred between the CPU's

int cuinitmgpurbuffers(struct params **p,    
real **d_gmpiviscr0,    
real **d_gmpiviscr1,    
real **d_gmpiviscr2,   
real **d_gmpiwr0, 
real **d_gmpiwmodr0,   
real **d_gmpiwr1, 
real **d_gmpiwmodr1,   
real **d_gmpiwr2, 
real **d_gmpiwmodr2)
{

  int szw,  szvisc0,szvisc1,szvisc2,szw0,szw1,szw2;
  #ifdef USE_SAC
  //real *dt;
  
  szw=4*(  ((*p)->n[1])  +  ((*p)->n[0])   );
  szw0=4*NDERV*(  ((*p)->n[1])     );
  szw1=4*NDERV*(  ((*p)->n[0])     );

  szvisc0=4*(  (((*p)->n[1])+2 )   );
  szvisc1=4*(    (((*p)->n[0]) +2 )  );

 //dt=(real *)calloc( NTEMP2*(((*p)->n[0])+2)* (((*p)->n[1])+2),sizeof(real));

  #endif
  #ifdef USE_SAC_3D
  
  szw=4*NDERV*(  ((*p)->n[1])*((*p)->n[2])  +  ((*p)->n[0])*((*p)->n[2])  +  ((*p)->n[0])*((*p)->n[1])  );
  szw0=4*NDERV*(  ((*p)->n[1])*((*p)->n[2])    );
  szw1=4*NDERV*(    ((*p)->n[0])*((*p)->n[2])   );
  szw2=4*NDERV*(    ((*p)->n[0])*((*p)->n[1])  );



  szvisc0=4*(  (((*p)->n[1])+2)*(((*p)->n[2])+2)  ); 
  szvisc1=4*(   (((*p)->n[0])+2)*(((*p)->n[2])+2)    );    
  szvisc2=4*(  (((*p)->n[1])+2)*(((*p)->n[2])+2)   );    

   
  //dt=(real *)calloc( NTEMP2*(((*p)->n[0])+2)* (((*p)->n[1])+2)* (((*p)->n[2])+2),sizeof(real));
  #endif

  	  cudaMalloc((void**)d_gmpiwmodr0, szw0*sizeof(real));
	  cudaMalloc((void**)d_gmpiwr0, szw0*sizeof(real));
	  cudaMalloc((void**)d_gmpiwmodr1, szw1*sizeof(real));
	  cudaMalloc((void**)d_gmpiwr1, szw1*sizeof(real));

  #ifdef USE_SAC_3D  
	  cudaMalloc((void**)d_gmpiwmodr2, szw2*sizeof(real));
	  cudaMalloc((void**)d_gmpiwr2, szw2*sizeof(real));
          cudaMalloc((void**)d_gmpiviscr2, szvisc2*sizeof(real));
  #else

          cudaMalloc((void**)d_gmpiviscr2, sizeof(real));
  #endif
          cudaMalloc((void**)d_gmpiviscr0, szvisc0*sizeof(real));
          cudaMalloc((void**)d_gmpiviscr1, szvisc1*sizeof(real));
  return 0;



}





//prepare data buffers used to copy data between gpu and cpu
//this will update only the ghost cells transferred between the CPU's


int cuinitmgpubuffers(struct params **p,real **w, real **wmod, real **temp2, real **gmpivisc0, real **gmpivisc1, real **gmpivisc2,   real **gmpiw0, real **gmpiwmod0,   real **gmpiw1, real **gmpiwmod1,   real **gmpiw2, real **gmpiwmod2, struct params **d_p,   real **d_w, real **d_wmod,real **d_wtemp2,    real **d_gmpivisc0,    real **d_gmpivisc1,    real **d_gmpivisc2,   real **d_gmpiw0, real **d_gmpiwmod0,   real **d_gmpiw1, real **d_gmpiwmod1,   real **d_gmpiw2, real **d_gmpiwmod2)
{

  int szw,  szvisc0,szvisc1,szvisc2,szw0,szw1,szw2;
  #ifdef USE_SAC
  //real *dt;
  
  szw=4*(  ((*p)->n[1])  +  ((*p)->n[0])   );
  szw0=4*NDERV*(  ((*p)->n[1])     );
  szw1=4*NDERV*(  ((*p)->n[0])     );

  szvisc0=4*(  (((*p)->n[1])+2 )   );
  szvisc1=4*(    (((*p)->n[0]) +2 )  );

 //dt=(real *)calloc( NTEMP2*(((*p)->n[0])+2)* (((*p)->n[1])+2),sizeof(real));

  #endif
  #ifdef USE_SAC_3D
  
  szw=4*NDERV*(  ((*p)->n[1])*((*p)->n[2])  +  ((*p)->n[0])*((*p)->n[2])  +  ((*p)->n[0])*((*p)->n[1])  );
  szw0=4*NDERV*(  ((*p)->n[1])*((*p)->n[2])    );
  szw1=4*NDERV*(    ((*p)->n[0])*((*p)->n[2])   );
  szw2=4*NDERV*(    ((*p)->n[0])*((*p)->n[1])  );



  szvisc0=4*(  (((*p)->n[1])+2)*(((*p)->n[2])+2)  ); 
  szvisc1=4*(   (((*p)->n[0])+2)*(((*p)->n[2])+2)    );    
  szvisc2=4*(  (((*p)->n[1])+2)*(((*p)->n[2])+2)   );    

   
  //dt=(real *)calloc( NTEMP2*(((*p)->n[0])+2)* (((*p)->n[1])+2)* (((*p)->n[2])+2),sizeof(real));
  #endif






  //temp2=&dt;
  //gmpiwmod=(real **)malloc(szw*sizeof(real));
  //gmpiw=(real **)malloc(szw*sizeof(real));

  //gmpiwmod0=(real **)malloc(szw0*sizeof(real));
  //gmpiw0=(real **)malloc(szw0*sizeof(real));
  //gmpiwmod1=(real **)malloc(szw1*sizeof(real));
  //gmpiw1=(real **)malloc(szw1*sizeof(real));

  #ifdef USE_SAC_3D
	//  gmpiwmod2=(real **)malloc(szw2*sizeof(real));
	//  gmpiw2=(real **)malloc(szw2*sizeof(real));
  #endif

  //gmpivisc=(real **)malloc(szvisc*sizeof(real));
	//  cudaMalloc((void**)d_gmpiwmod, szw*sizeof(real));
	//  cudaMalloc((void**)d_gmpiw, szw*sizeof(real));


  	  cudaMalloc((void**)d_gmpiwmod0, szw0*sizeof(real));
	  cudaMalloc((void**)d_gmpiw0, szw0*sizeof(real));
	  cudaMalloc((void**)d_gmpiwmod1, szw1*sizeof(real));
	  cudaMalloc((void**)d_gmpiw1, szw1*sizeof(real));

  #ifdef USE_SAC_3D  
	  cudaMalloc((void**)d_gmpiwmod2, szw2*sizeof(real));
	  cudaMalloc((void**)d_gmpiw2, szw2*sizeof(real));
          cudaMalloc((void**)d_gmpivisc2, szvisc2*sizeof(real));
  #else

          cudaMalloc((void**)d_gmpivisc2, sizeof(real));
  #endif
          cudaMalloc((void**)d_gmpivisc0, szvisc0*sizeof(real));
          cudaMalloc((void**)d_gmpivisc1, szvisc1*sizeof(real));
  return 0;
}

//copy gpu memory data to mpi send buffer for w and wmod
//just update the edges of w and wmod with values copied from gmpiw, gmpiwmod and gmpivisc
int cucopywtompiwmod(struct params **p,real **w, real **wmod,    real **gmpiw0, real **gmpiwmod0,    real **gmpiw1, real **gmpiwmod1,    real **gmpiw2, real **gmpiwmod2, struct params **d_p  ,real **d_w, real **d_wmod,   real **d_gmpiw0, real **d_gmpiwmod0,   real **d_gmpiw1, real **d_gmpiwmod1,   real **d_gmpiw2, real **d_gmpiwmod2, int order, int idir)
{
     int i1,i2,i3;
     int ii[NDIM];
     int var,dim,bound;
     int szw0,szw1,szw2;

     int szbuf;
     int dimp=(((*p)->n[0]))*(((*p)->n[1]));
     
     
     i3=0;
     #ifdef USE_SAC_3D  
       dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
     #endif 
     int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;

     szbuf=2*2*( ((*p)->n[0])+((*p)->n[1]));
     #ifdef USE_SAC_3D
     szbuf=2*2*( ((*p)->n[0])*((*p)->n[1])+ ((*p)->n[0])*((*p)->n[2]) + ((*p)->n[1])*((*p)->n[2])        );
     #endif


  #ifdef USE_SAC
  
  szw0=4*NVAR*(  ((*p)->n[1])     );
  szw1=4*NVAR*(  ((*p)->n[0])     );

  #endif
  #ifdef USE_SAC_3D
  
   szw0=4*NVAR*(  ((*p)->n[1])*((*p)->n[2])    );
  szw1=4*NVAR*(    ((*p)->n[0])*((*p)->n[2])   );
  szw2=4*NVAR*(    ((*p)->n[0])*((*p)->n[1])  );

  #endif

    //real **d_tgmpiw0;
    //real **tgmpiw0=(real **)malloc(szw0*sizeof(real));
    //cudaMalloc((void**)d_tgmpiw0, szw0*sizeof(real));
    // for(var=0; var<NVAR; var++)
    //   for(dim=0;dim<NDIM;dim++)
     gputompiwmod_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wmod,*d_gmpiw0,*d_gmpiwmod0,*d_gmpiw1,*d_gmpiwmod1,*d_gmpiw2,*d_gmpiwmod2,order,idir);

 
 #ifdef USE_GPUDIRECT
    // printf("using gpudirect in1\n");
     cudaThreadSynchronize();

#else
 
     
     cudaThreadSynchronize();

if(idir==0)
{
     cudaMemcpy(*gmpiwmod0, *d_gmpiwmod0, szw0*sizeof(real), cudaMemcpyDeviceToHost);

// cudaThreadSynchronize();

   //   cudaMemcpy(*gmpiw0, *d_gmpiw0, szw0*sizeof(real), cudaMemcpyDeviceToHost);




}

if(idir==1)
{
     cudaMemcpy(*gmpiwmod1, *d_gmpiwmod1, szw1*sizeof(real), cudaMemcpyDeviceToHost);
   //  cudaMemcpy(*gmpiw1, *d_gmpiw1, szw1*sizeof(real), cudaMemcpyDeviceToHost);
}
 

//struct params *tp;     







    

   #ifdef USE_SAC_3D
if(idir==2)
{
     cudaMemcpy(*gmpiwmod2, *d_gmpiwmod2, szw2*sizeof(real), cudaMemcpyDeviceToHost);
     //cudaMemcpy(*gmpiw2, *d_gmpiw2, szw2*sizeof(real), cudaMemcpyDeviceToHost);
}
   #endif 
   
   
   #endif

cudaThreadSynchronize();
}


//copy gpu memory data to mpi send buffer for w and wmod
//just update the edges of w and wmod with values copied from gmpiw, gmpiwmod and gmpivisc
int cucopywtompiw(struct params **p,real **w, real **wmod,    real **gmpiw0, real **gmpiwmod0,    real **gmpiw1, real **gmpiwmod1,    real **gmpiw2, real **gmpiwmod2, struct params **d_p  ,real **d_w, real **d_wmod,   real **d_gmpiw0, real **d_gmpiwmod0,   real **d_gmpiw1, real **d_gmpiwmod1,   real **d_gmpiw2, real **d_gmpiwmod2, int order, int idir)
{
     int i1,i2,i3;
     int ii[NDIM];
     int var,dim,bound;
     int szw0,szw1,szw2;

     int szbuf;
     int dimp=(((*p)->n[0]))*(((*p)->n[1]));
     
     
     i3=0;
     #ifdef USE_SAC_3D  
       dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
     #endif 
     int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;

     szbuf=2*2*( ((*p)->n[0])+((*p)->n[1]));
     #ifdef USE_SAC_3D
     szbuf=2*2*( ((*p)->n[0])*((*p)->n[1])+ ((*p)->n[0])*((*p)->n[2]) + ((*p)->n[1])*((*p)->n[2])        );
     #endif


  #ifdef USE_SAC
  
  szw0=4*NVAR*(  ((*p)->n[1])     );
  szw1=4*NVAR*(  ((*p)->n[0])     );

  #endif
  #ifdef USE_SAC_3D
  
   szw0=4*NVAR*(  ((*p)->n[1])*((*p)->n[2])    );
  szw1=4*NVAR*(    ((*p)->n[0])*((*p)->n[2])   );
  szw2=4*NVAR*(    ((*p)->n[0])*((*p)->n[1])  );

  #endif

    //real **d_tgmpiw0;
    //real **tgmpiw0=(real **)malloc(szw0*sizeof(real));
    //cudaMalloc((void**)d_tgmpiw0, szw0*sizeof(real));
    // for(var=0; var<NVAR; var++)
    //   for(dim=0;dim<NDIM;dim++)
     gputompiw_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wmod,*d_gmpiw0,*d_gmpiwmod0,*d_gmpiw1,*d_gmpiwmod1,*d_gmpiw2,*d_gmpiwmod2,order,idir);


#ifdef USE_GPUDIRECT
     
     cudaThreadSynchronize();

#else



     
     cudaThreadSynchronize();

if(idir==0)
{
     //cudaMemcpy(*gmpiwmod0, *d_gmpiwmod0, szw0*sizeof(real), cudaMemcpyDeviceToHost);


      cudaMemcpy(*gmpiw0, *d_gmpiw0, szw0*sizeof(real), cudaMemcpyDeviceToHost);
}

if(idir==1)
{
     //cudaMemcpy(*gmpiwmod1, *d_gmpiwmod1, szw1*sizeof(real), cudaMemcpyDeviceToHost);
     cudaMemcpy(*gmpiw1, *d_gmpiw1, szw1*sizeof(real), cudaMemcpyDeviceToHost);
}
      
    

   #ifdef USE_SAC_3D
if(idir==2)
{
    // cudaMemcpy(*gmpiwmod2, *d_gmpiwmod2, szw2*sizeof(real), cudaMemcpyDeviceToHost);
     cudaMemcpy(*gmpiw2, *d_gmpiw2, szw2*sizeof(real), cudaMemcpyDeviceToHost);
}
   #endif 
   
   #endif

cudaThreadSynchronize();


 /*if(((*p)->ipe)==3  && ((*p)->it)==2)
{


       printf("%d %d \n",szw0,szw1);

        for(i1=0;i1<(((*p)->n[0]));i1++ )
                  {
                       ii[0]=i1;
                       ii[1]=0;
                       bound=0;
                       var=0;
                            printf(" %d %d %d %lg %d \n",i1,i2,bound,(*gmpiwmod1)[encodempiw1(*p,i1,i2,i3,var,bound)],encodempiw1(*p,i1,i2,i3,var,bound));                                        
                     ;//  (*wmod)[fencode3_mpiu(*p,ii,var)]=(*gmpiwmod1)[encodempiw1(*p,i1,i2,i3,var,bound)];              
                     ;//  (*w)[fencode3_mpiu(*p,ii,var)]=(*gmpiw1)[encodempiw1(*p,i1,i2,i3,var,bound)];



                  }
         ;// printf("\n");
}*/









//encodempiw1(p,i,j,k,var,bound)]



   //free(tgmpiw0);
   //cudaFree(*d_tgmpiw0);
//gmpiw behaving OK but cannot display or access any of the gmpiwmod variables!
//printf("%f\n",(*gmpiwmod)[0]);
     
//encodempiw (struct params *dp,int ix, int iy, int iz, int field,int bound,int dim)
     //copy data to correct area in w and wmod
   /*  for(var=0; var<NVAR; var++)
       for(dim=0;dim<NDIM;dim++) 
         for(bound=0;bound<4;bound++)
         {
            switch(dim)
            {
                       case 0:
            #ifdef USE_SAC_3D
         i1=bound*(bound<2)+(((*p)->n[0])-(bound-1))*(bound>1);
         for(i2=0;i2<(((*p)->n[1]));i2++ )
                  for(i3=0;i3<(((*p)->n[2]));i3++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                       ii[2]=i3;                                                                     
                       (*wmod)[fencode3_mpiu(*p,ii,var)]=(*gmpiwmod0)[encodempiw0(*p,i1,i2,i3,var,bound)];              
                       (*w)[fencode3_mpiu(*p,ii,var)]=(*gmpiw0)[encodempiw0(*p,i1,i2,i3,var,bound)];
                  }
            #else
         ii[2]=0;
         i1=bound*(bound<2)+(((*p)->n[0])-(bound-1))*(bound>1);
         for(i2=0;i2<(((*p)->n[1]));i2++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                      


		
                       (*wmod)[fencode3_mpiu(*p,ii,var)]=(*gmpiwmod0)[encodempiw0(*p,i1,i2,i3,var,bound)];  
                      //if(encodempiw(*p,i1,i2,i3,var,bound,dim)<10)  
                      if(var==5 && ((*p)->ipe)==0)                        
                        printf(" %d %d %d %d actual %d  mpi data%d %g\n",i1,i2,bound,dim,var,encodempiw0(*p,i1,i2,i3,var,bound),(*gmpiwmod0)[encodempiw0(*p,i1,i2,i3,var,bound)]);

                     ;// if(encodempiw(*p,i1,i2,i3,var,bound,dim)<10239 )
                       (*w)[fencode3_mpiu(*p,ii,var)]=(*gmpiw0)[encodempiw0(*p,i1,i2,i3,var,bound)];
                                                                     
                      // *(wmod+encode3_mpiu(*p,ii,var))=*(gmpiwmod0+encodempiw0(*p,i1,i2,i3,var,bound));              
                      // (*w)[encode3_mpiu(*p,ii,var)]=(*gmpiw0)[encodempiw0(*p,i1,i2,i3,var,bound)];

                      /* if(var==4  && ((*p)->ipe)==1)
                       {
				(*wmod)[fencode3_mpiu(*p,ii,var)]=0.5;
				(*w)[fencode3_mpiu(*p,ii,var)]=0.5;
                       }*/


              /*    }            
            
            #endif
                       
                       break;   
                       case 1:
            #ifdef USE_SAC_3D
         i2=bound*(bound<2)+(((*p)->n[1])-(bound-1))*(bound>1);
         for(i1=0;i1<(((*p)->n[0]));i1++ )
                  for(i3=0;i3<(((*p)->n[2]));i3++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                       ii[2]=i3;                                                                     
                       (*wmod)[fencode3_mpiu(*p,ii,var)]=(*gmpiwmod1)[encodempiw1(*p,i1,i2,i3,var,bound)];              
                       (*w)[fencode3_mpiu(*p,ii,var)]=(*gmpiw1)[encodempiw1(*p,i1,i2,i3,var,bound)];
                  }

            #else
         ii[2]=0;
         i2=bound*(bound<2)+(   ((*p)->n[1])-(bound-1)   )*(bound>1);
         for(i1=0;i1<(((*p)->n[0]));i1++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                                                                     
                     ;//  (*wmod)[fencode3_mpiu(*p,ii,var)]=(*gmpiwmod1)[encodempiw1(*p,i1,i2,i3,var,bound)];              
                     ;//  (*w)[fencode3_mpiu(*p,ii,var)]=(*gmpiw1)[encodempiw1(*p,i1,i2,i3,var,bound)];



                  }
            
            
            #endif
                       
                       break; 
            #ifdef USE_SAC_3D
                       case 2:
         i3=bound*(bound<2)+( ((*p)->n[2])-(bound-1) )*(bound>1);
         for(i1=0;i1<(((*p)->n[0]));i1++ )
                  for(i2=0;i2<(((*p)->n[1]));i2++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                       ii[2]=i3;                                                                     
                       (*wmod)[fencode3_mpiu(*p,ii,var)]=(*gmpiwmod2)[encodempiw2(*p,i1,i2,i3,var,bound)];              
                       (*w)[fencode3_mpiu(*p,ii,var)]=(*gmpiw2)[encodempiw2(*p,i1,i2,i3,var,bound)];
                  }                            
                       break;                       
            #endif             
             }
                                     
         }    */

}




int cucopywdtompiwd(struct params **p,real **wd,    real **gmpiw0,    real **gmpiw1,    real **gmpiw2, struct params **d_p  ,real **d_wd,   real **d_gmpiw0,   real **d_gmpiw1,   real **d_gmpiw2, int order, int idir)
{
     int i1,i2,i3;
     int ii[NDIM];
     int var,dim,bound;
     int szw0,szw1,szw2;

     int szbuf;
     int dimp=(((*p)->n[0]))*(((*p)->n[1]));
     
     
     i3=0;
     #ifdef USE_SAC_3D  
       dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
     #endif 
     int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;

     szbuf=2*2*( ((*p)->n[0])+((*p)->n[1]));
     #ifdef USE_SAC_3D
     szbuf=2*2*( ((*p)->n[0])*((*p)->n[1])+ ((*p)->n[0])*((*p)->n[2]) + ((*p)->n[1])*((*p)->n[2])        );
     #endif


  #ifdef USE_SAC
  
  szw0=4*NDERV*(  ((*p)->n[1])     );
  szw1=4*NDERV*(  ((*p)->n[0])     );

  #endif
  #ifdef USE_SAC_3D
  
   szw0=4*NDERV*(  ((*p)->n[1])*((*p)->n[2])    );
  szw1=4*NDERV*(    ((*p)->n[0])*((*p)->n[2])   );
  szw2=4*NDERV*(    ((*p)->n[0])*((*p)->n[1])  );

  #endif

    //real **d_tgmpiw0;
    //real **tgmpiw0=(real **)malloc(szw0*sizeof(real));
    //cudaMalloc((void**)d_tgmpiw0, szw0*sizeof(real));
    // for(var=0; var<NVAR; var++)
    //   for(dim=0;dim<NDIM;dim++)
     gputompiwd_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_wd,*d_gmpiw0,*d_gmpiw1,*d_gmpiw2,order, idir);

;//#ifdef USE_GPUDIRECT
     
;//     cudaThreadSynchronize();

;//#else


     
     cudaThreadSynchronize();
if(idir==0)
      cudaMemcpy(*gmpiw0, *d_gmpiw0, szw0*sizeof(real), cudaMemcpyDeviceToHost);

if(idir==1)
     cudaMemcpy(*gmpiw1, *d_gmpiw1, szw1*sizeof(real), cudaMemcpyDeviceToHost);

      
    

   #ifdef USE_SAC_3D
if(idir==2)
      cudaMemcpy(*gmpiw2, *d_gmpiw2, szw2*sizeof(real), cudaMemcpyDeviceToHost);
   #endif 

cudaThreadSynchronize();
;//#endif

 /*if(((*p)->ipe)==3  && ((*p)->it)==2)
{


       printf("%d %d \n",szw0,szw1);

        for(i1=0;i1<(((*p)->n[0]));i1++ )
                  {
                       ii[0]=i1;
                       ii[1]=0;
                       bound=0;
                       var=0;
                            printf(" %d %d %d %lg %d \n",i1,i2,bound,(*gmpiwmod1)[encodempiw1(*p,i1,i2,i3,var,bound)],encodempiw1(*p,i1,i2,i3,var,bound));                                        
                     ;//  (*wmod)[fencode3_mpiu(*p,ii,var)]=(*gmpiwmod1)[encodempiw1(*p,i1,i2,i3,var,bound)];              
                     ;//  (*w)[fencode3_mpiu(*p,ii,var)]=(*gmpiw1)[encodempiw1(*p,i1,i2,i3,var,bound)];



                  }
         ;// printf("\n");
}*/





}




//copy mpi recv buffer to gpu memory     
int cucopywfrommpiw(struct params **p,real **w, real **wmod,    real **gmpiw0, real **gmpiwmod0,    real **gmpiw1, real **gmpiwmod1,    real **gmpiw2, real **gmpiwmod2, struct params **d_p  ,real **d_w, real **d_wmod,   real **d_gmpiw0, real **d_gmpiwmod0,   real **d_gmpiw1, real **d_gmpiwmod1,   real **d_gmpiw2, real **d_gmpiwmod2, int order, int idir)
{
       int i1,i2,i3;
     int ii[NDIM];
     int var,dim,bound;     
       int szbuf;
     int szw0,szw1,szw2;

  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   
 #ifdef USE_SAC_3D  
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif      
     szbuf=2*2*( ((*p)->n[0])+((*p)->n[1]));
     #ifdef USE_SAC_3D
     szbuf=2*2*( ((*p)->n[0])*((*p)->n[1])+ ((*p)->n[0])*((*p)->n[2]) + ((*p)->n[1])*((*p)->n[2])        );
     #endif
        int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;


  #ifdef USE_SAC
  
  szw0=4*NVAR*(  ((*p)->n[1])     );
  szw1=4*NVAR*(  ((*p)->n[0])     );

  #endif
  #ifdef USE_SAC_3D
  
   szw0=4*NVAR*(  ((*p)->n[1])*((*p)->n[2])    );
  szw1=4*NVAR*(    ((*p)->n[0])*((*p)->n[2])   );
  szw2=4*NVAR*(    ((*p)->n[0])*((*p)->n[1])  );

  #endif




      //copy data from w and wmod to correct gmpiw and gmpiwmod

//encodempiw (struct params *dp,int ix, int iy, int iz, int field,int bound,int dim)
     //copy data to correct area in w and wmod
   /*  for(var=0; var<NVAR; var++)
       for(dim=0;dim<NDIM;dim++) 
         for(bound=0;bound<4;bound++)
         {
            switch(dim)
            {
                       case 0:
            #ifdef USE_SAC_3D
         i1=bound*(bound<2)+(((*p)->n[0])-(bound-1))*(bound>1);
         for(i2=0;i2<(((*p)->n[1]));i2++ )
                  for(i3=0;i3<(((*p)->n[2]));i3++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                       ii[2]=i3;                                                                     
                       (*gmpiwmod0)[encodempiw0(*p,i1,i2,i3,var,bound)]=(*wmod)[fencode3_mpiu(*p,ii,var)];              
                       (*gmpiw0)[encodempiw0(*p,i1,i2,i3,var,bound)]=(*w)[fencode3_mpiu(*p,ii,var)];
                  }
            #else
         ii[2]=0;
         i1=bound*(bound<2)+(((*p)->n[0])-(bound-1))*(bound>1);
         for(i2=0;i2<(((*p)->n[1]));i2++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                       (*gmpiwmod0)[encodempiw0(*p,i1,i2,i3,var,bound)]=(*wmod)[fencode3_mpiu(*p,ii,var)];  

                 //   if(var==4 && ((*p)->ipe)==1) 
                                 
                //        printf(" %d %d %d %d actual %d  mpi data%d %g %g\n",i1,i2,bound,dim,var,encodempiw0(*p,i1,i2,i3,var,bound),(*gmpiwmod0)[encodempiw0(*p,i1,i2,i3,var,bound)],  (*wmod)[fencode3_mpiu(*p,ii,var)] );
            
                       (*gmpiw0)[encodempiw0(*p,i1,i2,i3,var,bound)]=(*w)[fencode3_mpiu(*p,ii,var)];

                      //  if(var==4  && ((*p)->ipe)==1)
                     //  {
			//	(*gmpiwmod0)[encodempiw0(*p,i1,i2,i3,var,bound)]=0.5;
			//	(*gmpiw0)[encodempiw0(*p,i1,i2,i3,var,bound)]=0.5;
                      // }

                  }            
            
            #endif
                       
                       break;   
                       case 1:
            #ifdef USE_SAC_3D
         i2=bound*(bound<2)+(((*p)->n[1])-(bound-1))*(bound>1);
         for(i1=0;i1<(((*p)->n[0]));i1++ )
                  for(i3=0;i3<(((*p)->n[2]));i3++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                       ii[2]=i3;  

                       (*gmpiwmod1)[encodempiw1(*p,i1,i2,i3,var,bound)]=(*wmod)[fencode3_mpiu(*p,ii,var)];              
                       (*gmpiw1)[encodempiw1(*p,i1,i2,i3,var,bound)]=(*w)[fencode3_mpiu(*p,ii,var)];

                  }

            #else
         ii[2]=0;
         i2=bound*(bound<2)+(   ((*p)->n[1])-(bound-1)   )*(bound>1);
         for(i1=0;i1<(((*p)->n[0]));i1++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                      (*gmpiwmod1)[encodempiw1(*p,i1,i2,i3,var,bound)]=(*wmod)[fencode3_mpiu(*p,ii,var)];              
                       (*gmpiw1)[encodempiw1(*p,i1,i2,i3,var,bound)]=(*w)[fencode3_mpiu(*p,ii,var)];      

                  }
            
            
            #endif
                       
                       break; 
            #ifdef USE_SAC_3D
                       case 2:
         i3=bound*(bound<2)+(((*p)->n[2])-(bound-1))*(bound>1);
         for(i1=0;i1<(((*p)->n[0]));i1++ )
                  for(i2=0;i2<(((*p)->n[1]));i2++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                       ii[2]=i3; 

                      (*gmpiwmod2)[encodempiw2(*p,i1,i2,i3,var,bound)]=(*wmod)[fencode3_mpiu(*p,ii,var)];              
                       (*gmpiw2)[encodempiw2(*p,i1,i2,i3,var,bound)]=(*w)[fencode3_mpiu(*p,ii,var)];      
                    }                            
                       break;                       
            #endif             
             }
                                     
         }    */

//encodempiw (struct params *dp,int ix, int iy, int iz, int field,int bound,int dim)
     //copy data to correct area in w and wmod
     /*for(var=0; var<NVAR; var++)
       for(dim=0;dim<NDIM;dim++) 
         for(bound=0;bound<4;bound++)
         {
            switch(dim)
            {
                       case 0:
            #ifdef USE_SAC_3D
         i1=bound*(bound<2)+(((*p)->n[0])-(bound-1))*(bound>1);
         for(i2=0;i2<(((*p)->n[1]));i2++ )
                  for(i3=0;i3<(((*p)->n[2]));i3++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                       ii[2]=i3;     

                      (*gmpiwmod0)[encodempiw0(*p,i1,i2,i3,var,bound)]=(*wmod)[fencode3_mpiu(*p,ii,var)];              
                       (*gmpiw0)[encodempiw0(*p,i1,i2,i3,var,bound)]=(*w)[fencode3_mpiu(*p,ii,var)];      
  
                  }
            #else
         ii[2]=0;
         i1=bound*(bound<2)+(((*p)->n[0])-(bound-1))*(bound>1);
         for(i2=0;i2<(((*p)->n[1]));i2++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;

                      (*gmpiwmod0)[encodempiw0(*p,i1,i2,i3,var,bound)]=(*wmod)[fencode3_mpiu(*p,ii,var)];              
                       (*gmpiw0)[encodempiw0(*p,i1,i2,i3,var,bound)]=(*w)[fencode3_mpiu(*p,ii,var)];      
                  }            
            
            #endif
                       
                       break;   
                       case 1:
            #ifdef USE_SAC_3D
         i2=bound*(bound<2)+(((*p)->n[1])-(bound-1))*(bound>1);
         for(i1=0;i1<(((*p)->n[0]));i1++ )
                  for(i3=0;i3<(((*p)->n[2]));i3++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                       ii[2]=i3; 

                      (*gmpiwmod1)[encodempiw1(*p,i1,i2,i3,var,bound)]=(*wmod)[fencode3_mpiu(*p,ii,var)];              
                       (*gmpiw1)[encodempiw1(*p,i1,i2,i3,var,bound)]=(*w)[fencode3_mpiu(*p,ii,var)];      
                   }

            #else
         ii[2]=0;
         i2=bound*(bound<2)+(   ((*p)->n[1])-(bound-1)   )*(bound>1);
         for(i1=0;i1<(((*p)->n[0]));i1++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;


                      (*gmpiwmod1)[encodempiw1(*p,i1,i2,i3,var,bound)]=(*wmod)[fencode3_mpiu(*p,ii,var)];              
                       (*gmpiw1)[encodempiw1(*p,i1,i2,i3,var,bound)]=(*w)[fencode3_mpiu(*p,ii,var)];      
                  }
            
            
            #endif
                       
                       break; 
            #ifdef USE_SAC_3D
                       case 2:
         i3=bound*(bound<2)+(((*p)->n[2])-(bound-1))*(bound>1);
         for(i1=0;i1<(((*p)->n[0]));i1++ )
                  for(i2=0;i2<(((*p)->n[1]));i2++ )
                  {
                       ii[0]=i1;
                       ii[1]=i2;
                       ii[2]=i3; 


                      (*gmpiwmod2)[encodempiw2(*p,i1,i2,i3,var,bound)]=(*wmod)[fencode3_mpiu(*p,ii,var)];              
                       (*gmpiw2)[encodempiw2(*p,i1,i2,i3,var,bound)]=(*w)[fencode3_mpiu(*p,ii,var)];      
                   }                            
                       break;                       
            #endif             
             }
                                     
         }  */ 

#ifdef USE_GPUDIRECT
      mpiwtogpu_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wmod,*d_gmpiw0,*d_gmpiwmod0,*d_gmpiw1,*d_gmpiwmod1,*d_gmpiw2,*d_gmpiwmod2,idir);    
     cudaThreadSynchronize();

#else

if(idir==0)
{
   	 cudaMemcpy(*d_gmpiw0, *gmpiw0, szw0*sizeof(real), cudaMemcpyHostToDevice);     
   	 //cudaMemcpy(*d_gmpiwmod0, *gmpiwmod0, szw0*sizeof(real), cudaMemcpyHostToDevice); 
}

if(idir==1)
{
   	 cudaMemcpy(*d_gmpiw1, *gmpiw1, szw1*sizeof(real), cudaMemcpyHostToDevice);     
   	 //cudaMemcpy(*d_gmpiwmod1, *gmpiwmod1, szw1*sizeof(real), cudaMemcpyHostToDevice);     
}
    
            #ifdef USE_SAC_3D
if(idir==2)
{
   	      cudaMemcpy(*d_gmpiw2, *gmpiw2, szw2*sizeof(real), cudaMemcpyHostToDevice);     
   	     // cudaMemcpy(*d_gmpiwmod2, *gmpiwmod0, szw2*sizeof(real), cudaMemcpyHostToDevice); 
}    
         #endif

    //printf("call mpiwtogpu\n");

     mpiwtogpu_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wmod,*d_gmpiw0,*d_gmpiwmod0,*d_gmpiw1,*d_gmpiwmod1,*d_gmpiw2,*d_gmpiwmod2,idir);
     cudaThreadSynchronize();
     
     
 #endif    
     
}

int cucopywmodfrommpiw(struct params **p,real **w, real **wmod,    real **gmpiw0, real **gmpiwmod0,    real **gmpiw1, real **gmpiwmod1,    real **gmpiw2, real **gmpiwmod2, struct params **d_p  ,real **d_w, real **d_wmod,   real **d_gmpiw0, real **d_gmpiwmod0,   real **d_gmpiw1, real **d_gmpiwmod1,   real **d_gmpiw2, real **d_gmpiwmod2, int order, int idir)
{
       int i1,i2,i3;
     int ii[NDIM];
     int var,dim,bound;     
       int szbuf;
     int szw0,szw1,szw2;

  int dimp=(((*p)->n[0]))*(((*p)->n[1]));



real *tgmpiwmod0=*gmpiwmod0;
real *tgmpiwmod1=*gmpiwmod1;


   
 #ifdef USE_SAC_3D  
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif      
     szbuf=2*2*( ((*p)->n[0])+((*p)->n[1]));
     #ifdef USE_SAC_3D
     szbuf=2*2*( ((*p)->n[0])*((*p)->n[1])+ ((*p)->n[0])*((*p)->n[2]) + ((*p)->n[1])*((*p)->n[2])        );
     #endif
        int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;


  #ifdef USE_SAC
  
  szw0=4*NVAR*(  ((*p)->n[1])     );
  szw1=4*NVAR*(  ((*p)->n[0])     );

  #endif
  #ifdef USE_SAC_3D
  
   szw0=4*NVAR*(  ((*p)->n[1])*((*p)->n[2])    );
  szw1=4*NVAR*(    ((*p)->n[0])*((*p)->n[2])   );
  szw2=4*NVAR*(    ((*p)->n[0])*((*p)->n[1])  );

  #endif


#ifndef USE_GPUDIRECT
 
//printf("gpudirect not define!\n");

      //copy data from w and wmod to correct gmpiw and gmpiwmod




if(idir==0)
{



		    i3=0;
                   //   int bound,iside,n=0;
                     //for(iside=0;iside<2;iside++)
		    //for(int i1=0;i1<=1;i1++)
		      //for(int i2=0;i2<(*p)->n[1];i2++)
		      //{
			//iside=0;
                        // bound=i1+2*(iside>0);
			 
			//if(((*p)->ipe==0) /*&&  (*p)->it != -1     && iside==1 && (100*(p->ipe)+10*dim+iside)==101*/ )
			//{
                          //  printf("tini %d %d %d %lg  \n",bound,i2,i1,tgmpiwmod0[encodempiw0 (*p,i1, i2, i3, 0,bound)]);
                            // printf(" %d %d %d   \n",bound,i2,i1);

			//}
                   //   n++;
                    //}





   	 cudaMemcpy(*d_gmpiwmod0, *gmpiwmod0, szw0*sizeof(real), cudaMemcpyHostToDevice); 
}

if(idir==1)
{
   	 cudaMemcpy(*d_gmpiwmod1, *gmpiwmod1, szw1*sizeof(real), cudaMemcpyHostToDevice);     
}
    
            #ifdef USE_SAC_3D
if(idir==2)
{
   	      cudaMemcpy(*d_gmpiwmod2, *gmpiwmod0, szw2*sizeof(real), cudaMemcpyHostToDevice); 
}    
         #endif


  if(idir==1 /*&& (*p)->ipe==0     &&  idir==0*/ )
    {
        ;//printf("ipe2 mpiw0 after bound \n");
        
        ;//for(int j=0; j<4;j++)
         //for(int i=0; i<((*p)->n[1]);i++) 
        ;// for(int i=0; i<10;i++)              
         ;//    printf("%d %d %lg %lg\n",i,j, (tgmpiwmod0[4*rhob*((*p)->n[0]) +i+j*((*p)->n[0])]), (tgmpiwmod1[4*rhob*((*p)->n[0]) +i+j*((*p)->n[0])]));
         ;//printf("\n");
     }
     
 #endif



    //printf("call mpiwtogpu\n");

     mpiwmodtogpu_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w,*d_wmod,*d_gmpiw0,*d_gmpiwmod0,*d_gmpiw1,*d_gmpiwmod1,*d_gmpiw2,*d_gmpiwmod2,idir,order);
     cudaThreadSynchronize();








}


//copy mpi recv buffer to gpu memory     
int cucopywdfrommpiwd(struct params **p,real **wd,     real **gmpiw0,     real **gmpiw1,     real **gmpiw2,  struct params **d_p  ,real **d_wd,    real **d_gmpiw0,   real **d_gmpiw1,   real **d_gmpiw2,  int order, int idir)
{
       int i1,i2,i3;
     int ii[NDIM];
     int var,dim,bound;     
       int szbuf;
     int szw0,szw1,szw2;

  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   
 #ifdef USE_SAC_3D  
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif      
     szbuf=2*2*( ((*p)->n[0])+((*p)->n[1]));
     #ifdef USE_SAC_3D
     szbuf=2*2*( ((*p)->n[0])*((*p)->n[1])+ ((*p)->n[0])*((*p)->n[2]) + ((*p)->n[1])*((*p)->n[2])        );
     #endif
        int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;


  #ifdef USE_SAC
  
  szw0=4*NDERV*(  ((*p)->n[1])     );
  szw1=4*NDERV*(  ((*p)->n[0])     );

  #endif
  #ifdef USE_SAC_3D
  
   szw0=4*NDERV*(  ((*p)->n[1])*((*p)->n[2])    );
  szw1=4*NDERV*(    ((*p)->n[0])*((*p)->n[2])   );
  szw2=4*NDERV*(    ((*p)->n[0])*((*p)->n[1])  );

  #endif


//#ifndef USE_GPUDIRECT


       if(idir==0)
   	 cudaMemcpy(*d_gmpiw0, *gmpiw0, szw0*sizeof(real), cudaMemcpyHostToDevice);     
 
	if(idir==1)
   	 cudaMemcpy(*d_gmpiw1, *gmpiw1, szw1*sizeof(real), cudaMemcpyHostToDevice);     

    
            #ifdef USE_SAC_3D
     if(idir==2)
   	      cudaMemcpy(*d_gmpiw2, *gmpiw2, szw2*sizeof(real), cudaMemcpyHostToDevice);     
         #endif
//#endif
    //printf("call mpiwtogpu\n");

     mpiwdtogpu_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,0,0,*d_wd,*d_gmpiw0,*d_gmpiw1,*d_gmpiw2,idir);
     cudaThreadSynchronize();
}



//copy gpu memory data to mpi send buffer for w and wmod
//just update the edges of w and wmod with values copied from gmpiw, gmpiwmod and gmpivisc
int cucopytompivisc(struct params **p,real **temp2, real **gmpivisc0, real **gmpivisc1, real **gmpivisc2,  struct params **d_p,real **d_wtemp2,    real **d_gmpivisc0,    real **d_gmpivisc1,    real **d_gmpivisc2)
{


     int szbuf,szbuf0,szbuf1,szbuf2;
     int dim,bound,var=0;
     int i1,i2,i3;

  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 
             int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;


     szbuf=2*2*( ((*p)->n[0])+((*p)->n[1]));
	  szbuf0=4*(  (((*p)->n[1])+2 )   );
	  szbuf1=4*(    (((*p)->n[0]) +2 )  );


     #ifdef USE_SAC_3D
     szbuf=2*2*( ((*p)->n[0])*((*p)->n[1])+ ((*p)->n[0])*((*p)->n[2]) + ((*p)->n[1])*((*p)->n[2])        );
  szbuf0=4*(  (((*p)->n[1])+2)*(((*p)->n[2])+2)  ); 
  szbuf1=4*(   (((*p)->n[0])+2)*(((*p)->n[2])+2)    );    
  szbuf2=4*(  (((*p)->n[1])+2)*(((*p)->n[2])+2)   );    


     #endif
     gputompivisc_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_wtemp2,*d_gmpivisc0,*d_gmpivisc1,*d_gmpivisc2);
     cudaThreadSynchronize();
     cudaMemcpy(*gmpivisc0, *d_gmpivisc0, szbuf0*sizeof(real), cudaMemcpyDeviceToHost);
     cudaMemcpy(*gmpivisc1, *d_gmpivisc1, szbuf1*sizeof(real), cudaMemcpyDeviceToHost);
     #ifdef USE_SAC_3D
     	cudaMemcpy(*gmpivisc2, *d_gmpivisc2, szbuf2*sizeof(real), cudaMemcpyDeviceToHost);
     #endif
     //copy data to correct area in temp2
//encodempiw (struct params *dp,int ix, int iy, int iz, int field,int bound,int dim)
     //copy data to correct area in w and wmod
     /*  for(dim=0;dim<NDIM;dim++) 
         for(bound=0;bound<2;bound++)
         {
            switch(dim)
            {
                       case 0:
            #ifdef USE_SAC_3D
         i1=bound*(((*p)->n[0])+1);
         for(i2=1;i2<(((*p)->n[1])+2);i2++ )
                  for(i3=1;i3<(((*p)->n[2])+2);i3++ )
                  {     
                        
          //i1=(p->n[0])+1;
         
          //temp2[encode3p2_sacmpi (p,i1, i2, i3, tmpnui)]=gmpitgtbufferr[0][i2+i3*((p->n[1])+2)];
          //temp2[encode3p2_sacmpi (p,0, i2, i3, tmpnui)]=gmpitgtbufferl[0][i2+i3*((p->n[1])+2)];
         
                       (*temp2)[encode3p2_mpiu(*p,i1,i2,i3,var)]=(*gmpivisc0)[encodempivisc0(*p,i1,i2,i3,bound,dim)];
                  }
            #else
         i3=0;
         i1=bound*(((*p)->n[0])+1);
                  for(i2=1;i2<(((*p)->n[1])+2);i2++ )
                  {
                       (*temp2)[encode3p2_mpiu(*p,i1,i2,i3,var)]=(*gmpivisc0)[encodempivisc0(*p,i1,i2,i3,bound,dim)];
                  }            
            
            #endif
                       
                       break;   
                       case 1:
            #ifdef USE_SAC_3D
         i2=bound*(((*p)->n[1])+1);
         for(i1=1;i1<(((*p)->n[0])+2);i1++ )
                  for(i3=1;i3<(((*p)->n[2])+2);i3++ )
                  {
                       (*temp2)[encode3p2_mpiu(*p,i1,i2,i3,var)]=(*gmpivisc1)[encodempivisc1(*p,i1,i2,i3,bound,dim)];
                  }

            #else
         i3=0;
         i2=bound*(((*p)->n[1])+1);
                  for(i1=1;i1<(((*p)->n[0])+2);i1++ )
                  {                                                       
                       (*temp2)[encode3p2_mpiu(*p,i1,i2,i3,var)]=(*gmpivisc1)[encodempivisc1(*p,i1,i2,i3,bound,dim)];
                  }
            
            
            #endif
                       
                       break; 
            #ifdef USE_SAC_3D
                       case 2:
                  i3=bound*(((*p)->n[2])+1);
        for(i1=1;i1<(((*p)->n[0])+2);i1++ )
                  for(i2=1;i2<(((*p)->n[1])+2);i2++ )
                  {                                                           
                       (*temp2)[encode3p2_mpiu(*p,i1,i2,i3,var)]=(*gmpivisc2)[encodempivisc2(*p,i1,i2,i3,bound,dim)];
                  }                            
                       break;                       
            #endif             
             }
                                     
         }    */

}

//copy mpi recv buffer to gpu memory     
int cucopyfrommpivisc(struct params **p,real **temp2,real **gmpivisc0,real **gmpivisc1,real **gmpivisc2,  struct params **d_p,real **d_wtemp2,    real **d_gmpivisc0,    real **d_gmpivisc1,    real **d_gmpivisc2)
{
      int dim,bound,var=0;
     int i1,i2,i3;      
 
     int szbuf,szbuf0,szbuf1,szbuf2;

  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   
 #ifdef USE_SAC_3D  
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 

        int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;

     
     szbuf=2*2*( ((*p)->n[0])+((*p)->n[1]));
     
 	  szbuf0=4*(  (((*p)->n[1])+2 )   );
	  szbuf1=4*(    (((*p)->n[0]) +2 )  );

     #ifdef USE_SAC_3D
     szbuf=2*2*( ((*p)->n[0])*((*p)->n[1])+ ((*p)->n[0])*((*p)->n[2]) + ((*p)->n[1])*((*p)->n[2])        );

     
  szbuf0=4*(  (((*p)->n[1])+2)*(((*p)->n[2])+2)  ); 
  szbuf1=4*(   (((*p)->n[0])+2)*(((*p)->n[2])+2)    );    
  szbuf2=4*(  (((*p)->n[1])+2)*(((*p)->n[2])+2)   ); 

     #endif

      //copy data from temp2 to gmpivisc
        /*     for(dim=0;dim<NDIM;dim++) 
         for(bound=0;bound<2;bound++)
         {
            switch(dim)
            {
                       case 0:
            #ifdef USE_SAC_3D
         i1=bound*(((*p)->n[0])+1);
         for(i2=1;i2<(((*p)->n[1])+2);i2++ )
                  for(i3=1;i3<(((*p)->n[2])+2);i3++ )
                  {     
                        
          //i1=(p->n[0])+1;
         
          //temp2[encode3p2_sacmpi (p,i1, i2, i3, tmpnui)]=gmpitgtbufferr[0][i2+i3*((p->n[1])+2)];
          //temp2[encode3p2_sacmpi (p,0, i2, i3, tmpnui)]=gmpitgtbufferl[0][i2+i3*((p->n[1])+2)];
         
                       (*gmpivisc0)[encodempivisc0(*p,i1,i2,i3,bound,dim)]=(*temp2)[encode3p2_mpiu(*p,i1,i2,i3,var)];
                  }
            #else
         i3=0;
         i1=bound*(((*p)->n[0])+1);
                  for(i2=1;i2<(((*p)->n[1])+2);i2++ )
                  {
                       (*gmpivisc0)[encodempivisc0(*p,i1,i2,i3,bound,dim)]=(*temp2)[encode3p2_mpiu(*p,i1,i2,i3,var)];
                  }            
            
            #endif
                       
                       break;   
                       case 1:
            #ifdef USE_SAC_3D
         i2=bound*(((*p)->n[1])+1);
         for(i1=1;i1<(((*p)->n[0])+2);i1++ )
                  for(i3=1;i3<(((*p)->n[2])+2);i3++ )
                  {
                       (*gmpivisc1)[encodempivisc1(*p,i1,i2,i3,bound,dim)]=(*temp2)[encode3p2_mpiu(*p,i1,i2,i3,var)];
                  }

            #else
         i3=0;
         i2=bound*(((*p)->n[1])+1);
                  for(i1=1;i1<(((*p)->n[0])+2);i1++ )
                  {
                                                                     
                       (*gmpivisc1)[encodempivisc1(*p,i1,i2,i3,bound,dim)]=(*temp2)[encode3p2_mpiu(*p,i1,i2,i3,var)];
                  }
            
            
            #endif
                       
                       break; 
            #ifdef USE_SAC_3D
                       case 2:
                  i3=bound*(((*p)->n[2])+1);
        for(i1=1;i1<(((*p)->n[0])+2);i1++ )
                  for(i2=1;i2<(((*p)->n[1])+2);i2++ )
                  {
                                                              
                       (*gmpivisc2)[encodempivisc2(*p,i1,i2,i3,bound,dim)]=(*temp2)[encode3p2_mpiu(*p,i1,i2,i3,var)];
                  }                            
                       break;                       
            #endif             
             }
                                     
         } */   


   	 cudaMemcpy(*d_gmpivisc0, *gmpivisc0, szbuf0*sizeof(real), cudaMemcpyHostToDevice);     
   	 cudaMemcpy(*d_gmpivisc1, *gmpivisc1, szbuf1*sizeof(real), cudaMemcpyHostToDevice);
       #ifdef USE_SAC_3D    
   	 cudaMemcpy(*d_gmpivisc2, *gmpivisc2, szbuf2*sizeof(real), cudaMemcpyHostToDevice); 
       #endif    

     mpivisctogpu_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_wtemp2,*d_gmpivisc0,*d_gmpivisc1,*d_gmpivisc2);
     cudaThreadSynchronize();
}


#endif



