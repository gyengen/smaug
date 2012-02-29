// iosac2.5d.cpp : Defines the entry point for the console application.
//


#include "../include/iosac2.5d.h"
#include "../include/step.h"
#include "../include/iobparams.h"
/*----------------------*/ 
real second()
{

   /*REAL secs;
   clock_t Time;
   Time = clock();

   secs = (real) Time / (real) CLOCKS_PER_SEC;
   return secs;*/
   real retval;
	static long zsec=0;
	static long zusec=0;
	real esec;
	
	struct timeval tp;
	struct timezone tzp;
	
	gettimeofday(&tp, &tzp);
	
	if(zsec==0) zsec=tp.tv_sec;
	if(zusec==0) zusec=tp.tv_usec;
	
	retval=(tp.tv_sec - zsec)+(tp.tv_usec-zusec)*0.000001;
	return retval;

}




int main(int argc, char* argv[])
{

int itype=-1;
int status=1;
int it=0; //test integer to be returned 

//getintparam_( int elist.id,char *sname,int *iv,  int elist.port, char *selist.server );
//int elist.id=0;
//int elist.port=8080;

int i1,i2,i3,j1;
int i,j,k,iv;


char *portfile=(char *)calloc(500,sizeof(char));
char *sdir=(char *)calloc(500,sizeof(char));
char *name=(char *)calloc(500,sizeof(char));
char *outfile=(char *)calloc(500,sizeof(char));
char *formfile=(char *)calloc(500,sizeof(char));


#include "../include/defs.h"
#include "../include/iosac2.5dparams.h"


struct bparams *d_bp;
struct bparams *bp=(struct bparams *)malloc(sizeof(struct bparams));


FILE *portf;

#ifdef USE_IOME
if(argc>1)
{ 
   sprintf(portfile,"%s0_port.txt",argv[1]) ;  
   //strcpy(name,argv[1]);
   portf=fopen(portfile,"r");
   fscanf(portf,"%d %s",&elist.port,elist.server);
   fclose(portf);

   printf("read file junk is %d %s\n",elist.port,elist.server);
}
#endif





//int cuprop(struct params **p, real **w, real **wnew, real **b,struct params **d_p, real **d_w, real **d_wnew, real **d_b, real **d_wmod, real **d_dwn1, real **d_dwn2, real **d_dwn3, real **d_dwn4, real **d_wd)




printf("rho %d mom1 %d mom2 %d\n",rho,mom1,mom2);


printf("calling cuinit\n");

#ifdef USE_MPI
     MPI::Init(argc, argv);
     mpiinit(p);
     sprintf(configfile,"%s",cfgout);
     ipe2iped(p);     
     mpineighbours(0,p);
     mpineighbours(1,p);
     
     #ifdef USE_SAC3D
          mpineighbours(2,p);
     #endif
     
#else
     sprintf(configfile,"%s",cfgout);
#endif

//printf("cfgfile %s\n",configfile);
 //   getintparam_( &elist.id,"i1",&it,  &elist.port, "localhost" );	
//	printf("Get integer %d\n",it);
    //Set input filename as first arg
	//if NULL use defaults
	char *method=NULL;
	//CIoSimulation *TestSimulation;
	//this should be executed by the iome start up application
	//exec('ioshallowwater.sce');

	//this application is started using the io  start scilab application
	//exec('paramssteeringtest1.sce');
	//stacksize('max');
	//stacksize(268435454)
	//open the file generated
	//sprintf(elist.portfile,"%s0_elist.port.txt",meta.name);
	//FILE *fd=fopen(elist.portfile,"r");
	//int elist.portelist.id;
	//fscanf(fd,"%d",&elist.portelist.id);
	//fclose(fd);
	//elist.elist.port=elist.portelist.id;

    #ifdef USE_IOME
        if(argc>2)
        {
          //simfile already read by 
          readsim(p,&meta,argv[2],elist);
          //if((p->readini)!=0)
          //   readconfig(meta.ini_file,*p,meta,w);
        }
        else
	  createsim(*p,meta,simfile,elist);

	sprintf(simfile,"%s.xml",meta.name);
        sprintf(newsimfile,"%s_update.xml",meta.name);
     #endif
	//NewSimulation(metadata.name,'test1.xsl',elist);

// Build empty u, v, b matrices
// Define h
printf("allocating w and wnew\n");

  #ifdef USE_SAC_3D
 w=(real *)calloc(ni*nj*nk*NVAR,sizeof(real ));
wd=(real *)calloc(ni*nj*nk*NDERV,sizeof(real ));
 wnew=(real *)calloc(ni*nj*nk*NVAR,sizeof(real ));
for(i=0;i<ni*nj*nk*NDERV;i++)
    wd[i]=0.0;
 #else
 w=(real *)calloc(ni*nj*NVAR,sizeof(real ));
wd=(real *)calloc(ni*nj*NDERV,sizeof(real ));
 wnew=(real *)calloc(ni*nj*NVAR,sizeof(real ));
 for(i=0;i<ni*nj*NDERV;i++)
    wd[i]=0.0;
#endif

//set initial time step to a large value
if(p->moddton==1.0)
{
	p->dt=1.0e50;
}

if((p->readini)==0)
 initconfig(p, &meta, w);
else
 readasciivacconfig(cfgfile,*p,meta,w,wd,hlines);

printf("after read\n");
p->it=0;

//writeasciivacconfig(cfgout,*p, meta , w,hlines,*state);
//writevacconfig(cfgout,0,*p, meta , w,*state);
  /*   for( j1=2;j1<5;j1++)
      {
        for( i1=0;i1<ni;i1++)
	{
        //j1=2;
	printf("%d %d %f %f %f\n",i1,j1,w[j1*ni+i1+(ni*nj*rho)],w[j1*ni+i1+(ni*nj*mom1)],w[j1*ni+i1+(ni*nj*mom2)]);
        }     
       //

      }
 printf("\n");*/  
  real *u,  *v,  *h;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;
  h=w+(ni)*(nj)*rho;
  u=w+(ni)*(nj)*mom1;
  v=w+(ni)*(nj)*mom2;

cuinit(&p,&bp,&w,&wnew,&wd,&state,&d_p,&d_bp,&d_w,&d_wnew,&d_wmod, &d_dwn1,  &d_wd, &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);
initgrid(&p,&w,&wnew,&state,&wd,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1,  &d_wd, &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);

#ifdef USE_MPI

  //initialise the mpi used memory locations
 cuinitmpibuffers(&p, &w, &wmod, &temp2, &gmpivisc,   &gmpiw, &gmpiwmod, &d_p, &d_w, &d_wmod,&d_wtemp2,  &d_gmpivisc, &d_gmpiw, &d_gmpiwmod);



  cucopywtompiw(&p,&w, &wmod,    &gmpiw, &gmpiwmod, &d_p,  &d_w, &d_wmod,   &d_gmpiw, &d_gmpiwmod, 0);
#endif




for(int ii=0; ii<NVAR; ii++)
for(int idir=0; idir<NDIM; idir++)
{
   p->it=-1;  //initialise fixed boundaries
   //printf("btype %d %d %d\n",idir,ii,p->boundtype[ii][idir]);
   if((p->boundtype[ii][idir])==5)  //period=0 mpi=1 mpiperiod=2  cont=3 contcd4=4 fixed=5 symm=6 asymm=7
   {

               cuboundary(&p, &bp, &d_p, &d_bp, &d_state, &d_w, 0,idir,ii);
 
   }
}

p->it=0;  

  cuboundary(&p,&bp,&d_p,&d_bp,&d_state,&d_w, 0,0,0);
#ifdef USE_MPI
   mpibound(NVAR, d_w ,d_p);
#endif
  cuboundary(&p,&bp,&d_p,&d_bp,&d_state,&d_wmod, 0,0,0);
#ifdef USE_MPI
   mpibound(NVAR, d_wmod ,d_p);
   cucopywfrommpiw(&p,&w, &wmod,    &gmpiw, &gmpiwmod, &d_p,  &d_w, &d_wmod,   &d_gmpiw, &d_gmpiwmod,0);
#endif

printf("after cuinit\n");
/*for( j1=0;j1<nj;j1++)
      {
        for( i1=0;i1<ni;i1++)
	{
              for(iv=0;iv<NVAR;iv++)               
                      printf("%f ", w[j1*ni+i1+(ni*nj*iv)]);
               printf("\n");

           }
         }*/


//For a steerable simulation generate and save a dxformfile that saves a single data step
//used for the steering dx module
//printf("here in runsim2a\n");

#ifdef USE_IOME
getmetadata_(elist.id,"directory",&sdir,elist.port,elist.server);
//sdir=metadata.directory

//name=metadata.name;

getmetadata_(elist.id,"name",&name,elist.port,elist.server);
//disp(sdir,name)
//printf("here in runsim3\n");
sprintf(outfile,"%s/%s.out",sdir,name);
#endif



//createlog(meta.log_file);



//while(finishsteering == 0)
//{
 
  //  if( steeringenabled==0)
  //    finishsteering=1;
 int n;  
// nt=24; 
real t1,t2,ttot;
int order=0;
int ordero=0;
int order1;
int orderb=0;
int ii,ii0,ii1;
real dtdiffvisc;
ttot=0;
real time=0.0;



   state->it=0;
   state->t=0;
   state->dt=p->dt;

//printf("dx dy%f %f",dx,dy);
for( n=1;n<=nt;n++)
//for( n=0;n<1;n++)
{
    p->it=n;
    if(((n-1)%(p->cfgsavefrequency))==0)
    {
      //writeconfig(name,n,*p, meta , w);
#ifndef USE_MPI
     // writevtkconfig(configfile,n,*p, meta , w);
#endif
      //writeasciivacconfig(cfgout,*p, meta , w,hlines,*state);

      writevacconfig(configfile,n,*p, meta , w,wd,*state);
       
    }
   order=0;
   t1=second();

   if(p->moddton==1.0)
   {
        //printf(" courant is %f \n",p->courant);
        //p->courant=0.1;
        p->maxcourant=0.0;
        courantmax=0.0;
        for(int dim=0; dim<=(NDIM-1); dim++)
        {
        cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
        cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
        cucomputemaxcourant(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
        //printf("maxcourant %d %16.10g  %16.10g  %16.10g\n",dim,p->maxcourant,p->cmax,p->dx[dim]);
        }
        
        //courantmax=0.0;
        /*for(int dim=0; dim<=(NDIM-1); dim++)
        {
           if((cmax[dim]/(p->dx[dim]))>(p->maxcourant))
             courantmax=cmax[dim]/(p->dx[dim]);
             printf("%d %16.10g\n",dim,courantmax);
           //  printf("cmax %g ",cmax[dim]);
        }*/
        
        //printf("old dt is %g ",p->dt);
        //if(courantmax>smalldouble) dt=min(dt,courantpar/courantmax)

        //if(((p->maxcourant)>1.0e-8) && (p->dt)>(((p->courant)/(p->maxcourant))   ))
        if( (p->dt)>1.0e-8 && (p->dt)>  ((  (p->courant)/(p->maxcourant)  ))   )
               p->dt=(p->courant)/(p->maxcourant);
        //printf("new dt is %g %g\n",(p->courant)/(p->maxcourant),p->dt);


 
        cugetdtvisc1(&p,&d_p,&d_wmod, &wd,&d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2);
          #ifdef USE_MPI
              mpiallreduce(&(p->maxviscoef), MPI_MAX);
          #endif
       /* for(int dim=0; dim<=(NDIM-1); dim++)
        {
        dtdiffvisc=0.25/(p->maxviscoef/((p->dx[dim])*(p->dx[dim])));
        printf("dim %d dtdiffvisc %20.10g  %20.10g %20.10g\n",dim,p->maxviscoef,dtdiffvisc,p->dx[dim]);
        if(dtdiffvisc>1.0e-8 && (p->dt)>dtdiffvisc )
         //   if( (p->dt)>dtdiffvisc )
                                      p->dt=dtdiffvisc;
        }*/

         // printf("dtdiffvisc %20.10g  %20.10g\n",p->maxviscoef,p->dtdiffvisc);
         if(1/(p->dtdiffvisc)>1.0e-8 && (p->dt)>(1/(p->dtdiffvisc)) )
         //   if( (p->dt)>dtdiffvisc )
                                      p->dt=(1/p->dtdiffvisc);
        //cugetdtvisc1(&p,&d_p,&d_wmod, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2);
        printf(" modified dt is %20.10g \n",p->dt);

        //include gravitational modification
        

   } 


if((p->rkon)==0)
{
  ordero=0;
 
  cucomputedervfields(&p,&d_p,&d_wmod, &d_wd,order);
  order=1;

//printf("\n");
 for(int dir=0;dir<NDIM; dir++)
 {


              //cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dir);
               //cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dir);

              // printf("cmax=%10.12f\n",p->cmax);  
  for(int f=rho; f<=(mom1+NDIM-1); f++)
  { 
      if((f==mom1 && dir==0)  ||  (f==mom2 && dir==1)  || (f==mom2 && dir==2) )
       cucomputept(&p,&d_p,&d_wmod, &d_wd,order,dir);
      cucentdiff1(&p,&d_p,&d_state,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt,f,dir);

     
  }

   //cucomputevels(&p,&d_p,&d_wmod, &d_wd,order,dir);


#ifndef ADIABHYDRO
   for(int f=energy; f<=(b1+(NDIM-1)); f++)
   {
     if(f==energy)
     {
         cucomputevels(&p,&d_p,&d_wmod, &d_wd,order,dir);
         cucomputepbg(&p,&d_p,&d_wmod, &d_wd,order,dir);
         cucomputept(&p,&d_p,&d_wmod, &d_wd,order,dir);
     }
      
      cucentdiff2(&p,&d_p,&d_state,&d_w,&d_wmod, &d_dwn1, &d_wd,order, ordero,p->dt,f,dir);
   }

#endif
  }
   
 
       cugrav(&p,&d_p,&d_state,&d_w,&d_wmod, &d_dwn1, &d_wd,order, ordero,p->dt);

   //cuboundary(&p,&d_p,&d_wmod, ordero);
   if(p->divbon==1)
	       cudivb(&p,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt);
   if(p->hyperdifmom==1)
   {
    dt=(p->dt);
             p->maxviscoef=0.0;
    
    #ifdef USE_SHOCKVISC
       cunushk1(&p,&d_p,&d_wmod, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2);
    #endif
    for(int dim=0; dim<=(NDIM-1); dim++)
     {
              cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
               cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);

              // printf("field %d dim %d cmax %20.16g\n",0,dim,p->cmax);
          #ifdef USE_MPI
              mpiallreduce(&(p->cmax), MPI_MAX);
          #endif
               //printf("cmax=%f\n",p->cmax);       //printf(" courant is %f \n",p->courant);
       cmax[dim]=p->cmax;
       cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);

       //cuhyperdifvisc1il(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);
       #ifdef USE_MPI
       
          cucopytompivisc(&p,&temp2, &gmpivisc,  &d_p,&d_wtemp2,    &d_gmpivisc);
          mpivisc(dim,p,temp2);
          cucopyfrommpivisc(&p,&temp2, &gmpivisc,  &d_p,&d_wtemp2,    &d_gmpivisc);
       #endif
       cuhyperdifvisc1r(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);
       cuhyperdifvisc1l(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);



      cuhyperdifrhosource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,rho,dim,p->dt);

     }


     for(int dim=0; dim<=(NDIM-1); dim++)
     {
              cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
               cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
           //   printf("field %d dim %d cmax %20.16g\n",3,dim,p->cmax);
          #ifdef USE_MPI
              mpiallreduce(&(p->cmax), MPI_MAX);
          #endif
               //printf("cmax=%f\n",p->cmax);       cuhyperdifvisc1(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim,0);
      //p->cmax=cmax[dim];
      cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);
      //cuhyperdifvisc1il(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);
       #ifdef USE_MPI
          cucopytompivisc(&p,&temp2, &gmpivisc,  &d_p,&d_wtemp2,    &d_gmpivisc);
          mpivisc(dim,p,temp2);
          cucopyfrommpivisc(&p,&temp2, &gmpivisc,  &d_p,&d_wtemp2,    &d_gmpivisc);
       #endif
      cuhyperdifvisc1r(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);
      cuhyperdifvisc1l(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);
       
       cuhyperdifesource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,energy,dim,dt);
   
     }

for(int dim=0; dim<=(NDIM-1); dim++)
       for(int f=0; f<=(NDIM-1); f++)
           	                 
	     {
               cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
               cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);

             //printf("field %d dim %d cmax %20.16g\n",1+f,dim,p->cmax);
          #ifdef USE_MPI
              mpiallreduce(&(p->cmax), MPI_MAX);
          #endif
               //printf("cmax=%f\n",p->cmax);
               //p->cmax=cmax[dim];
              //p->cmax=1.0;
               cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);
               //cuhyperdifvisc1il(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);
                      #ifdef USE_MPI
          cucopytompivisc(&p,&temp2, &gmpivisc,  &d_p,&d_wtemp2,    &d_gmpivisc);
          mpivisc(dim,p,temp2);
          cucopyfrommpivisc(&p,&temp2, &gmpivisc,  &d_p,&d_wtemp2,    &d_gmpivisc);

       #endif
               cuhyperdifvisc1r(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);
               cuhyperdifvisc1l(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);

                for(ii1=0;ii1<=1;ii1++)
                {

                          if (ii1 == 0)
                          {
                           ii=dim;
                           ii0=f;
                          }
                          else
                           {
                           ii=f;
                           ii0=dim;
                            }

                  if(ii==dim)
                    cuhyperdifmomsource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,p->dt);
                   else
                    cuhyperdifmomsourcene1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,p->dt);

                }
             }
            int jj,mm,kk;
             real sb;
             for(int dim=0; dim<=(NDIM-1); dim++)
	     for(int f=0; f<=(NDIM-1); f++) 
             if(f!=dim)           
	     {
               cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
               cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);

              // printf("field %d dim %d cmax %20.16g\n",4+f,dim,p->cmax);
          #ifdef USE_MPI
              mpiallreduce(&(p->cmax), MPI_MAX);
          #endif
               //printf("cmax=%f\n",p->cmax);
      //p->cmax=cmax[dim];

               cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);
              // cuhyperdifvisc1il(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);
                      #ifdef USE_MPI
          cucopytompivisc(&p,&temp2, &gmpivisc,  &d_p,&d_wtemp2,    &d_gmpivisc);
          mpivisc(dim,p,temp2);
          cucopyfrommpivisc(&p,&temp2, &gmpivisc,  &d_p,&d_wtemp2,    &d_gmpivisc);
         
       #endif
               cuhyperdifvisc1r(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);
               cuhyperdifvisc1l(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);


                for(ii1=0;ii1<=1;ii1++)
                {

                          if (ii1 == 0)
                          {
                           jj=dim;
                           mm=f;
                           sb=-1.0;
                           ii0=dim;
                          }
                          else
                           {
                           ii0=f;
                           mm=dim;
                           sb=1.0;
                           jj=f;
                           
                            }

                  if(mm==dim)
                     cuhyperdifbsource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,p->dt);
                  else
                     cuhyperdifbsourcene1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,p->dt);
 

                }
              } 


   }
   
           cusource(&p,&d_p,&d_state,&d_w,&d_wmod, &d_dwn1, &d_wd,order, ordero,p->dt);

   //cuadvance(&p,&d_p,&d_wmod,&d_w, order);
	#ifdef USE_MPI
	   cucopywtompiw(&p,&w, &wmod,    &gmpiw, &gmpiwmod, &d_p,  &d_w, &d_wmod,   &d_gmpiw, &d_gmpiwmod, order);
	   mpibound(NVAR, d_wmod ,d_p);
	   cucopywfrommpiw(&p,&w, &wmod,    &gmpiw, &gmpiwmod, &d_p,  &d_w, &d_wmod,   &d_gmpiw, &d_gmpiwmod,order);
	   
	#endif
   cuboundary(&p,&bp,&d_p,&d_bp,&d_state,&d_wmod, ordero,0,0);

}

   if((p->rkon)==1)
     for(order=0; order<4; order++) 
   {	   
           ordero=order+1;
           dt=(p->dt)/2.0;
           orderb=order+2;

           if(order==2)
           {
              dt=(p->dt);
              orderb=1;
            }


           if(order==3)
           {
              dt=(p->dt)/6.0;
              ordero=0;
              orderb=0;
           }


           cucomputedervfields(&p,&d_p,&d_wmod, &d_wd,order);
 for(int dir=0;dir<(NDIM-1); dir++)
 {
           cucomputevels(&p,&d_p,&d_wmod, &d_wd,order,dir);

           for(int f=rho; f<=mom1+(NDIM-1); f++)
	       cucentdiff1(&p,&d_p,&d_state,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,dt,f,dir);


#ifndef ADIABHYDRO
           for(int f=energy; f<=b1+(NDIM-1); f++)
	       cucentdiff2(&p,&d_p,&d_state,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt,f,dir);

#endif

}
	   //cuderivsource(&p,&w,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt);
	   if(p->divbon==1)
	       cudivb(&p,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd, order,ordero,p->dt);
           if(p->hyperdifmom==1)
           {
             p->maxviscoef=0.0;
         #ifdef USE_SHOCKVISC        
           cunushk1(&p,&d_p,&d_wmod,&d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2);
         #endif

	     for(int dim=0; dim<=(NDIM-1); dim++)
	     {
               cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
               cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
          #ifdef USE_MPI
              mpiallreduce(&(p->cmax), MPI_MAX);
          #endif
          
           cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);
 	       //cuhyperdifvisc1il(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);
       #ifdef USE_MPI
          cucopytompivisc(&p,&temp2, &gmpivisc,  &d_p,&d_wtemp2,    &d_gmpivisc);
          mpivisc(dim,p,temp2);
          cucopyfrommpivisc(&p,&temp2, &gmpivisc,  &d_p,&d_wtemp2,    &d_gmpivisc);
       #endif
	       cuhyperdifvisc1r(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);
	       cuhyperdifvisc1l(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);

 

	       cuhyperdifrhosource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,rho,dim,dt);


	     }

     for(int dim=0; dim<=(NDIM-1); dim++)
     {
               cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
               cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
          #ifdef USE_MPI
              mpiallreduce(&(p->cmax), MPI_MAX);
          #endif
       cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);
      //cuhyperdifvisc1il(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);
       #ifdef USE_MPI
          cucopytompivisc(&p,&temp2, &gmpivisc,  &d_p,&d_wtemp2,    &d_gmpivisc);
          mpivisc(dim,p,temp2);
          cucopyfrommpivisc(&p,&temp2, &gmpivisc,  &d_p,&d_wtemp2,    &d_gmpivisc);
       #endif
       cuhyperdifvisc1r(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);
      cuhyperdifvisc1l(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);


       cuhyperdifesource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,energy,dim,dt);

     }
             
for(int dim=0; dim<=(NDIM-1); dim++)
       for(int f=0; f<=(NDIM-1); f++)
           	                 
	     {
               cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
               //cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim);
          #ifdef USE_MPI
           ;//   mpiallreduce(&(p->cmax), MPI_MAX);
          #endif
               cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);
               //cuhyperdifvisc1il(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);
                      #ifdef USE_MPI
          cucopytompivisc(&p,&temp2, &gmpivisc,  &d_p,&d_wtemp2,    &d_gmpivisc);
          mpivisc(dim,p,temp2);
          cucopyfrommpivisc(&p,&temp2, &gmpivisc,  &d_p,&d_wtemp2,    &d_gmpivisc);

       #endif
               cuhyperdifvisc1r(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);
               cuhyperdifvisc1l(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);

                for(ii1=0;ii1<=1;ii1++)
                {

                          if (ii1 == 0)
                          {
                           ii=dim;
                           ii0=f;
                          }
                          else
                           {
                           ii=f;
                           ii0=dim;
                            }

                  if(ii==dim)
                    cuhyperdifmomsource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
                   else
                    cuhyperdifmomsourcene1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);

                }
             }

            int jj,mm,kk;
             real sb;
             for(int dim=0; dim<=(NDIM-1); dim++)
	     for(int f=0; f<=(NDIM-1); f++) 
             if(f!=dim)                     
	     {
               cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
               //cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim);
          #ifdef USE_MPI
              ;//mpiallreduce(&(p->cmax), MPI_MAX);
          #endif
               cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);
               //cuhyperdifvisc1il(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);
       #ifdef USE_MPI
          cucopytompivisc(&p,&temp2, &gmpivisc,  &d_p,&d_wtemp2,    &d_gmpivisc);
          mpivisc(dim,p,temp2);
          cucopyfrommpivisc(&p,&temp2, &gmpivisc,  &d_p,&d_wtemp2,    &d_gmpivisc);

       #endif
               cuhyperdifvisc1r(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);
               cuhyperdifvisc1l(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);


                for(ii1=0;ii1<=1;ii1++)
                {

                          if (ii1 == 0)
                          {
                           jj=dim;
                           mm=f;
                           sb=-1.0;
                           ii0=dim;
                          }
                          else
                           {
                           ii0=f;
                           mm=dim;
                           sb=1.0;
                           jj=f;
                           
                            }

                  if(mm==dim)
                     cuhyperdifbsource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
                   else
                     cuhyperdifbsourcene1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);


                }
              } 

           }
           //cuboundary(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1, &d_wd,ordero);
           cuadvance(&p,&d_p,&d_wmod,&d_w,order);
		#ifdef USE_MPI
	   cucopywtompiw(&p,&w, &wmod,    &gmpiw, &gmpiwmod, &d_p,  &d_w, &d_wmod,   &d_gmpiw, &d_gmpiwmod, order);
	   mpibound(NVAR, d_wmod ,d_p);
	   cucopywfrommpiw(&p,&w, &wmod,    &gmpiw, &gmpiwmod, &d_p,  &d_w, &d_wmod,   &d_gmpiw, &d_gmpiwmod,order);
	   
	#endif
           cuboundary(&p,&bp,&d_p,&d_bp,&d_state,&d_wmod, orderb,0,0);
	   

   }

   p->it=n+1;
   

  // cuupdate(&p,&w,&wmod,&wd,&temp2,&state,&d_p,&d_w,&d_wmod, &d_wtemp2 &d_state,n);
 cuupdate(&p,&w,&wmod,&temp2,&state,&d_p,&d_w,&d_wmod,&d_wtemp2,  &d_state,n);

printf("\n");
   //printf("nummaxthreads %d\n",p->mnthreads);

   t2=second()-t1;
   ttot+=t2;
   printf("step %d total time %f\n",n,ttot);

   state->it=n;
   state->t=time+(p->dt);
   time=state->t;
   state->dt=p->dt;



   //appendlog(meta.log_file,*p, *state);

 
    
    /*getintparam_(&elist.id,"steeringenabled",&steeringenabled,&elist.port,elist.server);
    if(steeringenabled==1)
    {
      //disp('getting updatea params');
      //for steering get the modified control params
      double dg;
      getintparam_(&elist.id,"finishsteering",&finishsteering,&elist.port,elist.server);//source y location  
        // Constants
      getdoubleparam_(elist.id,"g",&dg,elist.port,elist.server);

      g=dg;
     
    }*/
    
   /* for( j1=ngj;j1<nj-ngj;j1++)
        for( i1=ngi;i1<ni-ngi;i1++)
{

;//w[j1*ni+i1+(ni*nj*b1)]=wd[j1*ni+i1+(ni*nj*(hdnur))];
;//w[j1*ni+i1+(ni*nj*b2)]=wd[j1*ni+i1+(ni*nj*(hdnul))];

}*/
           
      //save file containing current data
     // sprintf(configfile,"tmp/%ss%d.out",name,n);-
     // printf("check dims %d %d \n",ni,nj);
     // FILE *fdt=fopen(configfile,"w");
     // fprintf(fdt,"%d\n",n);
     //for( j1=0;j1<nj;j1++)
     // {
      //  for( i1=0;i1<ni;i1++)
	//{
               // printf("%d %d ", i1,j1);
	//	fprintf(fdt,"%d %d %f %f %f %f %f %f %f %f\n",i1,j1,*(h+(j1*ni+i1)),(u[j1*ni+i1]),(v[j1*ni+i1]),w[j1*ni+i1+(ni*nj*mom3)],w[j1*ni+i1+(ni*nj*energy)],w[j1*ni+i1+(ni*nj*b1)],w[j1*ni+i1+(ni*nj*b2)],w[j1*ni+i1+(ni*nj*b3)]);
           //fprintf(fdt,"%d %f %f %f ",j1+i1*nj, u[j1+i1*nj],v[j1+i1*nj],h[j1+i1*nj]);
               // fprintf(fdt,"%f ",h[j1+i1*nj]);
       // }     
        //printf("\n");   
        //fprintf(fdt,"\n");
     // }
     // fclose(fdt);
   
 //disp('writing data');
    //if finishsteering==1
    //  fprintf(fd,"%d\n",n);
    //  for( j1=0;j1<nj;j1++)
    //  {
     //   for( i1=0;i1<ni;i1++)
	//{
     //     fprintf(fd,"%f %f %f ",u[j1*ni+i1],v[j1*ni+i1],h[j1*ni+i1]);
	   	//fprintf(fd,"%f ",h[j1*ni+i1]);
	//}
       // fprintf(fd,"\n");
    //  }


    }//end of testep
  printf("params %f %f %f %f\n",g,dx,dy,dt);  
    //end
  
   
    //end
    
    //disp('written data');
//}//end of steering test finished

  // force the final ouput file to be over written
 // if(steeringenabled==1)
 // {
  //  if(finishsteering==0)
 //   {
  //    fclose(fd);
  //    fd=fopen(outfile,"w");
  //  }
 //  }

//}
//}//disp('while finsish steering');
//}//end //while finishsteering loop
//cufinish(&p,&w,&wnew,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1,  &d_wd);
cufinish(&p,&w,&wnew,&state,&d_p,&d_bp,&d_w,&d_wnew,&d_wmod, &d_dwn1,  &d_wd, &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);

#ifdef USE_MPI
     mpifinalize(p);
     cufinishmpi(&p,&w, &wmod, &temp2,&gmpivisc,   &gmpiw, &gmpiwmod, &d_p,   &d_w, &d_wmod,&d_wtemp2,    &d_gmpivisc,   &d_gmpiw, &d_gmpiwmod);
#endif
free(hlines);
free(p);
free(bp);
free(sdir);
free(name);
free(outfile);
free(formfile);

//for the completed simulation
//int nsteps=nt;
//fclose(fd);






	//[consts,domain,source]=loadsim('test1_16_02_09.xml',elist);
	//chdir(metadata.directory);
        //readsimulation_(elist.elist.id,simfile,elist.elist.port,elist.elist.server);
	//runsim(consts,dom,src,meta,simfile,elist);
#ifdef USE_IOME
	writesimulation_(elist.id,newsimfile,elist.port,elist.server);
#endif


	return 0;
}

