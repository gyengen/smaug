// iosac2.5d.cpp : Main routine for GPU enabled SAC
#include "../include/iosac2.5d.h"
#include "../include/smaugcukernels.h"
#include "../include/iobparams.h"

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

char configinfile[300];
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



printf("calling cuinit\n");

   char ext[3];
   char tcfg[300];
   char stemp[300];
   char *pch1,*pch2;
   strcpy(stemp,cfgfile);
   pch1 = strtok (stemp,".");
   sprintf(tcfg,"%s",pch1);
   pch2 = strtok (NULL,".");

   sprintf(ext,"%s",pch2);

sprintf(configfile,"%s",cfgout);
#ifdef USE_MPI
     MPI::Init(argc, argv);
     mpiinit(p);
     ipe2iped(p);     
     mpineighbours(0,p);
     mpineighbours(1,p);
     #ifdef USE_SAC3D
          mpineighbours(2,p);
     #endif

     #ifdef USE_SAC3D
	      if(p->ipe>99)
		sprintf(configinfile,"%s_np%d%d%d_%d.%s",tcfg,p->pnpe[0],p->pnpe[1],p->pnpe[2],p->ipe,ext);
	      else if(p->ipe>9)
		sprintf(configinfile,"%s_np0%d0%d0%d_0%d.%s",tcfg,p->pnpe[0],p->pnpe[1],p->pnpe[2],p->ipe,ext);
	      else
		sprintf(configinfile,"%s_np00%d00%d00%d_00%d.%s",tcfg,p->pnpe[0],p->pnpe[1],p->pnpe[2],p->ipe,ext);  	     
     #else
	      if(p->ipe>99)
		sprintf(configinfile,"%s_np%d%d_%d.%s",tcfg,p->pnpe[0],p->pnpe[1],p->ipe,ext);
	      else if(p->ipe>9)
		sprintf(configinfile,"%s_np%d%d_%d.%s",tcfg,p->pnpe[0],p->pnpe[1],p->ipe,ext);
	      else
		sprintf(configinfile,"%s_np0%d0%d_00%d.%s",tcfg,p->pnpe[0],p->pnpe[1],p->ipe,ext);  	     	     
     #endif
     printf("config files\n%s \n %s\n",configinfile,configfile);

     
#else
     sprintf(configinfile,"%s",cfgfile);
#endif

char *method=NULL;







       /*********************************************************************************************************/
       /* Start of section initialising steering and auto metadata collection*/
       /*********************************************************************************************************/


        //printf("cfgfile %s\n",configfile);
        //   getintparam_( &elist.id,"i1",&it,  &elist.port, "localhost" );	
        //	printf("Get integer %d\n",it);
        //Set input filename as first arg
	//if NULL use defaults
	
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


       /*********************************************************************************************************/
       /* End of section initialising steering and auto metadata collection*/
       /*********************************************************************************************************/










       /*********************************************************************************************************/
       /* Start of section creating arrays on the host*/
       /*********************************************************************************************************/
          printf("allocating w and wnew\n");
	  #ifdef USE_SAC_3D
		 w=(real *)calloc(ni*nj*nk*NVAR,sizeof(real ));
		 wmod=(real *)calloc(ni*nj*nk*NVAR,sizeof(real ));

		 wd=(real *)calloc(ni*nj*nk*NDERV,sizeof(real ));
		 wnew=(real *)calloc(ni*nj*nk*NVAR,sizeof(real ));
		 for(i=0;i<ni*nj*nk*NDERV;i++)
		    wd[i]=0.0;
	 #else
		 w=(real *)calloc(ni*nj*NVAR,sizeof(real ));
		 wmod=(real *)calloc(ni*nj*NVAR,sizeof(real ));

		 wd=(real *)calloc(ni*nj*NDERV,sizeof(real ));
		 wnew=(real *)calloc(ni*nj*NVAR,sizeof(real ));
		 for(i=0;i<ni*nj*NDERV;i++)
		    wd[i]=0.0;
	 #endif



        #ifdef USE_MPI
          int szw,szw0,szw1,szw2,szvisc0,szvisc1,szvisc2;
	  #ifdef USE_SAC
		  szw=4*(  ((p)->n[1])  +  ((p)->n[0])   );
		  szw0=4*NVAR*(  ((p)->n[1])     );
		  szw1=4*NVAR*(  ((p)->n[0])     );

		  szvisc0=4*NVAR*(  (((p)->n[1])+2 )   );
		  szvisc1=4*NVAR*(    (((p)->n[0]) +2 )  );
	  #endif
	  #ifdef USE_SAC_3D	  
		  szw=4*NVAR*(  ((p)->n[1])*((p)->n[2])  +  ((p)->n[0])*((p)->n[2])  +  ((p)->n[0])*((p)->n[1])  );
		  szw0=4*NVAR*(  ((p)->n[1])*((p)->n[2])    );
		  szw1=4*NVAR*(    ((p)->n[0])*((p)->n[2])   );
		  szw2=4*NVAR*(    ((p)->n[0])*((p)->n[1])  );

		  szvisc0=4*NVAR*(  (((p)->n[1])+2)*(((p)->n[2])+2)  ); 
		  szvisc1=4*NVAR*(   (((p)->n[0])+2)*(((p)->n[2])+2)    );    
		  szvisc2=4*NVAR*(  (((p)->n[1])+2)*(((p)->n[2])+2)   );   
	  #endif


	  #ifdef USE_SAC
	  temp2=(real *)calloc(NTEMP2*(((p)->n[0])+2)* (((p)->n[1])+2),sizeof(real));
	  #endif
	  #ifdef USE_SAC_3D
	  temp2=(real *)calloc(NTEMP2*(((p)->n[0])+2)* (((p)->n[1])+2)* (((p)->n[2])+2),sizeof(real));
	  #endif


	  gmpiwmod0=(real *)malloc(szw0*sizeof(real));
	  gmpiw0=(real *)malloc(szw0*sizeof(real));
	  gmpiwmod1=(real *)malloc(szw1*sizeof(real));
	  gmpiw1=(real *)malloc(szw1*sizeof(real));

          gmpivisc0=(real *)malloc(szvisc0*sizeof(real));
          gmpivisc1=(real *)malloc(szvisc1*sizeof(real));

	  #ifdef USE_SAC_3D
		  gmpiwmod2=(real *)malloc(szw2*sizeof(real));
		  gmpiw2=(real *)malloc(szw2*sizeof(real));
                  gmpivisc2=(real *)malloc(szvisc2*sizeof(real));


 

	  #endif
        #endif
       /*********************************************************************************************************/
       /* End of section creating arrays on the host*/
       /*********************************************************************************************************/








	//set initial time step to a large value
	if(p->moddton==1.0)
	{
		p->dt=1.0e-8;
	}








       /*********************************************************************************************************/
       /* Start of section initialising the configuration 
          on the host and on GPU host memory*/
       /*********************************************************************************************************/
	if((p->readini)==0)
	 initconfig(p, &meta, w);
	else
	 readasciivacconfig(configinfile,*p,meta,w,wd,hlines);

	printf("after read\n");
	p->it=0;
 
        real *u,  *v,  *h;

	h=w+(ni)*(nj)*rho;
	u=w+(ni)*(nj)*mom1;
	v=w+(ni)*(nj)*mom2;

	cuinit(&p,&bp,&w,&wnew,&wd,&state,&d_p,&d_bp,&d_w,&d_wnew,&d_wmod, &d_dwn1,  &d_wd, &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);

        //same as the grid initialisation routine in SAC
        //ensures boundaries defined correctly
	initgrid(&p,&w,&wnew,&state,&wd,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1,  &d_wd, &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);

	#ifdef USE_MPI
	  //initialise the mpi used memory locations
	 cuinitmpibuffers(&p, &w, &wmod, &temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,   &gmpiw0, &gmpiwmod0,   &gmpiw1, &gmpiwmod1,   &gmpiw2, &gmpiwmod2, &d_p, &d_w, &d_wmod,&d_wtemp2,  &d_gmpivisc0,  &d_gmpivisc1,  &d_gmpivisc2, &d_gmpiw0, &d_gmpiwmod0, &d_gmpiw1, &d_gmpiwmod1, &d_gmpiw2, &d_gmpiwmod2);
	 cucopywtompiw(&p,&w, &wmod,    &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,  &d_w, &d_wmod,   &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2, 0);




	#endif
       /*********************************************************************************************************/
       /* End of section initialising the configuration */
       /*********************************************************************************************************/











       /*********************************************************************************************************/
       /* Start of section initialising the configuration */
       /*********************************************************************************************************/
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
	  mpibound(NVAR, w ,p);
	#endif
	  cuboundary(&p,&bp,&d_p,&d_bp,&d_state,&d_wmod, 0,0,0);
	#ifdef USE_MPI
	   mpibound(NVAR, wmod ,p);
	   cucopywfrommpiw(&p,&w, &wmod,      &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,  &d_w, &d_wmod,    &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2,0);
	#endif
       /*********************************************************************************************************/
       /* End of section initialising the configuration */
       /*********************************************************************************************************/
	printf("after cuinit\n");


       





       /*********************************************************************************************************/
       /* Start of section getting parameters from metadata file */
       /*********************************************************************************************************/

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
        /*********************************************************************************************************/
        /* End of section getting parameters from metadata file */
        /*********************************************************************************************************/







        /*********************************************************************************************************/
        /* Start of section to set steering control (n.b. commented out) */
        /*********************************************************************************************************/     
	//createlog(meta.log_file);
	//while(finishsteering == 0)
	//{
	 
	  //  if( steeringenabled==0)
	  //    finishsteering=1;
        /*********************************************************************************************************/
        /* End of section to enable steering control */
        /*********************************************************************************************************/     







	int n;
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







       /*********************************************************************************************************/
       /* Start looping over iterations*/
       /*********************************************************************************************************/
	for( n=1;n<=nt;n++)
	{
	    p->it=n;
	    if(((n-1)%(p->cfgsavefrequency))==0)
	    {
			//writeconfig(name,n,*p, meta , w);
		#ifndef USE_MPI
			// writevtkconfig(configfile,n,*p, meta , w);
		#endif
			//writeasciivacconfig(configfile,*p, meta , w,wd,hlines,*state);
		writevacconfig(configfile,n,*p, meta , w,wd,*state);	 
     
	    }
	    order=0;
	    t1=second();

            printf("finished write routine\n");


       /*********************************************************************************************************/
       /* Start single step  iteration*/
       /*********************************************************************************************************/
	   if(p->moddton==1.0)
	   {
		p->maxcourant=0.0;
		courantmax=0.0;
		for(int dim=0; dim<=(NDIM-1); dim++)
		{
			cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
			cucomputemaxcourant(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
		}
	
		if(     ((  (p->courant)/(p->maxcourant)  ))>1.0e-8  )
		       p->dt=(p->courant)/(p->maxcourant);
		printf("new dt is %g %g\n",(p->courant)/(p->maxcourant),p->dt);

		if(n>1)
		   cugetdtvisc1(&p,&d_p,&d_wmod, &wd,&d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2);
		#ifdef USE_MPI
		   mpiallreduce(&(p->maxviscoef), MPI_MAX);
		#endif

		printf("dtdiffvisc %20.10g  %20.10g\n",p->maxviscoef,p->dtdiffvisc);
		if(1/(p->dtdiffvisc)>1.0e-8 && (p->dt)>((p->dtdiffvisc)) )
			                      			p->dt=(p->dtdiffvisc);
		printf(" modified dt is %20.10g \n",p->dt);

		//include gravitational modification
	   } 
       /*********************************************************************************************************/
       /* End of single step  iteration*/
       /*********************************************************************************************************/






       /*********************************************************************************************************/
       /* Start single step  iteration*/
       /*********************************************************************************************************/
	if((p->rkon)==0)
	{
	  ordero=0;
	 
	  cucomputedervfields(&p,&d_p,&d_wmod, &d_wd,order);
	  order=1;

	 for(int dir=0;dir<NDIM; dir++)
	 {
		  for(int f=rho; f<=(mom1+NDIM-1); f++)
		  { 
		      if((f==mom1 && dir==0)  ||  (f==mom2 && dir==1)  || (f==mom2 && dir==2) )
		       cucomputept(&p,&d_p,&d_wmod, &d_wd,order,dir);
		      cucentdiff1(&p,&d_p,&d_state,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt,f,dir);	     
		  } //end looping over fields for cucentdiff1
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
		   }//end looping over fields for cucentdiff2
		   #endif
	  }//end loop over directions
	   
	 
	  cugrav(&p,&d_p,&d_state,&d_w,&d_wmod, &d_dwn1, &d_wd,order, ordero,p->dt);//gravitational contributions

	  if(p->divbon==1)
		       cudivb(&p,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt);


           /*********************************************************************************************************/
           /* Start  of hyperdiffusion contributions for single step*/
           /*********************************************************************************************************/
	   if(p->hyperdifmom==1)
	   {
	    dt=(p->dt);
		     p->maxviscoef=0.0;
	    
	    #ifdef USE_SHOCKVISC
	       cunushk1(&p,&d_p,&d_wmod, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2);
	    #endif
            //density hyperdiffusion term
	    for(int dim=0; dim<=(NDIM-1); dim++)
	    {
		      cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
		      cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
		      #ifdef USE_MPI
			      mpiallreduce(&(p->cmax), MPI_MAX);
		      #endif
		      cmax[dim]=p->cmax;
		      cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);
		      #ifdef USE_MPI
			  cucopytompivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
			   mpivisc(dim,p,temp2);
			  ;//cucopyfrommpivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
		      #endif
		      cuhyperdifvisc1r(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);
		      cuhyperdifvisc1l(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);

	              cuhyperdifrhosource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,rho,dim,p->dt);
	     } //end of rho hyperdif contributions for each direction


             //energy hyperdiffusion term
	     for(int dim=0; dim<=(NDIM-1); dim++)
	     {
	        cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
		cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
		   
		#ifdef USE_MPI
		      mpiallreduce(&(p->cmax), MPI_MAX);
		#endif
	        cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);

		#ifdef USE_MPI
			cucopytompivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
			mpivisc(dim,p,temp2);
			cucopyfrommpivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
		#endif

		cuhyperdifvisc1r(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);
		cuhyperdifvisc1l(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);	       
	        cuhyperdifesource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,energy,dim,dt);   
	     }

        //momentum hyperdiffusion term
	for(int dim=0; dim<=(NDIM-1); dim++)
	       for(int f=0; f<=(NDIM-1); f++)		   	                 
		{
			  cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			  cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);

			  #ifdef USE_MPI
			      mpiallreduce(&(p->cmax), MPI_MAX);
			  #endif
		          cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);
		          #ifdef USE_MPI
				  
                                  cucopytompivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
				  mpivisc(dim,p,temp2);
				  cucopyfrommpivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
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
		     }  //end of loop over dimensions and fields


                     //b field hyperdiffusion term
		     int jj,mm,kk;
		     real sb;
		     for(int dim=0; dim<=(NDIM-1); dim++)
		     for(int f=0; f<=(NDIM-1); f++) 
		     if(f!=dim)           
		     {
			       cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			       cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);

			       #ifdef USE_MPI
			      		mpiallreduce(&(p->cmax), MPI_MAX);
			       #endif
			       cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);

			       #ifdef USE_MPI
				  cucopytompivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
				  mpivisc(dim,p,temp2);
				  cucopyfrommpivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);		 
		               #endif
			       cuhyperdifvisc1r(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);
			       cuhyperdifvisc1l(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);

                                
				for(ii1=0;ii1<=1;ii1++)  //start of compute cross-product term of b field hyperdiffusion terms
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
				}//end of compute cross-product term of b field hyperdiffusion terms
		      }//end of loop over fields and dimensions

	   }//closes if(p->hyperdifmom==1)
           /*********************************************************************************************************/
           /* End of hyperdiffusion contributions for single step */
           /*********************************************************************************************************/

	  //source terms
          cusource(&p,&d_p,&d_state,&d_w,&d_wmod, &d_dwn1, &d_wd,order, ordero,p->dt);

	  cuboundary(&p,&bp,&d_p,&d_bp,&d_state,&d_wmod, ordero,0,0);

	} //end of if((p->rkon)==0)
       /*********************************************************************************************************/
       /* End single step  iteration*/
       /*********************************************************************************************************/




       /*********************************************************************************************************/
       /* Start runge-kutta  iteration*/
       /*********************************************************************************************************/

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

			for(int f=rho; f<=mom1+(NDIM-1); f++)//looping over fields for cucentdiff1
			cucentdiff1(&p,&d_p,&d_state,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,dt,f,dir);

			#ifndef ADIABHYDRO
				   for(int f=energy; f<=b1+(NDIM-1); f++)//looping over fields for cucentdiff2
				       cucentdiff2(&p,&d_p,&d_state,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt,f,dir);

			#endif
		    }

		   if(p->divbon==1)
		       cudivb(&p,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd, order,ordero,p->dt);

		   /*********************************************************************************************************/
		   /* End of hyperdiffusion contributions for single step */
		   /*********************************************************************************************************/
		   if(p->hyperdifmom==1)
		   {
		     p->maxviscoef=0.0;
		     #ifdef USE_SHOCKVISC        
			   cunushk1(&p,&d_p,&d_wmod,&d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2);
		     #endif
			
		     //density hyperdiffusion term
		     for(int dim=0; dim<=(NDIM-1); dim++)
		     {
			cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
			#ifdef USE_MPI
				mpiallreduce(&(p->cmax), MPI_MAX);
			#endif
		  
			cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);
			#ifdef USE_MPI
				cucopytompivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
				mpivisc(dim,p,temp2);
				cucopyfrommpivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
			#endif
			cuhyperdifvisc1r(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);
			cuhyperdifvisc1l(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);				
			cuhyperdifrhosource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,rho,dim,dt);
		     }//end of loop over  dimensions

		     //energy hyperdiffusion term
		     for(int dim=0; dim<=(NDIM-1); dim++)
		     {
			cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
			#ifdef USE_MPI
				mpiallreduce(&(p->cmax), MPI_MAX);
			#endif
			cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);
			#ifdef USE_MPI
			        cucopytompivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
				mpivisc(dim,p,temp2);
				cucopyfrommpivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
			#endif
			cuhyperdifvisc1r(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);
			cuhyperdifvisc1l(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);
			cuhyperdifesource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,energy,dim,dt);
		     }//end of loop over dimensions
		     
		     //momentum hyperdiffusion term
		     for(int dim=0; dim<=(NDIM-1); dim++)
		     for(int f=0; f<=(NDIM-1); f++)				   	                 
		     {
			cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			#ifdef USE_MPI
				;//   mpiallreduce(&(p->cmax), MPI_MAX);
			#endif
			cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);
			#ifdef USE_MPI
                                cucopytompivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
				mpivisc(dim,p,temp2);
				cucopyfrommpivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
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
		    }//end of loop over fields and dimensions

		     //b field hyperdiffusion term
		    int jj,mm,kk;
		    real sb;
		    for(int dim=0; dim<=(NDIM-1); dim++)
		    for(int f=0; f<=(NDIM-1); f++) 
		    if(f!=dim)                     
		    {
			cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			#ifdef USE_MPI
				;//mpiallreduce(&(p->cmax), MPI_MAX);
			#endif
			cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);
				//cuhyperdifvisc1il(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);
			#ifdef USE_MPI
				cucopytompivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
				mpivisc(dim,p,temp2);
				cucopyfrommpivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
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
		   } //closes if(f!=dim) and end of loop over fields and dimensions

		   }//closes if(p->hyperdifmom==1)

		   cuadvance(&p,&d_p,&d_wmod,&d_w,order);
		   #ifdef USE_MPI
			cucopywtompiw(&p,&w, &wmod,    &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,  &d_w, &d_wmod,    &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2, order);
			mpibound(NVAR, d_wmod ,d_p);
			cucopywfrommpiw(&p,&w, &wmod,   &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,  &d_w, &d_wmod,    &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2,order);

		   #endif
		   cuboundary(&p,&bp,&d_p,&d_bp,&d_state,&d_wmod, orderb,0,0);
		   

	   }//looping over orders
       /*********************************************************************************************************/
       /* End runge-kutta  iteration*/
       /*********************************************************************************************************/


	   p->it=n+1;

	 cuupdate(&p,&w,&wmod,&temp2,&state,&d_p,&d_w,&d_wmod,&d_wtemp2,  &d_state,n);

	  #ifdef USE_MPI
		   cucopywtompiw(&p,&w, &wmod,    &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,  &d_w, &d_wmod,   &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2, order);
		   mpibound(NVAR, w ,p);
		   cucopywfrommpiw(&p,&w, &wmod,    &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,  &d_w, &d_wmod,   &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2,order);		   
	  #endif

	printf("\n");

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

           // mpisync();
	    }//end of testep
       /*********************************************************************************************************/
       /* End of looping over iterations*/
       /*********************************************************************************************************/





	cufinish(&p,&w,&wnew,&state,&d_p,&d_bp,&d_w,&d_wnew,&d_wmod, &d_dwn1,  &d_wd, &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);

	#ifdef USE_MPI
	     
	     cufinishmpi(&p,&w, &wmod, &temp2,&gmpivisc0,&gmpivisc1,&gmpivisc2,   &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,   &d_w, &d_wmod,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2,   &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2);
             mpifinalize(p);


	#endif
	free(hlines);
	free(p);
	free(bp);
	free(sdir);
	free(name);
	free(outfile);
	free(formfile);


	#ifdef USE_IOME
		writesimulation_(elist.id,newsimfile,elist.port,elist.server);
	#endif


		return 0;
	}

