#include "../include/iotypes.h"

#ifdef USE_IOME
#include <iome/genericsimulationlib/IoGenericSimulationLib.h>

void createsim(params k,  meta metadata,char *simname, iome el)
{
int i;
int ni,nj;
double xmax,ymax;
double dx,dy,dt,tmax,wavespeed;
double courant;

//elist=list();  parameter used by iome to contain port and server address
//elist=list();
//printf("createsim1\n");
//printf("author %s %d %s\n",metadata.author,el.port,el.server);
//it   t   dt    rho m1 m2 e bx by
addmetadata_(el.id,"author",metadata.author,el.port,el.server);
addmetadata_(el.id,"directory",metadata.directory,el.port,el.server);
addmetadata_(el.id,"date",metadata.sdate,el.port,el.server);
addmetadata_(el.id,"platform",metadata.platform,el.port,el.server);
addmetadata_(el.id,"description",metadata.desc,el.port,el.server);
addmetadata_(el.id,"name",metadata.name,el.port,el.server);
addmetadata_(el.id,"ini_file",metadata.ini_file,el.port,el.server);
addmetadata_(el.id,"log_file",metadata.log_file,el.port,el.server);
addmetadata_(el.id,"output_file",metadata.out_file,el.port,el.server);

// Constants
//adddoubleparam_(el.id,"g",k.g,7,el.port,el.server);
//adddoubleparam_(el.id,"u0",k.u0,7,el.port,el.server);
//adddoubleparam_(el.id,"v0",k.v0,7,el.port,el.server);
//adddoubleparam_(el.id,"b",k.b0,7,el.port,el.server);
//adddoubleparam_(el.id,"h0",k.h0,7,el.port,el.server);

//Domain definition
// Define the x domain
//ni = 151; 
ni=k.n[0];
xmax = k.xmax[0];                      
k.dx[0] = xmax/(ni-1);
//x  = [0:dx:xmax];

// Define the y domain
//nj = 151;  
nj=k.n[1];
ymax = k.xmax[1];                      
k.dx[1] = ymax/(nj-1);
//y  = [0:dy:ymax];

tmax = k.tmax;

// Define the wavespeed
dt=k.dt;
//wavespeed = k.u0 + sqrt(k.g*(k.h0 - k.b0));

// Define time-domain
//dt = 0.68*dx/wavespeed;

//t = [0:dt:tdomain];
//t=[1:dt:tmax];
k.nt=(int)((tmax-1)/dt);
courant = wavespeed*dt/dx;

adddoubleparam_(el.id,"ni",k.n[0],7,el.port,el.server);
adddoubleparam_(el.id,"nj",k.n[1],7,el.port,el.server);
adddoubleparam_(el.id,"xmax",k.xmax[0],7,el.port,el.server);
adddoubleparam_(el.id,"ymax",k.xmax[1],7,el.port,el.server);
adddoubleparam_(el.id,"tmax",k.tmax,7,el.port,el.server);
addintparam_(el.id,"nt",k.nt,7,el.port,el.server);
addintparam_(el.id,"steeringenabled",k.steeringenabled,7,el.port,el.server);
addintparam_(el.id,"finishsteering",k.finishsteering,7,el.port,el.server);
addintparam_(el.id,"step",k.dt,7,el.port,el.server);







//simfile=sprintf('%s.xml',simname)



//endfunction
}

void readsim(params *k,  meta *md,char *simfile, iome el)
{
    //      readsimulation_(el.id,simfile,el.port,el.server);
double ni,nj;
double xmax,ymax;
double dx,dy,tmax;
int nt,steeringenabled, finishsteering;


getmetadata_(el.id,"author",&(md->author),el.port,el.server);
getmetadata_(el.id,"directory",&(md->directory),el.port,el.server);
getmetadata_(el.id,"date",&(md->sdate),el.port,el.server);
getmetadata_(el.id,"platform",&(md->platform),el.port,el.server);
getmetadata_(el.id,"description",&(md->desc),el.port,el.server);
getmetadata_(el.id,"name",&(md->name),el.port,el.server);
getmetadata_(el.id,"ini_file",&(md->ini_file),el.port,el.server);
getmetadata_(el.id,"log_file",&(md->log_file),el.port,el.server);
getmetadata_(el.id,"out_file",&(md->out_file),el.port,el.server);

getdoubleparam_(el.id,"ni",&ni,el.port,el.server);
getdoubleparam_(el.id,"nj",&nj,el.port,el.server);
getdoubleparam_(el.id,"xmax",&xmax,el.port,el.server);
getdoubleparam_(el.id,"ymax",&ymax,el.port,el.server);
getintparam_(&el.id,"nt",&nt,&el.port,el.server);
getintparam_(&el.id,"steeringenabled",&steeringenabled,&el.port,el.server);
getintparam_(&el.id,"finishsteering",&finishsteering,&el.port,el.server);




k->n[0]=ni;
k->n[1]=nj;
k->xmax[0]=xmax;
k->xmax[1]=ymax;
k->nt=nt;
k->steeringenabled=steeringenabled;
k->finishsteering=finishsteering;


  printf("read metadata\n");

}

#endif

void initconfig(params *k, meta *md, real *w)
{
	int i1,j1;
        int ni=k->n[0];
        int nj=k->n[1];
        for(i1=0; i1<(k->n[0]) ;i1++)
	  for(j1=0; j1<(k->n[1]) ;j1++)
          {
                    for(int f=rho; f<NVAR; f++)
                    {

                    switch(f)
		            {
		              case rho:
		            	w[j1*ni+i1+(ni*nj*f)]=1.0;
			      break;
		              case mom1:
		            	w[j1*ni+i1+(ni*nj*f)]=0.01;
			      break;
		              case mom2:
		            	w[j1*ni+i1+(ni*nj*f)]=0.01;
			      break;
		              //case mom3:
		            	//w[j1*ni+i1+(ni*nj*f)]=0.0;
			      //break;
		              case energy:
		            	w[j1*ni+i1+(ni*nj*f)]=0.0;
			      break;
		              case b1:
		            	w[j1*ni+i1+(ni*nj*f)]=0.0;
			      break;
		              case b2:
		            	w[j1*ni+i1+(ni*nj*f)]=0.0;
			      break;
		              //case b3:
		            //	w[j1*ni+i1+(ni*nj*f)]=0.0;
			     // break;
		            }; //end of switch to check for field

			}//end of loop over f

          }//end of loop over j and i


}
