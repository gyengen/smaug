/*
IOME LICENSE
IOME Version 1.1.1

IOME Development  Tools
Copyright (C) 2001-2004, Michael Kenneth Griffiths, All Rights Reserved.

--------------------------------------------------------------------------------
IOME public license.

The contents of this file are subject to the IOME Public License Version 1.3
(the "License"); you may not use this file except in compliance with the
License. You may obtain a copy of the License at
http://81.174.178.112/iome/licensing/iomelicense.html
Software distributed under the License is distributed on an "AS IS" basis,
WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License
for the specific language governing rights and limitations under the License.

The Initial Developer of the Original Code is Michael Kenneth Griffiths.
Copyright (C) 2000-2004 Michael Kenneth Griffiths. All Rights Reserved.
--------------------------------------------------------------------------------
GPL license.

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 59 Temple
Place, Suite 330, Boston, MA 02111-1307 USA

Author contact information:
mikeg@photon0.freeserve.co.uk
--------------------------------------------------------------------------------
*/







#ifdef USE_IOME
	#include <iome/simulation/IoInitialiser.h>
	//#include <iome/simulation/soapH.h>
	#include <iome/genericsimulationlib/IoGenericSimulationLib.h>

	#include <iome/simulation/IoSteerWS.nsmap>
#endif
#ifndef HYPERDIF_H_
#define HYPERDIF_H_
//	#include <iome/simulation/stdsoap2.h>
    	#include <unistd.h>
	#include <sys/stat.h>
	#include <sys/types.h>
	#include <sys/wait.h>
        #include <sys/time.h>

#ifdef USE_MPI
        #include "smaugmpi.h"
#endif





//#include "stdafx.h"
//#include "IoTestDEVSimulation.h"
//#include "IoTestAgentSimulation.h"
#include <time.h>
#include <pthread.h>	// use Pthreads
//#include "soapH.h"
//#include "IoSteerWS.nsmap"
#include <iostream>

#include "paramssteeringtest1.h"
#include "readwrite.h"
#include "dxroutines.h"
#include "initialisation.h"

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



using std::cout;
using std::endl;
//#include "IoTestGenSimulation.h"
//#include <afxwin.h>         // MFC core and standard components
//*/
#ifdef _DEBUG
#define new DEBUG_NEW
#endif
#ifdef CWDEBUG
	using namespace libcwd;
#endif

char simfile[300];
char newsimfile[300];
char portfile[300];



int getintparam_( int id,char *sname,int *iv,  int port, char *sserver );
int m_isimfinished=0;
char m_serverclient[300] = "localhost:8080";
char m_hostname[300] = "localhost";
int m_port=8080;
int port=m_port;
void readsim(params *k,  meta *md,char *simfile, iome el);

void createsim(params k, meta metadata,char *simname, iome el);
void gendxgen(char *dir,char *jobname,int nsteps,int n1,int n2);

#endif

