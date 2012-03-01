Guidelines for Using sac_cuda Application for MHD

Introduction

The Sheffield Advanced Code (SAC) is a novel fully non-linear MHD code,
designed for simulations of linear and non-linear wave propagation in
gravitationally strongly stratified magnetised plasma.
Ref.http://adsabs.harvard.edu/abs/2008A%26A...486..655S 


Requirements

CUDA-Enabled Tesla GPU Computing Product with at least compute capability 1.3.
See  
  http://developer.nvidia.com/cuda-gpus

CUDA toolkit
http://developer.nvidia.com/cuda-toolkit-40

Correctly installed and compiler on user path



Installation

Create a directory and from that directory,
using a subversion client checkout the latest distribution using the command:

svn checkout http://ccpforge.cse.rl.ac.uk/svn/sac/dev/sac_cuda

when prompted, Password for 'anonymous' just press return.

Download distribution from
http://wrgrid.group.shef.ac.uk/downloads/sac_cuda.tgz



Building and running a Model

From the distribution base directory change directory to the src folder

To make the Brio-Wu test (use the following commands)
make clean
make bw
make sac

Change back to the distribution base directory
Create directories for output data e.g.
mkdir vtk  out

Run the model
./iosac.sh a
 
The test models available are as follows.
The code used with make to make the model is shown in the second column
1d Brio-Wu                 bw
2d Orszag-Tang             ot
2d (or 3d)Hypernova blast  bach
1d Alfven Wave             alf

Model parameters and configurations can be edited by
edit the file iosac2.5dparams.h
move to the src folder and recompile the model using
make sac

Move back to the base directory and run the model.
Note that the iosac2.5dparams.h file will be overwritten by one of the models
in the models folder when the command
make amodel

is executed.




Guidelines for Users Developing Customised Models






Help Support



========================================

* Please make sure your PATH includes /usr/local/cuda/bin
* Please make sure your LD_LIBRARY_PATH
*   for 32-bit Linux distributions includes /usr/local/cuda/lib
*   for 64-bit Linux distributions includes /usr/local/cuda/lib64:/usr/local/cuda/lib
* OR
*   for 32-bit Linux distributions add /usr/local/cuda/lib
*   for 64-bit Linux distributions add /usr/local/cuda/lib64 and /usr/local/cuda/lib
* to /etc/ld.so.conf and run ldconfig as root

* Please read the release notes in /usr/local/cuda/doc/

* To uninstall CUDA, delete /usr/local/cuda
* Installation Complete




















