#!/bin/bash
#$ -j y
#$ -V 
#$ -l arch=intel*
#$ -l gpu=1
#$ -P cs-test
##$ -P gpu
#$ -N jexp3d
#$ -l mem=12G
#$ -l rmem=12G
##$ -l h_rt=168:00:00
#$ -l h_rt=672:00:00

module add libs/cuda/4.0.17


cd include
cp iosmaugparams.h.sedov3d iosmaugparams.h
cd ..

cd src
cp ../models/usersource.cu.default usersource.cu
cp ../models/boundary.cu.sedov3d boundary.cu
make clean
make -f Makefile_3d smaug
cd ..



export TIMECOUNTER=0
source timeused
bin/smaug a
source timeused

