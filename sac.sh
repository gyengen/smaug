#!/bin/bash

#remove ## markers from lines commented with ## to enable iome server

if [ $# -gt 0 ]; then
  IOME_SIMNAME=$1
else
  IOME_SIMNAME="mysim"
fi

if [ $# -gt 1 ]; then
  SIMFILE=$2
  echo "simfile is" $SIMFILE
fi


#iome the following lines comment using ##
##iogs initiome null $IOME_SIMNAME null >& iogs.err &
##sleep 3
##INPUT=`cat ${IOME_SIMNAME}0_port.txt`
##IOME_WSPORT=$(echo $INPUT | cut -d' ' -f1 )
##echo port is $IOME_WSPORT
##end of iome comment lines


if [ $# -gt 1 ]; then
#  iogs readsimulation simfile.xml 0 $IOME_WSPORT localhost
  bin/iosac $IOME_SIMNAME $SIMFILE
else
  bin/iosac $IOME_SIMNAME
fi

#iome the following lines comment using ##
##iogs exitiome 0 $IOME_WSPORT 
##end of iome comment lines

