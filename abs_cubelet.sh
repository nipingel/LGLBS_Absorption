#!/bin/bash
#
# abs_cubelet.sh
# execution script for imaging of sub cube around absorption source

## untar CASA
mv abs_cubelet.py tmp
cd tmp
cp /scratch/casa-6.5.0-15-py3.8.tar.xz .
xz -d casa-6.5.0-15-py3.8.tar.xz
tar -xvf casa-6.5.0-15-py3.8.tar

## change home directory so CASA will run
HOME=$PWD

## untar measurement set
ms_name="20A-346.sb41042474.eb41074466.59579.84285818287.speclines.ms.split_chan"$1

## copy measurement set to working directory
cp /projects/vla-processing/measurement_sets ./

## define input values
output_name="37W051_chan$"1

# make casa call to imaging script
casa-6.5.0-15-py3.8/bin/casa --logfile 37W051_chan$1.log -c abs_cubelet.py -v $ms_name -s $1 -o $output_name

tar -cvf $output_name"_output.tar" $output_name.*

mv $output_name"_output.tar" /projects/vla-processing/measurement_sets

## clean up
rm -rf $ms_name
