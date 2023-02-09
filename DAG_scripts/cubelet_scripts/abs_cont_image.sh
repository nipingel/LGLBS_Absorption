#!/bin/bash
#
# abs_cont_image.sh
# execution script for imaging of sub cube around absorption source

## untar CASA
mv abs_cont_image.py tmp
cd tmp
cp /scratch/casa-6.5.0-15-py3.8.tar.xz .
xz -d casa-6.5.0-15-py3.8.tar.xz
tar -xvf casa-6.5.0-15-py3.8.tar

## change home directory so CASA will run
HOME=$PWD

## set measurement set names
ms_name="M31_B_20A-346.sb41042474.eb41074466.59579.84285818287.continuum.ms.split.tar"
untar_name="20A-346.sb41042474.eb41074466.59579.84285818287.continuum.ms.split"
source_name="M31"

## copy measurement set to working directory & untar
tar -xvf /projects/vla-processing/measurement_sets/$source_name/$ms_name --directory .

## define input values
output_name="M31_1_cont"


# make casa call to imaging script
casa-6.5.0-15-py3.8/bin/casa --logfile $output_name".log" -c abs_cont_image.py -v $untar_name -o $output_name 

tar -cvf $output_name"_output.tar" $output_name*

mv $output_name"_output.tar" /projects/vla-processing/images/$source_name

## clean up
rm $ms_name
rm -rf $untar_name
