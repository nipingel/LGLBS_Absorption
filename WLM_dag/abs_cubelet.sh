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

## set inputs
ms_name=$1
ra_phase_center=$2
dec_phase_center=$3
start_velocity=$4
n_chan=$5
src_name=$6
field_id=$7

untar_name=(${ms_name//.tar/ })

## copy measurement set to working directory & untar
tar -xvf /projects/vla-processing/measurement_sets/${src_name}/raw_measurement_sets/${ms_name} --directory .

## define input values
output_name=$2'_'$3'_nouv'

# make casa call to imaging script
casa-6.5.0-15-py3.8/bin/casa --logfile $output_name".log" -c abs_cubelet.py -v ${untar_name} -o ${output_name} -r ${ra_phase_center} -d${dec_phase_center} -n ${n_chan} -f ${field_id} -s ${start_velocity}

tar -cvf $output_name"_output.tar" $output_name*

mv $output_name"_output.tar" /projects/vla-processing/images/${src_name}

## clean up
rm -rf $untar_name
