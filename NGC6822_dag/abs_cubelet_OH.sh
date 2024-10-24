#!/bin/bash
#
# abs_cubelet.sh
# execution script for imaging of sub cube around absorption source

## change home directory so CASA will run
HOME=$PWD

## set inputs
ms_name=$1
ra_phase_center=$2
dec_phase_center=$3
src_name=$4

full_path=/projects/vla-processing/measurement_sets/${src_name}/${ms_name}
output_path=/projects/vla-processing/images/${src_name}/Absorption/OH

## define input values
output_name=${output_path}"/"$2'_'$3

# make casa call to imaging script
casa --logfile $output_name".log" -c abs_cubelet_OH.py -v ${full_path} -o ${output_name} -r ${ra_phase_center} -d${dec_phase_center} 
