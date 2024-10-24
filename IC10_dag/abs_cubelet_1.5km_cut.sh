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
start_velocity=$4
n_chan=$5
src_name=$6

full_path=/projects/vla-processing/measurement_sets/${src_name}/${ms_name}
output_path=/projects/vla-processing/images/${src_name}/Absorption/pipeline_rerun

## define input values
output_name=${output_path}"/"${ra_phase_center}"_"${dec_phase_center}_"1.5km_cut"

# make casa call to imaging script
casa --logfile $output_name".log" -c abs_cubelet_1.5km_cut.py -v ${full_path} -o ${output_name} -r ${ra_phase_center} -d${dec_phase_center} -n ${n_chan} -s ${start_velocity}
