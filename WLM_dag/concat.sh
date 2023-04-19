#!/bin/bash
#
# statwt_indv.sh
# execution script reweight visibilities based on rms in emission-free channels

## change home directory so CASA will run
HOME=$PWD

## define variables
extension=$1
output_vis_name=$2
src_name=$3
full_path=/projects/vla-processing/measurement_sets/${src_name}/raw_measurement_sets

## make call to casa
/casa-6.5.0-15-py3.8/bin/casa --nologfile -c concat.py -e ${extension} -p ${full_path} -o ${output_vis_name}