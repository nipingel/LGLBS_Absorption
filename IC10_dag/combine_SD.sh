#!/bin/bash

# combine_SD.sh
# execution script to extract 30' region around absorption source & feather with GBT data

## unpack user arguments
## set inputs
src_name=$1
ra_phase_center=$2
dec_phase_center=$3
vla_cube_name=$4
sd_cube_name=$5
vla_cube_path=/projects/vla-processing/images/${src_name}/${vla_cube_name}

## set HOME variable
HOME=$PWD

## run script to pull out a 30' subregion around source & feather in single dish
casa --nologfile --log2term --nogui -c combine_SD.py -p ${vla_cube_path} -s ${sd_cube_name} -r ${ra_phase_center} -d${dec_phase_center}

mv VLA_ABCD_GBT_${ra_phase_center}_${dec_phase_center}_30arcmin.fits /projects/vla-processing/images/${src_name}/Absorption/pipeline_rerun
