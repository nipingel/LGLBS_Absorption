#!/bin/bash
#
# ms_transform.sh
# execution script to run ms_transform 

## change home directory so CASA will run
HOME=$PWD

## assign inputs & script variables
input_name=$1
src_name=$2
full_path=/projects/vla-processing/measurement_sets/${src_name}/raw_measurement_sets


# make casa call to imaging script
/casa-6.5.0-15-py3.8/bin/casa --nologfile -c ms_transform_combine_spw.py -p ${full_path} -n ${input_name}