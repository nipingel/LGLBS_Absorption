#!/bin/bash
#
# ms_transform.sh
# execution script to run ms_transform 

## change home directory so CASA will run
HOME=$PWD

## assign inputs & script variables
ms_name=$1
chan_width=$2
output_frame=$3
src_name=$4
untar_name=(${ms_name//.tar/ })
output_ms_name=${untar_name}".transformed"

## untar measurement set to current working directory
tar -xvf /projects/vla-processing/measurement_sets/${src_name}/raw_measurement_sets/${ms_name} --directory .

# make casa call to imaging script
/casa-6.5.0-15-py3.8/bin/casa --logfile ${output_ms_name}".log" -c ms_transform.py -p ${untar_name} -w ${chan_width} -o ${output_ms_name} -f ${output_frame}

tar -cvf ${output_ms_name}".tar" ${output_ms_name}
mv ${output_ms_name}".tar" /projects/vla-processing/measurement_sets/${src_name}/raw_measurement_sets

## clean up
rm -rf ${untar_name}