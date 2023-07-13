#!/bin/bash
#
# image_channel.sh
# execution script for imaging single channel of LGLBS sources 

## change home directory so CASA will run
HOME=$PWD

## assign variables
ms_name=$1
src_name=$2
ra_phase_center="19h44m56.6s" 
dec_phase_center="-14d47m21s"
output_name=${ms_name}"_nbcont"
full_path=/projects/vla-processing/measurement_sets/${src_name}/${ms_name}
# make casa call to imaging script
/casa-6.5.0-15-py3.8/bin/casa --logfile ${output_name}".log" -c narrowband_cont_image.py -v ${full_path} -o ${output_name} -r ${ra_phase_center} -d${dec_phase_center}

tar -cvf ${output_name}".tar" ${output_name}*
mv ${output_name}".tar" /projects/vla-processing/images/${src_name}
