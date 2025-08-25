#!/bin/bash
#
# untar_measurement_sets.sh
# execution script to untar input measurement set 

src_name=$1
file_name=$2

cd /projects/vla-processing/measurement_sets/${src_name}/
## untar 
tar -xvf ${file_name}
