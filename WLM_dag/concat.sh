#!/bin/bash
#
# statwt_indv.sh
# execution script reweight visibilities based on rms in emission-free channels

mv concat.py tmp
cd tmp
cp /scratch/casa-6.5.0-15-py3.8.tar.xz .
xz -d casa-6.5.0-15-py3.8.tar.xz
tar -xvf casa-6.5.0-15-py3.8.tar

## change home directory so CASA will run
HOME=$PWD

## define variables
output_vis_name=$1
src_name=$2

##  set measurement set name, copy and untar
cp /projects/vla-processing/measurement_sets/${src_name}/raw_measurement_sets*.transformed.tar .

## copy and untar measurement sets
for s in *.tar
do
	tar -xvf $s".tar" 
done

## make call to casa
casa-6.5.0-15-py3.8/bin/casa -c concat.py -o $output_vis_name

## pack up and copy back
tar -cvf ${output_vis_name}".tar" ${output_vis_name}

mv ${output_vis_name}".tar" /projects/vla-processing/measurement_sets/${src_name}/raw_measurement_sets

## clean up
#for s in ${ms_names}
#do 
#	rm -rf $s
#done
rm -rf ${output_vis_name}
