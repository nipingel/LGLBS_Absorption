#!/bin/bash
#
# source_finding.sh
# execution script to run source finding on narrowband continuum image to identify strong
# background sources with which to measure absorption against

## unpack user arguments
image_name=$1
src_name=$2
int_threshold=$3

## set HOME variable
HOME=$PWD

## activate python environment
source /miniconda3/etc/profile.d/conda.sh
conda activate astro_env

## copy narrowband image
cp /projects/vla-processing/images/${src_name}/Absorption/${image_name}.image.fits .

## now, run source finding
aegean ${image_name}.image.fits --out ${image_name}_source_catalog.csv

## parse log
python parse_log.py -n ${image_name}_source_catalog.csv -t ${int_threshold}

## phase-center-list.txt will be moved back to home directory after job completion, but clean up other files
rm ${image_name}.image.fits
