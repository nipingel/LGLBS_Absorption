#!/bin/bash
#
# extract_spectra.sh
# execution script to extract emission/absorption spectra and produce analysis products

## unpack user arguments
## set inputs
src_name=$1
ra_phase_center=$2
dec_phase_center=$3

## activate python environment
source /miniconda3/etc/profile.d/conda.sh
conda activate astro_env

## copy the data products to working area
cp /projects/vla-processing/images/${src_name}/Absorption/pipeline_rerun/${ra_phase_center}_${dec_phase_center}.image.pbcor.fits .
cp /projects/vla-processing/images/${src_name}/Absorption/pipeline_rerun/VLA_ABCD_GBT_${ra_phase_center}_${dec_phase_center}_30arcmin.fits .

# run analysis script that extracts spectra
python3 extract_spectra.py -n ${ra_phase_center}_${dec_phase_center}.image.pbcor.fits -c VLA_ABCD_GBT_${ra_phase_center}_${dec_phase_center}_30arcmin.fits

## tar output
tar -cvf ${ra_phase_center}_${dec_phase_center}_analysis_products.tar *.image.pbcor.mom0.fits *.csv *.pickle
mv *.tar /projects/vla-processing/images/${src_name}/Absorption/pipeline_rerun
rm *.pickle
rm *.fits
rm *.csv