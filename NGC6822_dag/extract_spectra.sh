#!/bin/bash
#
# extract_spectra.sh
# execution script to extract emission/absorption spectra and produce analysis products

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

## get python distribution
wget https://repo.anaconda.com/archive/Anaconda3-2023.07-1-Linux-x86_64.sh -O ~/anaconda.sh
bash ~/anaconda.sh -b -p $HOME/anaconda3                                                  
~/anaconda3/bin/conda clean -all -y                                                       
## create conda env for astro tools                                                                    
~/anaconda3/bin/conda create -c conda-forge -y -n astro_env astropy pip numpy scipy matplotlib pip
source ~/anaconda3/etc/profile.d/conda.sh
conda activate astro_env
pip install spectral_cube
pip install AegeanTools

## copy the data products to working area
cp /projects/vla-processing/images/${src_name}/${src_name}/Absorption/${ra_phase_center}_${dec_phase_center}.image.pbcor.fits .
cp /projects/vla-processing/images/${src_name}/${sd_cube_name} .

## run script to pull out a 30' subregion around source & feather in single dish
/casa-6.5.0-15-py3.8/bin/casa --nologfile --log2term --nogui -c combine_SD.py -p ${vla_cube_path} -s {sd_cube_name} -r ${ra_phase_center} -d${dec_phase_center}

# run analysis script that extracts spectra
python3 spectral_extraction -n ${ra_phase_center}_{$dec_phase_center}.image.pbcor.fits -c VLA_ABCD_GBT_${ra_phase_center}_${dec_phase_center}_30arcmin.fits -r ${ra_phase_center} -d${dec_phase_center}

## tar output
tar -cvf ${ra_phase_center}_${dec_phase_center}_analysis_products.tar *.csv *.pickle *.pdf
mv *.tar /projects/vla-processing/images/${src_name}/Absorption
