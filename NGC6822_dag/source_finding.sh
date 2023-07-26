#!/bin/bash
#
# source_finding.sh
# execution script to run source finding on narrowband continuum image to identify strong
# background sources with which to measure absorption against

## unpack user arguments
image_name=$1
src_name=$2
int_threshold=$3

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

## untar narrowband continuum products to working directory
tar -xvf /projects/vla-processing/images/${src_name}/Absorption/${image_name}.tar --directory .

## run exportfits to convert to FITS file
/casa-6.5.0-15-py3.8/bin/casa --nologfile --log2term --nogui -c exportfits.py -i ${image_name}.image -f ${image_name}.image.fits

## now, run source finding
aegean ${image_name}.image.fits --out ${image_name}_source_catalog.csv

## parse log
python parse_log.py -n ${image_name}_source_catalog.csv -t ${int_threshold}

## phase-center-list.txt will be moved back to home directory after job completion, but clean up other files
rm ${image_name}*
rm anaconda.sh 
rm -rf anaconda3