#!/bin/bash
#
# sub_cube_abs.sh
# execution script for imaging of sub cube around absorption source

## untar CASA
tar -xvf casa.tar

## untar measurement set
ms_name="M31_B_20A-346.sb41042474.eb41074466.59579.84285818287.speclines.ms.split.tar"
untar_name="20A-346.sb41042474.eb41074466.59579.84285818287.speclines.ms.split"

## copy measurement set to working directory
cp /staging/groups/astro_stanimirovic/$ms_name ./

## untar 
tar -xvf $ms_name

# make casa call to imaging script
casa-6.4.4-31-py3.8/bin/casa -c sub_cube_abs.py

tar -cvf output.tar output.*

mv output.tar /staging/groups/astro_stanimirovic

## clean up
rm -rf $ms_name
rm -rf $untar_name

