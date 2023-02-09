#!/bin/bash
#
# split_channels.sh
# execution script to split out range of channels from a staged and calibrated LGLBS measurement set 

## TODO include ms_name and untar_name variables as argument into this executable script
## copy measurement set to working directory
cp /projects/vla-processing/measurement_sets/$tarball_name ./

## untar measurement set
tarball_name="split_chans_"$1"_to_"$2".tar"

## untar 
tar -xvf $tarball_name
untar_str="20A-346.sb41042474.eb41074466.59579.84285818287.speclines.ms.split_chan"

mv $untar_str"*" /projects/vla-processing/measurement_sets

## clean up
rm -rf $untar_str"*"
rm -rf $tarball_name

