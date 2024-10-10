"""
07/26/2023
parse output log from aegean source finding software
User inputs:
-n --name - <required> name of input file
-t --threshold - <required> peak intensity threshold above which to keep sources
__author__="Nickolas Pingel"
__version__="1.0"
__email__="nmpingel@wisc.edu"
__status__="Production"
"""
# imports
import argparse
from astropy.coordinates import SkyCoord
import pandas as pd
import astropy.units as u
import numpy as np
import csv
import warnings
import sys
## avoid pandas warning
warnings.filterwarnings("ignore", category=FutureWarning) 

## parse user inputs
parser = argparse.ArgumentParser()
parser.add_argument('-n', '--name', help = '<required> name of input file', required = True)
parser.add_argument('-t', '--threshold', help = '<required> peak intensity threshold above which to keep sources', required = True)
args, unknown = parser.parse_known_args()

## unpack user arguments
file_name = args.name
threshold = float(args.threshold)

## function to produce a plain text file with just RA/Dec phase centers 
## used in the cubelet generation node
def write_phase_centers(f_name):
	## read in file
	data_array = np.loadtxt('%s' % f_name, dtype = 'str', skiprows=2)

	## select the peak intensity column
	peak_column = np.array(data_array[:, 9], dtype = 'float')

	## select good rows
	good_rows = np.where(peak_column > threshold)[0]

	## exit if no suitable sources
	if len(good_rows) == 0:
		exit()

	## otherwise, continue to write file by selecting ra/dec columns from input table
	ra_column = data_array[:, 3]
	dec_column = data_array[:, 4]

	## construct SkyCoord object
	c = SkyCoord(ra_column[good_rows], dec_column[good_rows], frame = 'fk5', unit = (u.hourangle, u.deg))

	## open file for writing
	f = open("phase-center-list.txt", "w")
	## loop through good rows and construct the string to write out to file
	for i in range(0, len(ra_column[good_rows])):
		ra_str = '%sh%sm%.2fs' % (int(c[i].ra.hms[0]), int(c[i].ra.hms[1]), c[i].ra.hms[2])
		dec_str = '%sd%sm%.2fs ' % (int(c[i].dec.dms[0]), int(np.abs(c[i].dec.dms[1])), np.abs(c[i].dec.dms[2]))
		f.write('%s, %s\n'  % (ra_str, dec_str))
	f.close()

## function to write out csv columns for sources that meet flux limit
def write_out_csv(f_name):
	## read in file into pandas
	df = pd.read_csv(f_name, delim_whitespace=True, skiprows=2)

	## select the peak intensity column
	peak_column = np.array(df['Peak'].values[2:], dtype = float)

	## select good rows
	good_rows = np.where(peak_column > threshold)

	## exit if no suitable sources
	if len(good_rows) == 0:
		sys.exit("No identified sources exceed given peak flux density threshold")

	## parse df, accounting for the first two bad rows
	parsed_df = df.iloc[good_rows[0]+2] 

	## rename columns
	parsed_df = parsed_df.rename(columns={"Peak":"Peak Flux Density [Jy/beam]", "err.2":"Peak err [Jy/beam]", "S_int":"S_int [Jy]", "err.3":"S_int err [Jy]"})
	parsed_df.to_csv(f_name.replace('.csv', '_parsed.csv'), 
		columns = ['RA', 'DEC', 'Peak Flux Density [Jy/beam]', 'Peak err [Jy/beam]', 'S_int [Jy]', 'S_int err [Jy]'], index = False)
	return

def main():

	## write out file with just phase centers to be used in next node
	write_phase_centers(file_name)

	## write out all columns of sources that meet flux threshold
	write_out_csv(file_name)

if __name__=='__main__':
	main()
	exit()