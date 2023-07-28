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
import astropy.units as u
import numpy as np

## parse user inputs
parser = argparse.ArgumentParser()
parser.add_argument('-n', '--name', help = '<required> name of input file', required = True)
parser.add_argument('-t', '--threshold', help = '<required> peak intensity threshold above which to keep sources', required = True)
args, unknown = parser.parse_known_args()

## unpack user arguments
file_name = args.name
threshold = float(args.threshold)


def main():
	## read in file
	data_array = np.loadtxt('%s' % file_name, dtype = 'str', skiprows=2)

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
if __name__=='__main__':
	main()
	exit()
