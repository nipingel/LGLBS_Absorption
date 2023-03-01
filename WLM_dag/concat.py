"""
11/14/2i022
concatenate input measurement sets into a single file
User inputs:
-o --output_name - name of output concatentated measurement set
__author__="Nickolas Pingel"
__version__="1.0"
__email__="nmpingel@wisc.edu"
__status__="Production"
"""
# imports
import argparse
import glob as glob

## parse user inputs
parser = argparse.ArgumentParser()
#parser.add_argument('-l', '--ms_list', nargs = '+', help='<required> list of measurement sets to concatentate', required = True)
parser.add_argument('-o', '--output', help='<required> name of output concatentated measurement set', required = True)
args, unknown = parser.parse_known_args()

## parse measurement set list & output
output_vis = args.output

## get list of input measurement sets
ms_list = glob.glob('./*.transformed')

def main():
	concat_params = {
		'vis': ms_list,
		'concatvis':output_vis,
		'freqtol':'0.4kHz',
		'dirtol': '0.1deg'}
	concat(**concat_params)
if __name__=='__main__':
	main()
	exit()
