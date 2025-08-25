"""
07/26/2023
small script to run exportfits on casa image file
User inputs:
-i --image - <required> name of image file
-f --fits - <required> output name of fits file 
__author__="Nickolas Pingel"
__version__="1.0"
__email__="nmpingel@wisc.edu"
__status__="Production"
"""
# imports
import argparse

## parse user inputs
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--image', help = '<required> name of image file', required = True)
parser.add_argument('-f', '--fits', help = '<required> name of output fits file', required = True)
args, unknown = parser.parse_known_args()

image_name = args.image
fits_name = args.fits
def main():
	## exportfits dictionary
	exportfits_params={
		'imagename': image_name,
		'fits_image': fits_name,
		'history': False,
		'dropdeg': True, 
		'dropstokes': True,
		}

	## run exportfits
	exportfits(**tclean_params)
if __name__=='__main__':
	main()
	exit()