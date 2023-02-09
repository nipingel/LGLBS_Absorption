"""
10/16/2022
Image continuum source around user provided phase center
User inputs:
-p --path - <required> path of ms file
-o --output_name - <required> name of output file
__author__="Nickolas Pingel"
__version__="1.0"
__email__="nmpingel@wisc.edu"
__status__="Production"
"""
# imports
import argparse

## parse user inputs
parser = argparse.ArgumentParser()
parser.add_argument('-v', '--vis_path', help = '<required> name of measurement set', required = True)
#parser.add_argument('-p', '--phase_center', help = '<required> phase center in form e.g.: J2000 00h40m13.8 +40d50m04.73', required = True)
parser.add_argument('-o', '--output_name', help = '<required> name of output file', required = True)

args, unknown = parser.parse_known_args()

## unpack user arguments
vis_path = args.vis_path
output_name = args.output_name

def main():
	## set default values for select tclean parameter
	phasecenter='J2000 00h42m18.757s +41d29m27.83'
	ref_freq='1.42040571183GHz'
	rest_freq='1.42040571183GHz'
	uvdist='>1.5Klambda'
	field_id = '49'
	tot_iter = 1000
	threshold_value = '0.1mJy'
	tclean(vis=vis_path, imagename=output_name, intent='OBSERVE_TARGET*', reffreq=ref_freq, restfreq=rest_freq, 
		phasecenter=phasecenter,imsize=50,uvrange=uvdist,weighting='natural',gridder='standard',
		pbcor=True, threshold = threshold_value, cell='1arcsec', specmode = 'mfs', usemask='pb', niter=tot_iter)
if __name__=='__main__':
	main()
	exit()
	
