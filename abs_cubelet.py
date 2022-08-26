"""
8/3/2022
Imaging cubelet around user provided phase center
User inputs:
-p --path - <required> path of ms file
-s --start_chan - <required> starting channel for splitting
-e --end_chan - <required> ending channel for splliting
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
#parser.add_argument('-s', '--start_chan', help = '<required> start channel to split from main measurement set', required = True, type = int)
parser.add_argument('-o', '--output_name', help = '<required> name of output file', required = True)

args, unknown = parser.parse_known_args()

## unpack user arguments
vis_path = args.vis_path
#phase_center = args.phase_center
#start_chan = args.start_chan
output_name = args.output_name

def main():
	## set default values for select tclean parameters
	phasecenter='J2000 00h40m13.8 +40d50m04.73'
	ref_freq='1.42040571183GHz'
	rest_freq='1.42040571183GHz'
	uvdist='>1.5Klambda'
	tclean(vis=vis_path, imagename=output_name, reffreq=ref_freq, restfreq=rest_freq, 
		phasecenter=phasecenter,imsize=50,uvrange=uvdist,weighting='natural',gridder='standard',
		pbcor=True,niter=1000,cell='1arcsec', specmode = 'cube')
if __name__=='__main__':
	main()
	exit()
	
