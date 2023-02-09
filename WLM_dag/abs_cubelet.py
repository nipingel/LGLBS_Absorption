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
parser.add_argument('-r', '--ra_phase_center', help = '<required> RA phase center in form e.g.: 00h40m13.8s', required = True)
parser.add_argument('-d', '--dec_phase_center', help = '<required> Dec phase center in form e.g.: +40d50m04.73', required = True)
parser.add_argument('-f', '--field_id', help = '<required> field ID', required = True)
parser.add_argument('-s', '--start_velocity', help = '<required> starting TOPO velocity (e.g., -500', required = True)
parser.add_argument('-n', '--n_chan', help = '<required> number of output channels', type = int, required = True)
parser.add_argument('-o', '--output_name', help = '<required> name of output file', required = True)
args, unknown = parser.parse_known_args()

## unpack user arguments
vis_path = args.vis_path
start_vel = args.start_velocity+'km/s'
n_chan = args.n_chan
f_id = args.field_id
output_name = args.output_name
ra_phase_center = args.ra_phase_center
dec_phase_center = args.dec_phase_center

def main():
	## set default values for select tclean parameters
	phasecenter='J2000 %s %s' % (ra_phase_center, dec_phase_center)
	ref_freq='1.42040571183GHz'
	rest_freq='1.42040571183GHz'
	uvdist='>1.5Klambda'
	field_id = f_id
	tot_iter = 750000
	threshold_value = '10mJy'
	tclean(vis=vis_path, imagename=output_name, field = field_id, reffreq=ref_freq, restfreq=rest_freq, 
		phasecenter=phasecenter,imsize=50,uvrange=uvdist,weighting='natural',gridder='standard',
		pbcor=True, threshold = threshold_value, cell='1arcsec', specmode = 'cube', start = start_vel, nchan = n_chan, usemask='pb', niter=tot_iter)
if __name__=='__main__':
	main()
	exit()
	
