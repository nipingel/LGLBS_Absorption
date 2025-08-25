"""
5/04/2023
Image LGLBS source around user provided phase center
User inputs:
-v --vis_path - <required> path of ms file
-r --ra - <required> ra phase center in form e.g.: 00h40m13.8 
-d --dec - <required> dec phase center in form e.g.: +40d50m04.73
-o --output_name - <required> name of output file
-t --threshold - <required> global rms threshold to stop cleaning
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
parser.add_argument('-r', '--ra', help = '<required> ra phase center in form e.g.: 00h40m13.8', required = True)
parser.add_argument('-d', '--dec', help = '<required> ra phase center in form e.g.: +40d50m04.73', required = True)
parser.add_argument('-o', '--output_name', help = '<required> name of output file', required = True)
parser.add_argument('-t', '--threshold', help = '<required> global rms threshold to stop cleaning', required = True)
args, unknown = parser.parse_known_args()

## unpack user arguments
vis_path = args.vis_path
ra_phase_center = args.ra
dec_phase_center = args.dec
output_name = args.output_name
major_threshold = args.threshold
def main():
	## define tclean variables below
	## image output properties
	im_size = 5120
	cell_size = '0.6arcsec'
	restore_beam = 'common'
	## automasking parameters ##
	use_mask = 'pb'
	sidelobe_threshold = 2.5
	noise_threshold = 3.5
	min_beam_frac = 0.3
	lownoisethreshold = 1.5
	negativethreshold = 0.0
	grow_iters=75
	dogrowprune=False
	verbose = True
	## deconvolution parameters
	deconvolver_mode = 'hogbom'
	tot_niter = 200000
	min_threshold = '%smJy' % major_threshold
	restart_parameter = False
	## tclean dictionary
	tclean_params={
		'vis':vis_path,
		'imagename':output_name,
		'phasecenter':'J2000 %s %s' % (ra_phase_center, dec_phase_center),
		'restfreq':'1.42040571183GHz',
		'selectdata': True,
		'datacolumn': 'data',
		'specmode':'mfs',
		'imsize':im_size,
		'cell':cell_size,
		'restoringbeam': restore_beam, 
		'pblimit':0.1, 
		'weighting':'briggs', 
		'robust':0.5, 
		'gridder':'wproject', 
		'pbcor':True, 
		'niter':tot_niter, 
		'deconvolver':deconvolver_mode, 
		'cyclefactor':0.8, 
		'minpsffraction':0.05, 
		'maxpsffraction':0.8, 
		'threshold':min_threshold, 
		'usemask':use_mask, 
		'pbmask':0.5, 
		'sidelobethreshold':sidelobe_threshold, 
		'noisethreshold':noise_threshold, 
		'minbeamfrac':min_beam_frac, 
		'lownoisethreshold':lownoisethreshold, 
		'negativethreshold':0.0, 
		'growiterations':grow_iters, 
		'dogrowprune':False, 
		'verbose':True, 
		'restart':restart_parameter
		}

	## run tclean
	tclean(**tclean_params)
	
	## run exportfits
	exportfits_params={
		'imagename':'%s.image' % output_name,
		'fitsimage':'%s.image.fits' % output_name,
		'dropstokes':True,
		'dropdeg':True
		}
	exportfits(**exportfits_params)
if __name__=='__main__':
	main()
	exit()
