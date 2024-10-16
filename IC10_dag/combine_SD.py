"""
8/1/2023
extract 30' subimage around provided ra/dec and feather with single dish data
User inputs:
-p --int_path - <required> path of ms file
-s --sd_name - <required> name of sd cube
-r --ra_phase_center - <required> RA phase center
-d --dec_phase_center -<required> Dec phase center
__author__="Nickolas Pingel"
__version__="1.0"
__email__="nmpingel@wisc.edu"
__status__="Production"
"""
# imports
import argparse

## parse user inputs
parser = argparse.ArgumentParser()
parser.add_argument('-p', '--int_path', help = '<required> path to vla cube', required = True)
parser.add_argument('-s', '--sd_name', help = '<required> name single dish cubemeasurement set', required = True)
parser.add_argument('-r', '--ra_phase_center', help = '<required> RA phase center in form e.g.: 00h40m13.8s', required = True)
parser.add_argument('-d', '--dec_phase_center', help = '<required> Dec phase center in form e.g.: +40d50m04.73', required = True)
args, unknown = parser.parse_known_args()

## unpack user arguments
int_path = args.int_path
sd_name = args.sd_name
ra_phase_center = args.ra_phase_center
dec_phase_center = args.dec_phase_center

def main():
	## read in sd_cube
	importfits_params = {
		'fitsimage': int_path,
		'imagename': int_path.replace('.fits', '.im')
	}
	#importfits(**importfits_params)

	## extract 30' subregion around source
	imsubimage_params = {
		'imagename': int_path,
		'outfile': 'VLA_ABCD_GBT_%s_%s_30arcmin.imsub' % (ra_phase_center, dec_phase_center),
		'region': 'circle[[%s, %s], 15arcmin]' % (ra_phase_center, dec_phase_center),
		'dropdeg': True,
	}
	imsubimage(**imsubimage_params)

	## regrid sd to vla cube
	regrid_params = {
		'imagename': sd_name.replace('.fits', '.im'),
		'template': 'VLA_ABCD_%s_%s_30arcmin.imsub' % (ra_phase_center, dec_phase_center), 
		'output': 'GBT_%s_%s_30arcmin.regrid' % (ra_phase_center, dec_phase_center),
		'asvelocity': True,
		'axes': [0, 1, 2], 
	}
	#imregrid(**regrid_params)

	## run imhead to put beam info in header
	imhead_params = {
		'imagename': 'VLA_ABCD_%s_%s_30arcmin.imsub' % (ra_phase_center, dec_phase_center),
		'mode': 'put',
		'hdkey': 'bmaj',
		'hdvalue': '7arcsec'
	}
	#imhead(**imhead_params)
	imhead_params = {
		'imagename': 'VLA_ABCD_%s_%s_30arcmin.imsub' % (ra_phase_center, dec_phase_center),
		'mode': 'put',
		'hdkey': 'bmin',
		'hdvalue': '5.2arcsec'
	}
	#imhead(**imhead_params)
	imhead_params = {
		'imagename': 'VLA_ABCD_%s_%s_30arcmin.imsub' % (ra_phase_center, dec_phase_center),
		'mode': 'put',
		'hdkey': 'bpa',
		'hdvalue': '-2.96deg'
	}
	#imhead(**imhead_params)
	imsmooth_params = {
		'imagename': 'VLA_ABCD_GBT_%s_%s_30arcmin.imsub' % (ra_phase_center, dec_phase_center),
		'outfile': 'VLA_ABCD_GBT_%s_%s_30arcmin.imsub.imsmooth' % (ra_phase_center, dec_phase_center),
		'kernel': 'commonbeam'
	}
	imsmooth(**imsmooth_params)
	## combine sd and vla cube
	feather_params = {
		'imagename':'VLA_ABCD_GBT_%s_%s_30arcmin.imsub' % (ra_phase_center, dec_phase_center),
		'highres': 'VLA_ABCD_%s_%s_30arcmin.imsub.imsmooth' % (ra_phase_center, dec_phase_center), 
		'lowres': 'GBT_%s_%s_30arcmin.regrid' % (ra_phase_center, dec_phase_center),
		'sdfactor': 1.0
	}
	#feather(**feather_params)

	## export the combined subcube
	exportfits_params = {
		'imagename': 'VLA_ABCD_GBT_%s_%s_30arcmin.imsub.imsmooth' % (ra_phase_center, dec_phase_center), 
		'fitsimage': 'VLA_ABCD_GBT_%s_%s_30arcmin.fits' % (ra_phase_center, dec_phase_center),
		'velocity': True, 
		'dropdeg': True,
		'dropstokes': True,
		'history': False
		}
	exportfits(**exportfits_params)
if __name__=='__main__':
	main()
	exit()
