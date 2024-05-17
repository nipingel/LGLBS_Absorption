"""
8/1/2023
run analysis on absorption cubelet to:
	extract emission spectrum
	extract absorption spectrum
	compute noise
	define spectral range of absorption features
	compute spin temperature
User inputs:
-n --cubelet_name - <required> name of cubelet file
-c --combined_name - <required> name of combined cube
__author__="Nickolas Pingel"
__version__="1.0"
__email__="nmpingel@wisc.edu"
__status__="Production"
"""
## imports
import numpy as np
import os
import argparse
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import glob as glob
import pickle
from spectral_cube import SpectralCube
from scipy.interpolate import interp1d
import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.gridspec as gridspec
matplotlib.rc('font', family='serif')
matplotlib.rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 10})
matplotlib.rcParams['figure.dpi']= 250
matplotlib.rc("savefig", dpi=250)

# These are the "Tableau 20" colors as RGB.
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
			 (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
			 (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
			 (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
			 (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib
# accepts.
for i in range(len(tableau20)):
	r, g, b = tableau20[i]
	tableau20[i] = (r / 255., g / 255., b / 255.)

## parse user inputs
parser = argparse.ArgumentParser()
parser.add_argument('-n', '--cubelet_name', help = '<required> name of cubelet file', required = True)
parser.add_argument('-c', '--combined_name', help = '<required> name of combined cube', required = True)
args, unknown = parser.parse_known_args()

## function to define spectral mask to avoid velocity range with emission
def define_mask(low_vel_list, high_vel_list, vel_axis):
	mask_axis = np.zeros_like(vel_axis, dtype='bool')
	for low_vel, high_vel in zip(low_vel_list, high_vel_list):
		inds = (((vel_axis) >= low_vel) *
		((vel_axis) <= high_vel))
		mask_axis[inds] = True
	return mask_axis

## function to return a spectral axis from an input cube
def return_spectral_axis(cube_name):
	cubelet = SpectralCube.read(cube_name)
	abs_vel_axis = cubelet.spectral_axis.value/1000.0 ## km/s
	return abs_vel_axis

## function to return beam info
def get_beam_info(cube_name):
	hdu = fits.open(cube_name)
	bmaj = hdu[0].header['BMAJ']
	bmin = hdu[0].header['BMIN']
	bpa = hdu[0].header['BPA']
	hdu.close()
	return bmaj, bmin, bpa

## function to run source finding on integrated intensity image
def source_find_mom0(cubelet_name, bmaj, bmin, bpa):
	## make integrated intensity image
	hdu = fits.open(cubelet_name)
	data = hdu[0].data
	mom0_image = np.sum(data, axis = 0)
	fits.writeto(cubelet_name.replace('.fits', '.mom0.fits'), mom0_image, header = hdu[0].header)
	hdu.close()

	## run source-finding
	os.system('aegean %s --beam %s %s %s --out %s_aegean_catalog.csv' % (cubelet_name.replace('.fits', '.mom0.fits'), bmaj, bmin, bpa, cubelet_name[:-5]))
	return mom0_image

## function to read in csv file from source finding on cubelet
def parse_csv(f):
	parse_array = np.loadtxt(f, dtype = 'str')
	## add additional axis if only single row
	if len(parse_array.shape) == 1:
		parse_array = np.expand_dims(parse_array, axis=0)
	## extract the useful columns
	ra_str_column = np.array(parse_array[:, 3], dtype = 'str')
	dec_str_column = np.array(parse_array[:, 4], dtype = 'str')
	ra_column = np.array(parse_array[:, 5], dtype = 'float')
	dec_column = np.array(parse_array[:, 7], dtype = 'float')
	a_column = np.array(parse_array[:, 13], dtype = 'float')
	b_column = np.array(parse_array[:, 15], dtype = 'float')
	pa_column = np.array(parse_array[:, 17], dtype = 'float')
	return ra_str_column, dec_str_column, ra_column, dec_column, a_column, b_column, pa_column

## function to convert ra/dec strings to hms, dms string for saving output
def convert_hms_dms(ra_str, dec_str):
	c = SkyCoord(ra_str, dec_str, frame = 'fk5', unit = (u.hourangle, u.deg))
	new_ra_str = '%dh%dm%.2fs' % (c.ra.hms[0], c.ra.hms[1], c.ra.hms[2])
	new_dec_str = '%dd%dm%.2fs' % (c.dec.dms[0], np.abs(c.dec.dms[1]), np.abs(c.dec.dms[2]))
	return new_ra_str, new_dec_str

## function to extract pixels within source ellipse
def extract_mean_pixel_spectrum(mask, cubelet_name, vel_axis, bmaj, bmin, bpa, ra, dec, a, b, pa):
	hdu = fits.open(cubelet_name)
	cubelet = hdu[0].data
	cubelet_wcs = WCS(hdu[0].header)
	pix_size = (hdu[0].header['CDELT1'])*3600
	pix_area = (hdu[0].header['CDELT1'])**2*3600**2
	hdu.close()

	## define radial pixel grid and place origin at centeral ra, dec pixels
	## get pixel indices of ellipse center
	central_c = SkyCoord(ra, dec, frame = 'fk5', unit = u.deg)
	central_pix = central_c.to_pixel(cubelet_wcs)
	## convert a/b to pixels
	a_pix = a/pix_size
	b_pix = b/pix_size


	## define 2D radial grid centered on source
	y = np.arange(-cubelet.shape[1]/2, cubelet.shape[1]/2)
	x = np.arange(-cubelet.shape[2]/2, cubelet.shape[2]/2)
	yy, xx = np.meshgrid(y, x, indexing='ij')
	y0 = central_pix[1]-cubelet.shape[1]/2
	x0 = central_pix[0]-cubelet.shape[2]/2
	pa_rad = np.pi/2+np.deg2rad(pa)
	z = (np.sin(pa_rad)**2/a_pix**2+np.cos(pa_rad)**2/b_pix**2)*(yy-y0)**2 + \
		2*np.cos(pa_rad)*np.sin(pa_rad)*(1/a_pix**2-1/b_pix**2)*(yy-y0)*(xx-x0) + \
		((xx-x0))**2*(np.sin(pa_rad)**2/b_pix**2+np.cos(pa_rad)**2/a_pix**2)

	## create a boolean array for pixels that fall within source ellipse
	source_image_mask = np.zeros([cubelet.shape[1], cubelet.shape[2]])
	source_image_mask[np.where(z<1)] = 1.0
	## extend spatial pixel mask along entire spectral axis
	source_cube_mask = np.repeat(source_image_mask[np.newaxis, :, :], len(vel_axis), axis=0)
	## compute mean brightness temperature along each sightline
	mean_tb_image = np.nanmean(cubelet, axis = 0)
	## compute weighted sum
	spectrum = np.nansum(mean_tb_image**2*cubelet, axis = (1,2), where = source_cube_mask.astype(bool))
	## compute continuum level
	c = np.mean(spectrum[mask])
	return spectrum/c

## function to generate the emission spectrum
def extract_emission_spectrum(em_cube_name, ra, dec):
	hdu = fits.open(em_cube_name)
	em_cube = hdu[0].data
	y_len, x_len = em_cube.shape[1], em_cube.shape[2]
	beam_maj = hdu[0].header['BMAJ']*3600
	beam_min = hdu[0].header['BMIN']*3600
	pix_size = np.abs(hdu[0].header['CDELT1'])*3600
	## convert to units of brightness temperature K
	em_cube_Tb = 1*em_cube #*6.07e5/(beam_maj*beam_min) ## Jy/beam->Tb
	em_wcs = WCS(hdu[0].header)
	central_c = SkyCoord(ra, dec, frame = 'fk5', unit = u.deg)
	c_pix = central_c.to_pixel(em_wcs)
	hdu.close()

	## define 2D radial grid centered on source
	y = np.arange(-y_len/2, y_len/2)
	x = np.arange(-x_len/2, x_len/2)
	yy, xx = np.meshgrid(y, x, indexing='ij')
	y0 = c_pix[1]-y_len/2
	x0 = c_pix[0]-y_len/2
	r_pixs = np.sqrt((yy-y0)**2 + (xx-x0)**2)
	## make a mask based on radius from source
	radial_mask = np.zeros([em_cube_Tb.shape[1], em_cube_Tb.shape[2]])
	## center an annulus on souce 2x beamwidths wide
	radial_mask[r_pixs < 2*beam_maj/pix_size] = 1
	## mask inner pixels corresponding to single synthesized beam
	radial_mask[r_pixs < beam_maj/pix_size] = 0
	## extend spatial pixel mask along entire spectral axis
	radial_cube_mask = np.repeat(radial_mask[np.newaxis, :, :], em_cube_Tb.shape[0], axis=0)
	s = np.nanmean(em_cube_Tb, axis = (1,2), where = radial_cube_mask.astype(bool))
	s_err = np.nanstd(em_cube_Tb, axis = (1,2), where = radial_cube_mask.astype(bool))
	return s, s_err

## function to compute mean emission spectrum
def compute_mean_em_spectrum(em_cube_name):
	hdu = fits.open(em_cube_name)
	em_cube = hdu[0].data
	beam_maj = hdu[0].header['BMAJ']*3600
	beam_min = hdu[0].header['BMIN']*3600
	hdu.close()

	## get emission spectral axis
	em_cube_sc = SpectralCube.read(em_cube_name)
	em_vel_axis = em_cube_sc.spectral_axis.value/1000.0 ## km/s

	## compute mean spectrum
	mean_em_spectrum = np.nanmean(em_cube, axis = (1, 2)) #*6.07e5/(beam_maj*beam_min) ## Jy/beam -> Tb
	return mean_em_spectrum, em_vel_axis

## function to find consecutive indices to identify absorption features
def consecutive(data, stepsize=1):
	return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)

## function to generate abs/emission plot
def generate_plot(vel_axis, abs_spectrum, em_spectrum, three_sig_env, one_sig_env, rms, detection_sort, ra, dec):
	fig = plt.figure(constrained_layout=True)
	gs = gridspec.GridSpec(ncols=1, nrows=2, figure=fig)
	ax1 = fig.add_subplot(gs[0])
	ax2 = fig.add_subplot(gs[1])

	## ABS SPECTRUM ##
	## set up tick mark parameters for absorption
	majorYLocator = MultipleLocator(0.1)
	majorYFormatter = FormatStrFormatter('%.2f')
	minorYLocator = MultipleLocator(0.1/5)
	majorXLocator = MultipleLocator(25)
	majorXFormatter = FormatStrFormatter('%d')
	minorXLocator = MultipleLocator(25/5)
	## set tick mark parameters
	#ax1.yaxis.set_major_locator(majorYLocator)
	#ax1.yaxis.set_major_formatter(majorYFormatter)
	#ax1.yaxis.set_minor_locator(minorYLocator)
	#ax1.xaxis.set_major_locator(majorXLocator)
	#ax1.xaxis.set_major_formatter(majorXFormatter)
	#ax1.xaxis.set_minor_locator(minorXLocator)
	#ax1.tick_params(axis='both', which='major')

	## do plotting
	ax1.plot(vel_axis, abs_spectrum)
	ax1.plot(vel_axis, three_sig_env, color = tableau20[2], linestyle = '--')
	ax1.fill_between(vel_axis, 1+one_sig_env, 1-one_sig_env, color = 'grey', alpha = 0.5)
	ax1.axhline(1.0, color = tableau20[6], linewidth = 0.75)
	## add markings to show spand of features above 3-sigma
	for l in detection_sort:
		if len(l) > 2:
			ax1.axvspan(vel_axis[l[0]], vel_axis[l[-1]], ymin = 0, ymax = 1.2, color = tableau20[8], alpha = 0.35)

	ax1.set_xlim(-120, 50)
	ax1.set_ylim()
	ax1.set_xlabel(r'Velocity [km s$^{-1}$]')
	ax1.set_ylabel(r'$e^{-\tau}$')
	ax1.set_title(r"$\alpha_{\rm J2000}$=%s, $\delta_{\rm J2000}$=%s" "\n" "Absorption" % (ra, dec))
	ax1.grid()

	## EMISSION SPECTRUM ##
	## set up tick mark parameters for absorption
	majorYLocator = MultipleLocator(10)
	majorYFormatter = FormatStrFormatter('%d')
	minorYLocator = MultipleLocator(10/5)
	majorXLocator = MultipleLocator(25)
	majorXFormatter = FormatStrFormatter('%d')
	minorXLocator = MultipleLocator(25/5)
	## set tick mark parameters
	#ax2.yaxis.set_major_locator(majorYLocator)
	#ax2.yaxis.set_major_formatter(majorYFormatter)
	#ax2.yaxis.set_minor_locator(minorYLocator)
	#ax2.xaxis.set_major_locator(majorXLocator)
	#ax2.xaxis.set_major_formatter(majorXFormatter)
	#ax2.xaxis.set_minor_locator(minorXLocator)

	## do plotting
	ax2.plot(vel_axis, em_spectrum)
	ax2.fill_between(vel_axis, em_spectrum-rms, em_spectrum+rms, color = 'grey', alpha = 0.5)
	ax2.axhline(0.0, color = tableau20[6], linewidth = 0.75)
	## add markings to show spand of features above 3-sigma
	for l in detection_sort:
		if len(l) > 2:
			ax2.axvspan(vel_axis[l[0]], vel_axis[l[-1]], ymin = 0, ymax = 1.2, color = tableau20[8], alpha = 0.35)
	ax2.set_xlim(-120, 50)
	ax2.set_xlabel(r'Velocity [km s$^{-1}$]')
	ax2.set_ylabel(r'$T_{\rm b}$ [K]')
	ax2.set_title('Emission')
	ax2.grid()
	plt.savefig('%s_%s_abs_em_spectrum.pdf' % (ra, dec))
	plt.show()
	plt.close()

## function to compute spin temperature
def compute_mean_spin_temperature(em_profile, em_profile_err, abs_profile, abs_profile_err, dv):
	tb_integral = np.sum(em_profile)*dv
	tb_integral_err = np.sqrt(np.sum(em_profile_err**2))*dv
	abs_integral = np.sum(1-np.exp(-1*abs_profile))*dv
	abs_integral_err = np.sqrt(np.sum(np.exp(-2*abs_profile)*abs_profile_err**2))*dv
	return tb_integral/abs_integral, tb_integral/abs_integral*np.sqrt((tb_integral_err/tb_integral)**2+(abs_integral_err/abs_integral)**2)

## function to set ra/dec str for the output file names
def set_ra_dec_output_name(cubelet_name):
	ra = cubelet_name.split('_')[0]
	dec = cubelet_name.split('_')[1].split('s')[0]
	return ra, dec

## function to pickle results
def pickle_results(vel_axis, abs_spectrum, em_spectrum, three_sig_env, one_sig_env, rms, detection_sort, ra, dec, spin_temp, spin_temp_err):
	pickle.dump([vel_axis, abs_spectrum, em_spectrum, three_sig_env, one_sig_env, rms, detection_sort, ra, dec, spin_temp, spin_temp_err], \
		open('%s_%s_analysis_variables.pickle' % (ra, dec), "wb" ))





## NEW IO CODE YAYAYAYAYAYY #######
def get_source_apperture(mask, cubelet_name, vel_axis, bmaj, bmin, bpa, ra, dec, a, b, pa):
	hdu = fits.open(cubelet_name, mode='readonly')
	cubelet = hdu[0].data
	cubelet_wcs = WCS(hdu[0].header)
	pix_size = (hdu[0].header['CDELT1'])*3600
	pix_area = (hdu[0].header['CDELT1'])**2*3600**2
	hdu.close()

	## define radial pixel grid and place origin at centeral ra, dec pixels
	## get pixel indices of ellipse center
	central_c = SkyCoord(ra, dec, frame = 'fk5', unit = u.deg)
	central_pix = central_c.to_pixel(cubelet_wcs)
	## convert a/b to pixels
	a_pix = a/pix_size
	b_pix = b/pix_size
	## define 2D radial grid centered on source
	y = np.arange(-cubelet.shape[1]/2, cubelet.shape[1]/2)
	x = np.arange(-cubelet.shape[2]/2, cubelet.shape[2]/2)
	yy, xx = np.meshgrid(y, x, indexing='ij')
	y0 = central_pix[1]-cubelet.shape[1]/2
	x0 = central_pix[0]-cubelet.shape[2]/2
	pa_rad = np.pi/2+np.deg2rad(pa)
	z = (np.sin(pa_rad)**2/a_pix**2+np.cos(pa_rad)**2/b_pix**2)*(yy-y0)**2 + \
		2*np.cos(pa_rad)*np.sin(pa_rad)*(1/a_pix**2-1/b_pix**2)*(yy-y0)*(xx-x0) + \
		((xx-x0))**2*(np.sin(pa_rad)**2/b_pix**2+np.cos(pa_rad)**2/a_pix**2)

	## create a boolean array for pixels that fall within source ellipse
	source_image_mask = np.zeros([cubelet.shape[1], cubelet.shape[2]])
	source_image_mask[np.where(z<1)] = 1.0
	
	return source_image_mask

def handle_overlaps(mask, cubelet_name, vel_axis, bmaj, bmin, bpa, ra, dec, a, b, pa):
        
    hdu = fits.open(cubelet_name)
    cubelet = hdu[0].data
    ellipse_masks = []
    hdu.close()

    for ra, dec, a, b, pa in zip(ra, dec, a, b, pa):
        ellipse_masks.append(get_source_apperture(mask, cubelet_name, vel_axis, bmaj, bmin, bpa, ra, dec, a, b, pa))

    union_mask = calculate_union_mask(ellipse_masks)
    try:
        if union_mask is None:
            print('No overlap found between the sources')
            return None
        else:
            source_cube_mask = np.repeat(union_mask[np.newaxis, :, :], len(vel_axis), axis=0)
            ## compute mean brightness temperature along each sightline
            mean_tb_image = np.nanmean(cubelet, axis = 0)
            ## compute weighted sum
            spectrum = np.nansum(mean_tb_image**2*cubelet, axis = (1,2), where = source_cube_mask.astype(bool))
            ## compute continuum level
            c = np.mean(spectrum[mask])

            return spectrum/c
    except:
        print('No overlap found between the sources')
        return None


def calculate_union_mask(ellipse_masks):
    # Check if there is more than one mask
    if len(ellipse_masks) <= 1:
        return ellipse_masks[0] if len(ellipse_masks) == 1 else None

    # Initialize the union mask with the first mask
    union_mask = np.copy(ellipse_masks[0])

    # Flag to track if there is an overlap
    overlap_flag = False

    # Loop through the remaining masks
    for mask in ellipse_masks[1:]:
        # Check if there is any overlap between the current mask and the union mask
        if np.any(np.logical_and(union_mask, mask)):
            # Update the union mask by taking the logical OR of the current mask with the union mask
            union_mask = np.logical_or(union_mask, mask)
            overlap_flag = True  # Set the flag to True if overlap is found

    # Return the union mask only if there was an overlap, otherwise return None
    return union_mask if overlap_flag else None

## NEW IO CODE YAYAYAYAYAYY #######
## unpack user arguments
cubelet_name = args.cubelet_name
combined_name = args.combined_name

## useful constants
tsys = 25.0
apeff = 0.35
rms_K = 7.5

def main():
	## useful constants
	low_vel_list = [-475, -230] 
	high_vel_list = [-445, -175]

	## get velocity axis for absorption cubelet and construct mask
	abs_vel_axis = return_spectral_axis(cubelet_name)

	## create mask
	mask = define_mask(low_vel_list, high_vel_list, abs_vel_axis)

	## get beam information
	bmaj, bmin, bpa = get_beam_info(cubelet_name)

	## produce integrated intensity image and run source finding
	mom0_image = source_find_mom0(cubelet_name, bmaj, bmin, bpa)

	## extract absorption spectrum based source catalog
	## get ellipse parameters
	comp_ra_str, comp_dec_str, comp_ra, comp_dec, comp_a, comp_b, comp_pa = parse_csv('%s_aegean_catalog.csv' % cubelet_name[:-5])
	for i in range(len(comp_ra)):
		## ensure ra/dec strings are consistent
		str_ra_name, str_dec_name = set_ra_dec_output_name(cubelet_name[:-5])
		if i > 0:
			str_ra_name+= '_%s' % i
			str_dec_name+= '_%s' % i

		str_ra, str_dec = convert_hms_dms(comp_ra_str[i], comp_dec_str[i])

		HI_abs_spectrum = extract_mean_pixel_spectrum(mask, cubelet_name, abs_vel_axis, bmaj, bmin, bpa, comp_ra[i], comp_dec[i], comp_a[i], comp_b[i], comp_pa[i])

		## extract emission spectrum
		em_spectrum, em_spectrum_err = extract_emission_spectrum(combined_name, comp_ra[i], comp_dec[i])
		## add em_spectrum in quadrature with rms; -> just use std of emission mean
		#rms_K_spectrum = np.empty(len(em_spectrum_err))
		#rms_K_spectrum.fill(rms_K)
		#em_spectrum_err_final = np.sqrt(em_spectrum_err**2+rms_K_spectrum**2)

		## compute mean emission spectrum
		mean_em_spectrum, em_vel_axis = compute_mean_em_spectrum(combined_name)

		## compute 1-sigma uncertainty in absorption profile
		sigma_cont = np.std(HI_abs_spectrum[mask])
		abs_noise_env = sigma_cont*((tsys+apeff*mean_em_spectrum)/tsys)

		## interpolate noise envelope to velocity axis of absorption
		abs_noise_func = interp1d(em_vel_axis, abs_noise_env, kind = 'linear', fill_value = 0.1/3, bounds_error = False)
		abs_noise_env = abs_noise_func(abs_vel_axis)
		three_sigma_env = 1-3*abs_noise_env
		## interpolate em_spectrum/em_spectrum_err to velocity axis of absorption
		em_spectrum_func = interp1d(em_vel_axis, em_spectrum, kind = 'linear', fill_value = 0.0, bounds_error = False)
		em_spectrum_interp = em_spectrum_func(abs_vel_axis)
		em_spectrum_err_func = interp1d(em_vel_axis, em_spectrum_err, kind = 'linear', fill_value = 0.0, bounds_error = False)
		em_spectrum_err_interp = em_spectrum_err_func(abs_vel_axis)

		## determine where absorption features are > 3-sigma dectection
		detection_inds = np.array(np.where((HI_abs_spectrum) <= three_sigma_env)[0])
		detection_sort = consecutive(detection_inds)

		## generate plot
		## generate_plot(abs_vel_axis, HI_abs_spectrum, em_spectrum_interp, three_sigma_env, abs_noise_env, rms_K, detection_sort, str_ra, str_dec)

		## loop through features to compute spin temperatures
		spin_temp_list = []
		spin_temp_err_list = []
		for l in detection_sort:
			if len(l) > 2:
				## compute spin temperature and add to lists
				spin_temp, spin_temp_err = compute_mean_spin_temperature(em_spectrum_interp[l[0]:l[-1]], em_spectrum_err_interp[l[0]:l[-1]], HI_abs_spectrum[l[0]:l[-1]], abs_noise_env[l[0]:l[-1]], 0.41)
				spin_temp_list.append(spin_temp)
				spin_temp_err_list.append(spin_temp_err)
		## save out results
		pickle_results(abs_vel_axis, HI_abs_spectrum, em_spectrum_interp, three_sigma_env, abs_noise_env, em_spectrum_err_interp, detection_sort, str_ra_name, str_dec_name, spin_temp_list, spin_temp_err_list)

  ## NEW IO CODE YAYAYAYAYAYY #######
	## add a step here that if there is an overlap adds an additional pickle result	
  	new_spectrum = []
  	trial = handle_overlaps(mask, cubelet_name, abs_vel_axis, bmaj, bmin, bpa, comp_ra, comp_dec, comp_a, comp_b, comp_pa)
  	if trial is not None:
  		i = 0
  		print("overlap found !")
  		new_spectrum = trial
  
  		em_spectrum, em_spectrum_err = extract_emission_spectrum(combined_name, comp_ra[i], comp_dec[i])
  	
  		## compute mean emission spectrum
  		mean_em_spectrum, em_vel_axis = compute_mean_em_spectrum(combined_name)
  		## compute 1-sigma uncertainty in absorption profile
  		sigma_cont = np.std(new_spectrum[mask])
  		abs_noise_env = sigma_cont*((tsys+apeff*mean_em_spectrum)/tsys)
  
  		## interpolate noise envelope to velocity axis of absorption
  		abs_noise_func = interp1d(em_vel_axis, abs_noise_env, kind = 'linear', fill_value = 0.1/3, bounds_error = False)
  		abs_noise_env = abs_noise_func(abs_vel_axis)
  		three_sigma_env = 1-3*abs_noise_env
  		## interpolate em_spectrum/em_spectrum_err to velocity axis of absorption
  		em_spectrum_func = interp1d(em_vel_axis, em_spectrum, kind = 'linear', fill_value = 0.0, bounds_error = False)
  		em_spectrum_interp = em_spectrum_func(abs_vel_axis)
  		em_spectrum_err_func = interp1d(em_vel_axis, em_spectrum_err, kind = 'linear', fill_value = 0.0, bounds_error = False)
  		em_spectrum_err_interp = em_spectrum_err_func(abs_vel_axis)
  
  		## determine where absorption features are > 3-sigma dectection
  		detection_inds = np.array(np.where((new_spectrum) <= three_sigma_env)[0])
  		detection_sort = consecutive(detection_inds)
  
  		## loop through features to compute spin temperatures
  		spin_temp_list = []
  		spin_temp_err_list = []
  		for l in detection_sort:
  			if len(l) > 2:
  				## compute spin temperature and add to lists
  				spin_temp, spin_temp_err = compute_mean_spin_temperature(em_spectrum_interp[l[0]:l[-1]], em_spectrum_err_interp[l[0]:l[-1]], HI_abs_spectrum[l[0]:l[-1]], abs_noise_env[l[0]:l[-1]], 0.41)
  				spin_temp_list.append(spin_temp)
  				spin_temp_err_list.append(spin_temp_err)
  		
  		pickle_results(abs_vel_axis, HI_abs_spectrum, em_spectrum_interp, three_sigma_env, abs_noise_env, em_spectrum_err_interp, detection_sort, str_ra, str_dec, spin_temp_list, spin_temp_err_list)
  	else:
  		print('No overlap found')
  	
  		
  	## NEW IO CODE YAYAYAYAYAYY #######

if __name__=='__main__':
	main()
	exit()
