## CASA python script to extract a sub cube around a source
## excludes short baselines 
## usage: casa --nologfile --log2term -c sub_cube_abs.py

def main():
	phasecenter='J2000 00h40m13.8 +40d50m04.73'
	ref_freq='1.42040571183GHz'
	rest_freq='1.42040571183GHz'
	uvdist='>1.5Klambda'
	output_name='test'
	vis_name='20A-346.sb41042474.eb41074466.59579.84285818287.speclines.ms.split'
	tclean(vis=vis_name, imagename=output_name, reffreq=ref_freq, restfreq=rest_freq, 
		phasecenter=phasecenter,imsize=50,uvrange=uvdist,weighting='natural',gridder='standard',
		pbcor=True,niter=1000,cell='1arcsec')
if __name__=='__main__':
	main()
	exit()
	
