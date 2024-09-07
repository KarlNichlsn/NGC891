from spectral_cube import SpectralCube
from astrodendro.dendrogram import Dendrogram
from astropy.io import fits
from scimes import SpectralCloudstering
from astrodendro.analysis import ppv_catalog
from astropy import units as u
from astropy.wcs import wcs
from astropy.io import fits
from pycprops.cloudalyze import cloudalyze
import numpy as np
import skimage.segmentation as seg

#locations in RTDC
dendroloc = '/reduction/nicholson/masked_cubes/' 
propsloc = '/reduction/nicholson/pycprops_out/'
#distance chosen is 8.5 Mpc
dist = 8.5*u.Mpc



#reading in the masked cube, the low and high snr (4,5) masking relative to sigma in the channels, minimum pixels (40)
cube = SpectralCube.read(dendroloc+'cube_lo_snr_4_hi_snr_5_lowminpix_40_masked.fits')

#reading in the dendrogram, same meanings as above
d = Dendrogram.load_from(dendroloc+ 'dendrogram_lo_snr_4_hi_snr_5_lowminpix_40_mindelta_1.fits')

metadata = {}
metadata['data_unit'] = u.K #unit of spectral cube is Kelvin 
metadata['wavelength'] = (299492758 * u.m /u.s / (230.530 * u.GHz)).to(u.m) #line wavelength
metadata['beam_major'] =  cube.beam.major.to(u.deg)
metadata['beam_minor'] =  cube.beam.major.to(u.deg)
metadata['spatial_scale'] = cube.header['CDELT2'] * u.deg #degrees
metadata['wcs'] = cube.wcs #world coords

#astrodendro position-position-velocity
cat = ppv_catalog(d, metadata)
rms1 = 1.0

#SCIMES
dclust = SpectralCloudstering(d, cat,
                              header=cube.header,
                              save_all_leaves=True,
                              user_iter=10, rms=rms1)

# Re-label the assignment cube with the clustered ids
# Shift to cprops expectation
asgn = dclust.clusters_asgn.data + 1

#skimage segmentation 
asgn,_, _ = seg.relabel_sequential(asgn)
hdu = fits.PrimaryHDU(asgn, header=cube.header) 

#writing out the clustered fits file
pycproptsout_name = 'pycpropsout_asgn_29thOct.fits'
hdu.writeto(propsloc+pycproptsout_name, overwrite=True)

#distance chosen is 8.5 Mpc
#dist = 8.5*u.Mpc

#pycprops cloud analysis, alpha_CO chosen as 6.7
catalog = cloudalyze(cube, asgn,
                     distance=dist,
                     alphaCO=6.7, extrapolate=True,
                     bootstrap=100, rmstorad=1.91,
                     verbose=True,channelcorr=0.0, noise=rms1)

#writing out the catalogue
catalgoue_name = 'pycprops_bootstrap100.fits'
catalog.write(propsloc+catalgoue_name, overwrite=True)
