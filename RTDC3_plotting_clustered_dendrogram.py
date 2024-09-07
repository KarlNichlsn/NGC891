import matplotlib.pyplot as plt
import numpy as np
from astrodendro.dendrogram import Dendrogram
from astrodendro.analysis import ppv_catalog
from astropy.io import fits
from astropy import units as u 
from spectral_cube import SpectralCube
from astropy.wcs import WCS
import scipy.ndimage as nd
from pycprops.pycprops import cloudalyze #changed 
from astropy.io import ascii as asc
from radio_beam import Beam
from astropy.io import fits
from astropy.table import Table
from scimes import SpectralCloudstering
from radio_beam import Beam
import scipy.ndimage as nd
from astropy.stats import sigma_clip
import skimage.segmentation as seg
# from astrodendro.analysis import PPVStatistic




#### reruns ppv catalogue and scimes - can just load them instead to reduce wait time ####


plt.rcParams.update({'font.size': 20})

filename_fullcube = "/reduction/erickoch/NGC891/postprocess/NGC891/NGC891_sub+com+ext+vex_co21_2p5kms_pbcorr_trimmed_k.fits"
final_direct = '/reduction/nicholson/final/' #the final directory where I will put everything 


#open the masked cube
hdu_list = fits.open(final_direct+'masked_cube_losnr4_hisnr5_lominpix40_highminpix5.fits')
image = hdu_list[0]
image_data = hdu_list[0].data 

header = fits.getheader(filename_fullcube) #header
masked_cube = SpectralCube.read(hdu_list) #masked cube spectral cube 
hdu_list.close()


#open dendrogram
d = Dendrogram.load_from(final_direct +'dendrogram_min3_mindelta1_minpix10.fits') 
p = d.plotter()

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

# Plot the whole tree
p.plot_tree(ax, color='black')

ax.set_xlabel("Structure")
ax.set_ylabel("Intensity")
ax.tick_params(axis='both', which='major', length=10, width=1.5)#,direction='inout')
ax.tick_params(axis='both', which='minor', length=4, width=1.2)#,direction='inout')

plt.savefig(final_direct +'14-Jan-Dendrogram.png',bbox_inches='tight')

#####################
### PPV CATALOGUE ###
#####################     


my_beam = Beam.from_fits_header(header) #don't think I need this, spectral cube works

##Metadata##
metadata = {}
metadata['data_unit'] = u.K
metadata['wavelength'] = (299492758 * u.m /u.s / (230.530 * u.GHz)).to(u.m)
metadata['beam_major'] =  masked_cube.beam.major.to(u.deg)
metadata['beam_minor'] =  masked_cube.beam.major.to(u.deg)
metadata['spatial_scale'] = masked_cube.header['CDELT2'] * u.deg
metadata['wcs'] = masked_cube.wcs


## Making the Catalogue & writing it out
cat = ppv_catalog(d, metadata=metadata)
# cat.write(final_direct+'ppvcatalog_masked_mps.fits',overwrite=True)


##############
### SCIMES ###
##############  
        

rms1 = 1.0
dclust = SpectralCloudstering(d, cat, header=header,save_all_leaves=True,user_iter=10,rms=rms1)#added save all leaves

#the clustered dendrogram - showing which parts of the dendrogram are within which clouds
#saving this dendrogram
dclust.showdendro(savefile = final_direct +'14-Jan-SCIMES_Clustered_Dendrogram.png')


