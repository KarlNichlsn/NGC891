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

plt.rcParams.update({'font.size': 20})

finaldirect = '/Users/karlnicholson/Desktop/Harvard-Smithsonian/KN_NGC891/minicube/'
subcube_filename = 'masked_cube_losnr4_hisnr5_lominpix40_highminpix5_subcube_x1130.1370_y1844_1888.fits'
subcube_region = 'subcube_region_cutout.reg'


#########################
### GET CUBE & HEADER ### 
#########################

hdu_list = fits.open(final_direct+subcube_filename)
image = hdu_list[0]
image_data = hdu_list[0].data 

header = fits.getheader(filename_fullcube) #header
masked_cube = SpectralCube.read(hdu_list) #masked cube spectral cube 
hdu_list.close()


###############
### Moments ###
###############


#making moment 0,1,2 for the masked and unmasked cube
masked_mom0 = masked_cube.with_spectral_unit(u.km/u.s).moment(order=0)
# mom0 = cube.with_spectral_unit(u.km/u.s).moment(order=0)

masked_mom1 = masked_cube.moment(order=1)

#writing the moments out 
masked_mom0.write(final_direct+'subcube_masked_mom0.fits',overwrite=True)
masked_mom1.write(final_direct+'subcube_masked_mom1.fits',overwrite=True)

"""wcsmom0 = WCS(moment_0.header)
fig = plt.figure(figsize=(18, 12))
ax = fig.add_subplot(1,1,1,projection=wcsmom0)

im = ax.imshow(moment_0.data, interpolation='nearest',cmap='gray',origin='lower')#,norm=norm)
# ax.arrow(1100,1200, 1100,1300)
x_start,y_start=1295,2095
ax.arrow(x_start, y_start, 0,20, color='red', width=0.5, head_width=0.5, head_length=1)
ax.arrow(x_start, y_start, 20,0, color='red', width=0.5, head_width=0.5, head_length=1)
ax.text(x_start+5, y_start+15,'N',size=1,color='red')
ax.text(x_start+10, y_start-10,'E',size=1,color='red')

x_start,y_start=1150,2250
distance = 5
angle_rad = np.radians(-67)
# Calculate the endpoint of the line
x_end = x_start + distance * np.cos(angle_rad)
y_end = y_start + distance * np.sin(angle_rad)

# Plot the line
ax.plot([x_start, x_end], [y_start, y_end], color='blue', linewidth=0.5)
ax.text(x_start-15, y_start-15,'46 parsec',size=0.5,color='blue',rotation=-67)
# ax.set_xlim(1100,1200)
# ax.set_ylim(2300,2400)
plt.savefig('/Users/karlnicholson/Downloads/'+'moment0plotforposter.png',bbox_inches='tight',dpi=1200)"""


####################
### Dendrograms ###
####################
  

minimum = 3
mindelta = 1

#going to just reopen the masked cube
# hdu_list = fits.open(final_direct+'masked_cube_losnr4_hisnr5_lominpix40_highminpix5.fits')
# image = hdu_list[0]
# image_data = hdu_list[0].data 

##added this:
# hdu_list.close()
#dk if it will be a problem though ????

d = Dendrogram.compute(image_data, min_value= minimum, min_delta= mindelta , min_npix=10)
# d.save_to(final_direct +'dendrogram_min3_mindelta1_minpix10.fits')  


#####################
### PPV CATALOGUE ###
#####################     


##Metadata##
metadata = {}
metadata['data_unit'] = u.K
metadata['wavelength'] = (299492758 * u.m /u.s / (230.530 * u.GHz)).to(u.m)
metadata['beam_major'] =  masked_cube.beam.major.to(u.deg)
metadata['beam_minor'] =  masked_cube.beam.major.to(u.deg)
metadata['spatial_scale'] = masked_cube.header['CDELT2'] * u.deg
metadata['wcs'] = masked_cube.wcs


## Making the PPV Catalogue 
cat = ppv_catalog(d, metadata=metadata)


##############
### SCIMES ###
##############  
        

rms1 = 1.0
dclust = SpectralCloudstering(d, cat, header=header,save_all_leaves=True,user_iter=10,rms=rms1)#added save all leaves


#######################
### DClust plotting: ##
#######################

moment_0 = Projection.from_hdu(fits.open(final_direct+'subcube_masked_mom0.fits'))

wcsmom0 = WCS(moment_0.header)
fig = plt.figure(figsize=(18, 12),dpi=1200)
ax = fig.add_subplot(1,1,1,projection=wcsmom0)

im = ax.imshow(moment_0.data, interpolation='nearest',cmap='gray',origin='lower')#,norm=norm)

clusts = dclust.leaves
colors = dclust.colors

count = 0
for c in clusts:
    #plots the contour of the clouds

    mask = d[c].get_mask()
    mask_coll = np.nanmax(mask,axis=0)

    ax.contour(mask_coll, colors=colors[count], linewidths=2, levels = [0])

    count = count+1
    
plt.savefig(final_direct+'moment0plotwithcontours.png',bbox_inches='tight',dpi=1200)















