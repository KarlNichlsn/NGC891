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

filename_fullcube = "/reduction/erickoch/NGC891/postprocess/NGC891/NGC891_sub+com+ext+vex_co21_2p5kms_pbcorr_trimmed_k.fits"
final_direct = '/reduction/nicholson/final/' #the final directory where I will put everything 
finaldirect = '/reduction/nicholson/final/' #the final directory where I will put everything 


#########################
### GET CUBE & HEADER ### 
#########################


hdu_list = fits.open(filename_fullcube)
cube = SpectralCube.read(hdu_list) #make a spectral cube from the full cube
hdu_list.close()
header = fits.getheader(filename_fullcube) #get the header

cube.allow_huge_operations=True 
mad_std_spectrum = cube.mad_std(axis=(1, 2)) #this was just to investigate, leaving it in anyway


######################
### SIGMA CLIPPING ### 
######################


#cube is sigmaclipped, then the mad stdev is taken without the 3+ sigma values 
new_sclip = sigma_clip(cube, axis=0,stdfunc='mad_std') #standard sigma clipping is 3 sigma 
sigma_clipped_cube = cube.with_mask(~new_sclip.mask) #mask the full cube with the sigma clip

mad_std_map_sclip = sigma_clipped_cube.mad_std(axis=0) #find the std of the sigma clipped cube, essentially an error map

hi_snr_val,lo_snr_val = 5,4 #for the masking

low_snr_mask = (cube > lo_snr_val * mad_std_map_sclip).include() #include the points that are 4 sigma above the sigmaclipped std
high_snr_mask = (cube > hi_snr_val * mad_std_map_sclip).include() #these are masks

structure = np.ones((3, 3, 3), dtype=bool)  #the structure to connect (nearest neighbours)


###############
### MASKING ### 
###############


low_snr_mask_labels, num_labels = nd.label(low_snr_mask,
                                           structure=structure)
high_snr_mask_labels, num_labels_high = nd.label(high_snr_mask,
                                           structure=structure)

num_pixels_in_high_snr_mask = nd.sum(high_snr_mask,
                                     labels=low_snr_mask_labels,
                                     index=range(1, num_labels + 1)) 
                                    

num_pixels_in_low_snr_mask = nd.sum(low_snr_mask,
                                    labels=low_snr_mask_labels,
                                    index=range(1, num_labels + 1)) # +1 offset for mask labels

signal_mask = np.zeros(low_snr_mask.shape, dtype=bool)

#choosing minimum pixels to include in the mask
low_min_pixels = 40
high_min_pixels = 5

for num, (high_pix_num, low_pix_num) in enumerate(zip(num_pixels_in_high_snr_mask, num_pixels_in_low_snr_mask)):
    if high_pix_num >= high_min_pixels and low_pix_num >= low_min_pixels:
                    # This region passes the criteria. Keep it in the mask.
        signal_mask[low_snr_mask_labels == num + 1] = True
        continue

signal_mask_labels, num_labels = nd.label(signal_mask, structure=structure)

#spatial expansion - to include the entirety of the signal
structure = np.ones((3, 3), dtype=bool)

#spectral expansion
structure_spec = np.zeros((3, 3), dtype=bool)# different expansion in spectral direction vs spatial
structure_spec[1, 1] = True

structure = np.dstack([structure_spec, structure, structure_spec])

#mask
signal_mask = nd.binary_dilation(signal_mask, structure=structure, iterations=1)

#the final masked cube
masked_cube = cube.with_mask(signal_mask) 

#writing out the masked cube:
masked_cube_name = 'masked_cube_losnr4_hisnr5_lominpix40_highminpix5.fits'
masked_cube.write(final_direct+masked_cube_name, overwrite=True)



###############
### Moments ###
###############

#### Creating and writing out statistical moments ####


#making moment 0,1,2 for the masked and unmasked cube
#masked_mom0 = masked_cube.with_spectral_unit(u.km/u.s).moment(order=0)
#mom0 = cube.with_spectral_unit(u.km/u.s).moment(order=0)

#masked_mom1 = masked_cube.moment(order=1)
#mom1 = cube.moment(order=1)

#masked_mom2 = masked_cube.moment(order=2)

#masked_linewidth = masked_cube.linewidth_sigma() #v sensitive to noise - won't bother making them w/out mask
#masked_fwhm = masked_cube.linewidth_fwhm() #difference between fwhm and sigma! 

#writing the moments out 
#masked_mom2.write(final_direct+'moments/'+'masked_mom2.fits',overwrite=True)
# masked_mom1.write(final_direct+'moments/'+'masked_mom1.fits',overwrite=True)
# masked_linewidth.write(final_direct+'moments/'+'masked_linewidth.fits',overwrite=True)
# masked_fwhm.write(final_direct+'moments/'+'masked_fwhm.fits',overwrite=True)

# mom0.write(final_direct+'moments/'+'mom0.fits',overwrite=True)
# mom1.write(final_direct+'moments/'+'mom1.fits',overwrite=True)




####################
### Dendrograms ###
####################
  

#dendrogram inputs
minimum = 3
mindelta = 1

#going to just reopen the masked cube
hdu_list = fits.open(final_direct+'masked_cube_losnr4_hisnr5_lominpix40_highminpix5.fits')
image = hdu_list[0]
image_data = hdu_list[0].data 

hdu_list.close() #may need to remove this for further down the code...

#compute dendrogram
d = Dendrogram.compute(image_data, min_value= minimum, min_delta= mindelta , min_npix=10) #minpix unnescessary as cube was masked with min pixels too

d.load(final_direct +'dendrogram_min3_mindelta1_minpix10.fits')  

#reopen denrogram
d = Dendrogram.load_from(final_direct +'dendrogram_min3_mindelta1_minpix10.fits') 

#### plot dendrogram ####
# p = d.plotter()
# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1)

# # Plot the whole tree
# p.plot_tree(ax, color='black')

# ax.set_xlabel("Structure")
# ax.set_ylabel("Intensity")

# plt.savefig(final_direct +'dendrogram_bboxtight.png',bbox_inches='tight')


#### Following are in different py file: ####


#####################
### PPV CATALOGUE ###
#####################     


##############
### SCIMES ###
##############  
        

################
### PYCPROPS ###
################



