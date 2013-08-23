from pylab import *;import scipy,matplotlib,pyfits;
import idlsave
import itertools
import scipy.ndimage
from astropy.nddata.convolution import make_kernel,convolve
from astropy.io import fits
import glob
import numpy as np


Targets = ('G010','G000','G012','W49','SgrB2','W51')
Targets = ('G000','G000')
Bands = ('0','1','2','3')

band_waves = {0: 1.98, 1: 1.3, 2: 1.0, 3: 0.86}
waves = np.array([1.98, 1.3, 1.0,0.86])
TestRadius = 1.1 # Pixel search area
InnerApp = 20
OuterApp = 23

for target in Targets:
    cube = np.zeros((4,86,86))
    for band in Bands:
        SearchString = '../fits/'+target+'_?_Band'+band+'_smooth.fits'
        CoAdds = glob.glob(SearchString)
        im = np.zeros((86,86))
        for files in CoAdds:
            im += fits.getdata(files)
        im = im/len(CoAdds)
        cube[(np.asarray(band)).astype('int8'),:,:] = im
# Pivot off Max Pixel near center in Band 1
    yy,xx = np.indices(im.shape)
    xcen,ycen=np.asarray(im.shape)/2
    rr = ((yy-ycen)**2 + (xx-xcen)**2)**0.5
    xcen,ycen = np.unravel_index((np.argmax((rr<TestRadius)*cube[1,:,:])),im.shape)
    rr = ((yy-ycen)**2 + (xx-xcen)**2)**0.5
    rr1,bandidx = np.meshgrid(rr,np.arange((4)))
    rr1.shape = cube.shape
    bandidx.shape = cube.shape
    counts = np.nansum(np.nansum((rr1 < TestRadius)*(cube>0),axis=1),axis=1)
    fluxes = np.nansum(np.nansum(cube*(rr1 < TestRadius)*(cube>0),axis=1),axis=1)
    medians = np.zeros(4)
    mads = np.zeros(4)
    for i in xrange(4):
        indices = (rr1>InnerApp)*(rr1<OuterApp)*(cube==cube)*(bandidx == i)
        medians[i] = np.median(cube[indices])
        mads[i] = np.median(np.abs(cube[indices]-medians[i]))/0.6745
    fluxes -= medians
    figure(figsize=(5,5))
    title(target)
    scatter(waves,fluxes)
    errorbar(waves,fluxes,yerr=mads)
    xlabel('$\lambda$ (mm)')
    ylabel('Who knows')
    savefig("%s_sed.pdf" % target,bbox_inches="tight")
    clf()
