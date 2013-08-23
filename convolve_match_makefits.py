import numpy as np
from astropy.io import fits
import scipy.ndimage
from astropy.nddata.convolution import make_kernel,convolve

def convolve_and_match(data, objectname, clobber=True, writefits=True, unsharpscale=80):
    """
    Convolve all bands to Band 0 resolution and resample all to Band 3
    pixelization.  The data values are appropriately scaled after smoothing
    such that they are in units of mJy/beam, where the beam is the Band 0 beam.
    Sanity checks on this front are warranted.

    Parameters
    ----------
    data : dict of IDLSAVE structs
    objectname : str
        Object name for saving purposes
    clobber : bool
        Overwrite FITS files if they exist?
    writefits : bool
        Write the fits files to disk?
    unsharpscale : float
        Unsharp mask angular scale in arcseconds

    Returns
    -------
    smoothdict : dict
        A dictionary of band number : smoothed & resampled map
    unsharpdict : dict
        A dictionary of band number : unsharp-masked map
    """
    smoothed = {}
    unsharped = {}

    # the band 3 map is used as the reference pixel scale
    band3 = data[3].mapstruct.map[0]
    yy1,xx1 = grid1 = np.indices(band3.shape)

    pixscale = float(data[3].mapstruct['OMEGA_PIX_AM']**0.5 / 60.)
    pixscale_as = pixscale*3600

    if writefits:

        header = fits.Header()
        header['CDELT1'] = pixscale
        header['CDELT2'] = pixscale
        ffile = fits.PrimaryHDU(data=band3, header=header)

    for ii in xrange(4):
        # grab the map from band i
        m = data[ii].mapstruct.map[0]
        # ratio of map sizes (assumes they're symmetric, sort of)
        ratio = m.shape[0]/float(band3.shape[0])

        # rescale the band i map to band3 scale
        newm = scipy.ndimage.map_coordinates(np.nan_to_num(m), grid1*ratio)
        # flag out the NANs again (we had to make them zeros for the previous step to work)
        bads = scipy.ndimage.map_coordinates(np.array(m!=m,dtype='float'), grid1*ratio)
        newm[bads>0.5] = np.nan

        # Determine the appropriate convolution kernel size
        # beamsize = np.array([60*(v.mapstruct.omega_beam_am/np.pi/2.)**0.5 * (8*np.log(2))**0.5 for v in data.values()]).ravel()
        # array([ 45.00000061,  31.00000042,  25.00000034,  23.00000031])
        # these check out...
        beamsize_delta = (np.abs(data[ii].mapstruct['OMEGA_BEAM_AM']-data[0].mapstruct['OMEGA_BEAM_AM'])/np.pi/2)**0.5

        # pixel scale in the *output* image
        am_per_pix = data[3].mapstruct['OMEGA_PIX_AM']**0.5
        # kernel width in pixels
        kernelwidth = beamsize_delta/am_per_pix

        if kernelwidth > 0:
            kernel = make_kernel.make_kernel(band3.shape, kernelwidth=kernelwidth)
            newm = convolve.convolve_fft(newm, kernel, interpolate_nan=True)

        # rescale the values to be mJy per a much larger beam; the values should therefore be larger
        newm *= data[0].mapstruct['OMEGA_BEAM_AM']/data[ii].mapstruct['OMEGA_BEAM_AM']

        # Now do an unsharp mask with a fairly large kernel to ensure the spatial filters are identical
        kernel = make_kernel.make_kernel(band3.shape, kernelwidth=unsharpscale/pixscale_as)
        smm = convolve.convolve_fft(newm,kernel, interpolate_nan=True)

        smoothed[ii] = newm
        unsharped[ii] = newm-smm

        if writefits:
            ffile.data = newm

            ffile.writeto("%s_Band%i_smooth.fits" % (objectname,ii), clobber=clobber)

            ffile.data = newm - smm

            ffile.writeto("%s_Band%i_smooth_unsharp.fits" % (objectname,ii), clobber=clobber)

    return smoothed,unsharped
