import numpy as np
import pylab as pl
from mad import MAD

band_waves = {0: 1.98, 1: 1.3, 2: 1.0, 3: 0.86}

def sed_from_dict(data, TestRadius=5.1, InnerApp=10, OuterApp=13):
    """
    Extract SEDs around the peak of the band1 image 

    Parameters
    ----------
    data : dict
        A dictionary containing maps sampled on the same pixel grid, e.g.
        as output from convolve_and_match
    TestRadius : float
        Pixel search radius to examine for peak location...
    InnerApp : float
        Inner aperture radius
    OuterApp : float
        Outer aperture radius

    Returns
    -------
    fluxes : length-4 array of fluxes (NOT background-subtracted)
    backgrounds : length-4 array of background median fluxes
    errors : length-4 array of background median absolute deviation values (scaled to sigma)
    """
    shape = data[0].shape
    cube = np.zeros((len(data),shape[0],shape[1]))
    for band in band_waves:
        cube[band,:,:] = data[band]

    yy,xx = np.indices(shape)

    # original radius, centered at center...
    xcen,ycen=np.asarray(shape)/2
    rr = ((yy-ycen)**2 + (xx-xcen)**2)**0.5

    # recenter on peak
    xcen,ycen = np.unravel_index((np.argmax((rr<TestRadius)*cube[1,:,:])),shape)
    rr = ((yy-ycen)**2 + (xx-xcen)**2)**0.5

    # not used counts = np.nansum(np.nansum((rr1 < TestRadius)*(cube>0),axis=1),axis=1)
    mask = (rr < InnerApp)
    fluxes = np.nansum(cube[:,mask],axis=1)
    outermask = (rr > InnerApp) * (rr < OuterApp)
    background = np.median(cube[:,outermask],axis=1)
    errors = MAD(cube[:,outermask],axis=1)

    return fluxes,background,errors

def plot_sed(fluxes, backgrounds, errors, **kwargs):
    """
    Trivial SED plotting
    """
    pl.errorbar(band_waves.values(),fluxes-backgrounds,yerr=errors,marker='s', **kwargs)
    pl.xlabel('$\lambda$ (mm)')
    pl.ylabel('mJy/beam')

def coadd_aligned_dicts(dict1,dict2):
    return {k:dict1[k]+dict2[k] for k in dict1}
