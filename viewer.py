import pylab as pl
import idlsave
import os
from astropy.io import ascii,fits
from astropy import coordinates as coords
import numpy as np

def load_header(obsstruct):
    tbl = ascii.read('targets.txt')
    source_name = obsstruct['SOURCE_NAME'][0]
    match = np.array([source_name == x.lower() for x in tbl['col1']])
    ra,dec = tbl['col4'][match][0],tbl['col5'][match][0]
    coord = coords.ICRSCoordinates(ra+" "+dec,unit=('hour','deg'))

    header = fits.Header()
    header['CRVAL1'] = coord.ra.degree
    header['CRVAL2'] = coord.dec.degree
    shape = obsstruct['MAP'][0].shape
    header['CRPIX1'] = shape[1]/2. + 1 # may be offset by 0.5-1 pixel!?
    header['CRPIX2'] = shape[2]/2. + 1 # may be offset by 0.5-1 pixel!?
    header['CDELT1'] = obsstruct['OMEGA_PIX_AM'][0]**0.5 / 60
    header['CDELT2'] = obsstruct['OMEGA_PIX_AM'][0]**0.5 / 60
    header['CUNIT1'] = 'deg'
    header['CUNIT2'] = 'deg'
    header['CTYPE1'] = 'RA---TAN'
    header['CTYPE2'] = 'DEC--TAN'
    header['BUNIT'] = 'mJy/beam'

    return header

def load_data(obsname):
    """
    load data file given an obs date/number 
    
    Example
    -------
    >>> data = load_data('gal_010.47+00.03/130820_ob3')
    """
    if 'coadd' in obsname:
        obsdir,obsname = os.path.split(obsname)
        filename = obsname+"_clean_music_20130815.sav"
        filename = filename[:15]+"%i"+filename[16:]
        filename = os.path.join(obsdir,filename)
    else:
        filename = obsname+"_band%ii_clean_music_20130815_map.sav"
    data = {k:idlsave.read(filename % k, verbose=False) for k in xrange(4)}
    return data

def viewer(data, cb=False, **kwargs):
    """
    Simple viewer.  Load the data with load_data, then quick-look with this.

    Examples
    --------
    >>> data = load_data('gal_010.47+00.03/130820_ob3')
    >>> viewer(data,vmin=-1000,vmax=8000)
    """
    for ii,(k,v) in enumerate(data.iteritems()):
        pl.subplot(2,2,ii+1)
        pl.imshow(v.mapstruct.map[0], **kwargs)
        if cb:
            pl.colorbar()

def dictviewer(data, cb=False, **kwargs):
    """
    Another simple viewer

    Examples
    --------
    >>> import convolve_match_makefits
    >>> sm,us = convolve_match_makefits.convolve_and_match(data,'G10',writefits=False)
    >>> dictviewer(sm)
    >>> dictviewer(us, cb=True)
    """
    for ii,(k,v) in enumerate(data.iteritems()):
        pl.subplot(2,2,ii+1)
        pl.imshow(v, **kwargs)
        if cb:
            pl.colorbar()
