import pylab as pl
import idlsave

def load_data(obsname):
    """
    load data file given an obs date/number 
    
    Example
    -------
    >>> data = load_data('gal_010.47+00.03/130820_ob3')
    """
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
