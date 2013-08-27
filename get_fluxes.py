from astropy.io import fits
from astropy import coordinates as coords
import aplpy
import pylab as pl
import idlsave
import itertools
import scipy.ndimage
from astropy.nddata.convolution import make_kernel,convolve
import numpy as np
import sys
import os
sys.path.append(os.path.split(os.getcwd())[0])

keys = {"Band%i" % ii:'Band%i_native' % ii for ii in xrange(4)}
keys['BGPS'] = 'MOSAIC_MAP_V2'
keys['ATLASGAL'] = 'ATLASGAL'
keys['SCUBA'] = 'scuba_850um'

def get_peaks(dirname, fnames):
    # get the highest-level directory, assume it is the target source ID
    obj = os.path.split(dirname)[-1]

    files = {}
    data = {}

    for fn in fnames:
        for k,v in keys.iteritems():
            if v in fn:
                if k in files and not 'coadd' in fn:
                    continue
                files[k] = os.path.join(dirname,fn)
                data[k] = fits.getdata(os.path.join(dirname,fn))

    if len(files) == 0:
        return

    peaks = {}
    for ii,(k,fn) in enumerate(files.iteritems()):
        d = data[k]
        OK = d==d
        if OK.sum() == 0:
            continue
        print k,fn,d[OK].max()
        peaks[k] = d[OK].max()

    return peaks

if __name__ == "__main__":
    pl.hot() # select colormap
    peaks = {}
    for dirpath, dirnames, filenames in os.walk('./'):
        if dirpath != './':
            if 'neptune' in dirpath or 'g34.3' in dirpath:
                continue
            peaks[dirpath] = get_peaks(dirpath,filenames)
