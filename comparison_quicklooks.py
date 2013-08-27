from astropy.io import fits
from astropy import coordinates as coords
import montage
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

def make_plots(dirname, fnames):
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

    if 'SCUBA' in files and 'ATLASGAL' in files:
        del files['SCUBA']
        del data['SCUBA']

    print "Making figures for ",dirname, files

    if 'Band0' in files:
        header0 = fits.getheader(files['Band0'])
        center = coords.ICRSCoordinates("%s %s" % (header0['CRVAL1'],header0['CRVAL2']),unit=('deg','deg'))
    elif 'BGPS' in files:
        header0 = fits.getheader(files['BGPS'])
        center = coords.GalacticCoordinates("%s %s" % (header0['CRVAL1'],header0['CRVAL2']),unit=('deg','deg')).fk5
    elif 'ATLASGAL' in files:
        header0 = fits.getheader(files['ATLASGAL'])
        center = coords.GalacticCoordinates("%s %s" % (header0['CRVAL1'],header0['CRVAL2']),unit=('deg','deg')).fk5
    else:
        raise ValueError("Really?  No BGPS, MUSIC, or ATLASGAL?  WTF.  BRB.  BBQ.")

    fig = pl.figure(1)
    pl.clf()
    for ii,(k,fn) in enumerate(files.iteritems()):
        print k,fn,
        try:
            F = aplpy.FITSFigure(fn, subplot=(2,3,ii+1), convention='calabretta', figure=fig, north=True)
            F.show_colorscale(vmin=-0.5,vmax=5,cmap=pl.cm.hot)
            F._ax1.set_title(k)
            F.tick_labels.set_xformat('d.dd')
            if (ii) % 3 == 0:
                F.tick_labels.set_yformat('d.dd')
            else:
                F.tick_labels.hide_y()
            F.add_colorbar()
        except montage.MontageError as e:
            print "Montage error: ",e
            continue
        try:
            F.recenter(center.ra.degree, center.dec.degree, 5/60.)
        except:
            continue
        #pl.subplot(2,3,ii+1)
        #pl.imshow(data[k], vmin=-0.5, vmax=5)
        #pl.title(k)
        #pl.colorbar()

    prefix = os.path.join(dirname,obj)
    pl.savefig(prefix+"_compare_MUSIC_BGPS.png",bbox_inches='tight')


if __name__ == "__main__":
    pl.hot() # select colormap
    for dirpath, dirnames, filenames in os.walk('./'):
        if dirpath != './':
            if 'neptune' in dirpath or 'g34.3' in dirpath:
                continue
            make_plots(dirpath,filenames)
