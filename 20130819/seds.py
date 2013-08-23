from astropy.io import fits
import pylab as pl
import idlsave
import itertools
import scipy.ndimage
from astropy.nddata.convolution import make_kernel,convolve
import numpy as np
import sys
import os
sys.path.append(os.path.split(os.getcwd())[0])
from convolve_match_makefits import convolve_and_match
from sed_from_dict import sed_from_dict,plot_sed
from viewer import load_data


def make_plots(dirname, fnames):
    # get the highest-level directory, assume it is the target source ID
    obj = os.path.split(dirname)[-1]

    for fn in fnames:
        # each band gets listed...
        if 'band0' not in fn:
            continue
        obs = fn[:10]

        print "Working on file ",os.path.join(obj,obs)

        data = load_data(os.path.join(obj,obs))
        sm,us = convolve_and_match(data,'W49',writefits=False)
        flux,bg,err = sed_from_dict(sm)

        pl.clf()
        plot_sed(flux,bg,err, label='Smooth')
        flux,bg,err = sed_from_dict(us)
        plot_sed(flux,bg,err, label='Unsharp')

        pl.title(obj)
        pl.legend(loc='best')
        pl.savefig(os.path.join(obj,obs)+"_SED.png",bbox_inches='tight')

if __name__ == "__main__":
    for dirpath, dirnames, filenames in os.walk('./'):
        if dirpath != './':
            if 'neptune' in dirpath or 'g34.3' in dirpath:
                continue
            make_plots(dirpath,filenames)
