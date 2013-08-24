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
from viewer import load_data,viewer,dictviewer

def make_plots(dirname, fnames):
    # get the highest-level directory, assume it is the target source ID
    obj = os.path.split(dirname)[-1]

    for fn in fnames:
        # each band gets listed...
        if 'band0' not in fn or '.sav' not in fn:
            continue
        # templates:
        # coadd_cleanband0in_clean_music_20130815_jk000.sav
        # 130820_ob1_band0i_clean_music_20130815_map.sav
        obs = "_".join(fn.split("_")[:2])
        #obs = fn[:11]

        print "Working on file ",os.path.join(dirname,obs)

        data = load_data(os.path.join(dirname,obs))
        sm,us = convolve_and_match(data,obj,writefits=False)
        vmin = max([sm[2].min(),-1000])
        vmax = max([sm[2].max(),5000])
        print obj, obs, vmin, vmax

        prefix = os.path.join(dirname,obj+obs)

        pl.clf()
        viewer(data, vmin=vmin, vmax=vmax, cb=True)
        pl.suptitle(obj)
        pl.savefig(prefix+"_quicklook.png",bbox_inches='tight')

        pl.clf()
        dictviewer(sm, vmin=vmin, vmax=vmax, cb=True)
        pl.suptitle(obj)
        pl.savefig(prefix+"_quicklook_smooth.png",bbox_inches='tight')

        pl.clf()
        dictviewer(us, vmin=vmin, vmax=vmax, cb=True)
        pl.suptitle(obj)
        pl.savefig(prefix+"_quicklook_unsharp.png",bbox_inches='tight')

if __name__ == "__main__":
    pl.hot() # select colormap
    for dirpath, dirnames, filenames in os.walk('./'):
        if dirpath != './':
            if 'neptune' in dirpath or 'g34.3' in dirpath:
                continue
            make_plots(dirpath,filenames)
