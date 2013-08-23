from astropy.io import fits
import idlsave
import itertools
import scipy.ndimage
from astropy.nddata.convolution import make_kernel,convolve
import numpy as np
import sys
import os
sys.path.append(os.path.split(os.getcwd())[0])
from convolve_match_makefits import convolve_and_match

data_w49_1   = {k:idlsave.read('w49/130819_ob5_band%ii_clean_music_20130815_map.sav' % k) for k in xrange(4)}
data_w49_2   = {k:idlsave.read('w49/130819_ob6_band%ii_clean_music_20130815_map.sav' % k) for k in xrange(4)}
data_w51_1   = {k:idlsave.read('w51_bgps/130819_ob1_band%ii_clean_music_20130815_map.sav' % k) for k in xrange(4)}
data_w51_2   = {k:idlsave.read('w51_bgps/130819_ob2_band%ii_clean_music_20130815_map.sav' % k) for k in xrange(4)}
data_sgrb2_1 = {k:idlsave.read('sgr_b2/130819_ob3_band%ii_clean_music_20130815_map.sav' % k) for k in xrange(4)}
data_sgrb2_2 = {k:idlsave.read('sgr_b2/130819_ob4_band%ii_clean_music_20130815_map.sav' % k) for k in xrange(4)}


datasets = dict(zip(['W49_1','W49_2','W51_1','W51_2','SgrB2_1','SgrB2_2'],
                    [data_w49_1, data_w49_2, data_w51_1, data_w51_2, data_sgrb2_1, data_sgrb2_2]))

for jj,(k,data) in enumerate(datasets.iteritems()):
    print(jj)
    convolve_and_match(data,k)
