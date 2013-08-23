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

data_g010_1 = {k:idlsave.read('gal_010.47+00.03/130820_ob3_band%ii_clean_music_20130815_map.sav' % k) for k in xrange(4)}
data_g010_2 = {k:idlsave.read('gal_010.47+00.03/130820_ob4_band%ii_clean_music_20130815_map.sav' % k) for k in xrange(4)}
data_g000_1 = {k:idlsave.read('gal_0.253+0.016/130820_ob1_band%ii_clean_music_20130815_map.sav' % k) for k in xrange(4)}
data_g000_2 = {k:idlsave.read('gal_0.253+0.016/130820_ob2_band%ii_clean_music_20130815_map.sav' % k) for k in xrange(4)}
data_g012_1 = {k:idlsave.read('gal_012.21-00.10/130820_ob5_band%ii_clean_music_20130815_map.sav' % k) for k in xrange(4)}
data_g012_2 = {k:idlsave.read('gal_012.21-00.10/130820_ob6_band%ii_clean_music_20130815_map.sav' % k) for k in xrange(4)}

datasets = dict(zip(['G010_1','G010_2','G000_1','G000_2','G012_1','G012_2'],
                    [data_g010_1, data_g010_2, data_g000_1, data_g000_2, data_g012_1, data_g012_2]))

for jj,(k,data) in enumerate(datasets.iteritems()):
    print(jj)
    convolve_and_match(data,k)
