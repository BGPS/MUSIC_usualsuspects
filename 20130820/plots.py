from pylab import *;import numpy,scipy,matplotlib,pyfits;
import idlsave
import itertools
import scipy.ndimage
from astropy.nddata.convolution import make_kernel,convolve

data_g010_1 = {k:idlsave.read('gal_010.47+00.03/130820_ob3_band%ii_clean_music_20130815_map.sav' % k) for k in xrange(4)}
data_g010_2 = {k:idlsave.read('gal_010.47+00.03/130820_ob4_band%ii_clean_music_20130815_map.sav' % k) for k in xrange(4)}
data_g000_1 = {k:idlsave.read('gal_0.253+0.016/130820_ob1_band%ii_clean_music_20130815_map.sav' % k) for k in xrange(4)}
data_g000_2 = {k:idlsave.read('gal_0.253+0.016/130820_ob2_band%ii_clean_music_20130815_map.sav' % k) for k in xrange(4)}
data_g012_1 = {k:idlsave.read('gal_012.21-00.10/130820_ob5_band%ii_clean_music_20130815_map.sav' % k) for k in xrange(4)}
data_g012_2 = {k:idlsave.read('gal_012.21-00.10/130820_ob6_band%ii_clean_music_20130815_map.sav' % k) for k in xrange(4)}

datasets = dict(zip(['G010_1','G010_2','G000_1','G000_2','G012_1','G012_2'],
                    [data_g010_1, data_g010_2, data_g000_1, data_g000_2, data_g012_1, data_g012_2]))

band_waves = {0: 1.98, 1: 1.3, 2: 1.0, 3: 0.86}

hot()

for jj,(k,data) in enumerate(datasets.iteritems()):
    figure(jj)
    clf()

    suptitle(k)
    for ii in xrange(4):
        subplot(2,2,ii+1)
        imshow(data[ii].mapstruct.map[0])
    
for jj,(k,data) in enumerate(datasets.iteritems()):
    figure(jj+6)
    clf()
#    suptitle(k)
    figure(jj+18)
    clf()

    for ii,(x,y) in enumerate(itertools.combinations(xrange(4),2)):
        figure(jj+6)
        x,y=y,x
        subplot(3,2,ii+1)
        m1 = data[x].mapstruct.map[0]
        m2 = data[y].mapstruct.map[0]
        yy1,xx1 = grid1 = np.indices(m1.shape)
        ratio = m2.shape[0]/float(m1.shape[0])
        newm2 = scipy.ndimage.map_coordinates(np.nan_to_num(m2), grid1*ratio)
        mask = (m1==m1)*(newm2==newm2)*(m1>0)*(newm2>0)
        #mask *= (m1>m1[mask].std()) * (newm2>newm2[mask].std())
        rr = ((xx1-m1.shape[1]/3.)**2 + (yy1-m1.shape[0]/3.)**2)**0.5
        #mask *= rr < 5

        beamsize_delta = (np.abs(data[y].mapstruct['OMEGA_BEAM_AM']-data[x].mapstruct['OMEGA_BEAM_AM'])/np.pi/2)**0.5
        am_per_pix = data[y].mapstruct['OMEGA_PIX_AM']**0.5
        kernelwidth = beamsize_delta/am_per_pix
        kernel = make_kernel.make_kernel(m1.shape, kernelwidth=kernelwidth)#, normalize_kernel=np.max)
        newm1 = convolve.convolve_fft(m1, kernel, interpolate_nan=True)
        newm1 *= data[y].mapstruct['OMEGA_BEAM_AM']/data[x].mapstruct['OMEGA_BEAM_AM']

        plot(newm1[mask],newm2[mask],'o', alpha=0.5)

        #kernel = make_kernel.make_kernel(m1.shape, kernelwidth=5)
        #smm1 = convolve.convolve_fft(newm1,kernel, interpolate_nan=True)
        #smm2 = convolve.convolve_fft(newm2,kernel, interpolate_nan=True)
        #plot((newm1-smm1)[mask],(newm2-smm2)[mask],'rs', alpha=0.5)

        xx = np.linspace(newm1[mask].min(),newm1[mask].max())
        plot(xx, xx*(band_waves[x]/band_waves[y])**3.5, 'r--')
        plot(xx, xx*(band_waves[x]/band_waves[y])**2, 'b:')
        subplots_adjust(hspace=0.3,wspace=0.3)
        xlabel("Band %i: %0.2f mm" % (x,band_waves[x]))
        ylabel("Band %i: %0.2f mm" % (y,band_waves[y]))

        figure(jj+18)
        subplot(3,4,ii*2+1)
        imshow(newm1)
        title("Band %i: %0.2f mm plot %i" % (x,band_waves[x],ii))
        subplot(3,4,ii*2+2)
        imshow(newm2)
        title("Band %i: %0.2f mm plot %i" % (y,band_waves[y],ii))

    figure(jj+6)
    savefig("CCD_%s.pdf" % k, bbox_inches='tight')

for jj,(k,data) in enumerate(datasets.iteritems()):
    figure(jj+12)
    clf()
    suptitle(k.replace("_"," "))

    band3 = data[3].mapstruct.map[0]
    yy1,xx1 = grid1 = np.indices(band3.shape)

    for ii in xrange(4):
        subplot(2,2,ii+1)
        m = data[ii].mapstruct.map[0]
        ratio = m.shape[0]/float(band3.shape[0])
        newm = scipy.ndimage.map_coordinates(np.nan_to_num(m), grid1*ratio)
        bads = scipy.ndimage.map_coordinates(np.array(m!=m,dtype='float'), grid1*ratio)
        newm[bads>0.5] = np.nan

        #kernel = make_kernel.make_kernel(band3.shape, kernelwidth=5)
        #smm = convolve.convolve_fft(newm,kernel, interpolate_nan=True)

        beamsize_delta = (np.abs(data[ii].mapstruct['OMEGA_BEAM_AM']-data[0].mapstruct['OMEGA_BEAM_AM'])/np.pi/2)**0.5
        am_per_pix = data[0].mapstruct['OMEGA_PIX_AM']**0.5
        kernelwidth = beamsize_delta/am_per_pix
        if kernelwidth > 0:
            kernel = make_kernel.make_kernel(band3.shape, kernelwidth=kernelwidth)
            newm = convolve.convolve_fft(newm, kernel, interpolate_nan=True)
            print ii,kernelwidth

        imshow(newm)

    savefig("resampled_%s.pdf" % k, bbox_inches='tight')

for jj,(k,data) in enumerate(datasets.iteritems()):
    figure(jj+12)
    clf()
    suptitle(k.replace("_"," "))

    band3 = data[3].mapstruct.map[0]
    yy1,xx1 = grid1 = np.indices(band3.shape)

    for ii in xrange(4):
        subplot(2,2,ii+1)
        m = data[ii].mapstruct.map[0]
        ratio = m.shape[0]/float(band3.shape[0])
        newm = scipy.ndimage.map_coordinates(np.nan_to_num(m), grid1*ratio)
        bads = scipy.ndimage.map_coordinates(np.array(m!=m,dtype='float'), grid1*ratio)
        newm[bads>0.5] = np.nan

        #kernel = make_kernel.make_kernel(band3.shape, kernelwidth=5)
        #smm = convolve.convolve_fft(newm,kernel, interpolate_nan=True)

        beamsize_delta = (np.abs(data[ii].mapstruct['OMEGA_BEAM_AM']-data[0].mapstruct['OMEGA_BEAM_AM'])/np.pi/2)**0.5
        am_per_pix = data[0].mapstruct['OMEGA_PIX_AM']**0.5
        kernelwidth = beamsize_delta/am_per_pix
        if kernelwidth > 0:
            kernel = make_kernel.make_kernel(band3.shape, kernelwidth=kernelwidth)
            newm = convolve.convolve_fft(newm, kernel, interpolate_nan=True)

        newm /= data[3].mapstruct['OMEGA_BEAM_AM']/data[ii].mapstruct['OMEGA_BEAM_AM']
        mask = (band3==band3)*(newm==newm)*(band3>0)*(newm>0)
        plot(band3[mask],newm[mask],'o', alpha=0.5)

        xx = np.linspace(band3[mask].min(),band3[mask].max())
        plot(xx, xx*(band_waves[x]/band_waves[y])**3.5, 'r--')
        plot(xx, xx*(band_waves[x]/band_waves[y])**2, 'b:')
        subplots_adjust(hspace=0.3,wspace=0.3)
        xlabel("Band %i: %0.2f mm" % (x,band_waves[0]))
        ylabel("Band %i: %0.2f mm" % (y,band_waves[y]))


    savefig("CCD_resampledband3_%s.pdf" % k, bbox_inches='tight')
