#!/usr/bin/env python                                                          #
#                                                                              #
# Autor: Michela Negro, University of Torino.                                  #
# On behalf of the Fermi-LAT Collaboration.                                    #
#                                                                              #
# This program is free software; you can redistribute it and/or modify         #
# it under the terms of the GNU GengReral Public License as published by       #
# the Free Software Foundation; either version 3 of the License, or            #
# (at your option) any later version.                                          #
#                                                                              #
#------------------------------------------------------------------------------#

import os
import re
import glob
import numpy as np
import pyfits as pf
import healpy as hp
import scipy.special as sp
from GRATools import GRATOOLS_OUT
from GRATools import GRATOOLS_CONFIG
from GRATools.utils.logging_ import logger
from GRATools.utils.matplotlib_ import pyplot as plt
from GRATools.utils.matplotlib_ import overlay_tag, save_current_figure
from GRATools.utils.gSpline import xInterpolatedBivariateSplineLinear
from GRATools.utils.gSpline import xInterpolatedUnivariateSplineLinear
from GRATools.utils.gFTools import get_energy_from_txt



def get_pl_vs_th(l, _costh):
    """return an array of Pl values correspoding to a given _costh array
    """
    _pl_th = sp.eval_legendre(l, _costh)
    return _pl_th

def get_l_max_autocorr(wbeam_spline, wpix_spline, energy, cut=0.5):
    """discarded...
    """
    _l = np.linspace(0, len(wbeam_spline.x), 2000)
    wb_en = wbeam_spline.hslice(energy)
    wl2 = (wb_en(_l))**2#*wpix_spline(_l))
    cut_arr = np.empty(len(wl2))
    cut_arr.fill(cut)
    mask = np.isclose(wl2, cut_arr, atol=1e-3)
    index = np.where(mask == True)
    l_max = _l[index]
    if len(l_max) < 1:
        l_max = len(wbeam_spline.x)
    else:
        l_max = l_max[-1]
    return l_max
    
    
def build_wpix(nside, l_max=1500):
    """Returns a spline of the pixel window function.
    
       nside : int
           nside of the map you are analizing
       l_max : float
           the maximum multipol at wich to calculate the pixel window function
    """
    wpix = hp.sphtfunc.pixwin(nside)[:l_max]
    fmt = dict(xname='$l$', xunits='', yname='W_{pixel}', \
                           yunits='')
    wpix_spline = xInterpolatedUnivariateSplineLinear(np.arange(l_max),
                                                      wpix, **fmt)
    return wpix_spline

def build_wbeam(psf, _l, out_file):
    """Calculates the Wbeam(l, E) and return a bivariate slpine.
    
       psf : numpy spline
          psf spline generated by get_psf
       _l : numpy array
          array of multipoles at wicht to compute the Wbeam
       out_file : str
          name of the outut txt file that will be created
    """
    out_txt = open(out_file, 'w')
    _en = psf.x
    _wb = []
    energy = str(list(_en)).replace('[','').replace(']','').replace(', ', ' ')
    out_txt.write('l\t%s\n'%energy)
    for e in _en:
        wb_e = np.array([])
        psf_th = psf.vslice(e)
        for l in _l:
            pl_th = get_pl_vs_th(l, np.cos(psf_th.x))
            fmt = dict(xname='th', xunits='rad', yname='convolution', \
                           yunits='MeV')
            _conv = xInterpolatedUnivariateSplineLinear(psf_th.x, \
                                    np.sin(psf_th.x)*psf_th.y*pl_th, **fmt)
            wb_e_l = min(1., 2*np.pi*(_conv.integral(np.amin(psf_th.x), \
                                                         np.amax(psf_th.x))))
            print 'Wbeam(%i, %.2f)'%(l, e), wb_e_l
            wb_e = np.append(wb_e, [wb_e_l])
        _wb.append(wb_e)
    _wb = np.array(_wb)
    fmt = dict(xname='$l$', xunits='', yname='Energy',
                   yunits='MeV', zname='W$_{beam}$(E,$l$)')
    wbeam = xInterpolatedBivariateSplineLinear(_l, _en, _wb.T, **fmt)
    for i, l in enumerate(_l):
        wb = str(list(_wb.T[i])).replace('[','').replace(']','').\
            replace(', ', ' ')
        out_txt.write('%i\t%s\n'%(l, wb))
    out_txt.close()
    return wbeam

def get_powerlaw_spline(index=2.3):
    """Returns a simple power law spline with a given spectral index.
       
       index : float
          spectral index of the powerlaw you want to generate
    """
    EMIN = 2
    EMAX = 6
    en = np.logspace(EMIN, EMAX, 2000)
    pl = en**(-index)
    fmt = dict(xname='$E$', xunits='MeV', yname='E$^{-%.1f}$'%index,
                   yunits='')
    spectrum = xInterpolatedUnivariateSplineLinear(en, pl, **fmt)
    return spectrum

def get_intensity_spectrum(label, mask_file=None):
    """discarded... 
       but it evaluates the counts spectra of a given map outside 
       a given mask.
    """
    _unmask = []
    if mask_file is not None:
        mask = hp.read_map(mask_file)
        _unmask = np.where(mask > 1e-30)[0]
    else:
        mask_file = os.path.join(GRATOOLS_CONFIG, 'fits/Mask_src1p5_gp30.fits')
        mask = hp.read_map(mask_file)
        _unmask = np.where(mask > 1e-30)[0]
    COUNTS_BIN = re.compile('\_\d+\.')
    counts_files = glob.glob(os.path.join(GRATOOLS_OUT, 
                                       'output_counts/%s_counts_*.fits'%label))
    _counts = []
    _bin = []
    for ff in counts_files:
        m = re.search(COUNTS_BIN, ff)
        if m:
            _bin.append(int(m.group(0).replace('_', '').replace('.', '')))
            cmap = hp.read_map(ff)
            cmap_unmask = cmap*mask
            _counts.append(np.sum(cmap_unmask))
    _counts = np.array(_counts)
    _bin = np.array(_bin)
    _bin, index_bin = np.unique(_bin, return_index=True)
    _counts = _counts[index_bin]
    emin, emax, _energy = get_energy_from_txt(os.path.join(GRATOOLS_OUT, 
                                                           'ebinning.txt'), 
                                            get_binning=False, mean='log')
    fmt = dict(xname='$E$', xunits='MeV', yname='counts',
                   yunits='')
    spectrum = xInterpolatedUnivariateSplineLinear(_energy[
            np.amin(_bin):np.amax(_bin)+1], _counts, **fmt)
    return spectrum

def get_wbeam(wb_file, show=False):
    """Retrive the bivariate spline of the Wbeam function recorded in a txt file.
       
       wb_file : str 
          .txt file generated by build_wbeam
       show : bool
          if True the Wbeam matrix is plotted
    """
    _ebin, _l, _z = wbeam_parse(wb_file)
    fmt = dict(xname='$l$', xunits='', yname='Energy',
                   yunits='MeV', zname='W$_{beam}$(E,$l$)')
    wbeam = xInterpolatedBivariateSplineLinear( _l, _ebin, _z, **fmt)
    if show == True:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        cax = ax.matshow(_z.T, origin='lower', aspect='auto', cmap='Spectral')
        en_tick = ['1$\cdot$10$^{2}$', '5.8$\cdot$10$^{2}$', 
                   '3.4$\cdot$10$^{3}$','1.9$\cdot$10$^{4}$', 
                   '1.1$\cdot$10$^{5}$']
        ax.set_yticklabels(['']+en_tick)
        plt.title('$W_{beam}(l, E)$')
        plt.xlabel('$l$')
        plt.ylabel('$Energy$')
        plt.colorbar(cax)
        plt.grid()
        plt.show()
    
    return wbeam

def get_integral_wbeam(wb_file, spectrum_spline, e_min, e_max):
    """Returns the integra Wbeam as a function of the multiple in a certain
       energy interval.
       
       wb_file : str 
          .txt file generated by build_wbeam
       spectrum_spline : numpy spline 
          spline of the spectrum returned by get_powerlaw_spline
       e_min : float
          lower limit of the energy interval 
       e_max : float
          upper limit of the energy interval 
    """
    _ebin, _l, _z = wbeam_parse(wb_file)
    spectrum_spline = xInterpolatedUnivariateSplineLinear(_ebin, 
                                                       spectrum_spline(_ebin))
    int_spec = spectrum_spline.integral(e_min, e_max)
    fmt = dict(xname='$l$', xunits='', yname='Energy',
                   yunits='MeV', zname='W$_{beam}$(E,$l$)')
    int_wbeam = xInterpolatedBivariateSplineLinear( _l, _ebin,
                                                   _z*spectrum_spline(_ebin),
                                                   **fmt)
    int_wbeam_l = []
    for l in _l:
        l_slice = int_wbeam.vslice(l)
        int_wbeam_l.append(l_slice.integral(e_min, e_max))
    int_wbeam_l = np.array(int_wbeam_l)
    
    fmt2 = dict(xname='$l$', xunits='', 
                yname='integral W$_{beam}$ ($%.1f - %.1f$)'%(e_min, e_max),
                   yunits='')
    wbeam = xInterpolatedUnivariateSplineLinear(_l, int_wbeam_l/int_spec,\
                                                  **fmt2)
    norm = 1/wbeam.y[0]
    return wbeam

def get_psf_ref(psf_file):
    """discarded ...
       but it gets the published curve of the psf as a func of the energy
    """
    f = open(psf_file, 'r')
    _e, _ang = [], []
    for line in f:
        try:
            e, ang = [float(item) for item in line.split()]
            _e.append(e)
            _ang.append(ang)
        except:
            pass 
    fmt = dict(xname='Energy', xunits='MeV', yname='Containment Angle',
               yunits='deg')
    psf = xInterpolatedUnivariateSplineLinear(np.array(_e), np.array(_ang),\
                                                  **fmt)
    f.close()
    return psf

def get_psf(psf_file, show=False):
    """Get the PSF from the fits file created by gtpsf
    """
    #show = True
    hdu_list = pf.open(psf_file)
    _th = np.radians(hdu_list['THETA'].data.field('Theta'))
    _en = hdu_list['PSF'].data.field('ENERGY')
    _psf = hdu_list['PSF'].data.field('PSF')
    fmt = dict(yname='Theta', yunits='rad', xname='Energy',
               xunits='MeV', zname='PSF')
    psf_th_e = xInterpolatedBivariateSplineLinear(_en, _th, _psf, **fmt)
    hdu_list.close()
    if show == True:
        for e in _en:
            psf_th = psf_th_e.vslice(e)
            fmt = dict(xname='cos(th)', xunits='', yname='psf[cos(th)]',\
                           yunits='')
            _psfxsin = xInterpolatedUnivariateSplineLinear(psf_th.x, \
                                       psf_th.y*np.sin(psf_th.x), **fmt)
            plt.plot(np.degrees(psf_th.x), psf_th.y, label='%.2f'%e)
            print 'INT(E=%.2f) ='%e, 2*np.pi*_psfxsin.integral(0, \
                                                      np.amax(psf_th.x[:-2]))
        plt.yscale('log')
        plt.legend()
        plt.show()
    return psf_th_e

def wbeam_parse(wbeam_file, l_max =1500):
    """Created to parse the txt file given in output by build_wbeam.
    
       wbeam_file :
       l_max =1500 :

    """
    f = open(wbeam_file, 'r')
    _e = []
    _l = np.arange(l_max)
    _z = []
    for i, line in enumerate(f):
        if 'l' in line:
            _e.extend(float(item) for item in line.split()[1:])
            #print str(i), line.split()[0]
        elif str(i-1) == line.split()[0]:
            _z.append(np.array([float(item) for item in line.split()[1:]]))
    return np.array(_e), _l, np.array(_z)

def pix2DecRa(nside, num_pix):
    """Converts a pixel index to DEC and RA position in the sky, and returns 
       RA and DEC.

       nside : int
          nside of the healpix map you are considering
       num_pix : int or numpy array
          number of the i-th pixel, or an array of the number of pixels you 
          want to compute the position in ra and dec
    """
    theta, phi = hp.pixelfunc.pix2ang(nside, num_pix)
    ra = np.degrees(np.pi*2-phi)
    dec = -np.degrees(theta-np.pi/2.)
    return dec, ra


def main():
    """Test module
    """
    
    mask_f = os.path.join(GRATOOLS_CONFIG, 'fits/Mask_src1p5_gp30.fits')
    out_wbeam_txt = 'output/Wbeam_P8R2_ULTRACLEANVETO_V6_56.txt'
    out_wbeam_txt = 'output/Wbeam_P8R2_SOURCE_V6_32.txt'
    wb_2d = get_wbeam(out_wbeam_txt, show=True)

    
    logger.info('Wait...computing integral Wbeam at some energy intervals!')
    gamma = 2.3
    _emin = np.array([158.49,301.00,524.81])#,1000.00,1737.80, 2754.23,
                      #4786.30, 8317.64])a
    _emax =  np.array([301.00,524.81,1000.00])#, 1737.80, 2754.23, 4786.30, 
                       #8317.64, 14454.40])
    _emin = np.array([1000.00,1737.80, 2754.23, 4786.30, 8317.64,14454.40,
                      22908.68,39810.71,69183.09, 109647.81,301995.16])
    _emax =  np.array([1737.80, 2754.23, 4786.30, 8317.64,14454.40,22908.68,
                       39810.71,69183.09,120226.44,301995.16,1000000.0])
    spec = get_powerlaw_spline(gamma)
    c = ['0.', '0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','0.95']
    plt.figure(facecolor='white')
    for i, (emin, emax) in enumerate(zip(_emin, _emax)):
        wb = get_integral_wbeam(out_wbeam_txt, spec, e_min=emin, e_max=emax)
        wb.plot(show=False, label='%.2f-%.2f GeV'%(emin/1000, emax/1000), 
                color='%s'%c[i])
        plt.ylabel('W$_{beam}$')
        plt.title('ULTRACLEANVETO PSF 1+2+3')
        plt.ylim(0,1)
        plt.legend(loc=3, fontsize=10)
    plt.show()
   


if __name__ == '__main__':
    main()
