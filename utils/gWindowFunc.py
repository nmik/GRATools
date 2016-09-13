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
import numpy as np
import pyfits as pf
import healpy as hp
import scipy.special as sp
from GRATools.utils.logging_ import logger
from GRATools.utils.matplotlib_ import pyplot as plt
from GRATools.utils.matplotlib_ import overlay_tag, save_current_figure
from GRATools.utils.gSpline import xInterpolatedBivariateSplineLinear
from GRATools.utils.gSpline import xInterpolatedUnivariateSplineLinear


def get_pl_vs_th(l, _costh):
    """return an array of Pl values correspoding to a given _costh array
    """
    _pl_th = sp.legendre(l)(_costh)
    fmt = dict(xname='th', xunits='rad', yname='Pl(cos(th))',
               yunits='')
    return _pl_th

def build_wbeam(psf, _l, out_file):
    """Calculate the Wbeam(l, E) and return a bivariate slpine.
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

def get_wbeam(wb_file):
    """Retrive the bivariate spline of the Wbeam function if a 
       file has been created
    """
    _ebin, _l, _e1, _e2, _e3, _e4, _e5, _e6, _e7, _e8, _e9, _e10 = \
        wbeam_parse(wb_file)
    _z = np.array([_e1, _e2, _e3, _e4, _e5, _e6, _e7, _e8, _e9, _e10])
    fmt = dict(xname='$l$', xunits='', yname='Energy',
                   yunits='MeV', zname='W$_{beam}$(E,$l$)')
    wbeam = xInterpolatedBivariateSplineLinear(_l, _ebin, _z.T, **fmt)
    return wbeam

def get_psf_ref(psf_file):
    """get the published curve of the psf as a func of the energy
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

def wbeam_parse(wbeam_file):
    """Created to parse the txt file given in output by build_wbeam
       ATT: it works only with 10 energy bins, if a different nuber of bins 
       is present in the file, this function have to be modified!!
    """
    f = open(wbeam_file, 'r')
    _e, _l = [], []
    _e1, _e2, _e3, _e4, _e5 = [], [], [], [], []
    _e6, _e7, _e8, _e9, _e10 = [], [], [], [], []
    for line in f:
        if 'l' in line:
            _e.extend(float(item) for item in line.split()[1:])
            
        try:
            l, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10 = [float(item) for item \
                                                     in line.split()]
            _l.append(l)
            _e1.append(e1)
            _e2.append(e2)
            _e3.append(e3)
            _e4.append(e4)
            _e5.append(e5)
            _e6.append(e6)
            _e7.append(e7)
            _e8.append(e8)
            _e9.append(e9)
            _e10.append(e10)
        except:
            pass
    _e = np.array(_e)
    _l = np.array(_l)
    _e1 = np.array(_e1)
    _e2 = np.array(_e2)
    _e3 = np.array(_e3)
    _e4 = np.array(_e4)
    _e5 = np.array(_e5)
    _e6 = np.array(_e6)
    _e7 = np.array(_e7)
    _e8 = np.array(_e8)
    _e9 = np.array(_e9)
    _e10 = np.array(_e10)
    f.close()
    return _e, _l, _e1, _e2, _e3, _e4, _e5, _e6, _e7, _e8, _e9, _e10

def IndexToDeclRa(NSIDE, index):
    """Converts a pixel index to DEC and RA position in the sky
    """
    theta, phi = hp.pixelfunc.pix2ang(NSIDE,index)
    ra = np.degrees(np.pi*2-phi)
    dec = -np.degrees(theta-np.pi/2.)
    return dec, ra


def main():
    """Test module
    """
    """
    plt.figure(figsize=(10, 7), dpi=80)
    _l = np.arange(0, 10, 2)
    for l in _l:
        pl_th = get_pl_vs_th(l, np.arange(-1, 1, 0.00001))
        plt.plot(np.arange(-1, 1, 0.00001), pl_th, '.', label='l = %i'%l)
    plt.legend()
    #plt.show()
    """
    NSIDE = 1
    NPIX = hp.nside2npix(NSIDE)
    iii = np.arange(NPIX)
    dec, ra = IndexToDeclRa(NSIDE, iii)
    index = np.where(abs(dec)>10)
    ra = ra[index]
    dec = dec[index]
    #ra = ra
    #dec = dec
    plt.figure(figsize=(10, 7), dpi=80)
    hp.mollview(iii, title="Mollview image RING")
    plt.figure(figsize=(10, 7), dpi=80)
    lab, plots = [], []
    for i in range(0, len(ra)):
        print ra[i], dec[i]
        out_wbeam_txt = 'output/%i_prova.txt'%i
        psf_file = 'output/%i_prova.fits'%i
        dict_gtpsf = {'expcube':'/data1/data/FT-files/output/output_gtltcube'+\
                          '/Allyrs_filtered_gti_ltcube.fits', 
                      'outfile': psf_file, 
                      'irfs': 'P8R2_ULTRACLEANVETO_V6',
                      'evtype': 56,
                      'ra': ra[i], 
                      'dec': dec[i], 
                      'emin': 500, 
                      'emax': 600000,
                      'nenergies': 10, 
                      'thetamax': 30, 
                      'ntheta': 300}
        from GRATools.utils.ScienceTools_ import gtpsf
        gtpsf(dict_gtpsf)
        _l = np.arange(0, 1000, 4)
        psf = get_psf(psf_file)
        if not os.path.exists(out_wbeam_txt):
            wb = build_wbeam(psf, _l, out_wbeam_txt)
        else:
            wb = get_wbeam(out_wbeam_txt)
        wb_1GeV = wb.hslice(1000)
        wl = wb_1GeV.plot(show=False, label='RA:%i, Dec:%i'%(ra[i],dec[i]))
        #lab.append('%i-%i'%(ra[i],dec[i]))
        #plots.append(wl)
        #wb_50GeV = wb.hslice(50000)
        #wb_50GeV.plot(show=False)
        #wb_100000GeV = wb.hslice(100000)
        #wb_100000GeV.plot(show=False)
        #wb_500000GeV = wb.hslice(500000)
        #wb_500000GeV.plot(show=False)
    plt.legend()
    plt.show()


if __name__ == '__main__':
    main()
