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
import scipy.special as sp
from GRATools.utils.matplotlib_ import pyplot as plt
from GRATools.utils.matplotlib_ import overlay_tag, save_current_figure
from GRATools.utils.gSpline import xInterpolatedBivariateSplineLinear
from GRATools.utils.gSpline import xInterpolatedUnivariateSplineLinear


def Pl(l, _x):
    """return an array of Pl values correspoding to a given x array
    """
    _pl = sp.legendre(l)(_x)
    return _pl

def get_pl_vs_th(l, _th):
    """
    """
    _pl = Pl(l, _th)
    fmt = dict(xname='cos(th)', xunits='', yname='Pl(cos(th))',
               yunits='')
    pl_th = xInterpolatedUnivariateSplineLinear(_th, _pl, **fmt)
    return pl_th

def get_wbeam(wb_file):
    f = open(wb_file)
    _ebin, _l, _e1, _e2, _e3, _e4, _e5, _e6, _e7, _e8 = wbeam_parse(wb_file)
    _z = np.array([_e2, _e3, _e4, _e5, _e6, _e7, _e8])
    fmt = dict(xname='$l$', xunits='', yname='Energy',
                   yunits='MeV', zname='W$_{beam}$(E,$l$)')
    wbeam = xInterpolatedBivariateSplineLinear(_l, _ebin[1:], _z.T, **fmt)
    return wbeam

def get_psf_ref(psf_file):
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
    return psf

def king_function(_x, _s, _g):
    """                                
    """
    K = (1./(2*np.pi*_s))*(1. - 1./_g)*(1. + (1/(2*_g))*(_x**2/_s**2))**(-_g)
    return K

def get_psf_sinth(psf_file):
    """
    """
    hdu_list = pf.open(psf_file)
    _th = np.radians(hdu_list['THETA'].data.field('Theta'))
    _en = hdu_list['PSF'].data.field('ENERGY')
    _psf = hdu_list['PSF'].data.field('PSF')*np.sin(_th)
    fmt = dict(yname='Theta', yunits='rad', xname='Energy',
               xunits='MeV', zname='PSF')
    psf_th_e = xInterpolatedBivariateSplineLinear(_en, _th, _psf, **fmt)
    return psf_th_e

def wbeam_parse(wbeam_file):
    f = open(wbeam_file, 'r')
    _e, _ebin, _l = [], [], []
    _e1, _e2, _e3, _e4, _e5, _e6, _e7, _e8 = [], [], [], [], [], [], [], []
    for line in f:
        if 'k' in line:
            _e.extend(float(item) for item in line.split()[1:])
        try:
            l, e1, e2, e3, e4, e5, e6, e7, e8 = [float(item) for item \
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
        except:
            pass
    _e = np.unique(np.array(_e))
    for emin, emax in zip(_e[:-1], _e[1:]):
        _ebin.append((emax-emin)/2)
    _ebin = np.array(_ebin)
    _l = np.array(_l)
    _e1 = np.array(_e1)
    _e2 = np.array(_e2)
    _e3 = np.array(_e3)
    _e4 = np.array(_e4)
    _e5 = np.array(_e5)
    _e6 = np.array(_e6)
    _e7 = np.array(_e7)
    _e8 = np.array(_e8)
    return _ebin, _l, _e1, _e2, _e3, _e4, _e5, _e6, _e7, _e8


def main():
    """Test module
    """
    psf_file = 'prova.fits'
    psf = get_psf_sinth(psf_file)
    #psf.plot(logz=True)
    #save_current_figure('prova.png', clear=False)
    psf_1GeV = psf.vslice(1000)
    print psf_1GeV.integral(0, 30)*2*np.pi
    psf_1GeV.plot()
    wb_file = 'config/ascii/Wbeam_p8_clean_v6.txt'
    wb = get_wbeam(wb_file)
    plt.figure(figsize=(10, 7), dpi=80)
    wb.plot(show=False)
    overlay_tag()
    save_current_figure('Wbeam_p8_clean_v6.png', clear=False)
    plt.figure(figsize=(10, 7), dpi=80)
    wb_1GeV = wb.hslice(1000)
    wb_1GeV.plot(show=False, logy=True)
    overlay_tag()
    save_current_figure('Wbeam_p8_clean_v6_at1GeV.png', clear=False)

if __name__ == '__main__':
    main()
