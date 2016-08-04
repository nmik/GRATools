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


"""Analysis module                                                              
"""


import os
import imp
import numpy as np
import healpy as hp
import pyfits as pf


__description__ = 'Makes the analysis'



"""Command-line switches.                                                       
"""

import argparse
from GRATools import GRATOOLS_OUT
from GRATools.utils.logging_ import logger, startmsg
from GRATools.utils.matplotlib_ import pyplot as plt
from GRATools.utils.matplotlib_ import overlay_tag, save_current_figure

GRATOOLS_OUT_FLUX = os.path.join(GRATOOLS_OUT, 'output_flux')

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('--config', type=str, required=True,
                    help='the input configuration file')

def get_var_from_file(filename):
    f = open(filename)
    global data
    data = imp.load_source('data', '', f)
    f.close()

def mkCl(**kwargs):
    """                                      
    """
    logger.info('Starting Cl analysis...')
    get_var_from_file(kwargs['config'])
    in_label = data.IN_LABEL
    out_label = data.OUT_LABEL
    binning_label = data.BINNING_LABEL
    wbeam_file = data.WBEAM_FILE
    psf_ref_file = data.PSF_REF_FILE 
    from GRATools.utils.gWindowFunc import get_wbeam
    from GRATools.utils.gWindowFunc import get_psf_ref
    wb = get_wbeam(wbeam_file)
    psf = get_psf_ref(psf_ref_file)
    #psf.plot(show=False, logx=True)
    cl_param_file = os.path.join(GRATOOLS_OUT, '%s_%s_parameters.txt' \
                                     %(in_label, binning_label))
    from GRATools.utils.gFTools import get_cl_param
    _emin, _emax, _emean, _f, _ferr, _cn, _fsky = get_cl_param(cl_param_file)
    """Here can be moved the part to draw the flux (now in bin/mkanalysis.py)
    """
    _cls = []
    _cls_err = []
    #plt.figure(figsize=(10, 7), dpi=80)
    for i, (emin, emax) in enumerate(zip(_emin, _emax)):
        psf_en = psf(_emean[i])
        logger.info('Enegry = %i -> PSF = %.2f'%(_emean[i], psf_en))
        l_max = np.pi/np.radians(psf_en)+30
        _l = np.arange(l_max)
        #l_max = 2000
        logger.info('Truncating l at %i' %l_max)
        logger.info('Considering bin %.2f - %.2f ...'%(emin, emax))
        wb_en = wb.hslice(_emean[i])(_l)
        #plt.figure(figsize=(10, 7), dpi=80) 
        #plt.plot(_l, wb_en**2, label='W$_{beam}$$^{2}$')
        flux_map_name = in_label+'_flux_%i-%i.fits'%(emin, emax)
        flux_map = hp.read_map(os.path.join(GRATOOLS_OUT_FLUX, flux_map_name))
        nside = hp.npix2nside(len(flux_map))
        wpix = hp.sphtfunc.pixwin(nside)[:l_max+1]
        #plt.plot(_l, wpix[:2001]**2, label='W$_{pix}$$^{2}$')
        #plt.plot(_l, (wb_en*wpix[:2001])**2, \
        #   label='(W$_{beam}\cdot$W$_{pix}$)$^{2}$')
        #plt.legend()
        #plt.xscale('log')
        #plt.yscale('log')
        #plt.xlabel('$l$')
        #overlay_tag(y=0.05)
        #save_current_figure('Window_funcs_example.png',clear=False)
        _cl = hp.sphtfunc.anafast(flux_map, lmax=l_max, alm=False, pol=False)
        plt.figure(figsize=(10, 7), dpi=80) 
        wl = wb_en*wpix
        _cl_err = _cl/10
        #cn_fit = np.average(_cl[1900:1990])
        print 'cn paper', _cn[i]
        #print 'cn fit', cn_fit
        plt.errorbar(_l[:l_max+1], _cl/_fsky[i]-_cn[i], fmt='o', markersize=3, \
                     elinewidth=1, xerr=np.array([0.5]*len(_l[:l_max+1])),\
                     yerr=_cl_err, color='red')
        #print (_cl/_fsky[i])[:100]
        _cl = (_cl/_fsky[i] - _cn[i])/(wl**2)
        #_cl = (_cl/_fsky[i] - _cn[i])/(wl**2)
        #_cls.append(_cl)
        _cl_err = np.sqrt(2./((2*_l+1)*_fsky[i]))*(_cl+(_cn[i]/wl**2))
        #_cls_err.append(_cl_err)
        print wl[:100]**2
        plt.errorbar(_l, _cl, fmt='o', markersize=3, \
                     elinewidth=1, xerr=np.array([0.5]*len(_l[:l_max+1])),\
                     yerr=_cl_err, color='green')
        plt.xlabel('$l$')
        plt.ylabel('$C_{l}$')
        plt.ylim(-1e-16, 1e-16)
        #plt.ylim(-1e-14, 1e-14)
        plt.xlim(50, l_max+1)
        #plt.xscale('log')
        #plt.yscale('log')
        plt.show()

if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    mkCl(**args.__dict__)
