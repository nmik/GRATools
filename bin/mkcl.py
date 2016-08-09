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
from GRATools import GRATOOLS_OUT, GRATOOLS_CONFIG
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
    get_var_from_file(kwargs['config'])

    logger.info('Calculating PSF with gtpsf...')
    dict_gtpsf = data.DICT_GTPSF
    logger.info('Calculating Wbeam Function...')
    out_wb_label = data.OUT_W_LABEL
    out_wb_txt = os.path.join(GRATOOLS_OUT, 'Wbeam_%s.txt'%out_wb_label)
    if not os.path.exists(out_wb_txt):
        from GRATools.utils.ScienceTools_ import gtpsf
        gtpsf(dict_gtpsf)
        from GRATools.utils.gWindowFunc import get_psf
        psf_file = data.PSF_FILE
        psf = get_psf(psf_file)
        _l = np.arange(0, 1000, 4)
        from GRATools.utils.gWindowFunc import build_wbeam
        wb = build_wbeam(psf, _l, out_wb_txt)
    else:
        from GRATools.utils.gWindowFunc import get_wbeam
        wb = get_wbeam(out_wb_txt)
    wb.plot(show=False)
    save_current_figure('Wbeam_%s.png'%out_wb_label, clear=True)

    logger.info('Starting Cl analysis...')
    in_label = data.IN_LABEL
    out_label = data.OUT_LABEL
    binning_label = data.BINNING_LABEL
    cl_param_file = os.path.join(GRATOOLS_OUT, '%s_%s_parameters.txt' \
                                     %(in_label, binning_label))
    from GRATools.utils.gFTools import get_cl_param
    _emin, _emax, _emean, _f, _ferr, _cn, _fsky = get_cl_param(cl_param_file)
    cl_txt = open(os.path.join(GRATOOLS_OUT, '%s_%s_cls.txt' \
                                   %(out_label, binning_label)), 'w')
    for i, (emin, emax) in enumerate(zip(_emin, _emax)):
        logger.info('Considering bin %.2f - %.2f ...'%(emin, emax))
        gamma = data.WEIGHT_SPEC_INDEX
        Im = (1/(1-gamma))*(emax**(1-gamma)-emin**(1-gamma))/(emax-emin)
        eweightedmean = np.power(1/Im, 1/gamma)
        cl_txt.write('ENERGY\t %.2f %.2f %.2f\n'%(emin, emax, eweightedmean))
        l_max= 1000
        _l = np.arange(l_max)
        wb_en = wb.hslice(eweightedmean)(_l)
        flux_map_name = in_label+'_flux_%i-%i.fits'%(emin, emax)
        flux_map = hp.read_map(os.path.join(GRATOOLS_OUT_FLUX, flux_map_name))
        nside = hp.npix2nside(len(flux_map))
        wpix = hp.sphtfunc.pixwin(nside)[:l_max]
        cn_fit = np.average(hp.sphtfunc.anafast(flux_map)[-100:-1])/99
        cn = _cn[i] #cn_fit/_fsky[i] #(_cn[i] + cn_fit) / 2
        _cl = hp.sphtfunc.anafast(flux_map, lmax=l_max-1, alm=False, pol=False)
        wl = wb_en*wpix
        _cl = (_cl/_fsky[i] - cn)/(wl**2)
        cl_txt.write('Cl\t%s\n'%str(list(_cl)).replace('[',''). \
                         replace(']','').replace(', ', ' '))
        _cl_err = np.sqrt(2./((2*_l+1)*_fsky[i]))*(_cl+(cn/wl**2))
        cl_txt.write('Cl_ERR\t%s\n\n'%str(list(_cl_err)).replace('[',''). \
                         replace(']','').replace(', ', ' '))
    cl_txt.close()
    logger.info('Created %s'%(os.path.join(GRATOOLS_OUT, '%s_%s_cls.txt' \
                                               %(out_label, binning_label))))



if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    mkCl(**args.__dict__)
