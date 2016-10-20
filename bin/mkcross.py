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
import ast
import argparse
import numpy as np
import healpy as hp
import pyfits as pf


__description__ = 'Makes the analysis'


"""Command-line switches.                                                       
"""


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
PARSER.add_argument('--show', type=ast.literal_eval, choices=[True, False],
                    default=False,
                    help='True if you want to see the maps')

def get_var_from_file(filename):
    f = open(filename)
    global data
    data = imp.load_source('data', '', f)
    f.close()

def mkCross(**kwargs):
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
        _l = np.arange(0, 1000)
        from GRATools.utils.gWindowFunc import build_wbeam
        wb = build_wbeam(psf, _l, out_wb_txt)
    else:
        from GRATools.utils.gWindowFunc import get_wbeam
        wb = get_wbeam(out_wb_txt)
    save_current_figure('Wbeam_%s.png'%out_wb_label, clear=True)

    logger.info('Starting Cl analysis...')
    in_label1 = data.IN_LABEL1
    in_label2 = data.IN_LABEL2
    out_label = data.OUT_LABEL
    binning_label = data.BINNING_LABEL
    mask_file = data.MASK_FILE
    mask = hp.read_map(mask_file)
    cl_param_file1 = os.path.join(GRATOOLS_OUT, '%s_%s_parameters.txt' \
                                     %(in_label1, binning_label))
    cl_param_file2 = os.path.join(GRATOOLS_OUT, '%s_%s_parameters.txt' \
                                      %(in_label2, binning_label))
    from GRATools.utils.gFTools import get_cl_param
    _emin, _emax, _emean, _f, _ferr, _cn, _fsky = get_cl_param(cl_param_file1)
    _emin2, _emax2, _emean2, _f2, _ferr2, _cn2, _fsky2 = get_cl_param(cl_param_file2)
    cross_txt = open(os.path.join(GRATOOLS_OUT, '%s_%s_cross.txt' \
                                   %(out_label, binning_label)), 'w')
    for i, (emin, emax) in enumerate(zip(_emin, _emax)):
        logger.info('Considering bin %.2f - %.2f ...'%(emin, emax))
        gamma = data.WEIGHT_SPEC_INDEX
        Im = (1/(1-gamma))*(emax**(1-gamma)-emin**(1-gamma))/(emax-emin)
        eweightedmean = np.power(1/Im, 1/gamma)
        cross_txt.write('ENERGY\t %.2f %.2f %.2f\n'%(emin, emax, eweightedmean))
        l_max= 1000
        _l = np.arange(l_max)
        wb_en = wb.hslice(eweightedmean)(_l)
        flux_map_name1 = in_label1+'_flux_%i-%i.fits'%(emin, emax)
        flux_map_name2 = in_label2+'_flux_%i-%i.fits'%(emin, emax)
        flux_map1 = hp.read_map(os.path.join(GRATOOLS_OUT_FLUX, flux_map_name1))
        flux_map2 = hp.read_map(os.path.join(GRATOOLS_OUT_FLUX, flux_map_name2))
        flux_map_masked1 = hp.ma(flux_map1)
        flux_map_masked1.mask = np.logical_not(mask)
        flux_map_masked2 = hp.ma(flux_map2)
        flux_map_masked2.mask = np.logical_not(mask)
        fsky1 = 1.-(len(np.where(flux_map_masked1.filled() == hp.UNSEEN)[0])/\
                       float(len(flux_map1)))
        fsky2 = 1.-(len(np.where(flux_map_masked2.filled() == hp.UNSEEN)[0])/\
                       float(len(flux_map2)))
        if kwargs['show'] == True:
            hp.mollview(flux_map_masked1.filled(), title='f$_{sky}$ = %.3f'%fsky,
                        min=1e-7, max=1e-4, norm='log')
            plt.show()
        print 'fsky = ', fsky1, fsky2
        fsky = np.sqrt(fsky1*fsky2)
        nside1 = hp.npix2nside(len(flux_map1))
        nside2 = hp.npix2nside(len(flux_map1))
        wpix1 = hp.sphtfunc.pixwin(nside1)[:l_max]
        wpix2 = hp.sphtfunc.pixwin(nside2)[:l_max]
        _cl_cross = hp.sphtfunc.anafast(flux_map_masked1.filled(),
                                        flux_map_masked2.filled(), lmax=l_max-1, \
                                             iter=5)
        wl2 = wb_en*wb_en*wpix1*wpix2
        _cl_cross = (_cl_cross/fsky)/(wl2)
        cross_txt.write('Cl\t%s\n'%str(list(_cl_cross)).replace('[',''). \
                         replace(']','').replace(', ', ' '))
        _cl_cross_err = np.sqrt(2./((2*_l+1)*fsky))*(_cl_cross)
        cross_txt.write('Cl_ERR\t%s\n\n'%str(list(_cl_cross_err)).replace('[',''). \
                         replace(']','').replace(', ', ' '))
    cross_txt.close()
    logger.info('Created %s'%(os.path.join(GRATOOLS_OUT, '%s_%s_cross.txt' \
                                               %(out_label, binning_label))))



if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    mkCross(**args.__dict__)
