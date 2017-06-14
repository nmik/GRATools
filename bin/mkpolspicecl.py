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


__description__ = 'Makes the analysis'


"""Command-line switches.                                                       
"""


from GRATools import GRATOOLS_OUT, GRATOOLS_CONFIG
from GRATools.utils.gPolSpice import *
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

def mkCl(**kwargs):
    """                                      
    """
    get_var_from_file(kwargs['config'])
    
    lmax = data.LMAX
    _l = np.arange(0, lmax)
    gamma = data.WEIGHT_SPEC_INDEX
    out_label = data.OUT_LABEL
    logger.info('Calculating energy spectrum...')
    from GRATools.utils.gWindowFunc import get_powerlaw_spline
    energy_spec = get_powerlaw_spline(gamma)

    logger.info('Calculating PSF with gtpsf...')
    dict_gtpsf = data.DICT_GTPSF
    logger.info('Calculating Wbeam Function...')
    out_wb_label = data.OUT_W_LABEL 
    out_wb_txt = os.path.join(GRATOOLS_OUT, 'Wbeam_%s.txt'%out_wb_label)
    in_label = data.IN_LABEL
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
        pass

    logger.info('Starting Cl analysis...')
    mask_label = data.MASK_LABEL
    in_label = in_label + '_' + mask_label
    binning_label = data.BINNING_LABEL
    mask_file = data.MASK_FILE
    cl_param_file = os.path.join(GRATOOLS_OUT, '%s_%s_parameters.txt' \
                                     %(in_label, binning_label))
    from GRATools.utils.gFTools import get_cl_param
    _emin, _emax, _emean, _f, _ferr, _cn, _fsky = get_cl_param(cl_param_file)
    cl_txt = open(os.path.join(GRATOOLS_OUT, '%s_%s_polspicecls.txt' \
                                   %(out_label, binning_label)), 'w')
    from GRATools.utils.gWindowFunc import get_integral_wbeam
    for i, (emin, emax) in enumerate(zip(_emin, _emax)):
        logger.info('Considering bin %.2f - %.2f ...'%(emin, emax))
        mask_f = mask_file
        if type(mask_file) == list:
            mask_f = mask_file[i]
        mask = hp.read_map(mask_f)
        wb1_en = get_integral_wbeam(out_wb_txt1, energy_spec, emin, emax)
        cl_txt.write('ENERGY\t %.2f %.2f %.2f\n'%(emin, emax, _emean[i]))
        l_max= lmax
        _l = np.arange(l_max)
        wb_en = wb_en(_l)
        flux_map_name = in_label+'_flux_%i-%i.fits'%(emin, emax)
        flux_map_f = os.path.join(GRATOOLS_OUT_FLUX, flux_map_name)
        flux_map_f_mdclean = remove_monopole_dipole(flux_map_f)
        flux_map = hp.read_map(flux_map_f_mdclean)
        fsky = _fsky[i]
        cn = _cn[i]
        if kwargs['show'] == True:
            hp.mollview(flux_map_masked.filled(), title='f$_{sky}$ = %.3f'%fsky,
                        min=1e-7, max=1e-4, norm='log')
            plt.show()
        logger.info('fsky = '%fsky)
        nside = hp.npix2nside(len(flux_map))
        wpix = hp.sphtfunc.pixwin(nside)[:l_max]
        out_name = '%s_%i-%i' %(out_label, emin, emax)
        out_folder =  os.path.join(GRATOOLS_OUT, 'output_pol')
        logger.info('cn poisson = %e'%cn)
        wl = wb_en*wpix
        if not os.path.exists(out_folder):
            os.makedirs(out_folder)
        pol_dict = data.POLCEPICE_DICT
        for key in pol_dict:
            if key == 'clfile':
                pol_dict[key] = os.path.join(out_folder,'%s_cl.txt'%out_name)
            if key == 'cl_outmap_file':
                pol_dict[key] = os.path.join(out_folder,'%s_clraw.txt'%out_name)
            if key == 'covfileout':
                 pol_dict[key] = os.path.join(out_folder,'%s_cov.fits'%out_name)
            if key == 'mapfile':
                pol_dict[key] = flux_map_f
            if key == 'maskfile':
                pol_dict[key] = mask_f
        config_file_name = 'pol_%s'%(out_name)
        if os.path.exists(os.path.join(out_folder,'%s_cl.txt'%out_name)) and \
                os.path.exists(os.path.join(out_folder,'%s_cov.fits'%out_name)):
            logger.info('ATT: Retriving power spectrum...')
            _l, _cl, _clerr = pol_cl_parse(os.path.join(out_folder,
                                                        '%s_cl.txt'%out_name),
                                           os.path.join(out_folder,
                                                        '%s_cov.fits'%out_name),
                                           raw_corr=(cn, wl),
                                           rebin=True)
            logger.info('... and covariance matrix.')
            _cov = pol_cov_parse(os.path.join(out_folder,
                                              '%s_cov.fits'%out_name),
                                 wl_array=wl,
                                 rebin=True, show=True)
        else:
            _l, _cl, _cl_err = pol_cl_calculation(pol_dict, 
                                                  config_file_name,
                                                  raw_corr=(cn, wl),
                                                  rebin=True,show=True) 
        cl_txt.write('multipole\t%s\n'%str(list(_l)).replace('[',''). \
                         replace(']','').replace(', ', ' '))
        cl_txt.write('Cl\t%s\n'%str(list(_cl)).replace('[',''). \
                         replace(']','').replace(', ', ' '))
        cl_txt.write('Cl_ERR\t%s\n\n'%str(list(_clerr)).replace('[',''). \
                         replace(']','').replace(', ', ' '))
    cl_txt.close()
    logger.info('Created %s'%(os.path.join(GRATOOLS_OUT,'%s_%s_polspicecls.txt'
                                           %(out_label, binning_label))))



if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    mkCl(**args.__dict__)
