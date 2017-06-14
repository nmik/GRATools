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

def mkCross(**kwargs):
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

    logger.info('Calculating PSF with gtpsf 1...')
    dict_gtpsf1 = data.DICT_GTPSF1
    logger.info('Calculating Wbeam Function1...')
    out_wb_label1 = data.OUT_W_LABEL1
    out_wb_txt1 = os.path.join(GRATOOLS_OUT, 'Wbeam_%s.txt'%out_wb_label1)
    in_label1 = data.IN_LABEL1
    if not os.path.exists(out_wb_txt1):
        from GRATools.utils.ScienceTools_ import gtpsf
        gtpsf(dict_gtpsf1)
        from GRATools.utils.gWindowFunc import get_psf
        psf_file1 = data.PSF_FILE1
        psf1 = get_psf(psf_file1)
        from GRATools.utils.gWindowFunc import build_wbeam
        wb1 = build_wbeam(psf1, _l, out_wb_txt1)
    else:
        pass
    logger.info('Calculating PSF with gtpsf 2...')
    dict_gtpsf2 = data.DICT_GTPSF2
    logger.info('Calculating Wbeam Function2...')
    out_wb_label2 = data.OUT_W_LABEL2
    out_wb_txt2 = os.path.join(GRATOOLS_OUT, 'Wbeam_%s.txt'%out_wb_label2)
    in_label2 = data.IN_LABEL2
    if not os.path.exists(out_wb_txt2):
        from GRATools.utils.ScienceTools_ import gtpsf
        gtpsf(dict_gtpsf2)
        from GRATools.utils.gWindowFunc import get_psf
        psf_file2 = data.PSF_FILE2
        psf2 = get_psf(psf_file2)
        from GRATools.utils.gWindowFunc import build_wbeam
        wb2 = build_wbeam(psf2, _l, out_wb_txt2)
    else:
        pass
    
    logger.info('Starting Cross analysis...')
    mask_label1 = data.MASK_LABEL1
    mask_label2 = data.MASK_LABEL2
    in_label1 = in_label1 + '_' + mask_label1
    in_label2 = in_label2 + '_' + mask_label2
    binning_label1 = data.BINNING_LABEL1
    binning_label2 = data.BINNING_LABEL2
    mask_file1 = data.MASK_FILE1
    mask_file2 = data.MASK_FILE2
    cl_param_file1 = os.path.join(GRATOOLS_OUT, '%s_%s_parameters.txt' \
                                     %(in_label1, binning_label1))
    cl_param_file2 = os.path.join(GRATOOLS_OUT, '%s_%s_parameters.txt' \
                                     %(in_label2, binning_label2))
    from GRATools.utils.gFTools import get_cl_param
    _emin, _emax, _emean, _f1, _ferr1, _cn1, _fsky1 = \
        get_cl_param(cl_param_file1)
    _emin, _emax, _emean, _f2, _ferr2, _cn2, _fsky2 = \
        get_cl_param(cl_param_file2)
    cl_txt = open(os.path.join(GRATOOLS_OUT, '%s_%s_polspicecross.txt' \
                                   %(out_label, binning_label1)), 'w')
    from GRATools.utils.gWindowFunc import get_integral_wbeam
    for i, (emin, emax) in enumerate(zip(_emin, _emax)):
        logger.info('Considering bin %.2f - %.2f ...'%(emin, emax))
        mask_f1 = mask_file1
        mask_f2 = mask_file2
        if type(mask_file1) == list:
            mask_f1 = mask_file1[i]
        if type(mask_file2) == list:
            mask_f2 = mask_file2[i]
        mask1 = hp.read_map(mask_f1)
        mask2 = hp.read_map(mask_f2)
        wb1_en = get_integral_wbeam(out_wb_txt1, energy_spec, emin, emax)
        wb2_en = get_integral_wbeam(out_wb_txt2, energy_spec, emin, emax) 
        cl_txt.write('ENERGY\t %.2f %.2f %.2f\n'%(emin, emax, _emean[i]))
        l_max = lmax
        _l = np.arange(l_max)
        wb1_en = wb1_en(_l)
        wb2_en = wb2_en(_l)
        flux_map_name1 = in_label1+'_flux_%i-%i.fits'%(emin, emax)
        flux_map_name2 = in_label2+'_flux_%i-%i.fits'%(emin, emax)
        flux_map_f1 = os.path.join(GRATOOLS_OUT_FLUX, flux_map_name1)
        flux_map_f2 = os.path.join(GRATOOLS_OUT_FLUX, flux_map_name2)
        flux_map_f_mdclean1 = remove_monopole_dipole(flux_map_f1)
        flux_map_f_mdclean2 = remove_monopole_dipole(flux_map_f2)
        flux_map1 = hp.read_map(flux_map_f_mdclean1)
        flux_map2 = hp.read_map(flux_map_f_mdclean2)
        fsky1 = _fsky1[i]
        fsky2 = _fsky2[i]
        cn = 0.
        if kwargs['show'] == True:
            hp.mollview(flux_map1*mask1, title='f$_{sky}$ = %.3f'%fsky,
                        min=1e-7, max=1e-4, norm='log')
            hp.mollview(flux_map2*mask2, title='f$_{sky}$ = %.3f'%fsky,
                        min=1e-7, max=1e-4, norm='log')
            plt.show()
        logger.info('fsky1 = %.3f'%fsky1)
        logger.info('fsky2 = %.3f'%fsky2)
        nside1 = hp.npix2nside(len(flux_map1))
        nside2 = hp.npix2nside(len(flux_map2))
        wpix1 = hp.sphtfunc.pixwin(nside1)[:l_max]
        wpix2 = hp.sphtfunc.pixwin(nside2)[:l_max]
        out_name = '%s_%i-%i' %(out_label, emin, emax)
        out_folder =  os.path.join(GRATOOLS_OUT, 'output_pol')
        logger.info('cn poisson = %e'%cn)
        wl = np.sqrt(wb1_en*wpix1*wb2_en*wpix2)
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
                pol_dict[key] = flux_map_f_mdclean1
            if key == 'maskfile':
                pol_dict[key] = mask_f1
            if key == 'mapfile2':
                pol_dict[key] = flux_map_f_mdclean2
            if key == 'maskfile2':
                pol_dict[key] = mask_f2
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
            _l, _cl, _clerr, _cov = pol_cl_calculation(pol_dict, 
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
                                           %(out_label, binning_label1))))



if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    mkCross(**args.__dict__)
