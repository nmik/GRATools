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
from itertools import combinations

__description__ = 'Makes the analysis'


"""Command-line switches.                                                       
"""


from GRATools import GRATOOLS_OUT, GRATOOLS_CONFIG
from GRATools.utils.gPolSpice import *
from GRATools.utils.logging_ import logger, startmsg
from GRATools.utils.matplotlib_ import pyplot as plt
from GRATools.utils.matplotlib_ import overlay_tag, save_current_figure
from GRATools.utils.gWindowFunc import get_integral_wbeam


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
    logger.info('Calculating energy spectrum...')
    from GRATools.utils.gWindowFunc import get_powerlaw_spline
    energy_spec = get_powerlaw_spline(gamma)

    wbeamo_list = data.WBEAM_LIST
    mask_list = data.MASK_LIST
    maps_list = data.MAPS_LIST
    ebin_list = data.EBINS
    index = np.arange(len(maps_list))
    index_comb = list(combinations(index, r=2))
    maps_comb = list(combinations(maps_list, r=2))
    mask_comb =  list(combinations(mask_list, r=2))
    wbeam_comb =  list(combinations(wbeamo_list, r=2))
    ebin_comb = list(combinations(ebin_list, r=2))
    outfile_label = data.OUTFILE_LABEL
    cl_txt = open(os.path.join(GRATOOLS_OUT, '%s_polspiceEcross.txt' \
                                   %outfile_label), 'w')
    for i, (map1, map2) in enumerate(maps_comb):
        logger.info('Cross correlation:')
        logger.info('>>> Map %i: %s' %(index_comb[i][0], map1))
        logger.info('>>> Map %i: %s' %(index_comb[i][1], map2))
        emin1, emax1 = ebin_comb[i][0][0], ebin_comb[i][0][1]
        emin2, emax2 = ebin_comb[i][1][0], ebin_comb[i][1][1]
        mask1_f = mask_comb[i][0]
        mask2_f = mask_comb[i][1]
        mask_f = ''
        mask1 = hp.read_map(mask1_f)
        mask2 = hp.read_map(mask2_f)
        nside = hp.npix2nside(len(mask1))
        if np.sum(mask1) >  np.sum(mask1):
            logger.info('taking Mask %i'%index_comb[i][1])
            mask_f = mask2_f
        else:
            logger.info('taking Mask %i'%index_comb[i][0])
            mask_f = mask1_f
        out_label = '%i-%i_%i-%i'%(emin1, emax1, emin2, emax2)
        cl_txt.write('ENERGY\t %.2f %.2f - %.2f %.2f\n'
                     %(emin1, emax1, emin2, emax2))
        wb_en1 = get_integral_wbeam(wbeam_comb[i][0], energy_spec, emin1, emax1)
        wb_en2= get_integral_wbeam(wbeam_comb[i][1], energy_spec, emin2, emax2)
        wb_en1 = wb_en1(_l)
        wb_en2 = wb_en2(_l)
        wpix1 = wpix2 = hp.sphtfunc.pixwin(nside)[:lmax]
        wl = np.sqrt(wb_en1*wpix1*wb_en2*wpix2)
        #plt.plot(_l, wl**2, label='wb1*wb2')
        #plt.plot(_l, wb_en1, label='wb1')
        #plt.plot(_l, wb_en2, label='wb2')
        #plt.legend()
        #plt.show()
        cn = 0
        
        out_folder =  os.path.join(GRATOOLS_OUT, 'output_pol')
        pol_dict = data.POLCEPICE_DICT
        for key in pol_dict:
            if key == 'clfile':
                pol_dict[key] = os.path.join(out_folder,'%s_Ecrosscl.txt'%out_label)
            if key == 'cl_outmap_file':
                pol_dict[key] = os.path.join(out_folder,'%s_Ecrossclraw.txt'%out_label)
            if key == 'covfileout':
                 pol_dict[key] = os.path.join(out_folder,'%s_Ecrosscov.fits'%out_label)
            if key == 'mapfile':
                pol_dict[key] = map1
            if key == 'mapfile2':
                pol_dict[key] = map2
            if key == 'maskfile':
                pol_dict[key] = mask_f
            if key == 'maskfile2':
                pol_dict[key] = mask_f
        config_file_name = 'pol_Ecross_%s'%(out_label)
        if os.path.exists(os.path.join(out_folder,'%s_Ecrosscl.txt'%out_label)) and \
                os.path.exists(os.path.join(out_folder,'%s_Ecrosscov.fits'%out_label)):
            logger.info('ATT: Retriving power spectrum...')
            _ll, _cl, _clerr = pol_cl_parse(os.path.join(out_folder,
                                            '%s_Ecrosscl.txt'%out_label),
                                             os.path.join(out_folder,
                                            '%s_Ecrosscov.fits'%out_label),
                                            raw_corr=(cn, wl),
                                            rebin=True)
            logger.info('... and covariance matrix.')
            _cov = pol_cov_parse(os.path.join(out_folder,
                                              '%s_Ecrosscov.fits'%out_label),
                                 wl_array=wl,
                                 rebin=True, show=True)
        else:
            _ll, _cl, _clerr, _cov = pol_cl_calculation(pol_dict, 
                                                  config_file_name,
                                                  raw_corr=(cn, wl),
                                                  rebin=True,show=True)
        cl_txt.write('multipole\t%s\n'%str(list(_ll)).replace('[',''). \
                         replace(']','').replace(', ', ' '))
        cl_txt.write('Cl\t%s\n'%str(list(_cl)).replace('[',''). \
                         replace(']','').replace(', ', ' '))
        cl_txt.write('Cl_ERR\t%s\n\n'%str(list(_clerr)).replace('[',''). \
                         replace(']','').replace(', ', ' '))
    cl_txt.close()
    logger.info('Created %s'%(os.path.join(GRATOOLS_OUT,'%s_polspiceEcross.txt'
                                           %outfile_label)))



if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    mkCross(**args.__dict__)
