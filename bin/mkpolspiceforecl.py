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
import re
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
from GRATools.utils.gSpline import xInterpolatedUnivariateSplineLinear


FORE_EN = re.compile('\_\d+\-\d+\.')

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

    logger.info('Starting Cl analysis...')
    e_min = data.E_MIN
    e_max = data.E_MAX
    in_label = 'fore'
    out_label = in_label
    mask_file = data.MASK_FILE

    cl_txt = open(os.path.join(GRATOOLS_OUT, '%s_polspicecls.txt'
                               %out_label), 'w')
    for i, (emin, emax) in enumerate(zip(e_min, e_max)):
        logger.info('Considering bin %.2f - %.2f ...'%(emin, emax))
        mask_f = mask_file
        if type(mask_file) == list:
            mask_f = mask_file[i]
        cl_txt.write('ENERGY\t %.2f %.2f\n'%(emin, emax))
        l_max = 1000
        _l = np.arange(l_max)
        flux_map_name = in_label+'_%i-%i.fits'%(emin, emax)
        flux_map_f = os.path.join(GRATOOLS_OUT, 'output_fore/'+flux_map_name)
        if kwargs['show'] == True:
            hp.mollview(flux_map_masked.filled(), title='fore map',
                        min=1e-7, max=1e-4, norm='log')
            plt.show()
        out_name = '%s_%i-%i' %(out_label, emin, emax)
        out_folder =  os.path.join(GRATOOLS_OUT, 'output_pol')
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
        _l, _cl, _cl_err = pol_cl_calculation(pol_dict, config_file_name) 
        cl_txt.write('Cl\t%s\n'%str(list(_cl)).replace('[',''). \
                         replace(']','').replace(', ', ' '))
        cl_txt.write('Cl_ERR\t%s\n\n'%str(list(_cl_err)).replace('[',''). \
                         replace(']','').replace(', ', ' '))
    cl_txt.close()
    logger.info('Created %s'%(os.path.join(GRATOOLS_OUT,
                                           '%s_polspicecls.txt'
                                           %out_label)))



if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    mkCl(**args.__dict__)
