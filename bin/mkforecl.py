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

GRATOOLS_OUT_FORE = os.path.join(GRATOOLS_OUT, 'output_fore')

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

def mkForeCl(**kwargs):
    """                                      
    """
    get_var_from_file(kwargs['config'])

    logger.info('Starting Cl analysis...')
    in_label = data.IN_LABEL
    out_label = data.OUT_LABEL
    binning_label = data.BINNING_LABEL
    from GRATools.utils.gFTools import get_cl_param
    cl_param_file = os.path.join(GRATOOLS_OUT, '%s_%s_parameters.txt'
                                 %(in_label, binning_label))
    _emin, _emax, _emean, _f, _ferr, _cn, _fsky = get_cl_param(cl_param_file)
    cl_txt = open(os.path.join(GRATOOLS_OUT, '%s_%s_forecls.txt' \
                                   %(out_label, binning_label)), 'w')
    for i, (emin, emax) in enumerate(zip(_emin, _emax)):
        mask_file = data.MASK_FILE
        if type(mask_file) == list:
            mask_file = mask_file[i]
        mask = hp.read_map(mask_file)
        logger.info('Considering bin %.2f - %.2f ...'%(emin, emax))
        cl_txt.write('ENERGY\t %.2f %.2f %.2f\n'%(emin, emax, _emean[i]))
        l_max= 1000
        _l = np.arange(l_max)
        flux_map_name = in_label+'_fore_%i-%i.fits'%(emin, emax)
        flux_map = hp.read_map(os.path.join(GRATOOLS_OUT_FORE, flux_map_name))
        flux_map_masked = hp.ma(flux_map)
        flux_map_masked.mask = np.logical_not(mask)
        fsky = 1.-(len(np.where(flux_map_masked.filled() == hp.UNSEEN)[0])/\
                       float(len(flux_map)))
        if kwargs['show'] == True:
            hp.mollview(flux_map_masked.filled(), title='f$_{sky}$ = %.3f'%fsky,
                        min=1e-7, max=1e-4, norm='log')
            plt.show()
        print 'fsky = ', fsky
        nside = hp.npix2nside(len(flux_map))
        wpix = hp.sphtfunc.pixwin(nside)[:l_max]
        _cl = hp.sphtfunc.anafast(flux_map_masked.filled(), lmax=l_max-1, \
                                      iter=5)
        _cl_fit = hp.sphtfunc.anafast(flux_map_masked.filled(), iter=4)
        cl_txt.write('Cl\t%s\n'%str(list(_cl)).replace('[',''). \
                         replace(']','').replace(', ', ' '))
    cl_txt.close()
    logger.info('Created %s'%(os.path.join(GRATOOLS_OUT, '%s_%s_forecls.txt' \
                                               %(out_label, binning_label))))



if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    mkForeCl(**args.__dict__)
