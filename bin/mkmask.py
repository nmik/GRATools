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

"""Produces masks
"""

import os
import imp
import numpy as np
import healpy as hp
import pyfits as pf


__description__ = 'Computes fluxes'


"""Command-line switches.
"""
import ast
import argparse
from GRATools import GRATOOLS_OUT
from GRATools.utils.logging_ import logger, startmsg

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('--config', type=str, required=True,
                    help='the input configuration file')
PARSER.add_argument('--srcmask', type=ast.literal_eval, choices=[True, False],
                    default=True,
                    help='sources mask activated')
PARSER.add_argument('--gpmask', type=ast.literal_eval, choices=[True, False],
                    default=True,
                    help='galactic plain mask activated')
PARSER.add_argument('--northmask', type=ast.literal_eval, choices=[True, False],
                    default=False,
                    help='galactic plain mask activated')
PARSER.add_argument('--southmask', type=ast.literal_eval, choices=[True, False],
                    default=False,
                    help='galactic plain mask activated')


def get_var_from_file(filename):
    f = open(filename)
    global data
    data = imp.load_source('data', '', f)
    f.close()

def mkMask(**kwargs):
    """
    """
    logger.info('Starting mask production...')
    get_var_from_file(kwargs['config'])
    bad_pix = []
    nside = data.NSIDE
    out_label = data.OUT_LABEL
    npix = hp.nside2npix(nside)
    mask = np.array([1]*npix)
    if kwargs['srcmask'] == True:
        from GRATools.utils.gMasks import mask_src
        src_mask_rad = data.SRC_MASK_RAD
        cat_file = data.SRC_CATALOG
        src_cat = pf.open(cat_file)
        bad_pix += mask_src(src_cat, src_mask_rad, nside)
        src_cat.close()
    if kwargs['gpmask'] == True:
        from GRATools.utils.gMasks import mask_gp
        gp_mask_lat = data.GP_MASK_LAT
        bad_pix += mask_gp(gp_mask_lat, nside)
    if kwargs['northmask'] == True:
        from GRATools.utils.gMasks import mask_hemi_north
        bad_pix += mask_hemi_north(nside) 
    if kwargs['southmask'] == True:
        from GRATools.utils.gMasks import mask_hemi_south
        bad_pix += mask_hemi_south(nside) 
    for bpix in np.unique(bad_pix):
        mask[bpix] = 0
    out_name = os.path.join(GRATOOLS_OUT, out_label+'.fits')
    hp.write_map(out_name, mask, coord='G')
    logger.info('Created %s' %out_name)
    
if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    mkMask(**args.__dict__)
