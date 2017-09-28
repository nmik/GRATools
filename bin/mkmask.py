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
from GRATools import GRATOOLS_CONFIG
from GRATools.utils.logging_ import logger, startmsg

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('--config', type=str, required=True,
                    help='the input configuration file')
PARSER.add_argument('--srcmask', type=ast.literal_eval, choices=[True, False],
                    default=False,
                    help='sources mask activated')
PARSER.add_argument('--extsrcmask', type=ast.literal_eval, 
                    choices=[True, False], default=False,
                    help='extended sources mask activated')
PARSER.add_argument('--srcmaskweight', type=ast.literal_eval, 
                    choices=[True, False],
                    default=False,
                    help='Flux-weighted sources mask activated')
PARSER.add_argument('--srcmaskweight_custom', type=ast.literal_eval, 
                    choices=[True, False],
                    default=False,
                    help='Flux-weighted sources mask activated (Custom cat)')
PARSER.add_argument('--reversesrcmask', type=ast.literal_eval,
                    choices=[True, False],
                    default=False,
                    help='True to reverse the sources mask')
PARSER.add_argument('--gpmask', type=ast.literal_eval, choices=[True, False],
                    default=False,
                    help='galactic plain mask activated')
PARSER.add_argument('--northmask', type=ast.literal_eval, choices=[True, False],
                    default=False,
                    help='northern hemisphere mask activated')
PARSER.add_argument('--southmask', type=ast.literal_eval, choices=[True, False],
                    default=False,
                    help='souththern hemisphere mask activated')
PARSER.add_argument('--eastmask', type=ast.literal_eval, choices=[True, False],
                    default=False,
                    help='eastthern hemisphere mask activated')
PARSER.add_argument('--westmask', type=ast.literal_eval, choices=[True, False],
                    default=False,
                    help='westhern hemispheremask activated')


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
    energy = data.ENERGY
    npix = hp.nside2npix(nside)
    mask = np.ones(npix)
    if kwargs['srcmask'] == True:
        from GRATools.utils.gMasks import mask_src
        src_mask_rad = data.SRC_MASK_RAD
        cat_file = data.SRC_CATALOG
        if kwargs['reversesrcmask'] == True:
            logger.info('Reversing source mask activated...')
            allpix = np.arange(npix)
            srcmask = np.array(mask_src(cat_file, src_mask_rad, nside))
            mask = np.in1d(allpix, srcmask)
            bad_pix += list(np.where(~mask)[0])
        else:
            bad_pix += mask_src(cat_file, src_mask_rad, nside)
    if kwargs['extsrcmask'] == True:
        from GRATools.utils.gMasks import mask_extsrc
        src_mask_rad = data.SRC_MASK_RAD
        cat_file = data.SRC_CATALOG
        if kwargs['reversesrcmask'] == True:
            logger.info('Reversing source mask activated...')
            allpix = np.arange(npix)
            srcmask = np.array(mask_extsrc(cat_file, src_mask_rad, nside))
            mask = np.in1d(allpix, srcmask)
            bad_pix += list(np.where(~mask)[0])
        else:
            bad_pix += mask_extsrc(cat_file, src_mask_rad, nside)
    if kwargs['srcmaskweight'] == True:
        from GRATools.utils.gMasks import mask_src_weighted
        src_mask_rad = data.SRC_MASK_RAD
        src_cat_file = data.SRC_CATALOG
        extsrc_cat_file =data.EXTSRC_CATALOG
        if kwargs['reversesrcmask'] == True:
            logger.info('Reversing source mask activated...')
            allpix = np.arange(npix)
            srcmask = np.array(mask_src_weighted(src_cat_file, extsrc_cat_file, 
                                                 energy, nside))
            mask = np.in1d(allpix, srcmask)
            bad_pix += list(np.where(~mask)[0])
        else:
            bad_pix += mask_src_weighted(src_cat_file, extsrc_cat_file,
                                         energy, nside)
    if kwargs['srcmaskweight_custom'] == True:
        from GRATools.utils.gMasks import mask_src_weighted_custom
        src_mask_rad = data.SRC_MASK_RAD
        cat_file = data.SRC_CATALOG
        if kwargs['reversesrcmask'] == True:
            logger.info('Reversing source mask activated...')
            allpix = np.arange(npix)
            srcmask = np.array(mask_src_weighted_custom(cat_file, 
                                                        energy, nside))
            mask = np.in1d(allpix, srcmask)
            bad_pix += list(np.where(~mask)[0])
        else:
            bad_pix += mask_src_weighted_custom(cat_file, energy, nside)
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
    if kwargs['eastmask'] == True:
        from GRATools.utils.gMasks import mask_hemi_east
        bad_pix += mask_hemi_east(nside) 
    if kwargs['westmask'] == True:
        from GRATools.utils.gMasks import mask_hemi_west
        bad_pix += mask_hemi_west(nside) 
    for bpix in np.unique(bad_pix):
        mask[bpix] = 0
    out_name = os.path.join(GRATOOLS_CONFIG, 'fits/'+out_label+'.fits')
    fsky = 1-(len(np.unique(bad_pix))/float(npix))
    logger.info('f$_{sky}$ = %.3f'%fsky)
    hp.write_map(out_name, mask, coord='G')
    logger.info('Created %s' %out_name)
    
if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    mkMask(**args.__dict__)
