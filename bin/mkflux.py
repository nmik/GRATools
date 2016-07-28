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


"""Flux analysis
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
PARSER.add_argument('--infile', type=str, required=True,
                    help='txt file with paths of ST output files')
PARSER.add_argument('--config', type=str, required=True,
                    help='the input configuration file')
PARSER.add_argument('--srcmask', type=ast.literal_eval, choices=[True, False],
                    default=True,
                    help='sources mask activated')
PARSER.add_argument('--gpmask', type=ast.literal_eval, choices=[True, False],
                    default=True,
                    help='galactic plain mask activated')

def get_var_from_file(filename):
    f = open(filename)
    global data
    data = imp.load_source('data', '', f)
    f.close()

def mkFlux(**kwargs):
    """
    """
    logger.info('Starting flux analysis...')
    get_var_from_file(kwargs['config'])
    txt = open(kwargs['infile'],'r')
    count_map, exposure_map, energy = [], [], []
    nbins = data.MICRO_NBINS
    logger.info('Retriving count and exposure maps...')
    for line in txt:
        if 'gtbin' in line:
            count_map = hp.read_map(line, field=range(0,nbins))
            from  GRATools.utils.gFTools import get_energy_from_fits
            energy = get_energy_from_fits(line)
        if 'gtexpcube2' in line:
            exposure_map = hp.read_map(line, field=range(0,nbins+1))
    flux_map = []
    logger.info('Computing the flux')
    nside = data.NSIDE
    sr = 4*np.pi/hp.nside2npix(nside)
    for i, cmap in enumerate(count_map):
        emap = np.sqrt(exposure_map[i]*exposure_map[i+1])
        flux_map.append(cmap/emap/sr/(energy[i]/1000))
    macro_bins = data.MACRO_BINS
    gamma = data.POWER_LOW_INDEX
    out_label = data.OUT_LABEL
    logger.info('Rebinning...')
    new_ebinning_txt = open(os.path.join(GRATOOLS_OUT,out_label+\
                                             '_new_ebinning_cn.txt'), 'w')
    new_ebinning_txt.write('#MACRO ENERGY BINNING\n\n')
    CNs = []
    for minb, maxb in macro_bins:
        logger.info('Merging fluxes from %.2f to %.2f MeV' \
                        %(energy[minb]/1000, energy[maxb]/1000))
        new_ebinning_txt.write('%.2f %.2f\n' \
                                   %(energy[minb]/1000, energy[maxb]/1000))
        macro_flux = flux_map[minb]
        macro_fluxerr = (energy[minb]/energy[0])**(-gamma)/(exposure_map[minb]*\
                                                                energy[minb])**2
        CN_bin = np.average(count_map[minb]/(exposure_map[minb])**2)/sr
        for b in range(minb+1, maxb):
            macro_flux = macro_flux + flux_map[b]
            macro_fluxerr = macro_fluxerr + \
                (energy[b]/energy[0])**(-gamma)/(exposure_map[b]*energy[b])**2
            CN_bin = CN_bin + np.average(count_map[b]/(exposure_map[b])**2)/sr
        logger.info('CN (white noise) term = %e'%CN_bin)
        macro_fluxerr = np.sqrt(count_map[0]*macro_fluxerr)/sr
        CNs.append(CN_bin)
        bad_pix = []
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
        for bpix in bad_pix:
            macro_flux[bpix] = hp.UNSEEN
            macro_fluxerr[bpix] = hp.UNSEEN
        hp.mask_bad(macro_flux)
        hp.mask_bad(macro_fluxerr)
        out_folder = os.path.join(GRATOOLS_OUT,'output_flux')
        if not os.path.exists(out_folder):
            os.makedirs(out_folder)
        out_name = os.path.join(out_folder,out_label+'_flux_%i-%i.fits'\
                                      %(minb, maxb))
        out_name_err = os.path.join(out_folder,out_label+'_fluxerr_%i-%i.fits'\
                                        %(minb, maxb))
        hp.write_map(out_name, macro_flux, coord='G')
        hp.write_map(out_name_err, macro_flux, coord='G')
        logger.info('Created %s' %out_name)
        logger.info('Created %s' %out_name_err)
    new_ebinning_txt.write('\n\n#MACRO BIN WHITE NOISE\n\n')
    for cn in CNs:
        new_ebinning_txt.write(str(cn)+'\n')
    new_ebinning_txt.close()
    logger.info('Created %s' %os.path.join(GRATOOLS_OUT,out_label+\
                                                   '_new_ebinning_cn.txt'))
    logger.info('done!')

if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    mkFlux(**args.__dict__)
