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

def mkRestyle(**kwargs):
    """
    """
    logger.info('Starting flux analysis...')
    get_var_from_file(kwargs['config'])
    macro_bins = data.MACRO_BINS
    gamma = data.POWER_LOW_INDEX
    out_label = data.OUT_LABEL
    binning_label = data.BINNING_LABEL
    in_labels_list = data.IN_LABELS_LIST
    new_txt_name = os.path.join(GRATOOLS_OUT, '%s_%s_parameters.txt' \
                                    %(out_label, binning_label))
    if os.path.exists(new_txt_name):
        new_txt_name = new_txt_name.replace('.txt','_2.txt')
    new_txt = open(new_txt_name,'w')
    new_txt.write('# \t E_MIN \t E_MAX \t E_MEAN \t F_MEAN \t FERR_MEAN \t CN \t FSKY \n')
    for minb, maxb in macro_bins:
        maxb = maxb + 1
        logger.info('Considering bins from %i to %i...' %(minb, maxb))
        logger.info('Retriving count and exposure maps...')
        E_MIN, E_MAX, E_MEAN = 0, 0, 0
        count_map, exposure_map = [], []
        emin, emax, emean = [], [], []
        for label in in_labels_list:
            txt_name = os.path.join(GRATOOLS_OUT, '%s_outfiles.txt' %label)
            txt = open(txt_name,'r')
            logger.info('Ref: %s'%label)
            for line in txt:
                if 'gtbin' in line:
                    count_map.append(np.asarray(hp.read_map(line, \
                                                    field=range(minb, maxb))))
                    from  GRATools.utils.gFTools import get_energy_from_fits
                    emin, emax, emean = get_energy_from_fits(line,
                                                             minbinnum=minb,
                                                             maxbinnum=maxb)
                    E_MIN, E_MAX = emin[0], emax[-1]
                    E_MEAN = (emax[0] + emin[-1])*0.5
                if 'gtexpcube2' in line:
                    exposure_map.append(np.asarray(hp.read_map(line, \
                                                    field=range(minb, maxb+1))))
            txt.close()
        logger.info('Summing in time...')
        all_counts, all_exps = count_map[0], exposure_map[0]
        for t in range(1, len(in_labels_list)):
            all_counts = all_counts + count_map[t]
            all_exps = all_exps + exposure_map[t]
        logger.info('Computing the flux for each micro energy bin...')
        flux_map, exp_mean_map = [], []
        nside = data.NSIDE
        sr = 4*np.pi/hp.nside2npix(nside)
        for i, cmap in enumerate(all_counts):
            emap = np.sqrt(all_exps[i]*all_exps[i+1])
            exp_mean_map.append(emap)
            flux_map.append(cmap/emap/sr)

        # now I have finelly gridded (in energy) summed in time fluxes
        logger.info('Rebinning...')
        logger.info('Merging fluxes from %.2f to %.2f MeV' %(E_MIN, E_MAX))
        macro_flux = flux_map[0]
        macro_fluxerr = (emean[0]/emean[0])**(-gamma)/(exp_mean_map[0])**2
        _mask = np.where(macro_flux != 1e50)[0]
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
        if len(bad_pix) != 0:
            _mask = np.unique(np.array(bad_pix))
        CN = np.mean(all_counts[0][~_mask]/(exp_mean_map[0][~_mask])**2)/sr
        for b in range(1, len(flux_map)):
            macro_flux = macro_flux + flux_map[b]
            macro_fluxerr = macro_fluxerr + \
                (emean[b]/emean[0])**(-gamma)/(exp_mean_map[b])**2
            CN = CN + np.mean(all_counts[b][~_mask]/ \
                                     (exp_mean_map[b][~_mask])**2)/sr
        logger.info('CN (white noise) term = %e'%CN)
        macro_fluxerr = np.sqrt(all_counts[0]*macro_fluxerr)/sr

        # now mask the rebinned flux and error maps
        for bpix in bad_pix:
            macro_flux[bpix] = hp.UNSEEN
            macro_fluxerr[bpix] = hp.UNSEEN
        hp.mask_bad(macro_flux)
        hp.mask_bad(macro_fluxerr)
        out_folder = os.path.join(GRATOOLS_OUT, 'output_flux')
        if not os.path.exists(out_folder):
            os.makedirs(out_folder)
        out_name = os.path.join(out_folder,out_label+'_flux_%i-%i.fits'\
                                      %(E_MIN, E_MAX))
        out_name_err = os.path.join(out_folder, out_label+'_fluxerr_%i-%i.fits'\
                                        %(E_MIN, E_MAX))
        logger.info('Created %s' %out_name)
        logger.info('Created %s' %out_name_err)
        hp.write_map(out_name, macro_flux, coord='G')
        hp.write_map(out_name_err, macro_flux, coord='G')
        _index = np.where(macro_flux != hp.UNSEEN)[0]
        F_MEAN = np.sum(macro_flux[_index])/len(macro_flux[_index])
        FERR_MEAN = np.sqrt(np.sum(macro_fluxerr[_index]**2))/\
                                   len(macro_flux[_index])
        FSKY = float(len(macro_flux[_index]))/float(len(macro_flux))
        logger.info('Fsky = %.3f'%FSKY)
        print 'F_MEAN, FERR_MEAN = ', F_MEAN, FERR_MEAN

        new_txt.write('%.2f \t %.2f \t %.2f \t %e \t %e \t %e \t %f \n' \
                          %(E_MIN, E_MAX, E_MEAN, F_MEAN, FERR_MEAN, CN, FSKY))
    new_txt.close()
    logger.info('Created %s' %os.path.join(GRATOOLS_OUT, '%s_%s_parameters.txt'\
                                             %(out_label, binning_label)))
    
    logger.info('done!')


if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    mkRestyle(**args.__dict__)
