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


"""This app is amed to prepare the data sample to the Anisotropy analysis.
   Several steps are computed:
       1) Sum in time: the data selection has been done in 1-year-wide time 
          ------------ bins, so now it is possible to merge those bins to get 
                       the total counts and exposures.
       2) Foreground subtraction: (only if foresub is True) for each micro 
          ----------------------- energy bin the foreground model fit and 
                                  subtraction is performed. The fit is 
                                  poissonian and given a galactic diffuse 
                                  emissin model map in a given energy range and
                                  a data map in the same energy range, it 
                                  returns a constant parameter and a 
                                  normalization to apply to the model map 
                                  before performing the subtration. The data 
                                  map can be either a flux map or a count map, 
                                  in which case also the exposure map has to 
                                  be given to the fit function.
       3) Poissonian noise: the poissonian noise term is computed in each micro
          ----------------- energy bin, then summed up those in the wider 
                            energy bins to obtain the total poissonian noise 
                            term of the rebinned maps (see step 4). This 
                            calculation takes into account only the pixels 
                            outside the mask we-ll apply to the maps for the 
                            APS analysis, so the mask for each macro energy bin
                            must be given.
       3) Flux map Computation: the flux maps in each micro energy bin is 
          --------------------- obtained by dividing the counts maps for the 
                                relative exposure. If the foreground 
                                subtraction is active, the normaized foreground
                                flux map is subtracted from the data map.
       4) Energy rebinning: at this point you have flux maps in fine energy 
          ----------------- bins. Prior to the anysotropy analysis we want 
                            wider energy bins.
       The output of this app is a txt file with all the relevant infos 
       obtained running this script. This file will have the suffix 
       '_parameters' and can be parse by the function  get_cl_param in
       GRATools/utils. 

   #######
   Example: 
   python bin/mkdatarestyle.py --config restyle_config.py --foresub True
   #######       
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
print GRATOOLS_OUT
from GRATools import GRATOOLS_CONFIG
from GRATools.utils.matplotlib_ import pyplot as plt
from GRATools.utils.logging_ import logger, startmsg
from GRATools.utils.gFTools import get_energy_from_fits
from GRATools.utils.gSpline import xInterpolatedUnivariateSplineLinear

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('--config', type=str, required=True,
                    help='the input configuration file')
PARSER.add_argument('--udgrade', type=int, default=512,
                    help='down/up-grade of the maps')
PARSER.add_argument('--foresub', type=str, choices=['gal', 'galsrc'],
                    default='gal',
                    help='galactic foreground subtraction activated')

def get_var_from_file(filename):
    f = open(filename)
    global data
    data = imp.load_source('data', '', f)
    f.close()

def mkRestyle(**kwargs):
    """
    """
    logger.info('Starting the restyling...')
    get_var_from_file(kwargs['config'])
    fore_files = data.FORE_FILES_LIST
    macro_bins = data.MACRO_BINS
    gamma = data.POWER_LOW_INDEX
    out_label = data.OUT_LABEL
    mask_label = data.MASK_LABEL
    binning_label = data.BINNING_LABEL
    in_labels_list = data.IN_LABELS_LIST
    new_txt_name = os.path.join(GRATOOLS_OUT, '%s_%s_%s_parameters.txt' \
                                    %(out_label, mask_label, binning_label))
    if os.path.exists(new_txt_name):
        new_txt_name = new_txt_name.replace('.txt','_2.txt')
    new_txt = open(new_txt_name,'w')
    new_txt.write(
        '# \t E_MIN \t E_MAX \t E_MEAN \t F_MEAN \t FERR_MEAN \t CN \t FSKY \n')
    fore_mean_list = []
    norm_list, norm_sx_list, norm_dx_list = [], [], []
    norm1_list, norm1_sx_list, norm1_dx_list = [], [], []
    norm2_list, norm2_sx_list, norm2_dx_list = [], [], []
    const_list, const_sx_list, const_dx_list = [], [], []
    for i, (minb, maxb) in enumerate(macro_bins):
        all_counts, all_exps = [], []
        flux_map = []
        micro_bins = np.arange(minb, maxb+1)
        print micro_bins
        logger.info('Considering bins from %i to %i...' %(minb, maxb-1))
        mask_file = data.MASK_FILE
        if type(mask_file) == list:
            mask_file = mask_file[i]
        mask = hp.read_map(mask_file)
        _unmask = np.where(mask != 0)[0]
        maxb = maxb + 1
        exists_counts_files, exists_exp_files = [], []
        for mb in micro_bins:
            micro_count_name = os.path.join(GRATOOLS_OUT, 
                                            'output_counts/%s_counts_%i.fits'
                                            %(out_label, mb))
            micro_exp_name = os.path.join(GRATOOLS_OUT,
                                          'output_counts/%s_exposure_%i.fits'
                                          %(out_label, mb))
            if os.path.exists(micro_count_name):
                exists_counts_files.append(micro_count_name)
                exists_exp_files.append(micro_exp_name)
        if len(exists_counts_files) == len(micro_bins):
            emin, emax, emean = [], [], []
            E_MIN, E_MAX, E_MEAN = 0, 0, 0
            txt_name = os.path.join(GRATOOLS_OUT, '%s_outfiles.txt' 
                                    %in_labels_list[0])
            txt = open(txt_name,'r')
            for line in txt:
                if 'gtbin' in line:
                    emin, emax, emean = get_energy_from_fits(line,
                                                             minbinnum=minb,
                                                             maxbinnum=maxb)
                    E_MIN, E_MAX = emin[0], emax[-1]
                    E_MEAN = np.sqrt(emax[0]*emin[-1])
            logger.info('Counts and exposure maps ready! Retriving them...')
            for j in range(0, len(micro_bins)):
                cc = hp.read_map(exists_counts_files[j])
                ee = hp.read_map(exists_exp_files[j])
                all_counts.append(cc)
                all_exps.append(ee)
        else:
            logger.info('Retriving count and exposure maps...')
            emin, emax, emean = [], [], []
            E_MIN, E_MAX, E_MEAN = 0, 0, 0
            count_map, exp_mean_map = [], []
            for label in in_labels_list:
                txt_name = os.path.join(GRATOOLS_OUT, '%s_outfiles.txt' %label)
                txt = open(txt_name,'r')
                logger.info('Ref: %s'%label)
                for line in txt:
                    if 'gtbin' in line:
                        cmap = hp.read_map(line, field=range(minb, maxb))
                        cmap_repix = hp.pixelfunc.ud_grade(cmap,
                                                           kwargs['udgrade'], 
                                                           pess=True,
                                                           power=-2)
                        count_map.append(np.asarray(cmap_repix))
                        emin, emax, emean = get_energy_from_fits(line,
                                                                 minbinnum=minb,
                                                                 maxbinnum=maxb)
                    
                        E_MIN, E_MAX = emin[0], emax[-1]
                        E_MEAN = np.sqrt(emax[0]*emin[-1])
                    if 'gtexpcube2' in line:
                        emap = hp.read_map(line, field=range(minb, maxb+1))
                        emap_repix = hp.pixelfunc.ud_grade(emap,
                                                           kwargs['udgrade'], 
                                                           pess=True)
                        emap_mean = []
                        for i in range(0,len(emap_repix)-1):
                            emap_mean.append(
                                np.sqrt(emap_repix[i]*emap_repix[i+1]))    
                        exp_mean_map.append(np.asarray(emap_mean))
                txt.close()
            logger.info('Summing in time...')
            all_counts, all_exps = count_map[0], exp_mean_map[0]
            for t in range(1, len(in_labels_list)):
                all_counts = all_counts + count_map[t]
                all_exps = all_exps + exp_mean_map[t]

            for i, cmap in enumerate(all_counts):
                micro_count_name = os.path.join(GRATOOLS_OUT, 
                                             'output_counts/%s_counts_%i.fits'
                                                %(out_label, micro_bins[i]))
                hp.write_map(micro_count_name, cmap)
                micro_exp_name = os.path.join(GRATOOLS_OUT, 
                                           'output_counts/%s_exposure_%i.fits'
                                              %(out_label, micro_bins[i]))
                hp.write_map(micro_exp_name, all_exps[i])

        logger.info('Computing the flux for each micro energy bin...')
        nside = kwargs['udgrade']
        npix = hp.nside2npix(nside)
        sr = 4*np.pi/npix
        iii = np.arange(npix)
        for i, cmap in enumerate(all_counts):
            flux_map.append(cmap/all_exps[i]/sr)
    
        # now I have finelly gridded (in energy) summed in time fluxes
        logger.info('Rebinning...')
        logger.info('Merging fluxes from %.2f to %.2f MeV' %(E_MIN, E_MAX))
        macro_fluxerr = (emean[0]/emean[0])**(-gamma)/(all_exps[0])**2
        macro_counts = all_counts[0]

        # implement foreground subtraction
        from GRATools.utils.gForeground import get_fore_integral_flux_map
        from GRATools.utils.gForeground import get_ref_igrb_spline
        if kwargs['foresub'] == 'gal':
            _norm_list, _norm_sx_list, _norm_dx_list = [], [], []
            _const_list, _const_sx_list, _const_dx_list = [], [], []
            from GRATools.utils.gForeground import fit_foreground_poisson
            out_fore_folder = os.path.join(GRATOOLS_OUT, 'output_fore')
            out_name_fore = os.path.join(out_fore_folder,'fore_%i-%i.fits'\
                                             %(E_MIN, E_MAX))
            if not os.path.exists(out_fore_folder):
                os.makedirs(out_fore_folder)
            all_fore = []
            all_c_guess = []
            for ii, (e1, e2) in enumerate(zip(emin, emax)):
                fore = get_fore_integral_flux_map(fore_files, e1, e2)
                all_fore.append(fore)
                all_c_guess.append(get_ref_igrb_spline()(emean[ii])*100)
            n0, c0, n0_sx, n0_dx, c0_sx, c0_dx = \
                fit_foreground_poisson(all_fore[0], 
                                       all_counts[0],
                                       mask_map=mask,
                                       exp=all_exps[0],
                                       n_guess=1., 
                                       c_guess=all_c_guess[0])
            _norm_list.append(n0)
            _norm_sx_list.append(n0_sx)
            _norm_dx_list.append(n0_dx)
            _const_list.append(c0)
            _const_sx_list.append(c0_sx)
            _const_dx_list.append(c0_dx)
            macro_flux = flux_map[0]-n0*all_fore[0]
            macro_fore = n0*all_fore[0]
            CN = np.mean(all_counts[0][_unmask]/(all_exps[0][_unmask])**2)/sr 
            for b in range(1, len(flux_map)):
                n, c, n_sx, n_dx, c_sx, c_dx = \
                    fit_foreground_poisson(all_fore[b], 
                                           all_counts[b], 
                                           mask_map=mask,
                                           exp=all_exps[b],
                                           n_guess=1., 
                                           c_guess=all_c_guess[b])
                _norm_list.append(n)
                _norm_sx_list.append(n_sx)
                _norm_dx_list.append(n_dx)
                _const_list.append(c)
                _const_sx_list.append(c_sx)
                _const_dx_list.append(c_dx)
                fluxerr = (emean[b]/emean[0])**(-gamma)/(all_exps[b])**2
                macro_fluxerr = macro_fluxerr + fluxerr
                macro_flux = macro_flux + flux_map[b]-n*all_fore[b]
                macro_counts = macro_counts + all_counts[b]
                macro_fore = macro_fore + n*all_fore[b]
                CN = CN + np.mean(all_counts[b][_unmask]/ \
                                      (all_exps[b][_unmask])**2)/sr
            logger.info('CN (white noise) term = %e'%CN)
            macro_fluxerr = (np.sqrt(all_counts[0]*macro_fluxerr)/sr)
            macro_fore_masked = hp.ma(macro_fore)
            macro_fore_masked.mask = np.logical_not(mask)
            hp.write_map(out_name_fore, macro_fore, coord='G')
            logger.info('Created %s' %out_name_fore)
            FORE_MEAN = np.mean(macro_fore[_unmask])
            fore_mean_list.append(FORE_MEAN)
            norm_list.append(np.mean(np.array(_norm_list)))
            norm_sx_list.append(np.amin(np.array(_norm_sx_list)))
            norm_dx_list.append(np.amax(np.array(_norm_dx_list)))
            const_list.append(np.sum(np.array(_const_list)))
            const_sx_list.append(np.sum(np.array(_const_sx_list)))
            const_dx_list.append(np.sum(np.array(_const_dx_list)))
        elif kwargs['foresub'] == 'galsrc':
            cat_file = data.CAT_FILE
            psf_file = data.PSF_FILE 
            _norm1_list, _norm1_sx_list, _norm1_dx_list = [], [], []
            _norm2_list, _norm2_sx_list, _norm2_dx_list = [], [], []
            _const_list, _const_sx_list, _const_dx_list = [], [], []
            from GRATools.utils.gSourceTemplate import build_src_template
            from GRATools.utils.gForeground import fit_fore_src_poisson
            out_fore_folder = os.path.join(GRATOOLS_OUT, 'output_fore')
            out_name_fore = os.path.join(out_fore_folder,'fore_%i-%i.fits'\
                                             %(E_MIN, E_MAX))
            if not os.path.exists(out_fore_folder):
                os.makedirs(out_fore_folder)
            all_fore = []
            all_srctempl = []
            all_c_guess = []
            for ii, (e1, e2) in enumerate(zip(emin, emax)):
                fore = get_fore_integral_flux_map(fore_files, e1, e2)
                all_fore.append(fore)
                cat = os.path.basename(cat_file).replace('.fit','')
                out_name_srctempl = os.path.join(GRATOOLS_OUT,
                                                 'output_src/src%s_%i-%i.fits'\
                                                     %(cat, e1, e2))       
                if not os.path.exists(out_name_srctempl):
                    srctempl_map = build_src_template(cat_file, psf_file,
                                                     emin=e1, emax=e2, b_cut=10)
                else:
                    logger.info('Retriving source template...')
                    srctempl_map = hp.read_map(out_name_srctempl)
                all_srctempl.append(srctempl_map)
                all_c_guess.append(get_ref_igrb_spline()(emean[ii])*100)
            n10, n20, c0, (n10_sx, n10_dx), (n20_sx, n20_dx), (c0_sx, c0_dx) = \
                fit_fore_src_poisson(all_fore[0], 
                                       all_counts[0],
                                       all_srctempl[0],
                                       mask_map=mask,
                                       exp=all_exps[0],
                                       n1_guess=1., 
                                       n2_guess=1., 
                                       c_guess=all_c_guess[0])
            _norm1_list.append(n10)
            _norm1_sx_list.append(n10_sx)
            _norm1_dx_list.append(n10_dx)
            _norm2_list.append(n20)
            _norm2_sx_list.append(n20_sx)
            _norm2_dx_list.append(n20_dx)
            _const_list.append(c0)
            _const_sx_list.append(c0_sx)
            _const_dx_list.append(c0_dx)
            macro_flux = flux_map[0]-n10*all_fore[0]-n20*all_srctempl[0]
            macro_fore = n10*all_fore[0]
            CN = np.mean(all_counts[0][_unmask]/(all_exps[0][_unmask])**2)/sr 
            for b in range(1, len(flux_map)):
                n1, n2, c, (n1_sx, n1_dx), (n2_sx, n2_dx), (c_sx, c_dx) = \
                    fit_fore_src_poisson(all_fore[b], 
                                         all_counts[b],
                                         all_srctempl[b],
                                         mask_map=mask,
                                         exp=all_exps[b],
                                         n1_guess=1., 
                                         n2_guess=1.,
                                         c_guess=all_c_guess[b])
                _norm1_list.append(n1)
                _norm1_sx_list.append(n1_sx)
                _norm1_dx_list.append(n1_dx)
                _norm2_list.append(n2)
                _norm2_sx_list.append(n2_sx)
                _norm2_dx_list.append(n2_dx)
                _const_list.append(c)
                _const_sx_list.append(c_sx)
                _const_dx_list.append(c_dx)
                fluxerr = (emean[b]/emean[0])**(-gamma)/(all_exps[b])**2
                macro_fluxerr = macro_fluxerr + fluxerr
                macro_flux = macro_flux+flux_map[b]-n1*all_fore[b]-n2*all_srctempl[b]
                macro_counts = macro_counts + all_counts[b]
                macro_fore = macro_fore + n1*all_fore[b]
                netcounts = all_counts[0]
                CN = CN + np.mean(all_counts[b][_unmask]/ \
                                      (all_exps[b][_unmask])**2)/sr
            logger.info('CN (white noise) term = %e'%CN)
            macro_fluxerr = (np.sqrt(all_counts[0]*macro_fluxerr)/sr)
            macro_fore_masked = hp.ma(macro_fore)
            macro_fore_masked.mask = np.logical_not(mask)
            hp.write_map(out_name_fore, macro_fore, coord='G')
            logger.info('Created %s' %out_name_fore)
            FORE_MEAN = np.mean(macro_fore[_unmask])
            fore_mean_list.append(FORE_MEAN)
            norm1_list.append(np.mean(np.array(_norm1_list)))
            norm1_sx_list.append(np.amin(np.array(_norm1_sx_list)))
            norm1_dx_list.append(np.amax(np.array(_norm1_dx_list)))
            norm2_list.append(np.mean(np.array(_norm2_list)))
            norm2_sx_list.append(np.amin(np.array(_norm2_sx_list)))
            norm2_dx_list.append(np.amax(np.array(_norm2_dx_list)))
            const_list.append(np.sum(np.array(_const_list)))
            const_sx_list.append(np.sum(np.array(_const_sx_list)))
            const_dx_list.append(np.sum(np.array(_const_dx_list)))
        else:
            CN = np.mean(all_counts[0][_unmask]/(all_exps[0][_unmask])**2)/sr
            macro_flux = flux_map[0]
            for b in range(1, len(flux_map)):
                fluxerr = (emean[b]/emean[0])**(-gamma)/(all_exps[b])**2
                macro_fluxerr = macro_fluxerr + fluxerr
                macro_flux = macro_flux + flux_map[b]
                macro_counts = macro_counts + all_counts[b]
                CN = CN + np.mean(all_counts[b][_unmask]/ \
                                      (all_exps[b][_unmask])**2)/sr
            logger.info('CN (white noise) term = %e'%CN)
            macro_fluxerr = (np.sqrt(all_counts[0]*macro_fluxerr)/sr)
        
        out_count_folder = os.path.join(GRATOOLS_OUT, 'output_counts')
        if not os.path.exists(out_count_folder):
            os.makedirs(out_count_folder)
        out_counts_name = os.path.join(out_count_folder,out_label+'_counts_%i-%i.fits'\
                                      %(E_MIN, E_MAX))
        logger.info('Created %s' %out_counts_name)
        hp.write_map(out_counts_name, macro_counts, coord='G')
        # now mask the rebinned flux maps        
        macro_flux_masked = hp.ma(macro_flux)
        macro_fluxerr_masked = hp.ma(macro_fluxerr)
        macro_flux_masked.mask = np.logical_not(mask)
        macro_fluxerr_masked.mask = np.logical_not(mask)
        out_folder = os.path.join(GRATOOLS_OUT, 'output_flux')
        if not os.path.exists(out_folder):
            os.makedirs(out_folder)
        out_name = os.path.join(out_folder,
                                out_label+'_%s_fluxmasked_%i-%i.fits'\
                                    %(mask_label, E_MIN, E_MAX))
        out_name_unmask = os.path.join(out_folder, 
                                       out_label+'_%s_flux_%i-%i.fits'\
                                           %(mask_label, E_MIN, E_MAX))
        hp.write_map(out_name, macro_flux_masked, coord='G')
        hp.write_map(out_name_unmask, macro_flux, coord='G')
        logger.info('Created %s' %out_name)
        logger.info('Created %s' %out_name_unmask)
        F_MEAN = np.sum(macro_flux[_unmask])/len(macro_flux[_unmask])
        FERR_MEAN = np.sqrt(np.sum(macro_fluxerr[_unmask]**2))/\
                                   len(macro_flux[_unmask])
        FSKY = float(len(macro_flux[_unmask]))/float(len(macro_flux))
        logger.info('Fsky = %.3f'%FSKY)
        print('F_MEAN, FERR_MEAN = ', F_MEAN, FERR_MEAN)
        new_txt.write('%.2f \t %.2f \t %.2f \t %e \t %e \t %e \t %f \n' \
                          %(E_MIN, E_MAX, E_MEAN, F_MEAN, FERR_MEAN, CN, FSKY))
    if kwargs['foresub'] == 'gal':
        new_txt.write('\n\n*** FOREGROUND PARAMETERS***\n\n')
        new_txt.write('MEAN_FORE_FLUX \t %s\n' %str(fore_mean_list))
        new_txt.write('NORM_FIT_PARAM \t %s\n' %str(norm_list))
        new_txt.write('NORM_FIT_PARAM  errsx \t %s\n' %str(norm_sx_list))
        new_txt.write('NORM_FIT_PARAM errdx \t %s\n' %str(norm_dx_list))
        new_txt.write('IGRB_FIT_PARAM \t %s\n' %str(const_list))
        new_txt.write('IGRB_FIT_PARAM errsx \t %s\n' %str(const_sx_list))
        new_txt.write('IGRB_FIT_PARAM errdx \t %s\n' %str(const_dx_list))
    if kwargs['foresub'] == 'galsrc':
        new_txt.write('\n\n*** FOREGROUND PARAMETERS***\n\n')
        new_txt.write('MEAN_FORE_FLUX \t %s\n' %str(fore_mean_list))
        new_txt.write('NORM_FORE_FIT PARAM \t %s\n' %str(norm1_list))
        new_txt.write('NORM_FORE_FIT PARAM  errsx \t %s\n' %str(norm1_sx_list))
        new_txt.write('NORM_FORE_FIT PARAM errdx \t %s\n' %str(norm1_dx_list))
        new_txt.write('NORM_SRC_FIT PARAM \t %s\n' %str(norm2_list))
        new_txt.write('NORM_SRC_FIT PARAM  errsx \t %s\n' %str(norm2_sx_list))
        new_txt.write('NORM_SRC_FIT PARAM errdx \t %s\n' %str(norm2_dx_list))
        new_txt.write('IGRB_FIT_PARAM \t %s\n' %str(const_list))
        new_txt.write('IGRB_FIT_PARAM errsx \t %s\n' %str(const_sx_list))
        new_txt.write('IGRB_FIT_PARAM errdx \t %s\n' %str(const_dx_list))
    new_txt.close()
    logger.info('Created %s' %os.path.join(GRATOOLS_OUT, 
                                     '%s_%s_%s_parameters.txt'\
                                      %(out_label, mask_label, binning_label)))
    logger.info('done!')


if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    mkRestyle(**args.__dict__)
