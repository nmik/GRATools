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
import numpy as np
import healpy as hp
import pyfits as pf


__description__ = 'Makes the analysis'



"""Command-line switches.
"""

import argparse
from GRATools import GRATOOLS_OUT
from GRATools.utils.logging_ import logger, startmsg
from GRATools.utils.matplotlib_ import pyplot as plt
from GRATools.utils.matplotlib_ import overlay_tag, save_current_figure

GRATOOLS_OUT_FLUX = os.path.join(GRATOOLS_OUT, 'output_flux')

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('--config', type=str, required=True,
                    help='the input configuration file')

def get_var_from_file(filename):
    f = open(filename)
    global data
    data = imp.load_source('data', '', f)
    f.close()

def mkAnalysis(**kwargs):
    """
    """
    logger.info('Starting analysis...')
    get_var_from_file(kwargs['config'])
    files_dict = data.FILES_DICT
    out_label = data.OUT_LABEL
    binning_label = data.BINNING_LABEL
    label_list = data.LABELS_LIST
    ebinning_file = os.path.join(GRATOOLS_OUT, '%s_%s_new_ebinning_cn.txt' \
                                     %(label_list[0], binning_label))
    from  GRATools.utils.gFTools import get_energy_from_txt
    from  GRATools.utils.gFTools import get_cn_from_txt
    _emean, _emin, _emax = get_energy_from_txt(ebinning_file, get_binning=True)
    final_maps, final_err_maps, = [], []
    _fmean = [0]*len(_emean)
    _fmeanerr = [0]*len(_emean)
    _cnmean = [0]*len(_emean)
    nyrs = len(label_list)
    fsky = 0
    for bin, labels in files_dict.iteritems():
        flux_bin, fluxerr_bin, = [], []
        _index = []
        bin_num = labels[0]
        en, emin, emax = _emean[bin_num], _emin[bin_num], _emax[bin_num]
        logger.info('considering bin: %i MeV...' %int(en))
        _cn = []
        for label in labels[1:]:
            cn_file = os.path.join(GRATOOLS_OUT,'%s_%s_new_ebinning_cn.txt' \
                                       %(label, binning_label))
            cn = get_cn_from_txt(cn_file)[bin_num]
            _cn.append(cn)
            flux_map_name = label+'_flux_'+ bin + '.fits'
            flux_map = hp.read_map(os.path.join(GRATOOLS_OUT_FLUX, \
                                                    flux_map_name))
            fluxerr_map_name = flux_map_name.replace('flux', 'fluxerr')
            fluxerr_map = hp.read_map(os.path.join(GRATOOLS_OUT_FLUX, \
                                                    fluxerr_map_name))
            _index = np.where(flux_map != hp.UNSEEN)[0]
            flux_bin.append(flux_map)
            fluxerr_bin.append(fluxerr_map)
        _cnmean[bin_num] = np.average(np.array(_cn))
        flux_bin_sum = flux_bin[0]
        fluxerr_bin_sum = fluxerr_bin[0]
        fluxerr_bin_sum[_index] = fluxerr_bin[0][_index]**2
        for i, m in enumerate(flux_bin[1:]):
            flux_bin_sum[_index] = flux_bin_sum[_index] + m[_index]
            fluxerr_bin_sum[_index] = fluxerr_bin_sum[_index] + \
                fluxerr_bin[i][_index]**2
        final_maps.append(flux_bin_sum)
        final_err_maps.append(np.sqrt(fluxerr_bin_sum))
        fsky = float(len(flux_bin_sum[_index]))/float(len(flux_bin_sum))
        mean_flux = np.sum(flux_bin_sum[_index])/len(flux_bin_sum[_index])/nyrs
        mean_fluxerr = np.sqrt(np.sum(flux_bin_sum[_index]**2)) \
            /len(flux_bin_sum[_index])/nyrs
        _fmean[bin_num] = mean_flux
        _fmeanerr[bin_num] = mean_fluxerr
        out_name = os.path.join(GRATOOLS_OUT_FLUX, '%s_flux_%i-%i.fits' \
                                    %(out_label, int(emin), int(emax)) )
        hp.write_map(out_name, flux_bin_sum, coord='G')
        logger.info('Created %s' %out_name)
    
    txt_outfile = open(os.path.join(GRATOOLS_OUT, out_label+'_report.txt'),'w')
    txt_outfile.write('E_MIN\t'+str([i for i in _emin]).replace('[','').\
                          replace(']','').replace(',','')+'\n')
    txt_outfile.write('E_MAX\t'+str([i for i in _emax]).replace('[','').\
                          replace(']','').replace(',','')+'\n')
    txt_outfile.write('E_MEAN\t'+str([i for i in _emean]).replace('[','').\
                          replace(']','').replace(',','')+'\n')
    txt_outfile.write('F_MEAN\t'+str([i for i in _fmean]).replace('[','').\
                          replace(']','').replace(',','')+'\n')
    txt_outfile.write('FERR_MEAN\t'+str([i for i in _fmeanerr]).replace('[','').\
                          replace(']','').replace(',','')+'\n')
    txt_outfile.write('CN_MEAN\t'+str([i for i in _cnmean]).replace('[','').\
                          replace(']','').replace(',','')+'\n')
    txt_outfile.write('FSKY\t'+str(fsky)+'\n')
    logger.info('Created %s' %os.path.join(GRATOOLS_OUT, out_label+'_report.txt'))
    
    plt.figure(figsize=(10, 7), dpi=80)
    _fmean = np.array(_fmean)
    _fmeanerr = np.array(_fmeanerr)
    from GRATools.utils.gDrawRef import ref_igrb_band
    from GRATools.utils.gDrawRef import ref_igrb_noFGsub
    leg, lab = ref_igrb_band()
    igrb, lab_igrb = ref_igrb_noFGsub()
    spec = plt.errorbar(_emean, _fmean*_emean*_emean, fmt='o', markersize=3, \
                     elinewidth=1, xerr=(_emax-_emin)/2, \
                     yerr=_fmeanerr*_emean*_emean)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel('Energy [MeV]')
    plt.ylabel('E$^{2}$ $\cdot$ Flux [MeV cm$^{-2}$ s$^{-1}$ sr$^{-1}$]')
    plt.title('Extra-Galactic Energy Spectrum')
    plt.legend([spec, igrb, leg], [out_label, lab_igrb, lab])
    overlay_tag(x=0.45, y=0.05)
    save_current_figure(out_label+'_ESpec.png')    


if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    mkAnalysis(**args.__dict__)
