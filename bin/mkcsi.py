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
from sets import Set
#from functools import reduce
import multiprocessing

__description__ = 'Makes the analysis'


"""Command-line switches.                                                       
"""


from GRATools import GRATOOLS_OUT, GRATOOLS_CONFIG
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

def csi_compute(param):
    """worker function"""
    get_var_from_file(os.path.join(GRATOOLS_CONFIG, 'Csi_config.py'))
    th_bins = data.TH_BINNING
    i, veci, dI, nside = param
    print i
    dIi = dI[i]
    dIij_list = [[] for l in range(0, len(th_bins)-1)]
    counts_list = [[] for l in range(0, len(th_bins)-1)]
    for th, (thmin, thmax) in enumerate(zip(th_bins[:-1], th_bins[1:])):
        pixintorad_min = hp.query_disc(nside, veci, thmin)
        pixintorad_max = hp.query_disc(nside, veci, thmax)
        pixintoring = np.setxor1d(pixintorad_max, pixintorad_min)
        dIj = dI[pixintoring]
        dIj = dIj[dIj > hp.UNSEEN]
        dIij = np.sum(dIi*dIj)
        #dIij = np.sum(np.array([dIi*j for j in dIj if j > hp.UNSEEN]))
        counts = len(dIj)
        #counts = len([j for j in dIj if j > hp.UNSEEN])
        dIij_list[th].append(dIij)
        counts_list[th].append(counts)
    return dIij_list, counts_list
        
def mkCsi(**kwargs):
    """                                      
    """
    get_var_from_file(kwargs['config'])
    logger.info('Calculating PSF with gtpsf...')
    dict_gtpsf = data.DICT_GTPSF
    from GRATools.utils.ScienceTools_ import gtpsf
    gtpsf(dict_gtpsf)
    from GRATools.utils.gWindowFunc import get_psf
    psf_file = data.PSF_FILE
    psf = get_psf(psf_file)
    p = multiprocessing.Pool(processes=8)
    logger.info('Starting Csi analysis...')
    in_label = data.IN_LABEL
    out_label = data.OUT_LABEL
    binning_label = data.BINNING_LABEL
    cl_param_file = os.path.join(GRATOOLS_OUT, '%s_%s_parameters.txt' \
                                     %(in_label, binning_label))
    from GRATools.utils.gFTools import get_cl_param
    _emin, _emax, _emean, _f, _ferr, _cn, _fsky = get_cl_param(cl_param_file)
    csi_txt = open(os.path.join(GRATOOLS_OUT, '%s_%s_csi.txt' \
                                   %(out_label, binning_label)), 'w')
    for i, (emin, emax) in enumerate(zip(_emin, _emax)):
        logger.info('Considering bin %.2f - %.2f ...'%(emin, emax))
        csi_txt.write('ENERGY\t %.2f %.2f %.2f\n'%(emin, emax, _emean[i]))
        flux_map_name = in_label+'_flux_%i-%i.fits'%(emin, emax)
        flux_map = hp.read_map(os.path.join(GRATOOLS_OUT_FLUX, flux_map_name))
        fsky = 1.-(len(np.where(flux_map == hp.UNSEEN)[0])/\
                       float(len(flux_map)))
        logger.info('fsky = %f'%fsky)
        npix = len(flux_map)
        nside = hp.npix2nside(npix)
        _unmask = np.where(flux_map != hp.UNSEEN)[0]
        npix_unmask = len(_unmask)
        dI = flux_map-_f[i]
        th_bins = data.TH_BINNING
        theta = []
        for thmin, thmax in zip(th_bins[:-1], th_bins[1:]):
            th_mean = np.sqrt(thmin*thmax)
            theta.append(th_mean)
        theta = np.array(theta)
        logger.info('Computing Csi...')
        diri = hp.pixelfunc.pix2ang(nside, _unmask)
        veci = hp.rotator.dir2vec(diri)
        xyz = np.array([(veci[0][i], veci[1][i], veci[2][i]) 
                        for i in range(0, len(veci[0]))])
        args = zip(_unmask, xyz, [dI]*npix_unmask,
                   [nside]*npix_unmask)
        a = np.array(p.map(csi_compute, args))
        SUMij_list = a[:,0]   
        SUMf_list = a[:, 1]
        SUMij_th = []
        SUMf_th = [] 
        for i, s in enumerate(SUMij_list[0]):
            SUMij_th.append(np.sum(SUMij_list[:, i]))
            SUMf_th.append(np.sum(SUMf_list[:, i]))
        csi = np.array(SUMij_th)/np.array(SUMf_th)
        csi_txt.write('THETA\t%s\n'%str(list(theta)).replace('[',''). \
                          replace(']','').replace(', ', ' ')) 
        csi_txt.write('CSI\t%s\n'%str(list(csi)).replace('[',''). \
                          replace(']','').replace(', ', ' '))
        plt.plot(theta, csi, 'o--')
        if kwargs['show'] == True:
            plt.show()
        save_current_figure('csi_th_%d-%d.png'%(emin, emax))
    csi_txt.close()
    p.close()
    p.join()
    logger.info('Created %s'%(os.path.join(GRATOOLS_OUT, '%s_%s_csi.txt' \
                                               %(out_label, binning_label))))

def main():
    nside = 256
    x = np.arange(hp.nside2npix(nside))
    f = np.arange(hp.nside2npix(nside))*1e-7
    m = np.where(x != 0)[0]
    p = multiprocessing.Pool(processes=8)
    n = len(x)
    args = zip(x, [f]*n, [m]*n, [nside]*n)
    a = np.array(p.map(csi_compute, args))
    SUMij_list = a[:,0]   
    SUMf_list = a[:, 1]
    SUMij_th = []
    SUMf_th = [] 
    for i, s in enumerate(SUMij_list[0]):
        SUMij_th.append(np.sum(SUMij_list[:, i]))
        SUMf_th.append(np.sum(SUMf_list[:, i]))
    print len(SUMij_th), SUMij_th
   

if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    mkCsi(**args.__dict__)
    #main()
