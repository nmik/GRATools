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


"""Analysis utilities 
"""

  
import os
import re
import numpy as np
import pyfits as pf
import healpy as hp
from GRATools import GRATOOLS_OUT
from GRATools import GRATOOLS_CONFIG
from GRATools.utils.logging_ import logger, abort

FORE_EN = re.compile('\_\d+\.')


def get_foreground_integral_flux_map(fore_files_list, e_min, e_max):
    """mean_energy: float
       the mean of the energy bin
    """
    fore_en = []
    for f in fore_files_list:
        m_en = FORE_EN.search(f).group()
        en = float(m_en.replace('_','').replace('.',''))
        fore_en.append(en)
    fore_en = np.array(fore_en)
    print fore_en
    out_name = fore_files_list[0].replace('_%i.fits'%fore_en[0], 
                                          '_%d-%d.fits'%(e_min, e_max))
    if os.path.exists(out_name):
        logger.info('ATT: file %s already exists and returned...'%out_name)
        fore_map = hp.read_map(out_name)
        return fore_map
    else: 
        logger.info('Computing the integral flux of the foreground model...')
        logger.info('...between %.2f - %.2f'%(e_min, e_max))
        if not list(np.where(fore_en>e_max)[0]):
            print np.where(fore_en>e_min)[0]
            if not list(np.where(fore_en>e_min)[0]):
                fore_map1 = hp.read_map(fore_files_list[-2])
                fore_map2 = hp.read_map(fore_files_list[-1])
                fore_e_min = fore_en[-2] 
                fore_e_max = fore_en[-1]
            else:
                fore_map2 = hp.read_map(fore_files_list[np.where(\
                            fore_en<e_max)[0][-1]])
                fore_e_max = fore_en[np.where(fore_en<e_max)[0][-1]] 
                fore_map1 = hp.read_map(fore_files_list[np.where(\
                            fore_en<e_min)[0][-1]])
                fore_e_min = fore_en[np.where(fore_en<e_min)[0][-1]]
        else:
            fore_map1 = hp.read_map(fore_files_list[np.where(\
                        fore_en<e_min)[0][-1]])
            fore_map2 = hp.read_map(fore_files_list[np.where(\
                        fore_en>e_max)[0][0]])
            fore_e_max = fore_en[np.where(fore_en>e_max)[0][0]]
            fore_e_min = fore_en[np.where(fore_en<e_min)[0][-1]]
        logger.info('getting foreground between %.2f - %.2f'\
                        %(fore_e_min, fore_e_max))
        A = (fore_map2 - fore_map1)/(fore_e_max - fore_e_min)
        B = fore_map1 - fore_e_min*A
        fore_integr = (A*np.sqrt(e_max*e_min) + B)
        hp.write_map(out_name, fore_integr)
        return fore_integr

def csi_parse(csi_file):
    """Parsing of the *_csi.txt files
    """
    logger.info('loading Csi values from %s'%csi_file)
    csi = []
    theta = []
    emin, emax, emean = [], [], []
    f = open(csi_file, 'r')
    for line in f:
        if 'ENERGY\t' in line:
            e1, e2, em = [float(item) for item in line.split()[1:]]
            emin.append(e1)
            emax.append(e2)
            emean.append(em)
        if 'CSI\t' in line:
            c = np.array([float(item) for item in line.split()[1:]])
            csi.append(c)
        if 'THETA\t' in line:
            th = np.array([float(item) for item in line.split()[1:]])
            theta.append(th)
    f.close()
    return np.array(emin), np.array(emax), np.array(emean), csi, theta

def cp_parse(cp_file):
    """Parsing of the *_cps.txt files
    """
    logger.info('loading Cp values from %s'%cp_file)
    ff = open(cp_file, 'r')
    _emin, _emax, _emean, _cp, _cperr = [], [], [], [], []
    for line in ff:
        try:
            emin, emax, emean, cp, cperr = [float(item) for item in \
                                                line.split()]
            _emin.append(emin)
            _emax.append(emax)
            _emean.append(emean)
            _cp.append(cp)
            _cperr.append(cperr)
        except:
            pass
    ff.close()
    return np.array(_emin), np.array(_emax),np.array(_emean), \
        np.array(_cp), np.array(_cperr)

def cl_parse(cl_file):
    """Parsing of the *_cls.txt files.
    """
    cls = []
    clerrs = []
    emin, emax, emean = [], [], []
    f = open(cl_file, 'r')
    for line in f:
        if 'ENERGY\t' in line:
            e1, e2, em = [float(item) for item in line.split()[1:]]
            emin.append(e1)
            emax.append(e2)
            emean.append(em)
        if 'Cl\t' in line:
            cl = np.array([float(item) for item in line.split()[1:]])
            cls.append(cl)
        if 'Cl_ERR' in line:
            cl_err = np.array([float(item) for item in line.split()[1:]])
            clerrs.append(cl_err)
    f.close()
    return np.array(emin), np.array(emax), np.array(emean), cls, clerrs

def get_cl_param(cl_param_file):
    """Parsing of *_parameters.txt files.
    """
    logger.info('loading parameters from %s'%cl_param_file)
    ff = open(cl_param_file, 'r')
    _emin, _emax, _emean, _f, _ferr, _cn, _fsky = [], [], [], [], [], \
        [], []
    for line in ff:
        try:
            emin, emax, emean, f, ferr, cn, fsky = [float(item) for item in \
                                                        line.split()]
            _emin.append(emin)
            _emax.append(emax)
            _emean.append(emean)
            _f.append(f)
            _ferr.append(ferr)
            _cn.append(cn)
            _fsky.append(fsky)
        except:
            pass
    ff.close()
    return np.array(_emin), np.array(_emax), np.array(_emean), np.array(_f), \
        np.array(_ferr), np.array(_cn), np.array(_fsky)
        
def get_crbkg(txt_file, str_evtype):
    """Not used
    """
    f = open(txt_file, 'r')
    _bkg, _en = [], []
    for line in f:
        if 'I_'+str_evtype in line:
            _bkg = [float(item) for item in line.split()[1:]]
        if 'E_'+str_evtype in line:
            _en = [float(item) for item in line.split()[1:]]
    f.close()
    return np.array(_en), np.array(_bkg)

def get_energy_from_txt(txt_file, get_binning=False, mean='log', ):
    """Returns a list with the center values of the energy bins

       txt_file: str
           txt file produced in output by mkflux.py application
       mean: str
           'log' or 'lin', dependind on the algorithm to use to 
           compute the mean
    """
    f = open(txt_file, 'r')
    _emin, _emax = [], []
    for line in f:
        try:
            emin, emax  = [float(item) for item in line.split()]
            _emin.append(emin)
            _emax.append(emax)
        except:
            pass
    f.close()
    emean = []        
    if mean == 'log':
        for emin, emax in zip(_emin, _emax):
            emean.append(np.sqrt(emin*emax))
    if mean == 'lin':
        for emin, emax in zip(_emin, _emax):
            emean.append(0.5*(emin+emax))
    if get_binning == True:
        return np.array(emean), np.array(_emin), np.array(_emax)
    return np.array(_emin), np.array(_emax), np.array(emean)

def get_energy_from_fits(fits_file, minbinnum=0, maxbinnum=100 ,mean='log'):
    """Returns a list with the center values of the energy bins

       fits_file: str
           fits file, usually we want to do to this at the level of
           gtbin output file
       mean: str
           'log' or 'lin', dependind on the algorithm to use to 
           compute the mean
    """
    f = pf.open(fits_file)
    ebounds = f[2].data
    _emin = ebounds['E_MIN'][minbinnum:maxbinnum]/1000
    _emax = ebounds['E_MAX'][minbinnum:maxbinnum]/1000
    emean = []        
    if mean == 'log':
        for emin, emax in zip(_emin, _emax):
            emean.append(np.sqrt(emin*emax))
    if mean == 'lin':
        for emin, emax in zip(_emin, _emax):
            emean.append(0.5*(emin+emax))
    f.close()
    return np.array(_emin), np.array(_emax), np.array(emean)

def ebinning_fits_file(ebinning_array):
    """Produces a fits file defining the enrgy binning to fed gtbin.

       ebinning_array: numpy array
           array in which the energy binnin is defined.
    """
    txt_file_name = os.path.join(GRATOOLS_OUT,'ebinning.txt')
    txt_file = open(txt_file_name,'w')
    fits_file = os.path.join(GRATOOLS_OUT,'ebinning.fits')
    for emin, emax in zip(ebinning_array[:-1], ebinning_array[1:]):
        txt_file.write('%.4f %.4f\n'%(emin, emax))
    txt_file.close()
    os.system('gtbindef bintype=E binfile=%s outfile=%s energyunits=MeV' \
                  %(txt_file_name, fits_file))
    logger.info('Created %s...'%fits_file)
    return fits_file
    

def mergeft(path_to_files, out_file_name, N1week, Nnweek):
    """creates a .txt file with the list of the FT files to merge.

       path_to_files: str
           path where datat files are stored
       out_file_name: str
           name of the txt output file (created in the same folder of data)
       N1week: int
           number of the starting week
       Nnweek: int
           number of the ending week
    """
    if N1week < 9:
        abort('Invalid number of weeks: the minimun must be > or = to 9')
    if Nnweek > 434:
        abort('Invalid number of weeks: the maximum must be < or = to 397')
    outtxtfile = os.path.join(path_to_files, out_file_name)
    if not os.path.exists(outtxtfile):
        out_file = open(outtxtfile, 'w')
        for i in range(N1week, Nnweek+1):
            if i == 9:
                out_file.write("%s/lat_photon_weekly_w00%i_p302_v001.fits \n" \
                                   %(path_to_files,i))
            if i >= 10 and i <= 99:
                out_file.write("%s/lat_photon_weekly_w0%i_p302_v001.fits \n" \
                                   %(path_to_files,i))
            if i > 99:
                out_file.write("%s/lat_photon_weekly_w%i_p302_v001.fits \n" \
                                   %(path_to_files,i))
        out_file.close()
    logger.info('Created %s...' %outtxtfile)
    return '@'+outtxtfile
    
def main():
    fore_files_list = [os.path.join(GRATOOLS_CONFIG, \
                                        'fits/gll_iem_v06_hp512_146647.fits'),
                       os.path.join(GRATOOLS_CONFIG, \
                                        'fits/gll_iem_v06_hp512_200561.fits'),
                       os.path.join(GRATOOLS_CONFIG, \
                                        'fits/gll_iem_v06_hp512_274296.fits'),
                       os.path.join(GRATOOLS_CONFIG, \
                                        'fits/gll_iem_v06_hp512_375138.fits'),
                       os.path.join(GRATOOLS_CONFIG, \
                                        'fits/gll_iem_v06_hp512_513056.fits'),
                       ]
    e_min, e_max = 524807.00, 575439.00
    get_foreground_integral_flux_map(fore_files_list, e_min, e_max)

if __name__ == '__main__':
    main()
