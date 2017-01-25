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
from scipy import optimize as opt
from GRATools import GRATOOLS_OUT
from GRATools import GRATOOLS_CONFIG
from GRATools import FT_DATA_FOLDER
from GRATools.utils.logging_ import logger, abort
from GRATools.utils.matplotlib_ import pyplot as plt
from GRATools.utils.gSpline import xInterpolatedUnivariateSplineLinear

FORE_EN = re.compile('\_\d+\.')

def flux2counts(flux_map, exposure_map):
    sr = 4*np.pi/len(flux_map)
    counts_map = flux_map*exposure_map*sr
    return counts_map

def fit_foreground(fore_map, data_map):
    """ATT: maps are intended to be healpix maps (namely numpy arrays)
    """
    nside_out = 64
    _notnull = np.where(data_map > 1e-30)[0]
    mask_f = os.path.join(GRATOOLS_CONFIG, 'fits/Mask64_src2_gp30.fits')
    mask = hp.read_map(mask_f)
    _unmask = np.where(mask > 1e-30)[0]
    logger.info('down grade...')
    fore_repix = np.array(hp.ud_grade(fore_map, nside_out=nside_out))
    data_repix =  np.array(hp.ud_grade(data_map, nside_out=nside_out))   
    A = np.vstack([fore_repix[_unmask], np.ones(len(fore_repix[_unmask]))]).T
    norm, const = np.linalg.lstsq(A, data_repix[_unmask])[0]
    #a, norm = best_fit(fore_map[_notnull], data_map[_notnull])
    logger.info('fit param (norm, const): %.3f, %e' %(norm, const))
    return norm, const

def get_foreground_integral_flux_map(fore_files_list, e_min, e_max):
    """fore_files_list: list of str
           Ordered list of the foreground files (one for each energy)
       e_min: float
           the min of the energy bin
       e_max: float 
           the max of the energy bin
    """
    input_file = os.path.join(FT_DATA_FOLDER, 'models/gll_iem_v06.fits')
    if not os.path.exists(input_file):
        abort("Map %s not found!"%input_file)
    frmaps = pf.open(input_file)
    fore_en = np.array([x[0] for x in frmaps['ENERGIES'].data])
    out_name = fore_files_list[0].replace('_%i.fits'%int(fore_en[0]), 
                                          '_%d-%d.fits'%(e_min, e_max))
    if os.path.exists(out_name):
        logger.info('ATT: file %s already exists and returned...'%out_name)
        fore_map = hp.read_map(out_name)
        return fore_map
    else: 
        logger.info('Computing the integral flux of the foreground model...')
        logger.info('...between %.2f - %.2f'%(e_min, e_max))
        if (fore_en < e_max).any() and (fore_en > e_min).any():
            if fore_en[(fore_en < e_max)&(fore_en > e_min)].size == 0:
                if (fore_en < e_max).all():
                    print 'no foreground model for e_max'
                    if (fore_en < e_min).all():
                        print 'no foreground model for e_min'
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
            else:
                e_int = fore_en[(fore_en < e_max)&(fore_en > e_min)]
                index = np.where(fore_en==e_int)[0]
                if abs(e_min-e_int) <= abs(e_max-e_int):
                    a = index-1
                    fore_e_min = fore_en[a]
                    fore_e_max = fore_en[index]
                    fore_map1 = hp.read_map(fore_files_list[a])
                    fore_map2 = hp.read_map(fore_files_list[index])
                else:
                    a = index+1
                    if a > 29:
                        a = index-1
                        fore_e_min = fore_en[a]
                        fore_e_max = fore_en[index]
                        fore_map1 = hp.read_map(fore_files_list[a])
                        fore_map2 = hp.read_map(fore_files_list[index])
                    else:
                        fore_e_min = fore_en[index]
                        fore_e_max = fore_en[a]
                        fore_map1 = hp.read_map(fore_files_list[index])
                        fore_map2 = hp.read_map(fore_files_list[a])
        elif (fore_en < e_max).all():
            print 'no foreground model for e_max'
            if (fore_en < e_min).all():
                print 'no foreground model for e_min'
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
        #fore_integr = (A/2*(e_max**2-e_min**2)+B*(e_max-e_min))
        #print 'e mean:', np.sqrt((e_max*e_min)), 'delta E:', e_max-e_min
        fore_integr1 = (A*e_max + B)
        fore_integr2 = (A*e_min + B)
        fore_integr = np.sqrt(fore_integr1*fore_integr2)*(e_max-e_min)
        hp.write_map(out_name, fore_integr)
        return fore_integr

def csi_parse(csi_file):
    """Parsing of the *_csi.txt files
    """
    logger.info('loading Csi values from %s'%csi_file)
    csi = []
    theta = []
    Rs = []
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
        if 'R\t' in line:
            r = np.array([float(item) for item in line.split()[1:]])
            Rs.append(r)
    f.close()
    return np.array(emin), np.array(emax), np.array(emean), np.array(csi), np.array(theta), np.array(Rs)

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

def clEcross_parse(cl_file):
    """Parsing of the *_cps.txt files
    """
    logger.info('loading Cp values from %s'%cl_file)
    ff = open(cl_file, 'r')
    _emin1, _emax1, _emean1, _cls, _clserr = [], [], [], [], []
    _emin2, _emax2, _emean2 = [], [], []
    for line in ff:
        if 'ENERGY1\t' in line:
            emin1, emax1, emean1 = [float(item) for item in line.split()[1:]]
            _emin1.append(emin1)
            _emax1.append(emax1)
            _emean1.append(emean1)
        if 'ENERGY2\t' in line:
            emin2, emax2, emean2 = [float(item) for item in line.split()[1:]]
            _emin2.append(emin2)
            _emax2.append(emax2)
            _emean2.append(emean2)
        if 'Cl\t' in line:
            cl = np.array([float(item) for item in line.split()[1:]])
            _cls.append(cl)
        if 'Cl_ERR\t' in line:
            clerr = np.array([float(item) for item in line.split()[1:]])
            _clserr.append(clerr)
    ff.close()
    return np.array(_emin1), np.array(_emax1),np.array(_emean1), \
        np.array(_emin2), np.array(_emax2),np.array(_emean2), \
        np.array(_cls), np.array(_clserr)

def clfore_parse(clfore_file):
    """Parsing of the *_forecls.txt files.
    """
    cls = []
    emin, emax = [], []
    f = open(clfore_file, 'r')
    for line in f:
        if 'ENERGY\t' in line:
            e1, e2 = [float(item) for item in line.split()[1:]]
            emin.append(e1)
            emax.append(e2)
        if 'Cl\t' in line:
            cl = np.array([float(item) for item in line.split()[1:]])
            cls.append(cl)
    f.close()
    return np.array(emin), np.array(emax), np.array(cls)

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
        
def get_crbkg(txt_file):
    """Get the CR residual bkg (spline) as a function of the energy
       from the txt files
    """
    logger.info('Getting CR residual bkg from file %s'%txt_file)
    f = open(txt_file, 'r')
    _bkg, _en = [], []
    for line in f:
        try:
            e, bkg = [float(item) for item in line.split()]
            _en.append(e)
            _bkg.append(bkg)
        except:
            pass 
    fmt = dict(xname='Energy', xunits='MeV', 
               yname='E$^{2}$ x CR Residual flux', 
               yunits='MeV cm$^{-2}$ s$^{-1}$ sr$^{-1}$')
    crbkg = xInterpolatedUnivariateSplineLinear(np.array(_en), np.array(_bkg),\
                                                  **fmt)
    f.close()
    return crbkg

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
    
def best_fit(X, Y):

    xbar = sum(X)/len(X)
    ybar = sum(Y)/len(Y)
    n = len(X) # or len(Y)

    numer = sum([xi*yi for xi,yi in zip(X, Y)]) - n * xbar * ybar
    denum = sum([xi**2 for xi in X]) - n * xbar**2

    b = numer / denum
    a = ybar - b * xbar

    print('best fit line:\ny = {:.2f} + {:.2f}x'.format(a, b))

    return a, b

def main():
    """
    crbkg_PSF1_file =  os.path.join(GRATOOLS_CONFIG, 
                                    'ascii/p8r2_bkg_flux_t8.txt')
    crbkg_PSF1 = get_crbkg(crbkg_PSF1_file)
    crbkg_PSF1.plot(show=False)
    x = crbkg_PSF1.x
    xx = x**2
    print x
    y1 =  crbkg_PSF1.y
    crbkg_PSF2_file =  os.path.join(GRATOOLS_CONFIG, 
                                    'ascii/p8r2_bkg_flux_t16.txt')
    crbkg_PSF2 = get_crbkg(crbkg_PSF2_file)
    crbkg_PSF2.plot(show=False)
    y2 =  crbkg_PSF2(x)
    crbkg_PSF3_file =  os.path.join(GRATOOLS_CONFIG, 
                                    'ascii/p8r2_bkg_flux_t32.txt')
    crbkg_PSF3 = get_crbkg(crbkg_PSF3_file)
    crbkg_PSF3.plot(show=False)
    y3 =  crbkg_PSF2(x)
    y = (y1/x+y2/x+y3/x)*x
    print y
    plt.yscale('log')
    plt.xscale('log')
    plt.show()
    
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
    get_foreground_integral_flux_map(fore_files_list, e_min, e_max)"""
    from scipy.optimize import curve_fit
    ClFORE_FILE = os.path.join(GRATOOLS_OUT, 'fore_polspicecls.txt')
    min1, emax1, clsfore = clfore_parse(ClFORE_FILE)
    l = np.arange(len(clsfore[0]))
    clsfore_spline = xInterpolatedUnivariateSplineLinear(l, clsfore[0])
    clsfore_spline.plot(show=False)
    plt.xscale('log')
    plt.yscale('log')
    plt.show()


if __name__ == '__main__':
    main()
