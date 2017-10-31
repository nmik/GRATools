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
from numba import jit
from scipy import optimize as opt
from itertools import product
from scipy.misc import factorial as fact
from GRATools import GRATOOLS_OUT
from GRATools import GRATOOLS_CONFIG
from GRATools import FT_DATA_FOLDER
from GRATools.utils.logging_ import logger, abort
from GRATools.utils.matplotlib_ import pyplot as plt
from GRATools.utils.gSpline import xInterpolatedUnivariateSplineLinear
from GRATools.utils.gSpline import xInterpolatedBivariateSplineLinear

FORE_EN = re.compile('\_\d+\.')


def get_ref_igrb_spline():
    """Returns a spline of the IGRB measurement as a funcion of the energy.
       this curve is used as initial guess of the constant parameter of the
       poissonian fit.
    """
    f = open(os.path.join(GRATOOLS_CONFIG,'ascii/IGRB_guess.txt'), 'r')
    x = np.array([float(l.split()[0]) for l in f])
    f = open(os.path.join(GRATOOLS_CONFIG,'ascii/IGRB_guess.txt'), 'r')
    y = np.array([float(l.split()[1]) for l in f])
    fmt = dict(xname='Energy' , xunits='MeV', yname='IGRB', 
               yunits='MeV$^{-1}$/cm$^{-2}$/s')
    igrb = xInterpolatedUnivariateSplineLinear(x, y, **fmt)
    return igrb

def flux2counts(flux_map, exposure_map):
    """Returns a map of counts given a flux map and an exposure map
       ATT: maps are intended to be healpix maps (namely numpy arrays)

       flux_map: numpy array
          healpy flux map
       exposure_map: numpy array
          healpy exposure map
    """
    sr = 4*np.pi/len(flux_map)
    counts_map = flux_map*exposure_map*sr
    return counts_map

def fit_foreground_lstsq(fore_map, data_map, mask_map):
    """Perform the gaussian fit (least square method) given a data map and
       a model map, with this expression: data = A + B*model; returns A and B,
       fit parametres.
       ATT: maps are intended to be healpix maps (namely numpy arrays)
       
       fore_map: numpy array
          healpy model map
       data_map: numpy array
          healpy data map
    """
    nside_out = 64
    _notnull = np.where(data_map > 1e-30)[0]
    #mask_f = os.path.join(GRATOOLS_CONFIG, 'fits/Mask64_src2_gp15.fits')
    #mask = hp.read_map(mask_f)
    _unmask = np.where(mask_map > 1e-30)[0]
    #logger.info('down grade...')
    #fore_repix = np.array(hp.ud_grade(fore_map, nside_out=nside_out))
    #data_repix = np.array(hp.ud_grade(data_map, nside_out=nside_out))
    A = np.vstack([fore_map[_unmask], np.ones(len(fore_map[_unmask]))]).T
    norm, const = np.linalg.lstsq(A, data_map[_unmask])[0]
    logger.info('fit param (norm, const): %.3f, %e' %(norm, const))
    return norm, const

@jit
def poisson_likelihood_2(norm1_guess, norm2_guess, const_guess, fore_map, 
                         data_map, srctempl_map, exp=None, sr=None):
    """Compute the log-likelihood as decribed here: 
       http://iopscience.iop.org/article/10.1088/0004-637X/750/1/3/pdf
       where the model to fit to data is given by norm*fore_map+const.

       norm1_guess : float
          initial guess for normalization parameter
       norm2_guess : float
          initial guess for normalization parameter
       const_guess : float
          initial guess for constant parameter
       fore_map : numpy array
          helapix map of foreground model
       data_map : numpy array
          helapix map of data. It could be either a count map or a flux map.
          If a counts map is given, an exposure map should be given too. See 
          next parameter.
       srctempl_map : numpy array
          helapix map of source template
       exp :  numpy array or None
          helapix map of the exposure. Should be given if the data map is in 
          counts (beacause foreground map is in flux units by default and it 
          needs to be turned to counts to be fitted). While, If data map is 
          in flux units, do not declare this parameter, which is None by 
          default.
       sr : float or None
          pixel area -> 4*pi/Npix
    """
    a = norm1_guess
    b = const_guess
    c = norm2_guess
    s = srctempl_map
    factorial_data = fact(data_map)
    lh = 0
    if exp is not None:
        for i, f in enumerate(fore_map):
            lh += (a*f+c*s[i]+b)*exp[i]*sr + np.log(factorial_data[i]) - \
                data_map[i]*np.log((a*f+c*s[i]+b)*exp[i]*sr)
    else:
        for i, f in enumerate(fore_map):
            lh += np.sum(((a*f+c*s[i]+b)+np.log(factorial_data[i]) -\
                              data_map[i]*np.log((a*f+c*s[i]+b))))
    return lh

@jit
def poisson_likelihood(norm_guess, const_guess, fore_map, data_map, exp=None, 
                       sr=None):
    """Compute the log-likelihood as decribed here: 
       http://iopscience.iop.org/article/10.1088/0004-637X/750/1/3/pdf
       where the model to fit to data is given by norm*fore_map+const.

       norm_guess : float
          initial guess for normalization parameter
       const_guess : float
          initial guess for constant parameter
       fore_map : numpy array
          helapix map of foreground model
       data_map : numpy array
          helapix map of data. It could be either a count map or a flux map.
          If a counts map is given, an exposure map should be given too. See 
          next parameter.
       exp :  numpy array or None
          helapix map of the exposure. Should be given if the data map is in 
          counts (beacause foreground map is in flux units by default and it 
          needs to be turned to counts to be fitted). While, If data map is 
          in flux units, do not declare this parameter, which is None by 
          default.
       sr : float or None
          pixel area -> 4*pi/Npix
    """
    a = norm_guess
    b = const_guess
    factorial_data = fact(data_map)
    lh = 0
    if exp is not None:
        for i, f in enumerate(fore_map):
            lh += (a*f+b)*exp[i]*sr + np.log(factorial_data[i]) - \
                data_map[i]*np.log((a*f+b)*exp[i]*sr)
    else:
        for i, f in enumerate(fore_map):
            lh += np.sum(((a*f+b)+np.log(factorial_data[i]) -\
                              data_map[i]*np.log((a*f+b))))
    return lh

def fit_fore_src_poisson(fore_map, data_map, srctempl_map, n1_guess=1., 
                         c_guess=0.1, n2_guess=1., exp=None, smooth=False, 
                         show=False):
    """Performs the poisonian fit, recursively computing the log likelihood 
       (using poisson_likelihood) for a grid of values of fit parameters around
       the guess. Returns the values of parameters which minimize the log 
       likelihood, togather to the 1-sigma error

       n1_guess : float
          initial guess for normalization parameter
       n1_guess : float
          initial guess for normalization parameter
       c_guess : float
          initial guess for constant parameter
       fore_map : numpy array
          helapix map of foreground model
       data_map : numpy array
          helapix map of data. It could be either a count map or a flux map.
          If a counts map is given, an exposure map should be given too. See 
          next parameter.
          exp :  numpy array or None
          helapix map of the exposure. Should be given if the data map is in 
          counts (beacause foreground map is in flux units by default and it 
          needs to be turned to counts to be fitted). While, If data map is 
          in flux units, do not declare this parameter, which is None by 
          default.
       smooth : bool
          not implemented yet...
       show : bool
          if true it shows some usefull plot to check if the fit is functioning
    """
    #show=True
    logger.info('Performing poissonian fit...')
    norm1_guess = n1_guess
    norm2_guess = n2_guess
    igrb_guess = c_guess
    nside_out = 64
    mask_f = os.path.join(GRATOOLS_CONFIG, 'fits/Mask64_src2_gp30.fits')
    mask = hp.read_map(mask_f)
    _unmask = np.where(mask > 1e-30)[0]
    logger.info('down grade...')
    fore_repix = np.array(hp.ud_grade(fore_map, nside_out=nside_out))
    data_repix =  np.array(hp.ud_grade(data_map, nside_out=nside_out, 
                                       power=-2))
    srct_repix = np.array(hp.ud_grade(srctempl_map, nside_out=nside_out))
    norm1_list = np.linspace(norm1_guess*0.3, norm1_guess*1.5, 21)
    norm2_list = np.linspace(norm2_guess*0.3, norm2_guess*1.5, 21)
    igrb_list = np.linspace(igrb_guess*0.01, igrb_guess*10., 101)
    logger.info('Minimization likelihood run1...')
    lh_list = []
    combinations = list(product(norm1_list, norm2_list, igrb_list))
    if exp is not None:
        exposure = exp
        exposure = np.array(hp.ud_grade(exposure, nside_out=nside_out))
        areapix = 4*np.pi/(len(data_repix))
        for i,j,k in combinations: 
            lh = poisson_likelihood_2(i, j, k, fore_repix[_unmask], 
                                    data_repix[_unmask], 
                                    srct_repix[_unmask],
                                    exp=exposure[_unmask],
                                    sr=areapix)
            lh_list.append(lh)
    else:
        for i,j,k in combinations:
            lh = poisson_likelihood_2(i, j, k, fore_repix[_unmask],
                                    data_repix[_unmask], srct_repix[_unmasked])
            lh_list.append(lh)
    lh_min = np.argmin(np.array(lh_list))
    (norm1_min, norm2_min, igrb_min) = combinations[lh_min]
    logger.info('Run1 results: n1=%.3f n2=%.3f c=%e'%(norm1_min, norm2_min, 
                                                      igrb_min))
    norm1_list = np.linspace(norm1_min*0.7, norm1_min*1.2, 21)
    norm2_list = np.linspace(norm2_min*0.7, norm2_min*1.2, 21)   
    igrb_list = np.linspace(igrb_min*0.5, igrb_min*1.5, 101)
    logger.info('Minimization likelihood run2...')
    lh_list = []
    combinations = np.array(list(product(norm1_list, norm2_list, igrb_list)))
    if exp is not None:
        exposure = exp
        exposure = np.array(hp.ud_grade(exposure, nside_out=nside_out))
        areapix = 4*np.pi/(len(data_repix))
        for i,j,k in product(norm1_list, norm2_list, igrb_list):
            lh = poisson_likelihood_2(i, j, k, fore_repix[_unmask],
                                    data_repix[_unmask],
                                    srct_repix[_unmask],
                                    exp=exposure[_unmask],
                                    sr=areapix)
            lh_list.append(lh)
    else:
        for i,j,k in product(norm1_list, norm2_list, igrb_list):
            lh = poisson_likelihood_2(i, j, k, fore_repix[_unmask],
                                    data_repix[_unmask],srct_repix[_unmask])
            lh_list.append(lh)
    lh_list = np.array(lh_list)
    lh_min = np.argmin(lh_list)
    (norm1_min, norm2_min, igrb_min) = combinations[lh_min]
    lh_delta = np.array(lh_list)[lh_min]+2.3
    index = np.where(np.array(lh_list) < lh_delta)[0]
    _norm1 = np.array([x[0] for x in combinations[index]])
    logger.info('Norm1 err: %.4f - %.4f'%(_norm1[0], _norm1[-1]))
    n1_err = (_norm1[0], _norm1[-1])
    _norm2 = np.array([x[0] for x in combinations[index]])
    logger.info('Norm2 err: %.4f - %.4f'%(_norm2[0], _norm2[-1]))
    n2_err = (_norm2[0], _norm2[-1])
    _igrb = np.array([x[1] for x in combinations[index]])
    logger.info('Igrb err: %.e - %.e'%(np.amin(_igrb), np.amax(_igrb)))
    igrb_err = (np.amin(_igrb), np.amax(_igrb))
    return norm1_min, norm2_min, igrb_min, n1_err, n2_err, igrb_err

def fit_foreground_poisson(fore_map, data_map, n_guess=1., c_guess=0.1,
                           exp=None, smooth=False, show=False):
    """Performs the poisonian fit, recursively computing the log likelihood 
       (using poisson_likelihood) for a grid of values of fit parameters around
       the guess. Returns the values of parameters which minimize the log 
       likelihood, togather to the 1-sigma error

       n_guess : float
          initial guess for normalization parameter
       c_guess : float
          initial guess for constant parameter
       fore_map : numpy array
          helapix map of foreground model
       data_map : numpy array
          helapix map of data. It could be either a count map or a flux map.
          If a counts map is given, an exposure map should be given too. See 
          next parameter.
          exp :  numpy array or None
          helapix map of the exposure. Should be given if the data map is in 
          counts (beacause foreground map is in flux units by default and it 
          needs to be turned to counts to be fitted). While, If data map is 
          in flux units, do not declare this parameter, which is None by 
          default.
       smooth : bool
          not implemented yet...
       show : bool
          if true it shows some usefull plot to check if the fit is functioning
    """
    #show=True
    logger.info('Performing poissonian fit...')
    norm_guess = n_guess
    igrb_guess = c_guess
    nside_out = 64
    mask_f = os.path.join(GRATOOLS_CONFIG, 'fits/Mask64_src2_gp30.fits')
    mask = hp.read_map(mask_f)
    _unmask = np.where(mask > 1e-30)[0]
    logger.info('down grade...')
    fore_repix = np.array(hp.ud_grade(fore_map, nside_out=nside_out))
    data_repix =  np.array(hp.ud_grade(data_map, nside_out=nside_out, 
                                       power=-2))
    norm_list = np.linspace(norm_guess*0.3, norm_guess*1.5, 50)
    igrb_list = np.linspace(igrb_guess*0.01, igrb_guess*10., 200)
    logger.info('Minimization likelihood run1...')
    lh_list = []
    combinations = list(product(norm_list, igrb_list))
    if exp is not None:
        exposure = exp
        exposure = np.array(hp.ud_grade(exposure, nside_out=nside_out))
        areapix = 4*np.pi/(len(data_repix))
        for i,j in product(norm_list, igrb_list): 
            lh = poisson_likelihood(i, j, fore_repix[_unmask], 
                                    data_repix[_unmask], 
                                    exp=exposure[_unmask],
                                    sr=areapix)
            lh_list.append(lh)
    else:
        for i,j in product(norm_list, igrb_list):
            lh = poisson_likelihood(i, j, fore_repix[_unmask],
                                    data_repix[_unmask])
            lh_list.append(lh)
    lh_min = np.argmin(np.array(lh_list))
    (norm_min, igrb_min) = combinations[lh_min]
    logger.info('Run1 results: n=%.3f c=%e'%(norm_min, igrb_min))
    norm_list = np.linspace(norm_min*0.7, norm_min*1.2, 51)   
    igrb_list = np.linspace(igrb_min*0.5, igrb_min*1.5, 101)
    logger.info('Minimization likelihood run2...')
    lh_list = []
    combinations = np.array(list(product(norm_list, igrb_list)))
    if exp is not None:
        exposure = exp
        exposure = np.array(hp.ud_grade(exposure, nside_out=nside_out))
        areapix = 4*np.pi/(len(data_repix))
        for i,j in product(norm_list, igrb_list):
            lh = poisson_likelihood(i, j, fore_repix[_unmask],
                                    data_repix[_unmask],
                                    exp=exposure[_unmask],
                                    sr=areapix)
            lh_list.append(lh)
    else:
        for i,j in product(norm_list, igrb_list):
            lh = poisson_likelihood(i, j, fore_repix[_unmask],
                                    data_repix[_unmask])
            lh_list.append(lh)
    lh_list = np.array(lh_list)
    lh_min = np.argmin(lh_list)
    (norm_min, igrb_min) = combinations[lh_min]
    lh_delta = np.array(lh_list)[lh_min]+2.3
    index = np.where(np.array(lh_list) < lh_delta)[0]
    _norm = np.array([x[0] for x in combinations[index]])
    logger.info('Norm err: %.4f - %.4f'%(_norm[0], _norm[-1]))
    _igrb = np.array([x[1] for x in combinations[index]])
    logger.info('Igrb err: %.e - %.e'%(np.amin(_igrb), np.amax(_igrb)))
    
    if show == True:
        n = np.array([x[0] for x in combinations])
        plt.figure(facecolor='white')
        plt.plot(n, lh_list, 'o', color='coral', alpha=0.3)
        plt.plot(norm_min, lh_list[lh_min] , 'r*')
        plt.plot([_norm[0], _norm[-1]], [lh_delta, lh_delta], 'r-')
        plt.xlabel('Normalization')
        plt.ylabel('-Log(Likelihood)')
        
        igrb = np.array([x[1] for x in combinations])
        plt.figure(facecolor='white')
        plt.plot(igrb, lh_list, 'o', color='coral', alpha=0.3)
        plt.plot(igrb_min, lh_list[lh_min] , 'r*')
        plt.plot([np.amin(_igrb), np.amax(_igrb)], [lh_delta, lh_delta], 'r-')
        plt.xlabel('Constant')
        plt.ylabel('-Log(Likelihood)')

        fig = plt.figure(facecolor='white')
        z = lh_list
        zmin = lh_list[lh_min]
        z.shape = (len(norm_list), len(igrb_list))
        ax = fig.add_subplot(111)
        cax = ax.matshow(z, origin='lower', cmap='Spectral',
                    aspect='auto')
        plt.xlabel('$C$ $[cm^{-2}s^{-1}sr^{-1}]$')
        plt.ylabel('$N$')
        #plt.title('$\Phi_{data}=N\cdot\Phi_{model}+C$')
        x_ticks = np.linspace(np.amin(igrb_list), np.amax(igrb_list), 6)
        formatting_function = np.vectorize(lambda f: format(f, '6.1E'))
        x_ticks = list(formatting_function(x_ticks))
        y_ticks = list(np.around(np.linspace(np.amin(norm_list), 
                                             np.amax(norm_list), 6), 
                                 decimals=3))
        ax.set_yticklabels(['']+y_ticks)
        ax.set_xticklabels(['']+x_ticks)
        ax.xaxis.set_ticks_position('bottom')
        cb = plt.colorbar(cax, format='$%.1e$')
        cb.set_label('-Log(Likelihood)', rotation=90)
        norm_min_ind = list(norm_list).index(norm_min)
        igrb_min_ind = list(igrb_list).index(igrb_min)
        _norm_ind = []
        _igrb_ind = []
        for i in range(0, len(index)):
            _norm_ind.append(list(norm_list).index(_norm[i]))
            _igrb_ind.append(list(igrb_list).index(_igrb[i]))
        _norm_ind = np.array(_norm_ind)
        _igrb_ind = np.array(_igrb_ind)
        plt.contourf(z, [zmin, zmin+2.3, zmin+4.61, zmin+5.99], 
                     colors='w', origin='lower', alpha=0.3)
        plt.scatter(igrb_min_ind, norm_min_ind, s=45, c='w', marker='+')
        plt.show()
    return norm_min, igrb_min, _norm[0], _norm[-1], np.amin(_igrb), \
        np.amax(_igrb)

def get_iso_integral_flux_map(iso_ascii_file, e_min, e_max, nside=512):
    """Returns a isotropic map of the integral IGRB between emin and emax from 
       'iso_P8R2_SOURCE_V6_v06.txt' file.

       iso_ascii_file: str
           Ascii file found here 
           https://fermi.gsfc.nasa.gov/ssc/data/access/lat/BackgroundModels.html
       emin: float
           minimum energy of the bin [MeV]
       emax: float
           maximum energy of the bin [MeV]
       nside: int
           healpix nside parameter
    """
    isofile = os.path.join(GRATOOLS_CONFIG, 'models', 
                           'iso_P8R2_ULTRACLEANVETO_V6_v06.txt')
    from GRATools.utils.gFTools import iso_parse
    e, difflux, diffluxerr = iso_parse(iso_ascii_file)
    index = np.where((e>e_min)*(e<e_max))
    erange = e[index]
    frange = difflux[index]
    f_e = xInterpolatedUnivariateSplineLinear(erange, frange)
    intf = f_e.integral(e_min, e_max)
    logger.info('Isotropic Bkg %e [cm-2s-1]'%intf)
    npix = hp.nside2npix(nside)
    intiso_map = np.full(npix, intf)
    return intiso_map

def find_outer_energies(en_val, en_arr):
    """Returns the first element on the right and the first on the left
       of a given value (en_val), among all values in an ordered array
       (en_arr).
       
       en_val : float
           mean energy
       en_arr : float
           array of the energies at which the foreground model is given.
    """
    en_sx_arr = en_arr[en_arr < en_val]
    en_dx_arr = en_arr[en_arr > en_val]
    if en_sx_arr.size == 0:
        en_sx = en_dx_arr[0]
        en_dx = en_dx_arr[1]
    elif en_dx_arr.size == 0:
        en_sx = en_sx_arr[-2]
        en_dx = en_sx_arr[-1]
    else:
        en_sx = en_sx_arr[-1]
        en_dx = en_dx_arr[0]
    return en_sx, en_dx

def get_fore_integral_flux_map(fore_files_list, e_min, e_max):
    """Returns the foreground map integrated between e_min and e_max
       A powerlaw is assumed fore the foregriunf energy spectrum, hence
       the interpolation between 2 given maps at given energies (given 
       by the model) is done in logarithmic scales. 
    
       fore_files_list: list of str
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
    fore_en = []#np.array([x[0] for x in frmaps['ENERGIES'].data])
    for ff in fore_files_list:
        m = re.search(FORE_EN, ff)
        en = int(m.group(0).replace('_', '').replace('.', ''))
        fore_en.append(en)
    fore_en = np.array(fore_en)
    out_name = fore_files_list[0].replace('_%i.fits'%fore_en[0], 
                                          '_%d-%d.fits'%(e_min, e_max))
    if os.path.exists(out_name):
        logger.info('ATT: file %s already exists and returned...'%out_name)
        fore_map = hp.read_map(out_name)
        return fore_map
    else: 
        logger.info('Computing the integral flux of the foreground model...')
        logger.info('...between %.2f - %.2f'%(e_min, e_max))
        fore_emin_sx, fore_emin_dx = find_outer_energies(e_min, fore_en)
        fore_emax_sx, fore_emax_dx = find_outer_energies(e_max, fore_en)
        fore_emin_sx_ind = np.where(fore_en == fore_emin_sx)[0]
        fore_emin_dx_ind = np.where(fore_en == fore_emin_dx)[0]
        fore_emax_sx_ind = np.where(fore_en == fore_emax_sx)[0]
        fore_emax_dx_ind = np.where(fore_en == fore_emax_dx)[0]
        fore_fmin_sx = hp.read_map(fore_files_list[fore_emin_sx_ind])
        fore_fmin_dx = hp.read_map(fore_files_list[fore_emin_dx_ind])
        fore_fmax_sx = hp.read_map(fore_files_list[fore_emax_sx_ind])
        fore_fmax_dx = hp.read_map(fore_files_list[fore_emax_dx_ind])
        m1 = (np.log10(fore_fmin_sx)-np.log10(fore_fmin_dx))/ \
            (np.log10(fore_emin_sx)-np.log10(fore_emin_dx))
        m2 = (np.log10(fore_fmax_sx)-np.log10(fore_fmax_dx))/ \
            (np.log10(fore_emax_sx)-np.log10(fore_emax_dx))
        logfore1 = m1*(np.log10(e_min)-np.log10(fore_emin_sx))+ \
            np.log10(fore_fmin_sx)
        logfore2 = m2*(np.log10(e_max)-np.log10(fore_emax_sx))+ \
            np.log10(fore_fmax_sx)
        fore1 = 10**(logfore1)
        fore2 = 10**(logfore2)
        fore_integ = np.sqrt(fore1*fore2)*(e_max - e_min)
        hp.write_map(out_name, fore_integ)
        logger.info('Created file %s'%out_name)
        return fore_integ

def main():
    """Test session
    """
    logger.info('No test module is available at the moment... bye bye!')
    return 0
    
    


if __name__ == '__main__':
    main()
