#!/usr/bin/env python                                                          #
#                                                                              #
# Autor: Michela Negro, University of Torino.                                  #
# On behalf of the Fermi-LAT Collaboration.                                    #
#                                                                              #
# This program is free software; you can redistribute it and/or modify         #
# it under the terms of the GNU General Public License as published by         #
# the Free Software Foundation; either version 3 of the License, or            #
# (at your option) any later version.                                          #
#                                                                              #
#------------------------------------------------------------------------------#


"""Functions to define map masks
"""

import os
import healpy as hp
import numpy as np
import pyfits as pf
import pickle as pc
from numba import jit
from GRATools import GRATOOLS_CONFIG, GRATOOLS_OUT
from GRATools.utils.logging_ import logger, abort
from GRATools.utils.matplotlib_ import pyplot as plt
from GRATools.utils.gWindowFunc import get_psf
from GRATools.utils.gSpline import xInterpolatedUnivariateSplineLinear

_ENERGY = np.logspace(1, 6,5000)

def LogParabola(norm, index1, index2, pivot_energy):
    """Returns a spline of a LogParabola spectrum with given parameters
    """
    _spec = norm*(_ENERGY/pivot_energy)**(-index1-index2*np.log( \
            _ENERGY/pivot_energy))
    _spec_spline = xInterpolatedUnivariateSplineLinear(_ENERGY, _spec)
    return _spec_spline

def PLSuperExpCutoff(norm, index, pivot_energy, index_cut, energy_cut):
    """Returns a spline of a Powerlaw with exponential cutoff spectrum with 
       given parameters
    """
    _spec = norm*(_ENERGY/pivot_energy)**(-index)*np.exp((pivot_energy/ \
             energy_cut)**(-index_cut)-(_ENERGY/pivot_energy)**(-index_cut))
    _spec_spline = xInterpolatedUnivariateSplineLinear(_ENERGY, _spec)
    return _spec_spline

def PLSuperExpCutoff2(norm, index, pivot_energy, index_cut, exp_fact):
    """Returns a spline of a Powerlaw with exponential cutoff spectrum with 
       given parameters
    """
    _spec = norm*(_ENERGY/pivot_energy)**(-index)*\
        np.exp(exp_fact*(pivot_energy**(-index_cut)-_ENERGY**(-index_cut)))
    _spec_spline = xInterpolatedUnivariateSplineLinear(_ENERGY, _spec)
    return _spec_spline

def PowerLaw(norm, index, pivot_energy):
    """Returns a spline of a Powerlaw spectrum with given parameters
    """
    _spec = norm*(_ENERGY/pivot_energy)**(-index)
    _spec_spline = xInterpolatedUnivariateSplineLinear(_ENERGY, _spec)
    return _spec_spline

def build_src_template(cat_file, psf_file, emin, emax, b_cut=20, nside=512):
    """Returns a map of sources inside the catalog smoothed by the PSF 
       at the gven energy.

       cat_file: str
           .fits fits file with the sorce catalog
       psf_file: str
           ascii file with 2 columns (Energy - PSF)
       emin: float
           minimum energy of the bin [MeV]
       emax: float
           maximum energy of the bin [MeV]
       b_cut: float
           cut on abs(b) --> only sources with abs(b)>b_cut will be
           taken
       nside: int
           healpix nside parameter
    """
    psf_th_e = get_psf(psf_file)
    emean = np.sqrt(emin*emax)
    psf_th = psf_th_e.vslice(emean)
    th = np.linspace(psf_th.x[0], psf_th.x[-1], 1000)
    sin_th = np.sin(th)
    psfXsin_th = sin_th*psf_th(th)
    psfXsin_th_spline = xInterpolatedUnivariateSplineLinear(th, psfXsin_th)
    psfXsin_th_spline = psfXsin_th_spline.scale(2*np.pi)
    integral_psf_norm = psfXsin_th_spline.integral(psf_th.x[0], psf_th.x[-1])
    logger.info('Integral PSF shaped curve = %f'%integral_psf_norm)
    ##########################################################
    src_cat = pf.open(cat_file)
    CAT = src_cat['LAT_Point_Source_Catalog']
    SOURCES = CAT.data
    GLAT = SOURCES.field('GLAT')
    logger.info('Considering |b| > %.1f'%b_cut)
    lat_index = np.where(abs(GLAT)>b_cut)[0]
    logger.info('Number of 3FGL sources: %i'%len(lat_index))
    GLAT =  GLAT[lat_index]
    GLON = SOURCES.field('GLON')[lat_index]
    FLUX1000 = SOURCES.field('Flux1000')[lat_index]
    spec_type = SOURCES.field('SpectrumType')[lat_index]
    E0 = SOURCES.field('Pivot_Energy')[lat_index]
    K = SOURCES.field('Flux_Density')[lat_index] #cm-2 s-1 MeV-1
    #alpha = GAMMA = SOURCES.field('Spectral_Index')[lat_index] #3FGL
    alpha = SOURCES.field('PL_Index')[lat_index]
    #beta = SOURCES.field('beta')[lat_index] #3FGL
    alpha_lp = SOURCES.field('LP_Index')[lat_index]
    beta = SOURCES.field('LP_beta')[lat_index]
    GAMMA_plec = SOURCES.field('PLEC_Index')[lat_index]
    b = SOURCES.field('PLEC_Exp_Index')[lat_index]
    a = SOURCES.field('PLEC_Expfactor')[lat_index]
    #Ecut = SOURCES.field('Cutoff')[lat_index]#3FGL, discarded in FL8Y!!
    src_name = SOURCES.field('Source_Name')[lat_index]
    src_cat.close()
    ##########################################################
    INT_FLUXES = [] #cm-2 s-1 sr-1
    for i, s in enumerate(spec_type):
        spec = 0
        if 'PowerLaw' in s:
            spec = PowerLaw(K[i], alpha[i], E0[i])
            integral_flux_bin = spec.integral(emin, emax)
            INT_FLUXES.append(integral_flux_bin)
        elif 'LogParabola' in s:
            spec = LogParabola(K[i], alpha_lp[i], beta[i], E0[i])
            integral_flux_bin = spec.integral(emin, emax)
            INT_FLUXES.append(integral_flux_bin)
        elif 'PLSuperExpCutoff2' in s:
            spec = PLSuperExpCutoff2(K[i], GAMMA_plec[i], E0[i], b[i], a[i])
            integral_flux_bin = spec.integral(emin, emax)
            INT_FLUXES.append(integral_flux_bin)
        elif 'PLSuperExpCutoff' in s:
            spec = PLSuperExpCutoff(K[i], alpha[i], E0[i], b[i], Ecut[i])
            integral_flux_bin = spec.integral(emin, emax)
            INT_FLUXES.append(integral_flux_bin)
    logger.info('Building...')
    src_map = src_builder(psfXsin_th_spline, nside, GLON, GLAT, INT_FLUXES)
    INT_FLUXES = np.array(INT_FLUXES)
    logger.info('Done!')
    out_srctempl_folder = os.path.join(GRATOOLS_OUT, 'output_src')
    if not os.path.exists(out_srctempl_folder):
        os.makedirs(out_srctempl_folder)
    cat = os.path.basename(cat_file).replace('.fit','')
    out_name_srctempl = os.path.join(out_srctempl_folder,'src%s_%i-%i.fits'\
                                             %(cat, emin, emax))
    hp.write_map(out_name_srctempl, src_map)
    logger.info('Created: %s'%out_name_srctempl)
    return src_map

@jit
def src_builder(psf_spline, nside, glon_, glat_, intflux_): 
    """Loop on every source and build the template of surces
    
       psf_spline: spline 
           spline of the PSF as a function of theta (given by Fermi IRFs)
       nside: int
           healpix nside parameter
       glon_: numpy array
           array of all the galactic longitutes of the sources
       glat_: numpy array
           array of all the galactic latitudes of the sources
       intflux_: numpy array
           array of all the total integral flux (in the considered energy
           bin) of the sources
    """
    x, y, z = hp.rotator.dir2vec(glon_, glat_, lonlat=True)
    src_pix = hp.pixelfunc.vec2pix(nside, x, y, z)
    isomap = np.zeros(hp.nside2npix(nside))
    areapix = 4*np.pi/len(isomap)
    for i, bn in enumerate(src_pix):
        pixdir1 = hp.pixelfunc.pix2ang(nside, bn)
        radintpix = hp.query_disc(nside, (x[i],y[i],z[i]), np.radians(9.5))
        pixdir2 = hp.pixelfunc.pix2ang(nside, radintpix)
        dist = hp.rotator.angdist(pixdir1, pixdir2)
        src_profile = psf_spline.scale(intflux_[i]/len(radintpix))
        pix_flux_values = src_profile(dist)
        isomap[radintpix] = isomap[radintpix] + pix_flux_values
    srcfluxtempl = isomap/areapix
    return srcfluxtempl

def main():
    """Simple test unit
    """
    sourcetemp = True
    
    if sourcetemp == True:
        abslatitude = 0
        out_srctempl_folder = os.path.join(GRATOOLS_OUT, 'output_src')
        cat_file = os.path.join(GRATOOLS_CONFIG,'catalogs/gll_psc_8year_v6.fit')
        cat = os.path.basename(cat_file).replace('.fit','')
        psf56_file = os.path.join(GRATOOLS_OUT, 
                                'psf_P8R2_ULTRACLEANVETO_V6_56.fits')
        mask_list = [os.path.join(GRATOOLS_CONFIG, 'fits/Mask_gp20.fits'),]
        
        EBINS = np.array([(1000.0,1096.00)])
                          
        emin = np.array([x[0] for x in EBINS]) 
        emax = np.array([x[1] for x in EBINS]) 
        emean = np.array([np.sqrt(1000.0*1096.00)])
        
        for i, (e_min, e_max) in enumerate(EBINS):
            out_name_srctempl = os.path.join(out_srctempl_folder,
                                             'src%s_%i-%i_test.fits'\
                                             %(cat, e_min, e_max))
            if os.path.exists(out_name_srctempl):
                src_map = hp.read_map(out_name_srctempl)
            else:
                src_map = build_src_template(cat_file, psf56_file, 
                                             emin=e_min, 
                                             emax=e_max, b_cut=abslatitude)
            hp.write_map(out_name_srctempl, src_map)
            mask_f = mask_list[i]
            mask = hp.read_map(mask_f)
            _unmask = np.where(mask > 1e-30)[0]
            tot_flux = np.mean(src_map[_unmask])
            areapix =  4*np.pi/len(src_map)
            logger.info('tot flux [cm-2s-1]: %e'%(tot_flux*areapix))
            logger.info('tot flux [cm-2s-1sr-1]: %e'%(tot_flux))
            src_map_masked = hp.ma(src_map)
            src_map_masked.mask = np.logical_not(mask)
            hp.mollview(src_map_masked, norm='log')
        plt.show()


if __name__ == '__main__':
    main()
