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

_ENERGY = np.logspace(1, 6,1000)

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
    FLUX30_100 = SOURCES.field('Flux30_100')[lat_index]
    FLUX100_300 = SOURCES.field('Flux100_300')[lat_index]
    FLUX300_1000 = SOURCES.field('Flux300_1000')[lat_index]
    E0 = SOURCES.field('Pivot_Energy')[lat_index]
    K = SOURCES.field('Flux_Density')[lat_index]
    alpha = GAMMA = SOURCES.field('Spectral_Index')[lat_index]
    spec_type = SOURCES.field('SpectrumType')[lat_index]
    beta = SOURCES.field('beta')[lat_index]
    b = SOURCES.field('Exp_Index')[lat_index]
    Ecut = SOURCES.field('Cutoff')[lat_index]
    src_name = SOURCES.field('Source_Name')[lat_index]
    src_cat.close()
    ##########################################################
    INT_FLUXES = []
    for i, s in enumerate(spec_type):
        spec = 0
        if s == 'PowerLaw':
            spec = PowerLaw(K[i], alpha[i], E0[i])
            integral_flux_bin = spec.integral(emin, emax)
            INT_FLUXES.append(integral_flux_bin)
        elif s == 'LogParabola':
            spec = LogParabola(K[i], alpha[i], beta[i], E0[i])
            integral_flux_bin = spec.integral(emin, emax)
            INT_FLUXES.append(integral_flux_bin)
        else:
            spec = PLSuperExpCutoff(K[i], alpha[i], E0[i], b[i], Ecut[i])
            integral_flux_bin = spec.integral(emin, emax)
            INT_FLUXES.append(integral_flux_bin)
    logger.info('Building...')
    src_map = src_builder(psf_th, nside, GLON, GLAT, INT_FLUXES)
    logger.info('Done!')
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
    for i, bn in enumerate(src_pix):
        pixdir1 = hp.pixelfunc.pix2ang(nside, bn)
        radintpix = hp.query_disc(nside, (x[i],y[i],z[i]), np.radians(10))
        pixdir2 = hp.pixelfunc.pix2ang(nside, radintpix)
        dist = hp.rotator.angdist(pixdir1, pixdir2)
        psf = psf_spline(dist)
        isomap[radintpix]= isomap[radintpix] + intflux_[i]*psf
    return isomap

def main():
    """Simple test unit
    """
    sourcetemp = True
    
    if sourcetemp == True:
        abslatitude = 20
        cat_file = os.path.join(GRATOOLS_CONFIG,'catalogs/gll_psc_v16.fit')
        psf_file = os.path.join(GRATOOLS_OUT, 
                                'psf_P8R2_ULTRACLEANVETO_V6_32.fits')
        isofile = os.path.join(GRATOOLS_CONFIG, 'models', 
                           'iso_P8R2_ULTRACLEANVETO_V6_v06.txt')
        e_min, e_max= 120226, 331131 #158, 350#1000, 1737
        src_templ_map = build_src_template(cat_file, psf_file, emin=e_min, 
                                           emax=e_max, b_cut=abslatitude)
        from GRATools.utils.gForeground import get_iso_integral_flux_map
        iso_map = get_iso_integral_flux_map(isofile, e_min, e_max)
        iso_src_map = iso_map + src_templ_map
        mask_f =  os.path.join(GRATOOLS_CONFIG, 
                              'fits/Mask_final_8317.fits')
        mask = hp.read_map(mask_f)
        iso_src_map_masked = hp.ma(iso_src_map)
        iso_src_map_masked.mask = np.logical_not(mask)
        titolo = 'Energy: %.1f - %.1f MeV'%(e_min, e_max)
        hp.mollview(iso_src_map_masked.filled(), title=titolo, 
                    coord='G')#,norm='log'), min=1.5e-7, max=4e-7)
        hp.graticule()
        plt.show()
        
        


if __name__ == '__main__':
    main()
