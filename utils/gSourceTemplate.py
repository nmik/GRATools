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
    integral_psf_th = psf_th.integral(psf_th.x[0], psf_th.x[-1])
    print 'integral', integral_psf_th
    psf_th = psf_th.scale(1/integral_psf_th)
    #psf_th.plot()
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
    K = SOURCES.field('Flux_Density')[lat_index] #cm-2 s-1 sr-1 MeV-1
    alpha = GAMMA = SOURCES.field('Spectral_Index')[lat_index]
    spec_type = SOURCES.field('SpectrumType')[lat_index]
    beta = SOURCES.field('beta')[lat_index]
    b = SOURCES.field('Exp_Index')[lat_index]
    Ecut = SOURCES.field('Cutoff')[lat_index]
    src_name = SOURCES.field('Source_Name')[lat_index]
    src_cat.close()
    ##########################################################
    INT_FLUXES = [] #cm-2 s-1 sr-1
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
    out_srctempl_folder = os.path.join(GRATOOLS_OUT, 'output_src')
    if not os.path.exists(out_srctempl_folder):
        os.makedirs(out_srctempl_folder)
    cat = os.path.basename(cat_file).replace('.fit','')
    out_name_srctempl = os.path.join(out_srctempl_folder,'src%s_%i-%i.fits'\
                                             %(cat, emin, emax))
    hp.write_map(out_name_srctempl, src_map)
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
    srcfluxtempl = isomap
    return srcfluxtempl

def main():
    """Simple test unit
    """
    sourcetemp = True
    
    if sourcetemp == True:
        abslatitude = 0
        out_srctempl_folder = os.path.join(GRATOOLS_OUT, 'output_src')
        cat_file = os.path.join(GRATOOLS_CONFIG,'catalogs/gll_psc_v16.fit')
        cat = os.path.basename(cat_file).replace('.fit','')
        psf32_file = os.path.join(GRATOOLS_OUT, 
                                'psf_P8R2_ULTRACLEANVETO_V6_32.fits')
        psf56_file = os.path.join(GRATOOLS_OUT, 
                                'psf_P8R2_ULTRACLEANVETO_V6_56.fits')
        mask_list = [os.path.join(GRATOOLS_CONFIG, 'fits/Mask_bat_420.fits'),
             os.path.join(GRATOOLS_CONFIG, 'fits/Mask_bat_420.fits'),
             os.path.join(GRATOOLS_CONFIG, 'fits/Mask_bat_524.fits'),
             os.path.join(GRATOOLS_CONFIG, 'fits/Mask_bat_1000.fits'),
             os.path.join(GRATOOLS_CONFIG, 'fits/Mask_bat_1737.fits'),
             os.path.join(GRATOOLS_CONFIG, 'fits/Mask_bat_2754.fits'),
             os.path.join(GRATOOLS_CONFIG, 'fits/Mask_bat_4786.fits'),
             os.path.join(GRATOOLS_CONFIG, 'fits/Mask_bat_8317.fits'),
             os.path.join(GRATOOLS_CONFIG, 'fits/Mask_bat3FHL_8317.fits'),
             os.path.join(GRATOOLS_CONFIG, 'fits/Mask_bat3FHL_8317.fits'),
             os.path.join(GRATOOLS_CONFIG, 'fits/Mask_bat3FHL_8317.fits'),
             os.path.join(GRATOOLS_CONFIG, 'fits/Mask_bat3FHL_8317.fits'),
             os.path.join(GRATOOLS_CONFIG, 'fits/Mask_bat3FHL_8317.fits'),
             os.path.join(GRATOOLS_CONFIG, 'fits/Mask_bat3FHL_8317.fits')]
        EBINS = np.array([(158.49,301.00), (301.00,524.81)
                          (524.82,1000.00), (1000.00,1737.80), 
                          (1737.80,2754.23), (2754.23,4786.30),
                          (4786.30,8317.64), (8317.64,14454.40), 
                          (14454.40,22908.68), (22908.68,39810.71), 
                          (39810.71,69183.09),(69183.09,120226.44), 
                          (120226.44,331131.12), (331131.12,1000000.00)])
        emin = np.array([x[0] for x in EBINS]) 
        emax = np.array([x[1] for x in EBINS]) 
        emean = np.array([218.78,398.11,724.44,1318.26, 2187.76,3630.78,
                          6309.57,10964.78,18197.01,30199.52,52480.75,91201.08,
                          199526.23,575439.94])
        cn = np.array([1.464618e-16,4.965573e-17,2.448300e-17,2.867573e-18,
                       1.129696e-18,6.847362e-19,3.189106e-19,1.397455e-19,
                       5.380968e-20,3.033270e-20,1.242414e-20,4.835532e-21,
                       2.443251e-21,4.048537e-22])
        
        cps = []
        for i, (e_min, e_max) in enumerate(EBINS):
            out_name_srctempl = os.path.join(out_srctempl_folder,
                                             'src%s_%i-%i.fits'\
                                             %(cat, e_min, e_max))
            if os.path.exists(out_name_srctempl):
                 src_map = hp.read_map(out_name_srctempl)
            else:
                if i < 3:
                    src_map = build_src_template(cat_file, psf32_file, 
                                                 emin=e_min, 
                                                 emax=e_max, b_cut=abslatitude)
                else:
                    src_map = build_src_template(cat_file, psf56_file, 
                                                 emin=e_min, 
                                                 emax=e_max, b_cut=abslatitude)
            mask_f = mask_list[i]
            mask = hp.read_map(mask_f)
            _unmask = np.where(mask > 1e-30)[0]
            src_map_masked = hp.ma(src_map)
            src_map_masked.mask = np.logical_not(mask)
            cp = np.sum((src_map[_unmask])**2)
            cps.append(cp)
            print cp-cn[i]
        #titolo = 'Energy: %.1f - %.1f MeV, E$^2$Flux$_{mean}$=%e'%(e_min, 
        #                                                           e_max, 
        #                                                           tot_flux)
        #hp.mollview(src_map_masked.filled(), title=titolo, 
        #            coord='G', norm='log')#,norm='log'), min=1.5e-7, max=4e-7)
        #hp.graticule()
        plt.figure()
        plt.plot(emean, (emean**4/(emax-emin)**2)*np.array(cps-cn), '.')
        plt.show()
            
        
        


if __name__ == '__main__':
    main()
