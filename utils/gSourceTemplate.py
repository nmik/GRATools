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


"""Functions to define map masks
"""

import os
import healpy as hp
import numpy as np
import pyfits as pf
from GRATools import GRATOOLS_CONFIG
from GRATools.utils.logging_ import logger, abort
from GRATools.utils.matplotlib_ import pyplot as plt
from GRATools.utils.gWindowFunc import get_psf_ref

def psf_en(psf_file, energy):
    psf_spline = get_psf_ref(psf_file)
    psf_energy = psf_spline(energy)
    return psf_energy
    

def buils_src_template(cat_file, psf_file, energy, NSIDE):
    """Returns a map of sources inside the catalog smoothed by the PSF 
       at the gven energy.

       cat_file: str
           .fits fits file with the sorce catalog
       psf_file: str
           ascii file with 2 columns (Energy - PSF)
       energy: float
           center of the energy (micro) bin [MeV]
       NSIDE: int
           healpix nside parameter
    """
    psf = psf_en(psf_file, energy)
    logger.info('ENERGY = %.2f MeV -> PSF = %.4f deg (%.f rad)'
                %(energy, psf, np.radians(psf)))
    src_cat = pf.open(cat_file)
    npix = hp.pixelfunc.nside2npix(NSIDE)
    src_map = np.full(npix, 0)
    CAT = src_cat['LAT_Point_Source_Catalog']
    SOURCES = CAT.data
    FLUX = SOURCES.field('Flux1000')
    GLON = SOURCES.field('GLON')
    GLAT = SOURCES.field('GLAT')
    src_cat.close()
    x, y, z = hp.rotator.dir2vec(GLON,GLAT,lonlat=True)
    src_pix= hp.pixelfunc.vec2pix(NSIDE, x, y, z)
    for i, bn in enumerate(src_pix):
        src_map[bn] = FLUX[i]
    src_map = hp.sphtfunc.smoothing(src_map, fwhm=np.radians(psf/1.18))
    return src_map


def main():
    """Simple test unit
    """
    sourcetemp = True
    
    if sourcetemp == True:
        cat_file = os.path.join(GRATOOLS_CONFIG,'catalogs/gll_psc_v16.fit')
        psf_file = os.path.join(GRATOOLS_CONFIG, 
                                'ascii/PSF_UCV_PSF1.txt')
        mask_f =  os.path.join(GRATOOLS_CONFIG, 
                              'fits/Mask_gp30.fits')
        mask = hp.read_map(mask_f)
        energy = 1340.69
        NSIDE = 512
        src_templ_map = buils_src_template(cat_file, psf_file, energy, NSIDE)
        src_templ_map_masked = hp.ma(src_templ_map)
        src_templ_map_masked.mask = np.logical_not(mask)
        titolo = 'Energy: %f MeV'%energy
        hp.mollview(src_templ_map_masked.filled(), title=titolo, 
                    coord='G', min=1e-15, max=6e-10, norm='log')
        hp.graticule()
        plt.show()
        
   


if __name__ == '__main__':
    main()
