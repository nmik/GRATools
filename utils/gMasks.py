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
from GRATools import FT_DATA_FOLDER
from GRATools import GRATOOLS_CONFIG
from GRATools.utils.logging_ import logger, abort
from GRATools.utils.matplotlib_ import pyplot as plt


def mask_src(cat_file, MASK_S_RAD, NSIDE):
    """Returns the 'bad pixels' defined by the position of a source and a 
       certain radius away from that point.

       SOURCE_CAT: str
           opened fits file with the sorce catalog
       MASK_S_RAD: float
           radius around each source definig bad pixels to mask
       NSIDE: int
           healpix nside parameter
    """
    logger.info('Mask for sources activated')
    src_cat = pf.open(cat_file)
    NPIX = hp.pixelfunc.nside2npix(NSIDE)
    CAT = src_cat['LAT_Point_Source_Catalog']
    BAD_PIX_SRC = []
    SOURCES = CAT.data
    RADrad = MASK_S_RAD*np.pi/180.
    for i in range (0,len(SOURCES)-1):
        GLON = SOURCES[i][3]
        GLAT = SOURCES[i][4]
        x, y, z = hp.rotator.dir2vec(GLON,GLAT,lonlat=True)
        b_pix= hp.pixelfunc.vec2pix(NSIDE, x, y, z)
        BAD_PIX_SRC.append(b_pix) 
    BAD_PIX_inrad = []
    for bn in BAD_PIX_SRC:
        pixVec = hp.pix2vec(NSIDE,bn)
        radintpix = hp.query_disc(NSIDE, pixVec, RADrad)
        BAD_PIX_inrad.extend(radintpix)  
    BAD_PIX_SRC.extend(BAD_PIX_inrad)
    src_cat.close()
    return BAD_PIX_SRC

def mask_gp(MASK_GP_LAT,NSIDE):
    """Returns the 'bad pixels' around the galactic plain .

       MASK_GP_LAT: float
           absolute value of galactic latitude definig bad pixels to mask
       NSIDE: int
           healpix nside parameter
    """
    logger.info('Mask for the galactic plane activated')
    NPIX = hp.pixelfunc.nside2npix(NSIDE)
    BAD_PIX_GP = []
    iii = range(NPIX)
    x,y,z = hp.pix2vec(NSIDE,iii)
    lon,lat = hp.rotator.vec2dir(x,y,z,lonlat=True)
    for i,b in enumerate(lat):
        if abs(b) <= MASK_GP_LAT:
            BAD_PIX_GP.append(iii[i])
    return BAD_PIX_GP

def mask_high_lat(MASK_LAT,NSIDE):
    """Returns the 'bad pixels' at high latitudes from the galactic plain .

       MASK_LAT: float
           absolute value of galactic latitude definig bad pixels to mask
       NSIDE: int
           healpix nside parameter
    """
    logger.info('Mask for high latitudes activated')
    NPIX = hp.pixelfunc.nside2npix(NSIDE)
    BAD_PIX_GP = []
    iii = range(NPIX)
    x,y,z = hp.pix2vec(NSIDE,iii)
    lon,lat = hp.rotator.vec2dir(x,y,z,lonlat=True)
    for i,b in enumerate(lat):
        if abs(b) >= MASK_LAT:
            BAD_PIX_GP.append(iii[i])
    return BAD_PIX_GP

def mask_hemi_north(NSIDE):
    """Returns the 'bad pixels' in the northen hemisphere.

       NSIDE: int
           healpix nside parameter
    """
    logger.info('Masking northen hemisphere...')
    NPIX = hp.pixelfunc.nside2npix(NSIDE)
    BAD_PIX_HEMI_N = []
    iii = range(NPIX)
    x,y,z = hp.pix2vec(NSIDE,iii)
    lon,lat = hp.rotator.vec2dir(x,y,z,lonlat=True)
    for i,b in enumerate(lat):
        if b >= 0:
            BAD_PIX_HEMI_N.append(iii[i])
    return BAD_PIX_HEMI_N


def mask_hemi_south(NSIDE):
    """Returns the 'bad pixels' in the southern hemisphere.

       NSIDE: int
           healpix nside parameter
    """
    logger.info('Masking southern hemisphere...')
    NPIX = hp.pixelfunc.nside2npix(NSIDE)
    BAD_PIX_HEMI_S = []
    iii = range(NPIX)
    x,y,z = hp.pix2vec(NSIDE,iii)
    lon,lat = hp.rotator.vec2dir(x,y,z,lonlat=True)
    for i,b in enumerate(lat):
        if b <= 0:
            BAD_PIX_HEMI_S.append(iii[i])
    return BAD_PIX_HEMI_S

def mask_hemi_east(NSIDE):
    """Returns the 'bad pixels' in the easthern hemisphere.

       NSIDE: int
           healpix nside parameter
    """
    logger.info('Masking easthern hemisphere...')
    NPIX = hp.pixelfunc.nside2npix(NSIDE)
    BAD_PIX_HEMI_E = []
    iii = range(NPIX)
    x,y,z = hp.pix2vec(NSIDE,iii)
    lon,lat = hp.rotator.vec2dir(x,y,z,lonlat=True)
    for i,b in enumerate(lon):
        if -180 <= b <= 0:
            BAD_PIX_HEMI_E.append(iii[i])
    return BAD_PIX_HEMI_E

def mask_hemi_west(NSIDE):
    """Returns the 'bad pixels' in the westhern hemisphere.

       NSIDE: int
           healpix nside parameter
    """
    logger.info('Masking westhern hemisphere...')
    NPIX = hp.pixelfunc.nside2npix(NSIDE)
    BAD_PIX_HEMI_W = []
    iii = range(NPIX)
    x,y,z = hp.pix2vec(NSIDE,iii)
    lon,lat = hp.rotator.vec2dir(x,y,z,lonlat=True)
    for i,b in enumerate(lon):
        if 0 <= b <= 180:
            BAD_PIX_HEMI_W.append(iii[i])
    return BAD_PIX_HEMI_W

def mask_src_weighted(cat_file, ENERGY, NSIDE):
    """Returns the 'bad pixels' defined by the position of a source and a 
       certain radius away from that point. The radii increase with the 
       brightness.

       SOURCE_CAT: str
           opened fits file with the sorce catalog
       NSIDE: int
           healpix nside parameter
    """
    from GRATools.utils.gWindowFunc import get_psf_ref
    psf_ref_file = os.path.join(GRATOOLS_CONFIG, 'ascii/PSF_UCV_PSF1.txt')
    src_cat = pf.open(cat_file)
    NPIX = hp.pixelfunc.nside2npix(NSIDE)
    CAT = src_cat['LAT_Point_Source_Catalog']
    CAT_EXTENDED = src_cat['ExtendedSources']
    BAD_PIX_SRC = []
    SOURCES = CAT.data
    EXT_SOURCES = CAT_EXTENDED.data
    src_cat.close()
    psf_ref = get_psf_ref(psf_ref_file)
    psf_en = psf_ref(ENERGY)
    psf_min, psf_max =  psf_ref.y[5], psf_ref.y[-1] 
    norm_min, norm_max = 1, 0.3
    norm = norm_min + psf_en*((norm_max - norm_min)/(psf_max - psf_min)) -\
        psf_min*((norm_max - norm_min)/(psf_max - psf_min))
    logger.info('Normalization of radii due to energy: %.3f'%norm)
    print 'Psf(%.2f)= %.2f'%(ENERGY, psf_en)
    FLUX = SOURCES.field('Flux1000')
    flux_min, flux_max = min(FLUX), max(FLUX)
    rad_min, rad_max = 2., 5.
    RADdeg = rad_min + FLUX*((rad_max - rad_min)/(flux_max - flux_min)) -\
        flux_min*((rad_max - rad_min)/(flux_max - flux_min))
    RADrad = np.radians(RADdeg)
    #*****
    #plt.title('Radius($\phi$)')
    #plt.plot(FLUX, RADdeg, 'ro', ms=3, alpha=0.75)
    #plt.plot((1e-12, 1e-5), (2, 2), '-', color='silver', linewidth=1.0)
    #plt.xlabel('$\phi$')
    #plt.ylabel('Radius [$\circ$]')
    #plt.yscale('log')
    #plt.xscale('log')
    #plt.ylim(0.7, 6)
    #plt.show()
    #*****
    logger.info('Masking the extended Sources')
    logger.info('-> 10deg around CenA and LMC')
    logger.info('-> 5deg around the remaining')
    for i, src in enumerate(EXT_SOURCES):
        NAME = EXT_SOURCES[i][0]
        GLON = EXT_SOURCES[i][4]
        GLAT = EXT_SOURCES[i][5]
        if NAME == 'LMC' or NAME == 'CenA Lobes':
            x, y, z = hp.rotator.dir2vec(GLON,GLAT,lonlat=True)
            b_pix= hp.pixelfunc.vec2pix(NSIDE, x, y, z)
            BAD_PIX_SRC.append(b_pix) 
            radintpix = hp.query_disc(NSIDE, (x, y, z), np.radians(10)*norm)
            BAD_PIX_SRC.extend(radintpix)
        else:
            x, y, z = hp.rotator.dir2vec(GLON,GLAT,lonlat=True)
            b_pix = hp.pixelfunc.vec2pix(NSIDE, x, y, z)
            BAD_PIX_SRC.append(b_pix) 
            radintpix = hp.query_disc(NSIDE, (x, y, z), np.radians(5)*norm)
            BAD_PIX_SRC.extend(radintpix)
    logger.info('Flux-weighted mask for sources activated')
    for i, src in enumerate(SOURCES):
        GLON = SOURCES[i][3]
        GLAT = SOURCES[i][4]
        x, y, z = hp.rotator.dir2vec(GLON,GLAT,lonlat=True)
        b_pix= hp.pixelfunc.vec2pix(NSIDE, x, y, z)
        BAD_PIX_SRC.append(b_pix) 
        radintpix = hp.query_disc(NSIDE, (x, y, z), RADrad[i]*norm)
        BAD_PIX_SRC.extend(radintpix)
    return BAD_PIX_SRC


def main():
    """Simple test unit
    """
    from GRATools.utils.gWindowFunc import get_psf_ref
    psf_ref_file = os.path.join(GRATOOLS_CONFIG, 'ascii/PSF_UCV_PSF1.txt')
    psf_ref = get_psf_ref(psf_ref_file)
    energy = np.array([743.73,1340.69,2208.67,3692.56,6416.93,11151.34,
                       18370.95,30713.33,53373.66,92752.78,152802.88,
                       255462.38,443942.75])
    psf_en = psf_ref(energy)
    psf_min, psf_max =  psf_ref.y[5], psf_ref.y[-1] 
    norm_min, norm_max = 1, 0.3
    norm = norm_min + psf_en*((norm_max - norm_min)/(psf_max - psf_min)) -\
        psf_min*((norm_max - norm_min)/(psf_max - psf_min))
    plt.title('Normalization Factor')
    plt.plot(energy, norm, 'ro--', ms=5, alpha=0.75)
    #plt.plot((1e-12, 1e-5), (2, 2), '-', color='silver', linewidth=1.0)
    plt.xlabel('Energy [MeV]')
    plt.ylabel('Normalization factor')
    plt.yscale('log')
    plt.xscale('log')
    plt.ylim(0.09, 2)
    plt.show()
    
    nside = 512
    SRC_CATALOG_FILE = os.path.join(FT_DATA_FOLDER,'catalogs/gll_psc_v16.fit')
    bad_pix = mask_src_weighted(SRC_CATALOG_FILE, 10000, nside)
    bad_pix += mask_gp(30, nside)
    npix = hp.nside2npix(nside)
    mask = np.ones(npix)
    for bpix in np.unique(bad_pix):
        mask[bpix] = 0
    fsky = 1-(len(np.unique(bad_pix))/float(npix))
    title = 'f$_{sky}$ = %.3f'%fsky
    hp.mollview(mask, title=title, coord='G')
    plt.show()


if __name__ == '__main__':
    main()
