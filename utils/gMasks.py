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
from GRATools import GRATOOLS_CONFIG
from GRATools.utils.logging_ import logger, abort
from GRATools.utils.matplotlib_ import pyplot as plt
from GRATools.utils.gWindowFunc import get_psf_ref


def mask_src(cat_file, MASK_S_RAD, NSIDE):
    """Returns the 'bad pixels' defined by the position of a source and a 
       certain radius away from that point.

       cat_file: str
           .fits file of the sorce catalog
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
        GLON = SOURCES.field('GLON')[i]
        GLAT = SOURCES.field('GLAT')[i]
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

def mask_extsrc(cat_file, MASK_S_RAD, NSIDE):
    """Returns the 'bad pixels' defined by the position of a source and a 
       certain radius away from that point.

       cat_file: str
           .fits file of the sorce catalog
       MASK_S_RAD: float
           radius around each source definig bad pixels to mask
       NSIDE: int
           healpix nside parameter
    """
    logger.info('Mask for extended sources activated')
    src_cat = pf.open(cat_file)
    NPIX = hp.pixelfunc.nside2npix(NSIDE)
    CAT_EXTENDED = src_cat['ExtendedSources']
    BAD_PIX_SRC = []
    EXT_SOURCES = CAT_EXTENDED.data
    src_cat.close()
    for i, src in enumerate(EXT_SOURCES):
        NAME = EXT_SOURCES.field('Source_Name')[i]
        GLON = EXT_SOURCES.field('GLON')[i]
        GLAT = EXT_SOURCES.field('GLAT')[i]
        if 'LMC' in NAME or 'CenA Lobes' in NAME:
            x, y, z = hp.rotator.dir2vec(GLON,GLAT,lonlat=True)
            b_pix= hp.pixelfunc.vec2pix(NSIDE, x, y, z)
            BAD_PIX_SRC.append(b_pix) 
            radintpix = hp.query_disc(NSIDE, (x, y, z), np.radians(10))
            BAD_PIX_SRC.extend(radintpix)
        else:
            x, y, z = hp.rotator.dir2vec(GLON,GLAT,lonlat=True)
            b_pix = hp.pixelfunc.vec2pix(NSIDE, x, y, z)
            BAD_PIX_SRC.append(b_pix) 
            radintpix = hp.query_disc(NSIDE, (x, y, z), np.radians(5))
            BAD_PIX_SRC.extend(radintpix)
    return BAD_PIX_SRC

def mask_gp(MASK_GP_LAT, NSIDE):
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

def mask_bat_gp(FLUX_CUT, NSIDE):
    """Returns the 'bad pixels' around the galactic plain. 
                      
       FLUX_CUT: float              
           absolute value of galactic latitude definig bad pixels to mask     
       NSIDE: int                                   
           healpix nside parameter                                
    """
    logger.info('Batman mask for the galactic plane activated')
    flux_map = hp.read_map(os.path.join(GRATOOLS_CONFIG,
                                 'fits/gll_iem_v06_hp512_1000-1096.fits'))
    NPIX = hp.pixelfunc.nside2npix(NSIDE)
    if len(flux_map) != NPIX:
        flux_map = hp.ud_grade(flux_map, nsideout=NSIDE)
    else:
        pass
    BAD_PIX_GP = []
    iii = range(NPIX)
    x,y,z = hp.pix2vec(NSIDE,iii)
    lon,lat = hp.rotator.vec2dir(x,y,z,lonlat=True)
    for i,b in enumerate(flux_map):
        if b >= FLUX_CUT:
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

def mask_src_weighted_custom(cat_file, ENERGY, NSIDE):
    """Returns the 'bad pixels' defined by the position of a source and a 
       certain radius away from that point. The radii increase with the 
       brightness and rescaled by a factor between 1 and 0.3 shaped as the PSF.

       cat_file: str
          .fits file with the sorce catalog
       ENERGY: float
          Mean energy of the map to be masked
       NSIDE: int
          healpix nside parameter
    """
    psf_ref_file = os.path.join(GRATOOLS_CONFIG, 'ascii/PSF_UCV_PSF1.txt')
    src_cat = pf.open(cat_file)
    NPIX = hp.pixelfunc.nside2npix(NSIDE)
    CAT = src_cat[1]
    BAD_PIX_SRC = []
    SOURCES = CAT.data
    src_cat.close()
    psf_ref = get_psf_ref(psf_ref_file)
    psf_en = psf_ref(ENERGY)
    psf_min, psf_max =  psf_ref.y[5], psf_ref.y[-1] 
    norm_min, norm_max = 1, 0.3
    norm = norm_min + psf_en*((norm_max - norm_min)/(psf_max - psf_min)) -\
        psf_min*((norm_max - norm_min)/(psf_max - psf_min))
    logger.info('Normalization of radii due to energy: %.3f'%norm)
    logger.info('Psf(%.2f)= %.2f'%(ENERGY, psf_en))
    FLUX = np.log10(SOURCES.field('eflux1000'))
    flux_min, flux_max = min(FLUX), max(FLUX)
    rad_min, rad_max = 1, 5.
    RADdeg = rad_min + FLUX*((rad_max - rad_min)/(flux_max - flux_min)) -\
        flux_min*((rad_max - rad_min)/(flux_max - flux_min))
    RADrad = np.radians(RADdeg)
    logger.info('Flux-weighted mask for sources activated')
    TS = SOURCES.field('ts') 
    indTS25 = TS > 25.
    GLON = SOURCES.field('GLON')[indTS25]
    GLAT = SOURCES.field('GLAT')[indTS25]
    logger.info('Num Src: %i'%len(TS))
    logger.info('Num Src TS>25: %i'%len(TS[indTS25]))
    for i, src in enumerate(SOURCES[indTS25]):
        x, y, z = hp.rotator.dir2vec(GLON[i],GLAT[i],lonlat=True)
        b_pix= hp.pixelfunc.vec2pix(NSIDE, x, y, z)
        BAD_PIX_SRC.append(b_pix)
        radintpix = hp.query_disc(NSIDE, (x, y, z), RADrad[i]*norm)
        BAD_PIX_SRC.extend(radintpix)
    return BAD_PIX_SRC
    

def mask_src_weighted(cat_src_file, cat_extsrc_file, ENERGY, NSIDE):
    """Returns the 'bad pixels' defined by the position of a source and a 
       certain radius away from that point. The radii increase with the 
       brightness and rescaled by a factor between 1 and 0.3 shaped as the PSF.

       cat_src_file: str
          .fits file with the source catalog
       cat_extsrc_file: str
          .fits file with the extended sources catalog
       ENERGY: float
          Mean energy of the map to be masked
       NSIDE: int
          healpix nside parameter
    """
    psf_ref_file = os.path.join(GRATOOLS_CONFIG, 'ascii/PSF_UCV_custom.txt')
    src_cat = pf.open(cat_src_file)
    extsrc_cat = pf.open(cat_extsrc_file)
    NPIX = hp.pixelfunc.nside2npix(NSIDE)
    CAT = src_cat['LAT_Point_Source_Catalog']
    CAT_EXTENDED = extsrc_cat['ExtendedSources']
    BAD_PIX_SRC = []
    SOURCES = CAT.data
    EXT_SOURCES = CAT_EXTENDED.data
    src_cat.close()
    psf_ref = get_psf_ref(psf_ref_file)
    if '4FGL' in cat_src_file:
        FLUX = np.log10(SOURCES.field('Flux_Density'))
    if 'psch' in cat_src_file:
        FLUX = np.log10(SOURCES.field('Flux'))
    else:
        FLUX = np.log10(SOURCES.field('Flux1000'))
    flux_min, flux_max = min(FLUX), max(FLUX)
    rad_min = 2*psf_ref(ENERGY)
    rad_max = 5*psf_ref(ENERGY)
    RADdeg = rad_min + FLUX*((rad_max - rad_min)/(flux_max - flux_min)) -\
        flux_min*((rad_max - rad_min)/(flux_max - flux_min))
    RADrad = np.radians(RADdeg)
    logger.info('Masking the extended Sources')
    logger.info('-> 10xPSF(%.2f GeV) around CenA and LMC'%(ENERGY/1000))
    logger.info('-> 5xPSF(%.2f GeV) around the remaining'%(ENERGY/1000))
    for i, src in enumerate(EXT_SOURCES):
        NAME = EXT_SOURCES[i][0]
        GLON = EXT_SOURCES.field('GLON')[i]
        GLAT = EXT_SOURCES.field('GLAT')[i]
        if 'LMC' in NAME or 'CenA Lobes' in NAME:
            rad = 10*psf_ref(ENERGY)
            x, y, z = hp.rotator.dir2vec(GLON,GLAT,lonlat=True)
            b_pix= hp.pixelfunc.vec2pix(NSIDE, x, y, z)
            BAD_PIX_SRC.append(b_pix) 
            radintpix = hp.query_disc(NSIDE, (x, y, z), np.radians(10))
            BAD_PIX_SRC.extend(radintpix)
        else:
            rad = 5*psf_ref(ENERGY)
            x, y, z = hp.rotator.dir2vec(GLON,GLAT,lonlat=True)
            b_pix = hp.pixelfunc.vec2pix(NSIDE, x, y, z)
            BAD_PIX_SRC.append(b_pix) 
            radintpix = hp.query_disc(NSIDE, (x, y, z), np.radians(5))
            BAD_PIX_SRC.extend(radintpix)
    logger.info('Flux-weighted mask for sources activated')
    for i, src in enumerate(SOURCES):
        GLON = SOURCES.field('GLON')[i]
        GLAT = SOURCES.field('GLAT')[i]
        x, y, z = hp.rotator.dir2vec(GLON,GLAT,lonlat=True)
        b_pix= hp.pixelfunc.vec2pix(NSIDE, x, y, z)
        BAD_PIX_SRC.append(b_pix) 
        radintpix = hp.query_disc(NSIDE, (x, y, z), RADrad[i])
        BAD_PIX_SRC.extend(radintpix)
    return BAD_PIX_SRC


def main():
    """Simple test unit
    """
    from scipy.interpolate import interp1d
    cat_file = os.path.join(GRATOOLS_CONFIG,'catalogs/gll_psc_v16.fit')
    src_cat = pf.open(cat_file)
    CAT = src_cat['LAT_Point_Source_Catalog']
    SOURCES = CAT.data
    src_cat.close()
    flux = np.sort(SOURCES.field('Flux1000'))
    flux_min, flux_max = min(flux), max(flux)
    """
    rad_min, rad_max = 0.5, 5.
    RADdeg = rad_min + flux*((rad_max - rad_min)/(flux_max - flux_min)) -\
        flux_min*((rad_max - rad_min)/(flux_max - flux_min))
    plt.figure(facecolor='white')
    plt.title('Disk Radius as a function of $\Phi_{src}$')
    plt.plot(flux, RADdeg, '*-', color='cornflowerblue', linewidth=1, 
             ms=5, alpha=0.2)
    plt.xlabel('Flux [photon/cm$^{2}$/s]')
    plt.ylabel('mask disk Radius [$^{\circ}$]')
    plt.yscale('log')
    plt.xscale('log')
    plt.ylim(1.5, 6)
    plt.grid(which='both', linestyle='--', linewidth=1)
    """
    """
    from GRATools.utils.gWindowFunc import get_psf_ref
    psf_ref_file = os.path.join(GRATOOLS_CONFIG, 'ascii/PSF_UCV_custom.txt')
    psf_ref = get_psf_ref(psf_ref_file)
    energy = np.logspace(1, 6, 10)
    energy = [301.99,524.81,1000.00,2754.23,
              8317.64,22908.68,331131.12]
    lab = ['PSF3','PSF3','PSF1+2+3', 'PSF1+2+3', 'PSF1+2+3','PSF1+2+3', 
           'PSF1+2+3']
    plt.figure(facecolor='white')  
    for i, en in enumerate(energy):
        r_min = 2*psf_ref(en)
        r_max = 5*psf_ref(en)
        print 'R min', r_min, 'R max', r_max
        print 'F min', flux_min, 'F max', flux_max
        RADdeg = r_min + np.log10(flux)*((r_max - r_min)/(np.log10(flux_max) - np.log10(flux_min))) -\
        np.log10(flux_min)*((r_max - r_min)/(np.log10(flux_max) - np.log10(flux_min)))
        #plt.figure(facecolor='white')
        plt.title('Disk Radius as a function of $\Phi_{src}$')
        plt.plot(np.log10(flux), RADdeg, 'o-', linewidth=1,
                 ms=3, alpha=0.4, label="%.1f GeV (%s)"%(en/1000., lab[i]))
        plt.xlabel('Flux [photon/cm$^{2}$/s]')
        plt.ylabel('Mask disk Radius [$^{\circ}$]')
        #plt.yscale('log')
        #plt.xscale('log')
        #plt.ylim(1.5, 6)
        plt.grid(which='both', linestyle='--', linewidth=1)
    plt.legend(fontsize=11, loc=2)
    """
    nside = 512
    SRC_CATALOG_FILE = os.path.join(GRATOOLS_CONFIG,'catalogs/gll_psc_v16.fit')
    bad_pix = mask_src_weighted(SRC_CATALOG_FILE, SRC_CATALOG_FILE,8317, nside)
    bad_pix += mask_bat_gp(2.5e-7, nside)
    npix = hp.nside2npix(nside)
    mask = np.ones(npix)
    for bpix in np.unique(bad_pix):
        mask[bpix] = 0
    fsky = 1-(len(np.unique(bad_pix))/float(npix))
    #title = 'Energy = 100 GeV, F$_{sky}$ = %.3f'%(fsky)
    title = 'F$_{sky}$ = %.3f'%(fsky)   
    hp.mollview(mask, title=title, coord='G', cmap='bone')
    hp.graticule(color='silver')
    plt.show()


if __name__ == '__main__':
    main()
