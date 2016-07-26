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


import healpy as hp
import numpy as np
from GRATools.utils.logging_ import logger, abort


def mask_src(SOURCE_CAT,MASK_S_RAD,NSIDE):
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
    NPIX = hp.pixelfunc.nside2npix(NSIDE)
    CAT = SOURCE_CAT['LAT_Point_Source_Catalog']
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
    return BAD_PIX_SRC

def mask_gp(MASK_GP_LAT,NSIDE):
    """Returns the 'bad pixels' around the galactic plain .

       MASK_GP_LAT: float
           absolute value of galactic latitude definig bad pixels to mask
       NSIDE: int
           healpix nside parameter
    """
    logger.info('Mask for the galactic plain activated')
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
