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
#------------------------------------------------------------------------------

"""Fits2Healpix and viceversa converter
"""

import os
import ast
import argparse
import numpy as np
import healpy as hp
import pyfits as pf
from matplotlib import cm
from GRATools import GRATOOLS_CONFIG
from GRATools.utils.logging_ import logger, abort, startmsg
from GRATools.utils.matplotlib_ import pyplot as plt
from GRATools.utils.matplotlib_ import overlay_tag, save_current_figure
from GRATools.utils.gSpline import xInterpolatedBivariateSplineLinear

__description__ = 'Converter from Cartesian to healpix format'


formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('--infile', type=str, required=True,
                    help='input fits file')
PARSER.add_argument('--nsideout', type=int, default=512,
                    help='')

def foreground_map_convert(**kwargs):
    """Viewer interface for healpix maps
    """
    input_file = kwargs['infile']
    nside_out = kwargs['nsideout']
    if not os.path.exists(input_file):
        abort("Map %s not found!"%input_file)
    frmaps = pf.open(input_file)
    maps_slices = frmaps[0].data
    energy = np.array([x[0] for x in frmaps['ENERGIES'].data])
    nside = 64#2048
    npix = hp.nside2npix(nside)
    iii = np.arange(npix)
    x,y,z = hp.pix2vec(nside, iii)
    lon_hp, lat_hp = hp.rotator.vec2dir(x,y,z,lonlat=True)
    hp_frmap = np.arange(npix, dtype=np.float64)
    lon_fits = np.arange(len(maps_slices[0][0]))
    nresx = 360./len(lon_fits)
    lon_fits_1 = (lon_fits[:1440]*nresx+180)
    lon_fits = np.append(lon_fits_1, lon_fits[1440:]*nresx-180)#+180
    lat_fits = np.arange(len(maps_slices[0]))
    lat_fits = lat_fits*nresx-90
    fr_e = []
    for i, en in enumerate(energy):
        logger.info('Running map convertion for energy %.2f...'%en)
        frmap = maps_slices[i]
        fmt = dict(xname='$l$', xunits='deg', yname='$b$',
                   yunits='deg', zname='Flux [cm$^{-2}$s$^{-1}$sr$^{-1}$]')
        lon, _indexx = np.unique(lon_fits, return_index=True)
        lat, _indexy = np.unique(lat_fits, return_index=True)
        frmap = frmap[:, _indexx]
        frspline = xInterpolatedBivariateSplineLinear(lon, lat, frmap.T, **fmt)
        for i, pix in enumerate(hp_frmap):
            hp_frmap[i] = frspline((lon_hp[i]+360)%360, lat_hp[i])
        out_name = os.path.basename(input_file).replace('.fits','_hp%i_%d.fits' 
                                                        %(nside_out, en))
        fr_e.append(hp_frmap[12426])
        out_path = os.path.join(GRATOOLS_CONFIG, 'fits', out_name)
        hp_frmap_out = hp.pixelfunc.ud_grade(hp_frmap, nside_out,  pess=True)
        hp.write_map(out_path, hp_frmap_out, coord='G')
        logger.info('Writed map %s'%out_path)
    frmaps.close()

if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    foreground_map_convert(**args.__dict__)


"""fr_e = np.array(fr_e)
plt.plot(energy, fr_e, 'o--')
plt.plot((524.81, 524.81),(-9e-9, 1e-6), '--', color='silver')
plt.plot((1000.00, 1000.00),(-9e-9, 1e-6), '--', color='silver')
plt.plot((1737.80, 1737.80),(-9e-9, 1e-6), '--', color='silver')
plt.plot((2754.23, 2754.23),(-9e-9, 1e-6), '--', color='silver')
plt.plot((4786.30, 4786.30),(-9e-9, 1e-6), '--', color='silver')
plt.plot((8317.64, 8317.64),(-9e-9, 1e-6), '--', color='silver')
plt.plot((14454.40, 14454.40),(-9e-9, 1e-6), '--', color='silver')
plt.plot((22908.68, 22908.68),(-9e-9, 1e-6), '--', color='silver')
plt.plot((39810.71, 39810.71),(-9e-9, 1e-6), '--', color='silver')
plt.plot((69183.09, 69183.09),(-9e-9, 1e-6), '--', color='silver')
plt.plot((120226.44, 120226.44),(-9e-9, 1e-6), '--', color='silver')
plt.plot((190546.06, 190546.06),(-9e-9, 1e-6), '--', color='silver')
plt.plot((331131.12, 331131.12),(-9e-9, 1e-6), '--', color='silver')
plt.plot((575439.94, 575439.94),(-9e-9, 1e-6), '--', color='silver')
plt.ylim(-9e-9, 1e-6)
plt.xlabel('Energy [MeV]')
plt.ylabel('Flux in a pixel')
#plt.yscale('log')
plt.xscale('log')
plt.show()
"""
