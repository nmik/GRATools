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

"""Healpix viewer
"""

import os
import ast
import argparse
import numpy as np
import healpy as hp
from matplotlib import cm
from GRATools.utils.logging_ import logger, abort, startmsg
from GRATools.utils.matplotlib_ import pyplot as plt
from GRATools.utils.matplotlib_ import overlay_tag, save_current_figure


__description__ = 'Viewer interface (.fits viewer)'


formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('--field', type=int, default=0,
                    help='the input configuration file')
PARSER.add_argument('--infile', type=str, required=True,
                    help='input fits file')
PARSER.add_argument('--udgrade', type=int, default=512,
                    help='')
PARSER.add_argument('--optimized', type=ast.literal_eval, choices=[True, False],
                    default=False,
                    help='set min and max of the map in mollview')
PARSER.add_argument('--counts', type=ast.literal_eval, choices=[True, False],
                    default=False,
                    help='set min and max of the map in mollview')


def maps_view(**kwargs):
    """Viewer interface for healpix maps
    """
    input_file = kwargs['infile']
    healpix_maps = hp.read_map(input_file, field=kwargs['field'])
    if not os.path.exists(input_file):
        abort("Map %s not found!"%input_file)
    t = os.path.basename(input_file)
    plt.figure(figsize=(10, 7), dpi=80)
    nside_out = kwargs['udgrade']
    logger.info('Returning a map with NSIDE=%i'%nside_out)
    if kwargs['field'] == 0:
        if kwargs['counts'] == True:
            healpix_maps = hp.pixelfunc.ud_grade(healpix_maps, nside_out, 
                                                 pess=True, power=-2)
        else:
            healpix_maps = hp.pixelfunc.ud_grade(healpix_maps, nside_out, 
                                                 pess=True)
        if kwargs['optimized'] == True:
            logger.info('Optimizing...')
            hp.mollview(healpix_maps, title=t.replace('.fits',''), \
                            coord='G', min=1e-7, max=1e-4, norm='log')
            hp.graticule()
            overlay_tag(color='silver', x=0.45)
            save_current_figure(t.replace('.fits','.png'))
        else:
            hp.mollview(healpix_maps, title=t.replace('.fits',''), \
                            coord='G')
            hp.graticule()
            overlay_tag(color='silver', x=0.45)
            save_current_figure(t.replace('.fits','.png'))
    else:
        for i, maps in enumerate(healpix_maps):
            healpix_maps = hp.pixelfunc.ud_grade(healpix_maps, nside_out, \
                                                     pess=True)
            if kwargs['optimized'] == True:
                logger.info('Optimizing...')
                hp.mollview(maps, title=t.replace('.fits','_%i'%i), \
                                coord='G', min=1e-7, max=1e-4, norm='log')
                hp.graticule()
                overlay_tag(color='silver', x=0.45)
                save_current_figure(t.replace('.fits','.png'))
            else:
                hp.mollview(healpix_maps, title=t.replace('.fits',''), \
                                coord='G')
                hp.graticule()
                overlay_tag(color='silver', x=0.05)
                save_current_figure(t.replace('.fits','_%i.png'%i))

if __name__ == '__main__':
    args = PARSER.parse_args()
    startmsg()
    maps_view(**args.__dict__)
