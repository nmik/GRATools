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

import os
import numpy as np
from GRATools import FT_DATA_FOLDER
from GRATools.utils.gFTools import ebinning_fits_file

OUT_LABEL = '7yr'

START_WEEK = 297 
END_WEEK = 344
EBINNING_ARRAY = np.logspace(2., 6., 101)
EBINNING_FILE = ebinning_fits_file(EBINNING_ARRAY)
FT2_FILE =  os.path.join(FT_DATA_FOLDER, \
                             'spacecraft/lat_spacecraft_merged.fits')
FILTER_CUT='DATA_QUAL==1&&LAT_CONFIG==1&&LAT_MODE==5&&IN_SAA!=T'+\
               '&&((ABS(ROCK_ANGLE)<52))'

GTSELECT_DICT = {'infile': 'DEFAULT',
                 'emin': 100,
                 'emax': 1000000,
                 'zmax': 90,
                 'evclass': 1024,
                 'evtype': 56,
                 'outfile': 'DEFAULT',
                 'clobber': 'yes'}

GTMKTIME_DICT = {'evfile': 'DEFAULT',
                 'scfile': FT2_FILE,
                 'filter': FILTER_CUT,
                 'roicut': 'no',
                 'outfile': 'DEFAULT',
                 'clobber': 'yes'}

GTBIN_DICT = {'evfile': 'DEFAULT',
              'algorithm': 'HEALPIX',
              'scfile': FT2_FILE,
              'hpx_ordering_scheme': 'RING',
              'hpx_order': 9,
              'coordsys': 'GAL',                                  
              'hpx_ebin': 'yes',
              'ebinalg': 'FILE',
              'ebinfile': EBINNING_FILE,
              'outfile': 'DEFAULT',
              'clobber': 'yes'}

GTLTCUBE_DICT = {'evfile':'DEFAULT',
                 'scfile': FT2_FILE,
                 'zmax': 90,                     
                 'dcostheta': 0.025,
                 'binsz': 1,
                 'outfile': 'DEFAULT',
                 'clobber': 'yes'}

GTEXPCUBE2_DICT = {'infile': 'DEFAULT',
                   'cmap': 'DEFAULT',
                   'irfs': 'CALDB',
                   'outfile': 'DEFAULT',
                   'clobber': 'yes'}




