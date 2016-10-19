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
from GRATools import GRATOOLS_CONFIG
from GRATools import GRATOOLS_OUT, FT_DATA_FOLDER

"""GTPSF VARIABLES
"""
LTCUBE = os.path.join(FT_DATA_FOLDER, \
                    'output/output_gtltcube/Allyrs_filtered_gti_ltcube.fits')
IRFS = 'P8R2_ULTRACLEANVETO_V6'
EVTYPE = 56
PSF_FILE = os.path.join(GRATOOLS_OUT, 'psf_csi-fit.fits')
DICT_GTPSF = {'expcube': LTCUBE,
              'outfile': PSF_FILE,
              'irfs': IRFS,
              'evtype': EVTYPE,
              'ra': 45,
              'dec': 45,
              'emin': 500,
              'emax': 600000,
              'nenergies': 10,
              'thetamax': 30,
              'ntheta': 600}

"""Csi VARIABLES
"""
IN_LABEL = 'Allyrs_UCV_t56_maskweighted' 
TH_BINNING = np.array([0, 0.001, 0.003, 0.008, 0.024, 0.070])
BINNING_LABEL = '13bins'
OUT_LABEL = IN_LABEL
