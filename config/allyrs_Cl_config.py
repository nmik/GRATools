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
from GRATools import GRATOOLS_OUT, FT_DATA_FOLDER

"""BEAM WINDOW FUNCTION VARIABLES
"""
LTCUBE = os.path.join(FT_DATA_FOLDER, \
                    'output/output_gtltcube/Allyrs_filtered_gti_ltcube.fits')
IRFS ='P8R2_ULTRACLEANVETO_V6'
EVTYPE = 56
OUT_W_LABEL = '%s_%i'%(IRFS, EVTYPE)
WEIGHT_SPEC_INDEX = 2.3
PSF_FILE = os.path.join(GRATOOLS_OUT, 'psf_%s.fits'%OUT_W_LABEL)
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
              'ntheta': 300}

"""Cl VARIABLES
"""
IN_LABEL = 'Allyrs_UCV_t56' #the label of the files to be considered
MASK_FILE = os.path.join(GRATOOLS_OUT, 'Mask_src1p5_gp30.fits')
BINNING_LABEL = '13bins'
OUT_LABEL = IN_LABEL+'_srcmask1p5'
