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
from GRATools import GRATOOLS_CONFIG
from GRATools import GRATOOLS_OUT, FT_DATA_FOLDER

"""BEAM WINDOW FUNCTION VARIABLES
"""
LTCUBE = os.path.join(FT_DATA_FOLDER, \
                    'output/output_gtltcube/Allyrs_filtered_gti_ltcube.fits')
IRFS = 'P8R2_ULTRACLEANVETO_V6'#'P8R2_SOURCE_V6'
EVTYPE = 56
OUT_W_LABEL = '%s_Rp_%i'%(IRFS, EVTYPE)
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
MASK_LIST = [os.path.join(GRATOOLS_CONFIG, 'fits/Mask_advanced_524-1000.fits'),
             os.path.join(GRATOOLS_CONFIG, 'fits/Mask_advanced_1000-1737.fits'),
             os.path.join(GRATOOLS_CONFIG, 'fits/Mask_advanced_1737-2754.fits'),
             os.path.join(GRATOOLS_CONFIG, 'fits/Mask_advanced_2754-4786.fits'),
             os.path.join(GRATOOLS_CONFIG, 'fits/Mask_advanced_4786-8317.fits'),
             os.path.join(GRATOOLS_CONFIG, 'fits/Mask_advanced_8317-14454.fits'),
             os.path.join(GRATOOLS_CONFIG, 'fits/Mask_advanced_14454-22908.fits'),
             os.path.join(GRATOOLS_CONFIG, 'fits/Mask_advanced_22908-39810.fits'),
             os.path.join(GRATOOLS_CONFIG, 'fits/Mask_advanced_39810-69183.fits'),
             os.path.join(GRATOOLS_CONFIG, 'fits/Mask_advanced_69183-120226.fits'),
             os.path.join(GRATOOLS_CONFIG, 'fits/Mask_advanced_120226-190546.fits'),
             os.path.join(GRATOOLS_CONFIG, 'fits/Mask_advanced_190546-331131.fits'),
             os.path.join(GRATOOLS_CONFIG, 'fits/Mask_advanced_331131-575439.fits')]
IN_LABEL = 'Allyrs_UCV_t56_Rp_maskweighted' #the label of the files to be considered
MASK_FILE = os.path.join(GRATOOLS_CONFIG, 'fits/Mask_no.fits')
#MASK_FILE = MASK_LIST
BINNING_LABEL = '13bins'
OUT_LABEL = IN_LABEL
