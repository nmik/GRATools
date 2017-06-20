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

IN_LABELS_LIST = ['1yr_UCV_t32', '2yr_UCV_t32',
                  '3yr_UCV_t32', '4yr_UCV_t32',
                  '5yr_UCV_t32', '6yr_UCV_t32',
                  '7yr_UCV_t32', '8yr_UCV_t32'] 

FORE_FILES_LIST = [os.path.join(GRATOOLS_CONFIG,'fits/gll_iem_v06_hp512_523.fits'),
                   os.path.join(GRATOOLS_CONFIG,'fits/gll_iem_v06_hp512_715.fits'),
                   os.path.join(GRATOOLS_CONFIG,'fits/gll_iem_v06_hp512_978.fits'),
                   os.path.join(GRATOOLS_CONFIG,'fits/gll_iem_v06_hp512_1338.fits'),
                   os.path.join(GRATOOLS_CONFIG,'fits/gll_iem_v06_hp512_1830.fits')]
OUT_LABEL = 'Allyrs_UCV_t32'

BINNING_LABEL = 'mybins'
MICRO_NBINS = 100
MACRO_BINS = [(18,24),(25,30)]

POWER_LOW_INDEX = 2.30

MASK_FILE = os.path.join(GRATOOLS_CONFIG, 'fits/Mask_src2_gp30.fits') 
MASK_LABEL = 'mymask'
# In MASK_FILE Can be also a list of Mask, if so they will be used in order 
# for each energy bin.


