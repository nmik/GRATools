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
from GRATools import FT_DATA_FOLDER


IN_LABELS_LIST = ['1yr_UCV_t2', '2yr_UCV_t2',
                  '3yr_UCV_t2', '4yr_UCV_t2',
                  '5yr_UCV_t2', '6yr_UCV_t2',
                  '7yr_UCV_t2', '8yr_UCV_t2'] 
OUT_LABEL = 'Allyrs_UCV_t2'
BINNING_LABEL = '13bins'
MICRO_NBINS = 100
MACRO_BINS = [(18,24),(25,30),(31,35),(36,41),(42,47),(48,53),(54,58),(59,64),(65,70),(71,76),(77,81),(82,87),(88,93)]
POWER_LOW_INDEX = 2.30
NSIDE = 512
SRC_CATALOG = os.path.join(FT_DATA_FOLDER,'catalogs/gll_psc_v16.fit')
SRC_MASK_RAD = 1.5 #[deg]
GP_MASK_LAT = 30.

