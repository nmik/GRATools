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


OUT_LABEL = '1yr_UCV_t56'
MICRO_NBINS = 100
MACRO_BINS = [(18,25),(25,31),(31,36),(36,42),(42,48),(48,54),(54,59),(59,65),(65,71),(71,77),(77,82),(82,88),(88,94)]

NSIDE = 512
SRC_CATALOG = os.path.join(FT_DATA_FOLDER,'catalogs/gll_psc_v16.fit')
SRC_MASK_RAD = 1.5 #[deg]
GP_MASK_LAT = 30.

