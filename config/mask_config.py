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

OUT_LABEL = 'Mask_src2_gp30' 

NSIDE = 512
SRC_CATALOG = os.path.join(FT_DATA_FOLDER,'catalogs/gll_psc_v16.fit')
SRC_MASK_RAD = 2 #[deg]
GP_MASK_LAT = 30.
ENERGY =  600. #If --srcweighted False it is not used

