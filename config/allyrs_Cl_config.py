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
from GRATools import GRATOOLS_OUT


IN_LABEL = 'Allyrs_UCV_t56' #the label of the files to be considered
BINNING_LABEL = '13bins'
OUT_LABEL = IN_LABEL
WBEAM_FILE = 'config/ascii/Wbeam_p8_clean_v6.txt'
PSF_REF_FILE = 'config/ascii/PSF_UCV_PSF1.txt'
