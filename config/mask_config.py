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

#OUT_LABEL = 'Mask_src2_gp30_maskwest'
OUT_LABEL = 'Mask_advanced_524-1000'
#OUT_LABEL = 'Mask_advanced_1000-1737'
#OUT_LABEL = 'Mask_advanced_1737-2754'
#OUT_LABEL = 'Mask_advanced_2754-4786'
#OUT_LABEL = 'Mask_advanced_4786-8317'
#OUT_LABEL = 'Mask_advanced_8317-14454'
#OUT_LABEL = 'Mask_advanced_14454-22908'
#OUT_LABEL = 'Mask_advanced_22908-39810'
#OUT_LABEL = 'Mask_advanced_39810-69183'
#OUT_LABEL = 'Mask_advanced_69183-120226'
#OUT_LABEL = 'Mask_advanced_120226-190546'
#OUT_LABEL = 'Mask_advanced_190546-331131'
#OUT_LABEL = 'Mask_advanced_331131-575439'

NSIDE = 512
SRC_CATALOG = os.path.join(FT_DATA_FOLDER,'catalogs/gll_psc_v16.fit')
SRC_MASK_RAD = 2 #[deg]
GP_MASK_LAT = 30.
ENERGY = 743.73

"""my energy bins:

E_MIN            E_MAX           E_MEAN         
524.81           1000.00         743.73          
1000.00          1737.80         1340.69  
1737.80          2754.23         2208.67  
2754.23          4786.30         3692.56  
4786.30          8317.64         6416.93  
8317.64          14454.40        11151.34 
14454.40         22908.68        18370.95 
22908.68         39810.71        30713.33 
39810.71         69183.09        53373.66 
69183.09         120226.44       92752.78 
120226.44        190546.06       152802.88
190546.06        331131.12       255462.38
331131.12        575439.94       443942.75
"""
