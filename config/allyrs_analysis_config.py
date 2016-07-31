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


OUT_LABEL = 'Allyrs_UCV_t56'
BINNING_LABEL = '13bins'
LABELS_LIST = ['1yr_UCV_t56', '2yr_UCV_t56',
               '3yr_UCV_t56', '4yr_UCV_t56',
               '5yr_UCV_t56', '6yr_UCV_t56',
               '7yr_UCV_t56', '8yr_UCV_t56']
FILES_DICT = {'524-1000': [0]+LABELS_LIST,
              '1000-1737': [1]+LABELS_LIST,
              '1737-2754': [2]+LABELS_LIST,
              '2754-4786': [3]+LABELS_LIST,
              '4786-8317': [4]+LABELS_LIST,
              '8317-14454': [5]+LABELS_LIST,
              '14454-22908': [6]+LABELS_LIST,
              '22908-39810': [7]+LABELS_LIST,
              '39810-69183': [8]+LABELS_LIST,
              '69183-120226': [9]+LABELS_LIST,
              '120226-190546': [10]+LABELS_LIST,
              '190546-331131': [11]+LABELS_LIST,
              '331131-575439': [12]+LABELS_LIST
              }
