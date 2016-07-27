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
EBINNING_FILE = os.path.join(GRATOOLS_OUT, '1yr_UCV_t56new_ebinning.txt')
LABELS_LIST = ['1yr_UCV_t56', '2yr_UCV_t56',
               '3yr_UCV_t56', '4yr_UCV_t56',
               '5yr_UCV_t56', '6yr_UCV_t56',
               '7yr_UCV_t56', '8yr_UCV_t56']
FILES_DICT = {'18-25': [0]+LABELS_LIST,
              '25-31': [1]+LABELS_LIST,
              '31-36': [2]+LABELS_LIST,
              '36-42': [3]+LABELS_LIST,
              '42-48': [4]+LABELS_LIST,
              '48-54': [5]+LABELS_LIST,
              '54-59': [6]+LABELS_LIST,
              '59-65': [7]+LABELS_LIST,
              '65-71': [8]+LABELS_LIST,
              '71-77': [9]+LABELS_LIST,
              '77-82': [10]+LABELS_LIST,
              '82-88': [11]+LABELS_LIST,
              '88-94': [12]+LABELS_LIST
              }
