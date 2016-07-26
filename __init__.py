#!/usr/bin/env python                                                          #
#                                                                              #
# Autor: Michela Negro, University of Torino.                                  #
#                                                                              #
# This program is free software; you can redistribute it and/or modify         #
# it under the terms of the GNU GengReral Public License as published by       #
# the Free Software Foundation; either version 3 of the License, or            #
# (at your option) any later version.                                          #
#                                                                              #
#------------------------------------------------------------------------------#

"""advlab: Framework for 'Advanced Laboratory' PhD Course
"""

import os

PACKAGE_NAME = 'GRATools'

"""Basic folder structure of the package.
"""
GRATOOLS_ROOT = os.path.abspath(os.path.dirname(__file__))
GRATOOLS_BIN = os.path.join(GRATOOLS_ROOT, 'bin')
GRATOOLS_CONFIG = os.path.join(GRATOOLS_ROOT, 'config')
GRATOOLS_UTILS = os.path.join(GRATOOLS_ROOT, 'utils')
GRATOOLS_DOC = os.path.join(GRATOOLS_ROOT, 'doc')
GRATOOLS_DATA = os.path.join(GRATOOLS_ROOT, 'data')

""" This is where we put the actual (FT1 and FT2) data sets.  
"""

from GRATools.utils.logging_ import logger
try:
    FT_DATA_FOLDER = os.environ['FT_DATA']
    logger.info('Base data folder set to $FT_DATA = %s...' % FT_DATA_FOLDER)
except KeyError:
    FT_DATA_FOLDER = '/data1/data/FT-files'
    logger.info('$FT_DATA not set, base data folder set to %s...' %\
                FT_DATA_FOLDER)

""" This is the output directory.
"""
try:
    GRATOOLS_OUT = os.environ['GRATOOLS_OUT']
    GRATOOLS_OUT_FIG = os.environ['GRATOOLS_OUT_FIG']
except:
    GRATOOLS_OUT = os.path.join(GRATOOLS_ROOT, 'output')
    GRATOOLS_OUT_FIG = os.path.join(GRATOOLS_ROOT, 'output/figures')

if __name__ == '__main__':
    print('GRATOOLS_ROOT: %s' % GRATOOLS_ROOT)
