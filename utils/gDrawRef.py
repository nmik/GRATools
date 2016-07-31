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


"""Drawing utilities 
"""

import os
import numpy as np

from GRATools import GRATOOLS_CONFIG
from GRATools.utils.logging_ import logger, abort
from GRATools.utils.matplotlib_ import pyplot as plt

def ref_igrb_noFGsub():
    f = open(os.path.join(GRATOOLS_CONFIG,'ascii/IGRB_P7REP_noFGsub.txt'), 'r')
    x = [float(l.split()[0]) for l in f]
    f = open(os.path.join(GRATOOLS_CONFIG,'ascii/IGRB_P7REP_noFGsub.txt'), 'r')
    y = [float(l.split()[1]) for l in f]
    igrb, = plt.plot(x, y, '--', color='red')
    label = 'LAT data (P7REP IGRB)'
    return igrb, label

def ref_igrb_band():
    fsub = open(os.path.join(GRATOOLS_CONFIG,'ascii/IGRB_sub_band.txt'), 'r')
    fup = open(os.path.join(GRATOOLS_CONFIG,'ascii/IGRB_up_band.txt'), 'r')
    sub_line_x = [float(l.split()[0]) for l in fsub]
    up_line_x = [float(l.split()[0]) for l in fup]
    fsub = open(os.path.join(GRATOOLS_CONFIG,'ascii/IGRB_sub_band.txt'), 'r')
    fup = open(os.path.join(GRATOOLS_CONFIG,'ascii/IGRB_up_band.txt'), 'r')
    sub_line_y = [float(l.split()[1]) for l in fsub]
    up_line_y = [float(l.split()[1]) for l in fup]
    igrb = plt.fill(up_line_x + sub_line_x[::-1], up_line_y + sub_line_y[::-1], \
                 color='silver', alpha='0.5')
    label='Galactic foreground modeling uncertainty'
    legend = plt.Rectangle((0, 0), 1, 1, fc='silver')
    #plt.xscale('log')
    #plt.yscale('log')
    return legend, label

def main():
    """Test module
    """
    plt.figure()
    leg, lab = ref_igrb_band()
    igrb, lab2 = ref_igrb_noFGsub()
    plt.xscale('log')                                                                                                                                  
    plt.yscale('log')
    plt.legend([leg, igrb], [lab, lab2])
    plt.show()

if __name__ == '__main__':
    main()
    
