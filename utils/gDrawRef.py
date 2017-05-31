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
    """Draw the reference plot of the Pass7 total flux
    """
    f = open(os.path.join(GRATOOLS_CONFIG,'ascii/IGRB_P7REP_noFGsub.txt'), 'r')
    x = [float(l.split()[0]) for l in f]
    f = open(os.path.join(GRATOOLS_CONFIG,'ascii/IGRB_P7REP_noFGsub.txt'), 'r')
    y = [float(l.split()[1]) for l in f]
    igrb, = plt.plot(x, y, '--', color='red')
    label = 'LAT data (P7REP)'
    return igrb, label

def ref_foreground_spec():
    """Draw the reference plot of the total flux of the foreground model
    """
    f = open(os.path.join(GRATOOLS_CONFIG,'ascii/fore_ref_spec.txt'), 'r')
    _x = [float(l.split()[0]) for l in f]
    f = open(os.path.join(GRATOOLS_CONFIG,'ascii/fore_ref_spec.txt'), 'r')
    _y = [float(l.split()[1]) for l in f]
    fore = plt.plot(_x, _y, color='silver')
    label='Published model (2016)'
    legend = plt.Rectangle((0, 0), 1, 1, fc='silver')
    return legend, label

def ref_igrb_band():
    """Draw the systematics band of the Pass7 IGRB flux
    """
    fsub = open(os.path.join(GRATOOLS_CONFIG,'ascii/IGRB_sub_band.txt'), 'r')
    fup = open(os.path.join(GRATOOLS_CONFIG,'ascii/IGRB_up_band.txt'), 'r')
    sub_line_x = [float(l.split()[0]) for l in fsub]
    up_line_x = [float(l.split()[0]) for l in fup]
    fsub = open(os.path.join(GRATOOLS_CONFIG,'ascii/IGRB_sub_band.txt'), 'r')
    fup = open(os.path.join(GRATOOLS_CONFIG,'ascii/IGRB_up_band.txt'), 'r')
    sub_line_y = [float(l.split()[1]) for l in fsub]
    up_line_y = [float(l.split()[1]) for l in fup]
    igrb = plt.fill(up_line_x + sub_line_x[::-1], up_line_y + sub_line_y[::-1],
                    color='0.8',  zorder=0)#, alpha='0.5')
    label='Pass7 IGRB measure (Ackermann et al. 2014)'
    legend = plt.Rectangle((0, 0), 1, 1, fc='0.8')
    return legend, label

def ref_cp_band():
    """Draw the statistics band of the previous Cp measure
    """
    fsub = open(os.path.join(GRATOOLS_CONFIG,'ascii/cp_sub_band.txt'), 'r')
    fup = open(os.path.join(GRATOOLS_CONFIG,'ascii/cp_up_band.txt'), 'r')
    sub_line_x = [float(l.split()[0]) for l in fsub]
    up_line_x = [float(l.split()[0]) for l in fup]
    fsub = open(os.path.join(GRATOOLS_CONFIG,'ascii/cp_sub_band.txt'), 'r')
    fup = open(os.path.join(GRATOOLS_CONFIG,'ascii/cp_up_band.txt'), 'r')
    sub_line_y = [float(l.split()[1]) for l in fsub]
    up_line_y = [float(l.split()[1]) for l in fup]
    cp = plt.fill(up_line_x + sub_line_x[::-1], up_line_y + sub_line_y[::-1], \
                 color='0.8')#, alpha='0.5')
    label= 'M. Fornasa et al. 2016'
    legend = plt.Rectangle((0, 0), 1, 1, fc='0.8')
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
    
