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
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

from GRATools import GRATOOLS_CONFIG
from GRATools.utils.logging_ import logger, abort
from GRATools.utils.matplotlib_ import pyplot as plt
from GRATools.utils.gSpline import xInterpolatedUnivariateSplineLinear

def make_error_boxes(ax, lx, hx, ydata, yerror, facecolor='k',
                     edgecolor='None', alpha=0.3):
    errorboxes = []
    for x1, x2, y, ye in zip(lx, hx, ydata, yerror):
        rect = Rectangle((x1, y-ye), abs(x2-x1), 2*ye)
        errorboxes.append(rect)
    pc = PatchCollection(errorboxes, facecolor=facecolor, alpha=alpha,
                         edgecolor=edgecolor)
    ax.add_collection(pc)
    art =  ax.errorbar(np.sqrt(lx*hx), ydata, xerr=x1*0, yerr=yerror*0,
                          fmt='None', ecolor='k', alpha=0.3)
    return art

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
    fmt = dict(xname='', xunits='', yname='', yunits='')
    up_line_y2 = np.interp(sub_line_x, up_line_x, up_line_y)
    lab='Pass7 IGRB (Ackermann et al. 2015)'
    igrb = plt.fill_between(sub_line_x, sub_line_y, up_line_y2, label=lab,
                    color='0.8',  zorder=0)
    label='Pass7 IGRB (Ackermann et al. 2015)'
    return igrb, lab

def ref_cp_band():
    """Draw the statistics band of the previous Cp measure
    """
    refcp = open(os.path.join(GRATOOLS_CONFIG,'ascii/PASS7_Cp.txt'), 'r')
    rows = (row.strip().split() for row in refcp)
    (lx, hx, y, yerr) = zip(*rows)
    lx = np.asarray(lx[1:], dtype='float')
    hx = np.asarray(hx[1:], dtype='float')
    y = np.asarray(y[1:], dtype='float')
    yerr = np.asarray(yerr[1:], dtype='float')
    fig, ax = plt.subplots(1)
    refcp_plot = make_error_boxes(ax, lx, hx, y*np.sqrt(lx*hx)**4/(hx-lx)**2, 
                                  yerr*np.sqrt(lx*hx)**4/(hx-lx)**2)
    label= 'M. Fornasa et al. 2016'
    legend = plt.Rectangle((0, 0), 1, 1, fc='k', alpha=0.3)
    return legend, label


def main():
    """Test module
    """
    leg, lab = ref_cp_band()
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
    
    plt.figure()
    leg, lab = ref_igrb_band()
    igrb, lab2 = ref_igrb_noFGsub()
    plt.xscale('log')       
    plt.yscale('log')
    plt.legend([leg, igrb], [lab, lab2])
    plt.show()

if __name__ == '__main__':
    main()
    
