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
import numpy as np
from GRATools import GRATOOLS_OUT, GRATOOLS_CONFIG
from GRATools.utils.logging_ import logger, startmsg
from GRATools.utils.matplotlib_ import pyplot as plt
from GRATools.utils.matplotlib_ import overlay_tag, save_current_figure

Cl_FILES = [os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t0_13bins_cls.txt'),
            #os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t0_Rm_13bins_cls.txt'),
            #os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t0_Rp_13bins_cls.txt')]
            os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t1_13bins_cls.txt'),
            os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t2_13bins_cls.txt'),
            os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t4_13bins_cls.txt'),
            os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t32_13bins_cls.txt')]
OUT_LABEL = 'Cl_t0-t1-t2-t4-t32'
#OUT_LABEL = 'Cl_t0_Rm-Rp'

rebinning = np.unique(np.int64(np.logspace(0, 3, 31)))
#print np.unique(rebinning)
psf_ref_file = os.path.join(GRATOOLS_CONFIG, 'ascii/PSF_UCV_PSF1.txt')

def cl_parse(cl_file):
    cls = []
    clerrs = []
    emin, emax, emean = [], [], []
    f = open(cl_file, 'r')
    for line in f:
        if 'ENERGY\t' in line:
            e1, e2, em = [float(item) for item in line.split()[1:]]
            emin.append(e1)
            emax.append(e2)
            emean.append(em)
        if 'Cl\t' in line:
            cl = np.array([float(item) for item in line.split()[1:]])
            cls.append(cl)
        if 'Cl_ERR' in line:
            cl_err = np.array([float(item) for item in line.split()[1:]])
            clerrs.append(cl_err)
    f.close()
    return np.array(emin), np.array(emax), np.array(emean), cls, clerrs

emin, emax, emean = [], [], []
cls_tocompare, clerrs_tocompare = [], []
for f in Cl_FILES:
    emin, emax, emean, cls, clerrs = cl_parse(f)
    cls_tocompare.append(cls)
    clerrs_tocompare.append(clerrs)


from GRATools.utils.gWindowFunc import get_psf_ref
psf_ref = get_psf_ref(psf_ref_file)
ymin, ymax = -1e-15, 1e-15
for i in range(0, len(cls_tocompare[0])):
    psf_en = psf_ref(emean[i])
    l_max = 2*(np.pi/np.radians(psf_en))
    plt.figure(figsize=(10, 7), dpi=80)
    for j, f in enumerate(Cl_FILES):
        _l = np.arange(1, len(cls_tocompare[j][i]))
        _l_rebin, _cls_rebin,  _clerrs_rebin = [], [], []
        xerrL, xerrR = [], []
        for bmin, bmax in zip(rebinning[:-1], rebinning[1:]):
            _l_rebin.append(np.sqrt(bmin*bmax))
            xerrL.append(abs(np.sqrt(bmin*bmax)-bmin))
            xerrR.append(abs(np.sqrt(bmin*bmax)-bmax))
            _index = np.where(np.logical_and(_l>=bmin, _l<bmax)) 
            clmean = np.average(cls_tocompare[j][i][_index])
            clmeanerr = np.sqrt(np.sum(clerrs_tocompare[j][i][_index]**2))
            _cls_rebin.append(clmean)
            _clerrs_rebin.append(clmeanerr)
        _l_rebin = np.array(_l_rebin)
        _cls_rebin = np.array(_cls_rebin)
        _clerrs_rebin = np.array(_clerrs_rebin)
        lab = os.path.basename(f).replace('_13bins_cls.txt', '')
        plt.errorbar(_l_rebin, _cls_rebin, fmt='o', markersize=3, \
                         elinewidth=1, xerr=[xerrL, xerrR], \
                         yerr=_clerrs_rebin, label=lab )
        plt.plot([l_max, l_max], [-5e-15, 5e-15], '--', color='silver')
    plt.xlim(50, _l[-1])
    plt.ylim(ymin, ymax)
    plt.xscale('log')
    plt.xlabel('$l$')
    plt.ylabel('$C_{sig,l}$')
    plt.title('%.2f - %.2f MeV'%(emin[i], emax[i]))
    plt.legend()
    save_current_figure(OUT_LABEL+'_%i-%i.png'%(emin[i], emax[i]))
    ymin = ymin + abs(ymin/(1.5+0.2*i))
    ymax = ymax - abs(ymax/(1.5+0.2*i))
    #plt.show()

