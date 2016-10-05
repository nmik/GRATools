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
from GRATools.utils.gFTools import cl_parse

Cl_FILES = [#os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t0_srcmask2_13bins_cls.txt'),
            #os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t56_Rm_srcmask2_13bins_cls.txt'),
            #os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t56_Rp_srcmask2_13bins_cls.txt'),
            #os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t1_srcmask2_13bins_cls.txt'),
            #os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t2_srcmask2_13bins_cls.txt'),
            #os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t4_srcmask2_13bins_cls.txt'),
            #os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t32_srcmask2_13bins_cls.txt')
            os.path.join(GRATOOLS_OUT, 
                         'Allyrs_UCV_t56_srcmask1p5_13bins_cls.txt'),
            #os.path.join(GRATOOLS_OUT, 
            #             'Allyrs_UCV_t56_srcmask2_raw_13bins_cls.txt'),
            #os.path.join(GRATOOLS_OUT,
            #             'Allyrs_UCV_t56_srcmask2_rawCn_13bins_cls.txt'),
            os.path.join(GRATOOLS_OUT,
                         'Allyrs_UCV_t56_srcmask2_13bins_cls.txt'),
            #os.path.join(GRATOOLS_OUT, 
            #             'Allyrs_UCV_t56_srcmask2_maskeast_13bins_cls.txt'),
            #os.path.join(GRATOOLS_OUT, 
            #             'Allyrs_UCV_t56_srcmask2_maskwest_13bins_cls.txt')
            #os.path.join(GRATOOLS_OUT, 
            #             'Allyrs_UCV_t56_srcmask2_cnfit_13bins_cls.txt'),
            #os.path.join(GRATOOLS_OUT, 
            #             'Allyrs_UCV_t56_srcmask2_clraw_13bins_cls.txt')
            ]
#OUT_LABEL = 'Cl_types_srcmask2'
#OUT_LABEL = 'Cl_t56_Rm-Rp'
#OUT_LABEL = 'Cl_t56_srcmask2_raw-Cn-Wbeam'
OUT_LABEL = 'Cl_t56_srcmask2-1p5'
#OUT_LABEL = 'Cl_t56_north-south'
#OUT_LABEL = 'Cl_t56_east-west'
rebinning = np.unique(np.int64(np.logspace(0, 3, 31)))
psf_ref_file = os.path.join(GRATOOLS_CONFIG, 'ascii/PSF_UCV_PSF1.txt')
_l_min = [100, 100, 100, 100, 100, 100, 49, 49, 49, 49, 49, 49, 49]
_l_max = [300, 600, 700, 1000, 1000, 1000, 1000, 1000, 1000, 
          1000, 1000, 1000, 1000]

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
    l_max = _l_max[i]#min(500, 1.9*(np.pi/np.radians(psf_en)))
    l_min = _l_min[i]#min(60, max(50-i*5,10))
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
            clmeanerr = np.sqrt(np.sum(clerrs_tocompare[j][i][_index]**2))/\
                np.sqrt(len(cls_tocompare[j][i][_index]))
            _cls_rebin.append(clmean)
            _clerrs_rebin.append(clmeanerr)
        _l_rebin = np.array(_l_rebin)
        l_range_fit = np.where(np.logical_and(_l_rebin>=l_min, _l_rebin<l_max))
        _cls_rebin = np.array(_cls_rebin)
        _clerrs_rebin = np.array(_clerrs_rebin)
        cp =  np.polyfit(_l_rebin[l_range_fit], _cls_rebin[l_range_fit], 0)
        lab = os.path.basename(f).replace('_13bins_cls.txt', '')
        plt.errorbar(_l_rebin, _cls_rebin, fmt='o', markersize=3, \
                         elinewidth=1, xerr=[xerrL, xerrR], \
                         yerr=_clerrs_rebin, label=lab )
        plt.plot([1, 1000], [cp[0], cp[0]], '-', linewidth=1, label='Fit '+lab)
    plt.plot([l_max, l_max], [-5e-15, 5e-15], '--', color='silver')
    plt.plot([l_min, l_min], [-5e-15, 5e-15], '--', color='silver')
    plt.xlim(1, _l[-1])
    plt.ylim(ymin, ymax)
    plt.xscale('log')
    #plt.yscale('log', nonposy='clip')
    plt.xlabel('$l$')
    plt.ylabel('$C_{sig,l}$')
    plt.title('%.2f - %.2f MeV'%(emin[i], emax[i]))
    plt.legend(loc=4, fontsize=10)
    save_current_figure(OUT_LABEL+'_%i-%i.png'%(emin[i], emax[i]))
    ymin = ymin + abs(ymin/(1.5+0.01*i))
    ymax = ymax - abs(ymax/(1.5+0.01*i))

    

