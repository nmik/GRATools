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
from GRATools import GRATOOLS_OUT
from GRATools.utils.logging_ import logger, startmsg
from GRATools.utils.matplotlib_ import pyplot as plt
from GRATools.utils.matplotlib_ import overlay_tag, save_current_figure
from GRATools.utils.gWindowFunc import get_psf_ref
from GRATools.utils.gFTools import cp_parse

#FLUX_REF = os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t0_srcmask2_13bins_cps.txt')
FLUX_REF = os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t56_srcmask2_masknorth_13bins_cps.txt')
#FLUX_REF_LABEL = 'UCV-Notype srcmask2  '
FLUX_REF_LABEL = 'UCV-t56 srcmask2 masknorth'
FLUX_FILES = [#os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t1_srcmask2_13bins_cps.txt'),
              #os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t2_srcmask2_13bins_cps.txt'),
              #os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t4_srcmask2_13bins_cps.txt'),
              #os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t32_srcmask2_13bins_cps.txt'),
              #os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t56_srcmask2_13bins_cps.txt'),
              #os.path.join(GRATOOLS_OUT, 
              #             'Allyrs_UCV_t56_srcmask1p5_13bins_cps.txt'),
              #os.path.join(GRATOOLS_OUT, 
              #             'Allyrs_UCV_t56_srcmask2_masknorth_13bins_cps.txt'),
              os.path.join(GRATOOLS_OUT, 
                           'Allyrs_UCV_t56_srcmask2_masksouth_13bins_cps.txt')
              ]
#FLUX_LABELS = ['UCV (t56) srcmask1p5']                     
FLUX_LABELS = ['UCV (t56) srcmask2 masksouth']#, 'UCV (t56) srcmask2 masknorth']
#FLUX_LABELS = ['UCV-FRONT srcmask2',
#               'UCV-BACK srcmask2',
#               'UCV-PSF0 srcmask2', 
#               'UCV-PSF3 srcmask2',
#               'UCV-PSF1+2+3 srcmask2']
#OUT_LABEL = 'Cp_t56_masksrc2vs1p5' 
OUT_LABEL = 'Cp_t56_northVSsouth'
#OUT_LABEL = 'Cp_types_srmask2'

plt.figure(figsize=(10, 7), dpi=80)
emin, emax, emean, cps_ref, cperrs_ref = cp_parse(FLUX_REF)
plt.plot((emin[0], emax[-1]), (0, 0), '--', color='gray')
plt.plot((0.1, 1000), (0.1, 0.1), '-', color='silver', linewidth=1.0)
plt.plot((0.1, 1000), (-0.1, -0.1), '-', color='silver', linewidth=1.0)
spec, spec_label = [], []
plt.title('Flux Residuals - Ref.: %s'%FLUX_REF_LABEL)
for f in FLUX_FILES:
    emin, emax, emean, cps, cperrs = cp_parse(f)
    _res = (cps_ref - cps)/cps_ref
    _res_err = np.sqrt(((cps/cps_ref**2)*cperrs_ref)**2+(cperrs/cps_ref)**2)
    cp_plot = plt.errorbar(emean, _res, fmt='o', markersize=3,
                           elinewidth=1, 
                           xerr=[(emean-emin), (emax-emean)], 
                           yerr=_res_err)
    spec.append(cp_plot)
plt.xscale("log")
plt.ylim(-1, 1)
plt.xlabel('Energy [GeV]')
plt.ylabel('(C$_{P,ref}$ - C$_{P}$) / C$_{P,ref}$')
plt.legend(spec, FLUX_LABELS, loc=3)
save_current_figure(OUT_LABEL+'_residuals.png')
