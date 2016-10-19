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

FLUX_REF = os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t56_srcmask2_13bins_parameters.txt')
#FLUX_REF_LABEL = 'UCV No evtype'
FLUX_REF_LABEL = 'UCV (t56) srcmask2'
#OUT_LABEL = 'Flux_Rm-Rp' 
OUT_LABEL = 'Flux_types_srcmask2'
FLUX_FILES = [#os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t56_srcmask2_13bins_parameters.txt'),
              #os.path.join(GRATOOLS_OUT, 'Allyrs_SRC_t32_srcmask2_13bins_parameters.txt'),
              #os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t56_srcmask1p5_13bins_parameters.txt'),
              #os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t56_maskweighted_13bins_parameters.txt'),
              #os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t0_13bins_parameters.txt'),
              #os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t56_Rm_13bins_parameters.txt'),
              #os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t56_Rp_13bins_parameters.txt')
              os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t4_srcmask2_13bins_parameters.txt'),
              os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t32_srcmask2_13bins_parameters.txt'),
              os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t1_srcmask2_13bins_parameters.txt'),
              os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t2_srcmask2_13bins_parameters.txt')
              ]
FLUX_LABELS = ['UCV (t56) PSF0', 'UCV (t56) PSF3', 'UCV (t56) FRONT', 'UCV (t56) BACK']
#OUT_LABEL = 'Flux_evtypes'

plt.figure(figsize=(10, 7), dpi=80)
from GRATools.utils.gFTools import get_cl_param
_emin, _emax, _emean, _f_ref, _ferr_ref, _cn_ref, fsky_ref = get_cl_param(FLUX_REF)
plt.plot((_emin[0], _emax[-1]), (0, 0), '--', color='gray')
plt.plot((100, 1000000), (0.1, 0.1), '-', color='silver', linewidth=1.0)
plt.plot((100, 1000000), (-0.1, -0.1), '-', color='silver', linewidth=1.0)
flux, flux_label = [], []
plt.title('Flux Residuals - Ref.: %s'%FLUX_REF_LABEL)
for f in FLUX_FILES:
    _emin, _emax, _emean, _f, _ferr, _cn, fsky = get_cl_param(f)
    #print _ferr
    _res = (_f_ref - _f)/_f_ref
    _res_err = np.sqrt(((_f/_f_ref**2)*_ferr_ref)**2+(_ferr/_f_ref)**2)
    spec = plt.errorbar(_emean, _res, fmt='o', markersize=3, \
                     elinewidth=1, xerr=(_emax-_emin)/2, yerr=_res_err)
    flux.append(spec)
plt.xscale("log")
#plt.yscale("log")
plt.ylim(-1, 1)
plt.xlabel('Energy [MeV]')
plt.ylabel('($\Phi_{fef}$ - $\Phi$) / $\Phi_{fef}$')
plt.legend(flux, FLUX_LABELS, loc=3)
#overlay_tag()
save_current_figure(OUT_LABEL+'_ESpecRes.png')
