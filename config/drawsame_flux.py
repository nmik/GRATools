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
from GRATools.utils.gFTools import get_crbkg
from GRATools.utils.matplotlib_ import pyplot as plt
from GRATools.utils.matplotlib_ import overlay_tag, save_current_figure

CRBKG_FILE = os.path.join(GRATOOLS_CONFIG,'ascii/CRbkg.txt')
FLUX_FILES = [os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t0_13bins_parameters.txt'),
              #os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t0_Rm_13bins_parameters.txt'),
              #os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t0_Rp_13bins_parameters.txt')]
              os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t4_13bins_parameters.txt'),
              os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t32_13bins_parameters.txt'),
              os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t1_13bins_parameters.txt'),
              os.path.join(GRATOOLS_OUT, 'Allyrs_UCV_t2_13bins_parameters.txt')]
#FLUX_LABELS = ['UCV No evtype', 'UCV No evtype Rock -', 'UCV No evtype Rock +']
FLUX_LABELS = ['UCV No evtype', 'UCV PSF0', 'UCV PSF3', 'UCV FRONT', 'UCV BACK'] 
OUT_LABEL = 'Flux_evtypes'

plt.figure(figsize=(10, 7), dpi=80)
from GRATools.utils.gDrawRef import ref_igrb_band
from GRATools.utils.gDrawRef import ref_igrb_noFGsub
leg, lab = ref_igrb_band()
igrb, lab_igrb = ref_igrb_noFGsub()
flux = []
for f in FLUX_FILES:
    from GRATools.utils.gFTools import get_cl_param
    _emin, _emax, _emean, _f, _ferr, _cn, fsky = get_cl_param(f)
    
    spec = plt.errorbar(_emean, _f*_emean, fmt='o', markersize=3, \
                     elinewidth=1, xerr=(_emax-_emin)/2, yerr=_ferr*_emean)
    label = os.path.basename(f).replace('_parameters.txt', '')
    flux.append(spec)
plt.xscale("log")
plt.yscale("log")
plt.xlabel('Energy [MeV]')
plt.ylabel('E$^{2}$ $\cdot$ Flux [MeV cm$^{-2}$ s$^{-1}$ sr$^{-1}$]')
plt.title(' Energy Spectrum')
plt.legend([igrb, leg]+flux, [lab_igrb, lab]+FLUX_LABELS, loc=3)
#overlay_tag()
save_current_figure(OUT_LABEL+'_ESpec.png')
