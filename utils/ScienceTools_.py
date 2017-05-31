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


"""Remake of some of the Fermi Science Tools
"""


import os
import time
import gt_apps as my_apps
from GRATools.utils.logging_ import logger
from GRATools import GRATOOLS_OUT
from GRATools import FT_DATA_FOLDER

FT_DATA_OUT = os.path.join(FT_DATA_FOLDER,'output')


def gtselect(label, filter_dict):
    """gtselect from Science Tools.

       label: str
          To automatically set the name of the output file
       filter_dict: python dict
          To define all the parameters
    """
    logger.info('Running gtselect...')
    LABEL = label
    if not os.path.exists(os.path.join(FT_DATA_OUT,'output_gtselect')):
        os.makedirs(os.path.join(FT_DATA_OUT, 'output_gtselect'))
    OUTPATH = os.path.join(FT_DATA_OUT, 'output_gtselect')
    OUTFILE = os.path.join(OUTPATH, LABEL + '_filtered.fits')
    if os.path.exists(OUTFILE):
        logger.info('ATT: Already created %s'%OUTFILE)
        return OUTFILE
    for key in filter_dict:
        my_apps.filter[key] = filter_dict[key]
    my_apps.filter['outfile'] = OUTFILE
    my_apps.filter.run()
    logger.info('Created %s'%OUTFILE)
    logger.info('gtselect --> CPU time spent: %.2f'%time.clock())
    return OUTFILE

def gtmktime(label, maketime_dict):
    """gtmktime from Science Tools.

       label: str
          To automatically set the name of the output file
       filter_dict: python dict
          To define all the parameters
    """
    logger.info('Running gtmktime...')
    LABEL = label
    if not os.path.exists(os.path.join(FT_DATA_OUT, 'output_gtmktime')):
        os.makedirs(os.path.join(FT_DATA_OUT, 'output_gtmktime'))
    OUTPATH = os.path.join(FT_DATA_OUT, 'output_gtmktime')
    OUTFILE = os.path.join(OUTPATH, LABEL + '_filtered_gti.fits')
    if os.path.exists(OUTFILE):
        logger.info('ATT: Already created %s'%OUTFILE)
        return OUTFILE
    for key in maketime_dict:
        if key == 'outfile':
            if maketime_dict[key] == 'DEFAULT':
                my_apps.maketime['outfile'] = OUTFILE
            else:
                my_apps.maketime[key] = maketime_dict[key]
            continue
        my_apps.maketime[key] = maketime_dict[key]
    my_apps.maketime.run()
    logger.info('Created %s'%OUTFILE)
    logger.info('gtmktime --> CPU time spent: %.2f'%time.clock())
    return OUTFILE

def gtbin(label, evtbin_dict):
    """gtbin from Science Tools. 

       label: str
          To automatically set the name of the output file
       filter_dict: python dict
          To define all the parameters
    """
    logger.info('Running gtbin...')
    LABEL = label
    if not os.path.exists(os.path.join(FT_DATA_OUT, 'output_gtbin')):
        os.makedirs(os.path.join(FT_DATA_OUT, 'output_gtbin'))
    OUTPATH = os.path.join(FT_DATA_OUT, 'output_gtbin')
    OUTFILE = os.path.join(OUTPATH, LABEL + '_filtered_gti_bin.fits')
    if os.path.exists(OUTFILE):
        logger.info('ATT: Already created %s'%OUTFILE)
        return OUTFILE
    for key in evtbin_dict:
        if key == 'outfile':
            if evtbin_dict[key] == 'DEFAULT':
                my_apps.evtbin['outfile'] = OUTFILE
            else:
                my_apps.evtbin[key] = evtbin_dict[key]
            continue
        my_apps.evtbin[key] = evtbin_dict[key]
    my_apps.evtbin.run()
    logger.info('Created %s'%OUTFILE)
    logger.info('gtbin --> CPU time spent: %.2f'%time.clock())
    return OUTFILE

def gtltcube(label, expcube_dict):
    """gtltcube from Science Tools.

       label: str
          To automatically set the name of the output file
       filter_dict: python dict
          To define all the parameters
    """
    logger.info('Running gtltcube...')
    LABEL = label
    if not os.path.exists(os.path.join(FT_DATA_OUT, 'output_gtltcube')):
        os.makedirs(os.path.join(FT_DATA_OUT, 'output_gtltcube'))
    OUTPATH = os.path.join(FT_DATA_OUT, 'output_gtltcube')
    OUTFILE = os.path.join(OUTPATH, LABEL + '_filtered_gti_ltcube.fits')
    if os.path.exists(OUTFILE):
        logger.info('ATT: Already created %s'%OUTFILE)
        return OUTFILE
    for key in expcube_dict:
        if key == 'outfile':
            if expcube_dict[key] == 'DEFAULT':
                my_apps.expCube['outfile'] = OUTFILE
            else:
                my_apps.expCube[key] = expcube_dict[key]
            continue
        my_apps.expCube[key] = expcube_dict[key]
    my_apps.expCube.run()
    logger.info('Created %s'%OUTFILE)
    logger.info('gtltcube --> CPU time spent: %.2f'%time.clock())
    return OUTFILE

def gtexpcube2(label, expcube2_dict):
    """gtexpcube2 from Science Tools.

       label: str
          To automatically set the name of the output file
       filter_dict: python dict
          To define all the parameters
    """
    logger.info('Running gtexpcube2...')
    LABEL = label
    if not os.path.exists(os.path.join(FT_DATA_OUT, 'output_gtexpcube2')):
        os.makedirs(os.path.join(FT_DATA_OUT, 'output_gtexpcube2'))
    OUTPATH = os.path.join(FT_DATA_OUT, 'output_gtexpcube2')
    OUTFILE = os.path.join(OUTPATH, LABEL + '_expcube.fits')
    if os.path.exists(OUTFILE):
        logger.info('ATT: Already created %s'%OUTFILE)
        return OUTFILE
    for key in expcube2_dict:
        if key == 'outfile':
            if expcube2_dict[key] == 'DEFAULT':
                my_apps.gtexpcube2['outfile'] = OUTFILE
            else:
                my_apps.gtexpcube2[key] = expcube2_dict[key]
            continue
        my_apps.gtexpcube2[key] = expcube2_dict[key]
    my_apps.gtexpcube2.run()
    logger.info('Created %s'%OUTFILE)
    logger.info('gtexpcube2 --> CPU time spent: %.2f'%time.clock())
    return OUTFILE

def gtpsf(gtpsf_dict):
    """gtpsf from Science Tools

       gtpsf_dict: python dict
          To define all the parameters
    """
    expcube = gtpsf_dict['expcube']
    outfile = gtpsf_dict['outfile']
    irfs = gtpsf_dict['irfs']
    evtype = gtpsf_dict['evtype']
    ra = gtpsf_dict['ra']
    dec = gtpsf_dict['dec']
    emin = gtpsf_dict['emin']
    emax = gtpsf_dict['emax']
    nenergies = gtpsf_dict['nenergies']
    thetamax = gtpsf_dict['thetamax']
    ntheta = gtpsf_dict['ntheta']
    os.system('gtpsf expcube=%s outfile=%s irfs=%s evtype=%i ra=%f dec=%f emin=%e emax=%e nenergies=%i thetamax=%i ntheta=%i' \
                  %(expcube, outfile, irfs, evtype, ra, dec, emin, emax, \
                        nenergies, thetamax, ntheta))

def main():
    """Test section.
    """
    import numpy as np
    from GRATools import FT_DATA_FOLDER
    from GRATools.utils.gFTools import mergeft, ebinning_fits_file
    PH = 'photon'
    SC = 'spacecraft'
    PH_FOLDER = os.path.join(FT_DATA_FOLDER, PH)
    SC_FOLDER = os.path.join(FT_DATA_FOLDER, SC)
    FT1_FILE = mergeft(PH_FOLDER, 'FT1_w9-12.txt', 9, 12)
    FILTER_CUT='DATA_QUAL==1&&LAT_CONFIG==1&&LAT_MODE==5&&IN_SAA!=T'+\
               '&&((ABS(ROCK_ANGLE)<52))'
    OUT_FILE_LABEL = 'test'
    FT2_FILE = os.path.join(SC_FOLDER, 'lat_spacecraft_merged.fits')
    EBINNING_FILE = ebinning_fits_file(np.logspace(3, 3.6, 11))
    data_filter = {'infile': FT1_FILE,
                   'emin': 1000,
                   'emax': 5000,
                   'zmax': 90,
                   'evclass': 1024,
                   'evtype': 32,
                   'clobber': 'yes'}
    out_gtselect = gtselect(OUT_FILE_LABEL, data_filter)
    data_mktime = {'evfile': out_gtselect,
                   'scfile': FT2_FILE,
                   'filter': FILTER_CUT,
                   'roicut': 'no',
                   'clobber': 'yes'}
    out_gtmktime = gtmktime(OUT_FILE_LABEL, data_mktime)
    data_bin = {'evfile': out_gtmktime,
                'algorithm': 'HEALPIX',
                'scfile': FT2_FILE,
                'hpx_ordering_scheme': 'RING',
                'hpx_order': 7,
                'coordsys': 'GAL',                                  
                'hpx_ebin': 'yes',
                'ebinalg': 'FILE',
                'ebinfile': EBINNING_FILE,
                'clobber': 'yes'}
    out_gtbin = gtbin(OUT_FILE_LABEL, data_bin)


if __name__ == '__main__':
    main()
