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


"""Analysis utilities 
"""

  
import os
import numpy as np
import pyfits as pf
from GRATools import GRATOOLS_OUT
from GRATools.utils.logging_ import logger, abort

def get_energy_from_fits(fits_file, mean='log'):
    """Returns a list with the center values of the energy bins

       fits_file: str
           fits file, usually we want to do to this at the level of
           gtbin output file
       mean: str
           'log' or 'lin', dependind on the algorithm to use to 
           compute the mean
    """
    f = pf.open(fits_file)
    ebounds = f[2].data
    emin_ = ebounds['E_MIN']
    emax_ = ebounds['E_MAX']
    emean = []        
    if mean == 'log':
        for emin, emax in zip(emin_, emax_):
            emean.append(np.sqrt(emin*emax))
    if mean == 'lin':
        for emin, emax in zip(emin_, emax_):
            emean.append(0.5*(emin+emax))
    return emean

def ebinning_fits_file(ebinning_array):
    """Produces a fits file defining the enrgy binning to fed gtbin.

       ebinning_array: numpy array
           array in which the energy binnin is defined.
    """
    txt_file_name = os.path.join(GRATOOLS_OUT,'ebinning.txt')
    txt_file = open(txt_file_name,'w')
    fits_file = os.path.join(GRATOOLS_OUT,'ebinning.fits')
    for emin, emax in zip(ebinning_array[:-1], ebinning_array[1:]):
        txt_file.write('%.4f %.4f\n'%(emin, emax))
    txt_file.close()
    os.system('gtbindef bintype=E binfile=%s outfile=%s energyunits=MeV' \
                  %(txt_file_name, fits_file))
    logger.info('Created %s...'%fits_file)
    return fits_file
    

def mergeft(path_to_files, out_file_name, N1week, Nnweek):
    """creates a .txt file with the list of the FT files to merge.

       path_to_files: str
           path where datat files are stored
       out_file_name: str
           name of the txt output file (created in the same folder of data)
       N1week: int
           number of the starting week
       Nnweek: int
           number of the ending week
    """
    if N1week < 9:
        abort('Invalid number of weeks: the minimun must be > or = to 9')
    if Nnweek > 397:
        abort('Invalid number of weeks: the maximum must be < or = to 397')
    outtxtfile = os.path.join(path_to_files, out_file_name)
    if not os.path.exists(outtxtfile):
        out_file = open(outtxtfile, 'w')
        for i in range(N1week, Nnweek+1):
            if i == 9:
                out_file.write("%s/lat_photon_weekly_w00%i_p302_v001.fits \n" \
                                   %(path_to_files,i))
            if i >= 10 and i <= 99:
                out_file.write("%s/lat_photon_weekly_w0%i_p302_v001.fits \n" \
                                   %(path_to_files,i))
            if i > 99:
                out_file.write("%s/lat_photon_weekly_w%i_p302_v001.fits \n" \
                                   %(path_to_files,i))
        out_file.close()
    logger.info('Created %s...' %outtxtfile)
    return '@'+outtxtfile
    
