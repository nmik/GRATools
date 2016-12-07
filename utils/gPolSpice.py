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
import healpy as hp
from scipy import optimize as opt
from GRATools import GRATOOLS_OUT
from GRATools import GRATOOLS_CONFIG
from GRATools import FT_DATA_FOLDER
from GRATools.utils.logging_ import logger, abort
from GRATools.utils.matplotlib_ import pyplot as plt
from GRATools.utils.gSpline import xInterpolatedUnivariateSplineLinear
from GRATools.utils.gSpline import xInterpolatedBivariateSplineLinear

def pol_create_config(pol_dict, config_file_name):
    """Creates PolSpice config file
    """
    pol_config = os.path.join(GRATOOLS_CONFIG, 'ascii/%s.txt'%config_file_name)
    pol_config_file = open(pol_config, 'w')
    for key in pol_dict:
        pol_config_file.write('%s = %s \n'%(key, str(pol_dict[key])))
    logger.info('Created config/ascii/%s.txt'%config_file_name)
    return pol_config 

def pol_run(config_file, spice_version='v03-02-00', spice_path='/opt'):
    """Runs PolSpice
    """
    os.system('%s/PolSpice_%s/src/spice -optinfile %s'
              %(spice_path, spice_version, config_file))

def pol_cl_parse(pol_cl_out_file):
    """
    """
    f = open(pol_cl_out_file, 'r')
    _l, _cl = [], []
    for line in f:
        try:
            l, cl = [float(item) for item in line.split()]
            _l.append(l)
            _cl.append(cl)
        except:
            pass
    return np.array(_l), np.array(_cl)

def pol_cov_parse(pol_cov_out_file):
    """deve salvare anche una figura della matrice
    """
    hdu = pf.open(pol_cov_out_file)
    _cov = hdu[0].data
    _i = np.arange(len(_cov[0]))
    _j =np.arange(len(_cov[0][0]))
    fmt = dict(yname='l$_{i}$', yunits='', xname='l$_{j}$',
               xunits='', zname='VAR')
    _cov_spline = xInterpolatedBivariateSplineLinear(_i, _j, _cov[0], **fmt)
    hdu.close()
    _clerr = np.array([np.sqrt(_cov[0][i][i]) for i in _i])
    return _clerr


def pol_cl_calculation(pol_dict, config_file_name):
    """
    """
    pol_cl_out_file = pol_dict['clfile']
    pol_cov_out_file = pol_dict['covfileout']
    config_file = pol_create_config(pol_dict, config_file_name)
    pol_run(config_file)
    _l, _cl = pol_cl_parse(pol_cl_out_file)
    _clerr = pol_cov_parse(pol_cov_out_file)
    return _l, _cl, _clerr
    
def main():
    """
    """
    cl_f = 'output/output_pol/mask_check_cl.txt'
    cov_f = 'output/output_pol/mask_check_cov.fits'
    l, cl = pol_cl_parse(cl_f)
    clerr = pol_cov_parse(cov_f)
    #plt.errorbar(l, cl, fmt='o', markersize=3, elinewidth=1, yerr=clerr)
    plt.plot(l, cl, 'o', markersize=3)
    plt.xscale('log')
    plt.yscale('log')
    #plt.ylim(-1e-15, 1e-15)
    plt.show()


if __name__ == '__main__':
    main()
