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

def iso_parse(isofile):
    """Parsing of the ascii file containing the isotropic emission.
        
        isofile : str
            ascii file containing the isotropic emission.
            It must contain 3 columns: the first with the energy, the second 
            with the flux, the third with the error of the flux.
    """
    f = open(isofile)
    _en, _flux, _fluxerr = [], [], []
    for line in f:
        try:
            en, flux, emeanfluxerr = [float(item) for item in \
                                                line.split()]
            _en.append(en)
            _flux.append(flux)
            _fluxerr.append(emeanfluxerr)
        except:
            pass
    f.close()
    return np.array(_en), np.array(_flux), np.array(_fluxerr)


def csi_parse(csi_file):
    """Parsing of the *_csi.txt files.
    
       csi_file : str
           This file is created by bin/mkcsi.py app.
    """
    logger.info('loading Csi values from %s'%csi_file)
    csi = []
    theta = []
    Rs = []
    emin, emax, emean = [], [], []
    f = open(csi_file, 'r')
    for line in f:
        if 'ENERGY\t' in line:
            e1, e2, em = [float(item) for item in line.split()[1:]]
            emin.append(e1)
            emax.append(e2)
            emean.append(em)
        if 'CSI\t' in line:
            c = np.array([float(item) for item in line.split()[1:]])
            csi.append(c)
        if 'THETA\t' in line:
            th = np.array([float(item) for item in line.split()[1:]])
            theta.append(th)
        if 'R\t' in line:
            r = np.array([float(item) for item in line.split()[1:]])
            Rs.append(r)
    f.close()
    return np.array(emin), np.array(emax), np.array(emean), np.array(csi), np.array(theta), np.array(Rs)

def cpEcross_parse(cp_file):
    """Parsing of the *_cps.txt files

       cp_file : str
           This file is created by bin/mkcpsEcross.py app.
    """
    logger.info('loading Cp values from %s'%cp_file)
    ff = open(cp_file, 'r')
    _emin, _emax = [], []
    _emin1, _emax1 = [], [] 
    _emin2, _emax2 = [], [] 
    _cp, _cperr = [], []
    for line in ff:
        if 'ENERGY' in line:
            emin, emax = [float(item) for item in \
                              line.split()[1:]]
            _emin.append(emin)
            _emax.append(emax)
    ff.close()
    ff = open(cp_file, 'r')
    for line in ff:
        try:
            emin1, emax1, emin2, emax2, cp, cperr = [float(item) \
                                                    for item in line.split()]
            _emin1.append(emin1)
            _emax1.append(emax1)
            _emin2.append(emin2)
            _emax2.append(emax2)
            _cp.append(cp)
            _cperr.append(cperr)
        except:
            pass
    ff.close()
    return  np.array(_emin), np.array(_emax), np.array(_emin1), \
        np.array(_emax1), np.array(_emin2), \
        np.array(_emax2), np.array(_cp), np.array(_cperr)

def cp_parse(cp_file):
    """Parsing of the *_cps.txt files

       cp_file : str
           This file is created by bin/mkcps.py app.
    """
    logger.info('loading Cp values from %s'%cp_file)
    ff = open(cp_file, 'r')
    _emin, _emax, _emean, _cp, _cperr = [], [], [], [], []
    for line in ff:
        try:
            emin, emax, emean, cp, cperr = [float(item) for item in \
                                                line.split()]
            _emin.append(emin)
            _emax.append(emax)
            _emean.append(emean)
            _cp.append(cp)
            _cperr.append(cperr)
        except:
            pass
    ff.close()
    return np.array(_emin), np.array(_emax),np.array(_emean), \
        np.array(_cp), np.array(_cperr)

def clEcross_pol_parse(Ecross_file):
    """Parsing of the *_polspiceEcross.txt files.

       Ecross_file : str
           This file is created by bin/mkpolspiceEcross.py app.
    """
    logger.info('loading Cl values from %s'%Ecross_file)
    f = open(Ecross_file, 'r')
    _emin1, _emax1, _emin2, _emax2 = [], [], [], []
    _ls, _cls, _clserr = [], [], []
    for line in f:
        if 'ENERGY' in line:
            emin1, emax1, emin2, emax2 = [float(item) for item in \
                                              line.split()[1:]]
            _emin1.append(emin1)
            _emax1.append(emax1)
            _emin2.append(emin2)
            _emax2.append(emax2)
        if 'multipole\t' in line:
            l = np.array([float(item) for item in line.split()[1:]])
            _ls.append(l)
        if 'Cl\t' in line:
            cl = np.array([float(item) for item in line.split()[1:]])
            _cls.append(cl)
        if 'Cl_ERR\t' in line:
            clerr = np.array([float(item) for item in line.split()[1:]])
            _clserr.append(clerr)  
    f.close()
    return np.array(_emin1), np.array(_emax1), np.array(_emin2), \
        np.array(_emax2), np.array(_ls), np.array(_cls), np.array(_clserr)

def clfore_parse(clfore_file):
    """Parsing of the *_forecls.txt files.
    
       clfore_file : str
           This file is created by bin/mkpolspiceforecl.py app.
    """
    ls = []
    cls = []
    emin, emax = [], []
    f = open(clfore_file, 'r')
    for line in f:
        if 'ENERGY\t' in line:
            e1, e2 = [float(item) for item in line.split()[1:]]
            emin.append(e1)
            emax.append(e2)
        if 'multipole\t' in line:
            l = np.array([float(item) for item in line.split()[1:]])
            ls.append(l)
        if 'Cl\t' in line:
            cl = np.array([float(item) for item in line.split()[1:]])
            cls.append(cl)
    f.close()
    return np.array(emin), np.array(emax), np.array(ls), np.array(cls)

def cl_pol_parse(cl_file):
    """Parsing of the *_cls.txt files.

       cl_file : str
           This file is created by bin/mkpolspicecl.py app.
    """
    ls = []
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
        if 'multipole\t' in line:
            l = np.array([float(item) for item in line.split()[1:]])
            ls.append(l)
        if 'Cl\t' in line:
            cl = np.array([float(item) for item in line.split()[1:]])
            cls.append(cl)
        if 'Cl_ERR' in line:
            cl_err = np.array([float(item) for item in line.split()[1:]])
            clerrs.append(cl_err)
    f.close()
    return np.array(emin),np.array(emax),np.array(emean),np.array(ls), \
        np.array(cls),  np.array(clerrs)


def get_cl_param(cl_param_file):
    """Parsing of *_parameters.txt files.

       cl_param_file : str
           This file is created by bin/mkdatarestyle.py app.
    """
    logger.info('loading parameters from %s'%cl_param_file)
    ff = open(cl_param_file, 'r')
    _emin, _emax, _emean, _f, _ferr, _cn, _fsky = [], [], [], [], [], \
        [], []
    for line in ff:
        try:
            emin, emax, emean, f, ferr, cn, fsky = [float(item) for item in \
                                                        line.split()]
            _emin.append(emin)
            _emax.append(emax)
            _emean.append(emean)
            _f.append(f)
            _ferr.append(ferr)
            _cn.append(cn)
            _fsky.append(fsky)
        except:
            pass
    ff.close()
    return np.array(_emin), np.array(_emax), np.array(_emean), np.array(_f), \
        np.array(_ferr), np.array(_cn), np.array(_fsky)
        
def get_crbkg(txt_file):
    """Get the CR residual bkg (spline) as a function of the energy
       from the txt files.

       txt_file : str
           It must contain 2 columns: the first one with the energy, the second
           one with the cosmic ray residual background flux.
    """
    logger.info('Getting CR residual bkg from file %s'%txt_file)
    f = open(txt_file, 'r')
    _bkg, _en = [], []
    for line in f:
        try:
            e, bkg = [float(item) for item in line.split()]
            _en.append(e)
            _bkg.append(bkg)
        except:
            pass 
    fmt = dict(xname='Energy', xunits='MeV', 
               yname='E$^{2}$ x CR Residual flux', 
               yunits='MeV cm$^{-2}$ s$^{-1}$ sr$^{-1}$')
    crbkg = xInterpolatedUnivariateSplineLinear(np.array(_en), np.array(_bkg),\
                                                  **fmt)
    f.close()
    return crbkg

def get_energy_from_txt(txt_file, get_binning=False, mean='log'):
    """Returns a list with the center values of the energy bins

       txt_file: str
           txt file produced in output by mkflux.py application
       mean: str
           'log' or 'lin', dependind on the algorithm to use to 
           compute the mean
    """
    f = open(txt_file, 'r')
    _emin, _emax = [], []
    for line in f:
        try:
            emin, emax  = [float(item) for item in line.split()]
            _emin.append(emin)
            _emax.append(emax)
        except:
            pass
    f.close()
    emean = []        
    if mean == 'log':
        for emin, emax in zip(_emin, _emax):
            emean.append(np.sqrt(emin*emax))
    if mean == 'lin':
        for emin, emax in zip(_emin, _emax):
            emean.append(0.5*(emin+emax))
    if get_binning == True:
        return np.array(emean), np.array(_emin), np.array(_emax)
    return np.array(_emin), np.array(_emax), np.array(emean)

def get_energy_from_fits(fits_file, minbinnum=0, maxbinnum=100 ,mean='log'):
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
    _emin = ebounds['E_MIN'][minbinnum:maxbinnum]/1000
    _emax = ebounds['E_MAX'][minbinnum:maxbinnum]/1000
    emean = []        
    if mean == 'log':
        for emin, emax in zip(_emin, _emax):
            emean.append(np.sqrt(emin*emax))
    if mean == 'lin':
        for emin, emax in zip(_emin, _emax):
            emean.append(0.5*(emin+emax))
    f.close()
    return np.array(_emin), np.array(_emax), np.array(emean)

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

def mergeft_P305month(path_to_files, out_file_name, N1month, Nnmonth):
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
    if N1month < 1:
        abort('Invalid number of weeks: the minimun must be > or = to 1')
    if Nnmonth > 108:
        abort('Invalid number of weeks: the maximum must be < or = to 108')
    outtxtfile = os.path.join(path_to_files, out_file_name)
    if not os.path.exists(outtxtfile):
        out_file = open(outtxtfile, 'w')
        for i in range(N1month, Nnmonth+1):
            if i > 1 and i < 10:
                out_file.write("%s/P305_Source_00%i_zmax105.fits \n" \
                                   %(path_to_files,i))
            if i >= 10 and i <= 99:
                out_file.write("%s/P305_Source_0%i_zmax105.fits \n" \
                                   %(path_to_files,i))
            if i > 99:
                out_file.write("%s/P305_Source_%i_zmax105.fits \n" \
                                   %(path_to_files,i))
        out_file.close()
    logger.info('Created %s...' %outtxtfile)
    return '@'+outtxtfile

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
    if Nnweek > 486:
        abort('Invalid number of weeks: the maximum must be < or = to 486')
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
    
def best_fit(X, Y):

    xbar = sum(X)/len(X)
    ybar = sum(Y)/len(Y)
    n = len(X) # or len(Y)

    numer = sum([xi*yi for xi,yi in zip(X, Y)]) - n * xbar * ybar
    denum = sum([xi**2 for xi in X]) - n * xbar**2

    b = numer / denum
    a = ybar - b * xbar

    print('best fit line:\ny = {:.2f} + {:.2f}x'.format(a, b))

    return a, b

def main():
    """Test Session
    """


if __name__ == '__main__':
    main()
