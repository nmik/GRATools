To use this package:
--------------------

0) - Install Fermi Science Tools
   - Install python 2.7 (or later)
   - Install HEASARC FTools (https://heasarc.gsfc.nasa.gov/ftools/)
   - Install matplotlib
   - Install Healpy (https://healpy.readthedocs.io/en/latest/)
   - Install PolSpice (http://www2.iap.fr/users/hivon/software/PolSpice/README.html)

1) Set the environment:
  
    export PYTHONPATH=:/path/to/this/package/${PYTHONPATH}
    export PATH=/path/to/this/package/GRATools/apps:${PATH}

2) Change the directory where the data files are stored in __init__.py.

   line 36:
   	FT_DATA_FOLDER = '/path/to/data-files'

   In this directory should there be the following folders:
        photon/      -> where FT1 files are stored
	spacecraft/  -> where FT2 files are stored
	catalogs/    -> where source catalog files are stored
	output/      -> where ST outputs will be stored


How to run the analysis:
------------------------
This will show all the possible settings of a given function in bin/ :

     python bin/ANYFUNCTION.py --h 


examples:

	python bin/mkdataselection.py --config config/1yr_st_aniso_config.py

	python bin/mkdatarestyle.py --config config/restyle_config.py

	python bin/mapview.py --infile output/output_flux/Allyrs_UCV_t56_flux_549-1000.fits --udgrade 128
