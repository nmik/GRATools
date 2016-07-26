To use this package:
--------------------

0) - Install Fermi Science Tools
   - Install python 2.7 (or later)
   - Install HEASARC FTools
   - Install matplotlib
   - Install pyRoot
   - Install pyfits
   - Install healpy

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
example:

python bin/mkdataselection.py --config config/1yr_st_aniso_config.py

python bin/mkflux.py --infile output/1yr_UCV_t56_outfiles.txt --config config/1yr_flux_config.py

python bin/mkanalysis.py --config config/allyrs_analysis.py
