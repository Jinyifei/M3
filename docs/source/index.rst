.. M3 documentation master file, created by
   sphinx-quickstart on Sat Dec 14 18:10:49 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


M3
================

*3D Astrophysical Photoionization Code*
---------------------------------------
.. figure:: m3logo6.jpg
   :width: 700
   
   Background Image Credits: NASA, ESA, CSA, STSci, Webb ERO Production Team

**M3**, Messenger Monte-Carlo MAPPINGS, is a 3D Monte-Carlo radiative transfer photoionization code.
It is a descendant of the MAPPINGS family. 
The MAPPINGS photoionization project began in 1976, led by Michael Dopita.
Since the 1980s, MAPPINGS have been re-written and maintained by Ralph Sutherland.
M3 is the latest version of MAPPINGS, incorporating the Monte-Carlo radiative transfer technique established by Lucy (1999).

M3 is a **uniform grid-based**, **MPI** code written in **Fortran77**, can be run on super computer clusters. It is designed for the high-spectral resolution (maximal 10240 frequency bins) and the customized spatial resolution.
To run M3, MPICH and Fortran compilers are required. 



Reference
*********

Please cite the following two publications if you use, modify **M3** for your own research. *Jin et al. 2022* describes the **M3** code. *Sutherland et al. 1993* describes the microphysics used in the **MAPPINGS** serial codes.

[1] Jin, Y.; Kewley, L.; Sutherland, R. 2022, 
https://iopscience.iop.org/article/10.3847/1538-4357/ac48f3/pdf

[2] Sutherland, R.; Dopita, M. 1993, 
https://articles.adsabs.harvard.edu/pdf/1993ApJS...88..253S

.. toctree::
   :titlesonly:
   :hidden:   
   :maxdepth: 10
   
   
   install
   data
   infile
   outfile
   example
