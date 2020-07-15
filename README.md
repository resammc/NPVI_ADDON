
# NPVI add-on for igaFEAP
[![DOI](https://zenodo.org/badge/279654786.svg)](https://zenodo.org/badge/latestdoi/279654786)

*Latest version: v1.4 (14 July 2020)*

It is a common practice in the isogeometric analysis post-processing to use FE-like mesh formats to plot the results. In the first step, an FE-mesh corresponding to the IGA-mesh is created. After that simulation results such as displacements, stresses, principal stresses, etc. are projected onto the new mesh. The code to perform these actions is already included in the NURBFEAP/igaFEAP for the linear case. We have used the same approach and fitted the algorithms for our purposes. A similar routine to `uparaview.f` is used to write the output file.

The routine includes the following features:
- linear interpolation of the geometry and the solution field (using linear FE elements)
- quadratic interpolation of the geometry and the solution field (using quadratic FE elements)
- uniform and non-uniform subdivisions of the NURBS elements
- standard FEAP output of the IGA mesh (in form of Serendipity and Lagrangian finite elements)
- **NEW** post-processing of NURBS meshes using the so-called Bézier elements in Paraview (only for FEAP 8.5 and FEAP 8.6)
- **NEW** compress VTU files on-demand (currently by executing an external Python script)

# FEAP
More info:

[FEAP - Official Webpage](http://projects.ce.berkeley.edu/feap/)

[NPVI Page in FEAP Forum](http://feap.berkeley.edu/forum/index.php?topic=1542.0)

# Contact
Reporting bugs and well as suggesting new features are always welcome. Contact me at resam.makvandi@ovgu.de

# Version history
## Ver. 1.0
Release date: 07 April 2017
## Ver. 1.1
Release date: 03 July 2017
- fixed a problem in defining the xl array 
## Ver. 1.2
Release date: 27 October 2017
- fixed an error when using different types of elements in a single problem
## Ver. 1.3 (major update)
Release date: 03 March 2018
- added support for FEAP 8.5
- there was a problem with stress/strain projection rotuines (igafeap) which is fixed now
- most of the routines are re-written (*_new.f files) based on igafeap plot routines. The performance (and stability) has increased significantly.
- more bug fixes
## Ver. 1.4 (current version)
Release date: 14 July 2020
- added support for FEAP 8.6
- post-processing of NURBS meshes using the so-called Bézier elements in Paraview ('npvi,[ ],4')
- added 'npvi,comp' to generate a Python script to compress the created VTU files
- bug fixes
