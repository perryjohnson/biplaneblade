bladedesign
===========

tool to set up structural models of [biplane] wind turbine blades

Usage: `create_blades.py`
-------------------------
start an IPython (qt)console with the pylab flag:  
`$ ipython qtconsole --pylab`  
Then, from the prompt, run this script:  
`|> %run create_blades`  
Once you are finished looking at the blades, you can clean up extra files:  
`|> %run clean`  

Usage: `create_meshes.py`
-------------------------
start an IPython (qt)console with the pylab flag:  
`$ ipython qtconsole --pylab`  
Then, from the prompt, run this script:  
`|> %run create_meshes`  

Usage: `create_stn16_mesh.py`
-------------------------
start an IPython (qt)console with the pylab flag:  
`$ ipython qtconsole --pylab`  
Then, from the prompt, run this script:  
`|> %run create_stn16_mesh`  

List of top-level scripts
-------------------------
* `create_blades.py` - Create the Sandia blade and 4 biplane blades
* `create_meshes.py` - Write initial TrueGrid files for Sandia blade stations
* `create_stn16_mesh.py` - Write initial TrueGrid file for Sandia station #16

List of cleanup scripts
-----------------------
* `clean.py` - Delete leftover files and paths in each `blade_path`
* `clear_logs.py` - Delete the contents of all the log files

Plan forward (as of March 6, 2014)
----------------------------------
Instead of automatically generating meshes with Python, let's only use Python to automatically generate boundary curves for each part in a blade station. (I think I already have this functionality.) Then, manually mesh each station in TrueGrid. (Start with Station #1 and work your way out to Station #34.) Save the final `*.tg` file in an archive of all the blade input files. Export the grid in ABAQUS format and use Python to tranlate it to VABS format. Use VABS to generate the mass and stiffness matrix. This will be tedious, but straightforward. In the end, you should have an archive of all the blade input files that looks like:

`sandia_blade/`
--`blade_definition.csv`
--`materials.csv`
--`airfoils/`
----`Cylinder.txt`
----...
----`NACA_64-618.txt`
--`stn01/`
----`mesh_start.tg` (initial TrueGrid input file with part boundary curves)
----`mesh_final.tg` (final TrueGrid input file with grids inside curves)
----`mesh_abq.txt` (exported grid file in ABAQUS format)
----`mesh_vabs.dat` (VABS input file of geometry and materials)
----`mesh_vabs.dat.K` (VABS output of mass and stiffness matrices)
--`stn02/`
--...
--`stn34/`
