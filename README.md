bladedesign
===========

tool to set up structural models of [biplane] wind turbine blades

Current workflow (as of April 10, 2014)
---------------------------------------
1. run `path_to_blade/create_stnXX_mesh.py` - write initial TrueGrid input file with boundary curves: `mesh_stnXX_start.tg`
2. manually edit `mesh_stnXX_start.tg` to create block meshes fitted to boundary curves; save as `mesh_stnXX_finish.tg`
3. run TrueGrid on `mesh_stnXX_finish.tg` to write ABAQUS output file: `mesh_stnXX.abq`
4. run `path_to_blade/layer_plane_angles_stnXX.py` - write updated grid object to VABS input file: `mesh_stnXX.vabs`
5. run `path_to_blade/run_vabs.py` - use VABS to calculate mass and stiffness matrices
6. `mesh_stnXX.vabs.K` - mass and stiffness matrices are in this file!
7. run `path_to_blade/plot_MK.py` - plot VABS data vs. Sandia published data
8. run `path_to_blade/write_DYMORE_input_file.py` - write VABS output for a DYMORE input file to `sandia_blade_OD.dat` and `sandia_blade_MK.dat`
9. manually copy contents of `sandia_blade_OD.dat` and `sandia_blade_MK.dat` into `sandia_blade.dat`.
10. run `rundymore.bat` to load the structural model
11. run `plot_DYMORE_results.py` to postprocess results in `FIGURES` directory
12. run `clean.bat` to erase all DYMORE results


Plan forward (as of April 10, 2014)
-----------------------------------
Create a beam model of the biplane blade (flapwise symmetric, no stagger configuration). Reuse cross-section meshes at stations 1-9 and 25-34 from the Sandia blade. Make new meshes for stations 10-24, then create new DYMORE input files for the biplane blade.

In the end, you should have an archive of all the blade input files that looks like:

* `biplane_blade/`
  * `blade_definition.csv`
  * `materials.csv`
  * `airfoils/`
    * `Cylinder.txt`
    * ...
    * `NACA_64-618.txt`
  * `stn01/`
    * `mesh_stn01_start.tg` (initial TrueGrid input file with part boundary curves)
    * `mesh_stn01_final.tg` (final TrueGrid input file with grids inside curves)
    * `mesh_stn01.abq` (exported grid file in ABAQUS format)
    * `mesh_stn01.vabs` (VABS input file of geometry and materials)
    * `mesh_stn01.vabs.K` (VABS output of mass and stiffness matrices)
  * `stn02/`
  * ...
  * `stn34/`
* `biplane_blade_lib/`
  * `plot_selected_stations.py`
  * `prep_stn01_mesh.py`
  * ...
  * `prep_stn34_mesh.py`
  * `layer_plane_angles_stn01.py`
  * ...
  * `layer_plane_angles_stn34.py`
  * `run_vabs.py`
  * `plot_MK.py`
  * `write_DYMORE_input_file.py`
  * `plot_DYMORE_results.py`
* `sandia_blade/`
  * `blade_definition.csv`
  * `materials.csv`
  * `airfoils/`
    * `Cylinder.txt`
    * ...
    * `NACA_64-618.txt`
  * `stn01/`
    * `mesh_stn01_start.tg` (initial TrueGrid input file with part boundary curves)
    * `mesh_stn01_final.tg` (final TrueGrid input file with grids inside curves)
    * `mesh_stn01.abq` (exported grid file in ABAQUS format)
    * `mesh_stn01.vabs` (VABS input file of geometry and materials)
    * `mesh_stn01.vabs.K` (VABS output of mass and stiffness matrices)
  * `stn02/`
  * ...
  * `stn34/`
* `sandia_blade_lib/`
  * `plot_selected_stations.py`
  * `prep_stn01_mesh.py`
  * ...
  * `prep_stn34_mesh.py`
  * `layer_plane_angles_stn01.py`
  * ...
  * `layer_plane_angles_stn34.py`
  * `run_vabs.py`
  * `plot_MK.py`
  * `write_DYMORE_input_file.py`
  * `plot_DYMORE_results.py`
