"""A script to create the Sandia blade and 4 biplane blades.

Five directories that contain blade definitions are located in the same
directory as this 'create_blades.py' script:
  sandia_blade/
    airfoils/
    blade_definition.csv
  biplane_blade/
    airfoils/
    blade_definition.csv
Each 'airfoils/' sub-directory contains a several text files for different 
airfoils, each of which list the corresponding airfoil coordinates.

Usage
-----
start an IPython (qt)console with the pylab flag:
$ ipython qtconsole --pylab
or
$ ipython --pylab
Then, from the prompt, run this script:
|> %run create_blades
Once you are finished looking at the blades, you can clean up extra files:
|> %run clean
(See the 'clean.py' script in this directory for details.)

Author: Perry Roth-Johnson
Last updated: April 10, 2014

"""


import lib.blade as bl
reload(bl)


biplane_flap_sym_no_stagger_flag = True
sandia_flag = False

# --- biplane blade, flapwise symmetric, no stagger----------------------------
if biplane_flap_sym_no_stagger_flag:
    b1 = bl.BiplaneBlade(
        'biplane blade, flapwise symmetric, no stagger, rj/R=0.452, g/c=1.25',
        'biplane_blade')
    b1.copy_all_airfoil_coords()

    # pre-process the airfoil coordinates
    for station in b1.list_of_stations:
        station.airfoil.create_polygon()
        # station.structure.create_all_layers()
        # station.structure.write_all_part_polygons()

    # make a 3D visualization of the entire blade with Mayavi's mlab
    b1.plot_blade(stn_nums=True, twist=True, export=False)

# --- sandia blade ------------------------------------------------------------
if sandia_flag:
    m = bl.MonoplaneBlade('Sandia blade SNL100-00', 'sandia_blade',
        rotate_airfoil_coords=False)

    # pre-process the airfoil coordinates
    for station in m.list_of_stations:
        station.airfoil.create_polygon()
        station.structure.create_all_layers()
        station.structure.write_all_part_polygons()

    # create some airfoil plots in Matplotlib
    # m.plot_selected_cross_sections(plot_edges=False, plot_parts=True)

    # calculate and plot blade quantities
    # m.plot_chord_schedule()
    # m.plot_twist_schedule()
    # m.plot_mass_schedule()
    # m.plot_percent_areas()
    # m.plot_percent_masses()
    # m.calculate_blade_mass()

    # make a 3D visualization of the entire blade with Mayavi's mlab
    for station in m.list_of_stations:
        station.find_SW_cs_coords()
    m.plot_blade(stn_nums=True, twist=True, export=False)


