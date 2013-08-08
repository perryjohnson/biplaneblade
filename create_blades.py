"""A script to create the Sandia blade and 4 biplane blades.

Five directories that contain blade definitions are located in the same
directory as this 'create_blades.py' script:
  sandia_blade/
    airfoils/
    blade_definition.csv
  biplane_flap-sym_no-stagger/
    airfoils/
    blade_definition.csv
  biplane_flap-sym_stagger/
    airfoils/
    blade_definition.csv
  biplane_flap-asym_no-stagger/
    airfoils/
    blade_definition.csv
  biplane_flap-asym_stagger/
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
Last updated: August 8, 2013

"""


import scripts.blade as bl
reload(bl)


sandia_flag = True
biplane_flap_sym_no_stagger_flag = True
biplane_flap_sym_stagger_flag = False
biplane_flap_asym_no_stagger_flag = False
biplane_flap_asym_stagger_flag = False

# --- sandia blade ------------------------------------------------------------
if sandia_flag:
    m = bl.MonoplaneBlade('Sandia blade SNL100-00', 'sandia_blade')
    m.copy_all_airfoil_coords()

    # make some plots of the chord and twist schedules
    # m.plot_chord_schedule()
    # m.plot_twist_schedule()

    # pre-process the airfoil coordinates
    for station in m.list_of_stations:
        station.airfoil.read_coords()
        station.airfoil.scale_and_translate_coords()
        station.airfoil.split_at_LE_and_TE()
        # station.part_edges()
        # station.find_all_part_cs_coords()

    # create some airfoil plots in Matplotlib
    # station = m.list_of_stations[10]
    for station in m.list_of_stations:
        (fig, axes) = station.create_plot()
        station.airfoil.plot_coords(fig, axes, split_flag=True)
        # station.plot_part_edges(axes)
        station.show_plot()
        # station.save_plot(fig)

    # make a 3D visualization of the entire blade with Mayavi's mlab
    # m.plot_blade(SW=False)


# --- biplane blade, flapwise symmetric, no stagger----------------------------
if biplane_flap_sym_no_stagger_flag:
    b1 = bl.BiplaneBlade(
        'biplane blade, flapwise symmetric, no stagger, rj/R=0.452, g/c=1.25',
        'biplane_flap-sym_no-stagger')
    b1.copy_all_airfoil_coords()

    # pre-process the airfoil coordinates
    for station in b1.list_of_stations:
        station.airfoil.read_coords()
        station.airfoil.scale_and_translate_coords()
        station.airfoil.split_at_LE_and_TE()
        # station.part_edges()
        # station.find_all_part_cs_coords()

    # create some airfoil plots in Matplotlib
    # station = b1.list_of_stations[10]
    for station in b1.list_of_stations:
        (fig, axes) = station.create_plot()
        station.airfoil.plot_coords(fig, axes, split_flag=True)
        # station.plot_part_edges(axes)
        station.show_plot()
        # station.save_plot(fig)

    # make a 3D visualization of the entire blade with Mayavi's mlab
    # b1.plot_blade(LE=True, TE=True, twist=True, SW=False, pitch_axis=True)


# --- biplane blade, flapwise symmetric, stagger-------------------------------
if biplane_flap_sym_stagger_flag:
    b2 = bl.BiplaneBlade(
        'biplane blade, flapwise symmetric, stagger, rj/R=0.452, g/c=1.25', 
        'biplane_flap-sym_stagger')
    b2.copy_all_airfoil_coords()

    # pre-process the airfoil coordinates
    for station in b2.list_of_stations:
        station.airfoil.read_coords()
        station.airfoil.scale_coords()
        # station.airfoil.split_at_LE_and_TE()
    #     # station.part_edges()
    #     # station.find_all_part_cs_coords()

    # # create some airfoil plots in Matplotlib
    # # station = b2.list_of_stations[10]
    # for station in b2.list_of_stations:
    #     (fig, axes) = station.create_plot()
    #     station.plot_airfoil_coords(fig, axes, split_flag=False)
    # #     station.plot_part_edges(axes)
    #     station.show_plot()
    #     # station.save_plot(fig)

    # make a 3D visualization of the entire blade with Mayavi's mlab
    # b2.plot_blade(LE=True, TE=True, twist=True, SW=False, pitch_axis=True)


# --- biplane blade, flapwise asymmetric, no stagger----------------------------
if biplane_flap_asym_no_stagger_flag:
    b3 = bl.BiplaneBlade(
        'biplane blade, flapwise asymmetric, no stagger, rj/R=0.452, g/c=1.25', 
        'biplane_flap-asym_no-stagger')
    b3.copy_all_airfoil_coords()

    # pre-process the airfoil coordinates
    for station in b3.list_of_stations:
        station.airfoil.read_coords()
        station.airfoil.scale_coords()
        # station.airfoil.split_at_LE_and_TE()
    #     # station.part_edges()
    #     # station.find_all_part_cs_coords()

    # # create some airfoil plots in Matplotlib
    # # station = b3.list_of_stations[10]
    # for station in b3.list_of_stations:
    #     (fig, axes) = station.create_plot()
    #     station.plot_airfoil_coords(fig, axes, split_flag=False)
    # #     station.plot_part_edges(axes)
    #     station.show_plot()
    #     # station.save_plot(fig)

    # make a 3D visualization of the entire blade with Mayavi's mlab
    # b3.plot_blade(LE=True, TE=True, twist=True, SW=False, pitch_axis=True)


# --- biplane blade, flapwise asymmetric, stagger-------------------------------
if biplane_flap_asym_stagger_flag:
    b4 = bl.BiplaneBlade(
        'biplane blade, flapwise asymmetric, stagger, rj/R=0.452, g/c=1.25', 
        'biplane_flap-asym_stagger')
    b4.copy_all_airfoil_coords()

    # pre-process the airfoil coordinates
    for station in b4.list_of_stations:
        station.airfoil.read_coords()
        station.airfoil.scale_coords()
        # station.airfoil.split_at_LE_and_TE()
    #     # station.part_edges()
    #     # station.find_all_part_cs_coords()

    # # create some airfoil plots in Matplotlib
    # # station = b4.list_of_stations[10]
    # for station in b4.list_of_stations:
    #     (fig, axes) = station.create_plot()
    #     station.plot_airfoil_coords(fig, axes, split_flag=False)
    # #     station.plot_part_edges(axes)
    #     station.show_plot()
    #     # station.save_plot(fig)

    # make a 3D visualization of the entire blade with Mayavi's mlab
    # b4.plot_blade(LE=True, TE=True, twist=True, SW=False, pitch_axis=True)