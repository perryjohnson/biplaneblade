"""A script to create the Sandia blade and a biplane blade from their
definition CSV files.

Author: Perry Roth-Johnson
Last updated: August 6, 2013

"""


import scripts.blade as bl
reload(bl)


sandia_flag = False
biplane_flap_sym_no_stagger_flag = True
biplane_flap_sym_stagger_flag = True
biplane_flap_asym_no_stagger_flag = True
biplane_flap_asym_stagger_flag = True

# --- sandia blade ------------------------------------------------------------
if sandia_flag:
    m = bl.MonoplaneBlade('Sandia blade SNL100-00', 'sandia_blade')
    m.copy_all_airfoil_coords()

    # make some plots of the chord and twist schedules
    # m.plot_chord_schedule()
    # m.plot_twist_schedule()

    # pre-process the airfoil coordinates
    for station in m.list_of_stations:
        station.read_airfoil_coords()
        station.scale_airfoil_coords()
        station.split_airfoil_at_LE_and_TE()
        station.part_edges()
        station.find_all_part_cs_coords()

    # create some airfoil plots in Matplotlib
    # station = m.list_of_stations[15]
    # for station in m.list_of_stations:
    #     (fig, axes) = station.create_plot()
    #     station.plot_airfoil_coords(fig, axes, upper_lower_flag=True)
    #     station.plot_part_edges(axes)
    #     station.show_plot()
        # station.save_plot(fig)

    # make a 3D visualization of the entire blade with Mayavi's mlab
    m.plot_blade()


# --- biplane blade, flapwise symmetric, no stagger----------------------------
if biplane_flap_sym_no_stagger_flag:
    b1 = bl.BiplaneBlade('biplane blade, flapwise symmetric, no stagger, rj/R=0.452, g/c=1.25', 'biplane_flap-sym_no-stagger')
    b1.copy_all_airfoil_coords()

    # pre-process the airfoil coordinates
    for station in b1.list_of_stations:
        station.read_airfoil_coords()
        station.scale_airfoil_coords()
    #     # station.split_airfoil_at_LE_and_TE()
    #     # station.part_edges()
    #     # station.find_all_part_cs_coords()

    # # create some airfoil plots in Matplotlib
    # # station = b1.list_of_stations[10]
    # for station in b1.list_of_stations:
    #     (fig, axes) = station.create_plot()
    #     station.plot_airfoil_coords(fig, axes, upper_lower_flag=False)
    # #     station.plot_part_edges(axes)
    #     station.show_plot()
    #     # station.save_plot(fig)

    # make a 3D visualization of the entire blade with Mayavi's mlab
    b1.plot_blade(LE=True, TE=True, twist=False, SW=False, pitch_axis=True)


# --- biplane blade, flapwise symmetric, stagger-------------------------------
if biplane_flap_sym_stagger_flag:
    b2 = bl.BiplaneBlade('biplane blade, flapwise symmetric, stagger, rj/R=0.452, g/c=1.25', 'biplane_flap-sym_stagger')
    b2.copy_all_airfoil_coords()

    # pre-process the airfoil coordinates
    for station in b2.list_of_stations:
        station.read_airfoil_coords()
        station.scale_airfoil_coords()
    #     # station.split_airfoil_at_LE_and_TE()
    #     # station.part_edges()
    #     # station.find_all_part_cs_coords()

    # # create some airfoil plots in Matplotlib
    # # station = b2.list_of_stations[10]
    # for station in b2.list_of_stations:
    #     (fig, axes) = station.create_plot()
    #     station.plot_airfoil_coords(fig, axes, upper_lower_flag=False)
    # #     station.plot_part_edges(axes)
    #     station.show_plot()
    #     # station.save_plot(fig)

    # make a 3D visualization of the entire blade with Mayavi's mlab
    b2.plot_blade(LE=True, TE=True, twist=False, SW=False, pitch_axis=True)


# --- biplane blade, flapwise asymmetric, no stagger----------------------------
if biplane_flap_asym_no_stagger_flag:
    b3 = bl.BiplaneBlade('biplane blade, flapwise asymmetric, no stagger, rj/R=0.452, g/c=1.25', 'biplane_flap-asym_no-stagger')
    b3.copy_all_airfoil_coords()

    # pre-process the airfoil coordinates
    for station in b3.list_of_stations:
        station.read_airfoil_coords()
        station.scale_airfoil_coords()
    #     # station.split_airfoil_at_LE_and_TE()
    #     # station.part_edges()
    #     # station.find_all_part_cs_coords()

    # # create some airfoil plots in Matplotlib
    # # station = b3.list_of_stations[10]
    # for station in b3.list_of_stations:
    #     (fig, axes) = station.create_plot()
    #     station.plot_airfoil_coords(fig, axes, upper_lower_flag=False)
    # #     station.plot_part_edges(axes)
    #     station.show_plot()
    #     # station.save_plot(fig)

    # make a 3D visualization of the entire blade with Mayavi's mlab
    b3.plot_blade(LE=True, TE=True, twist=False, SW=False, pitch_axis=True)


# --- biplane blade, flapwise asymmetric, stagger-------------------------------
if biplane_flap_asym_stagger_flag:
    b4 = bl.BiplaneBlade('biplane blade, flapwise asymmetric, stagger, rj/R=0.452, g/c=1.25', 'biplane_flap-asym_stagger')
    b4.copy_all_airfoil_coords()

    # pre-process the airfoil coordinates
    for station in b4.list_of_stations:
        station.read_airfoil_coords()
        station.scale_airfoil_coords()
    #     # station.split_airfoil_at_LE_and_TE()
    #     # station.part_edges()
    #     # station.find_all_part_cs_coords()

    # # create some airfoil plots in Matplotlib
    # # station = b4.list_of_stations[10]
    # for station in b4.list_of_stations:
    #     (fig, axes) = station.create_plot()
    #     station.plot_airfoil_coords(fig, axes, upper_lower_flag=False)
    # #     station.plot_part_edges(axes)
    #     station.show_plot()
    #     # station.save_plot(fig)

    # make a 3D visualization of the entire blade with Mayavi's mlab
    b4.plot_blade(LE=True, TE=True, twist=False, SW=False, pitch_axis=True)