"""A script to create the Sandia blade and a biplane blade from their
definition CSV files.

Author: Perry Roth-Johnson
Last updated: August 6, 2013

"""


import scripts.blade as bl
reload(bl)


sandia_flag = False
biplane_flag = True

# --- sandia blade ------------------------------------------------------------
if sandia_flag:
    b = bl.MonoplaneBlade('Sandia blade SNL100-00', 'sandia_blade')
    b.copy_all_airfoil_coords()

    # make some plots of the chord and twist schedules
    # b.plot_chord_schedule()
    # b.plot_twist_schedule()

    # pre-process the airfoil coordinates
    for station in b.list_of_stations:
        station.read_airfoil_coords()
        station.scale_airfoil_coords()
        station.split_airfoil_at_LE_and_TE()
        station.part_edges()
        station.find_all_part_cs_coords()

    # create some airfoil plots in Matplotlib
    # station = b.list_of_stations[15]
    # for station in b.list_of_stations:
    #     (fig, axes) = station.create_plot()
    #     station.plot_airfoil_coords(fig, axes, upper_lower_flag=True)
    #     station.plot_part_edges(axes)
    #     station.show_plot()
        # station.save_plot(fig)

    # make a 3D visualization of the entire blade with Mayavi's mlab
    b.plot_blade()


# --- biplane blade -----------------------------------------------------------
if biplane_flag:
    a = bl.BiplaneBlade('biplane blade, rj/R=0.452, g/c=1.25', 'biplane_blade')
    a.copy_all_airfoil_coords()

    # pre-process the airfoil coordinates
    for station in a.list_of_stations:
        station.read_airfoil_coords()
        station.scale_airfoil_coords()
    #     # station.split_airfoil_at_LE_and_TE()
    #     # station.part_edges()
    #     # station.find_all_part_cs_coords()

    # # create some airfoil plots in Matplotlib
    # # station = a.list_of_stations[10]
    # for station in a.list_of_stations:
    #     (fig, axes) = station.create_plot()
    #     station.plot_airfoil_coords(fig, axes, upper_lower_flag=False)
    # #     station.plot_part_edges(axes)
    #     station.show_plot()
    #     # station.save_plot(fig)

    # make a 3D visualization of the entire blade with Mayavi's mlab
    a.plot_blade(LE=False, TE=False, twist=False, SW=False, pitch_axis=True)
