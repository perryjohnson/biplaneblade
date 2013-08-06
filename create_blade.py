"""A script to create the Sandia blade from its definition CSV file.

Author: Perry Roth-Johnson
Last updated: August 5, 2013

"""


import scripts.blade as bl
reload(bl)


b = bl.Blade('Sandia blade SNL100-00', 'sandia_blade')
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
