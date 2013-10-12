"""A script to merge all the structural parts into one polygon at each station.

Author: Perry Roth-Johnson
Last updated: October 11, 2013

"""


import lib.blade as bl
reload(bl)

# --- sandia blade ------------------------------------------------------------
m = bl.MonoplaneBlade('Sandia blade SNL100-00', 'sandia_blade',
    rotate_airfoil_coords=False)

# pre-process the airfoil coordinates
for station in m.list_of_stations:
    # station.find_SW_cs_coords()
    station.airfoil.create_polygon()
    try:
        station.structure.create_all_layers()
    except Warning:
        print "  ...skipping Station #{0}************".format(station.station_num)
        print "length of list_of_polygons =", len(station.structure.list_of_polygons)
        pass
    # station.write_all_part_polygons()
    station.structure.calculate_area()
    # print "Station #{0}, area = {1:6.4f} m^2".format(station.station_num, station.structure.area)
    station.structure.calculate_all_percent_areas()

# m.plot_selected_cross_sections(plot_edges=False, plot_parts=True,
    # selected_stations=[7,11,14,16,18,19,23,26,30,31,32,33])

# stn31 = m.list_of_stations[30]
# stn32 = m.list_of_stations[31]
# stn33 = m.list_of_stations[32]
# stn31.plot_parts()
# stn32.plot_parts()
# stn33.plot_parts()
