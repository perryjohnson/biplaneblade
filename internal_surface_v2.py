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

m.plot_selected_cross_sections(plot_edges=False, plot_parts=True,
    selected_stations=[7,11,14,16,18,19,23,26,30,31,32,33])

# just have to fix polygon merging with lower spar cap and TE reinforcement in stn #32
# --PRJ, October 11, 2013