"""A script to merge all the structural parts into one polygon.

Author: Perry Roth-Johnson
Last updated: October 7, 2013

"""


import lib.blade as bl
reload(bl)
from shapely.geometry import Polygon
# from shapely.validation import explain_validity
import matplotlib.pyplot as plt

# --- sandia blade ------------------------------------------------------------
m = bl.MonoplaneBlade('Sandia blade SNL100-00', 'sandia_blade')
m.copy_all_airfoil_coords()

# pre-process the airfoil coordinates
for station in m.list_of_stations:
    station.airfoil.read_coords()
    station.airfoil.scale_and_translate_coords()
    station.airfoil.split_at_LE_and_TE()
    station.airfoil.make_polygon()
    station.find_part_edges()
    # station.find_all_part_cs_coords()
    station.find_all_part_polygons()

# check that parts are getting merged properly in stations that are plotting weird
for station in m.list_of_stations[:-1]:
    p = station.merge_all_parts()
    print "Station #{0} has {1} interior loops.".format(station.station_num, len(p.interiors))
    for n, interior in enumerate(p.interiors):
        c = interior.centroid.wkt
        print "  Interior #{0}, centroid at {1}".format(n,c)
# stn34 = m.list_of_stations[-1]
# p = stn34.merge_all_parts(plot_flag=True)
# stn34.plot_parts()

# m.plot_selected_cross_sections(plot_edges=False, plot_parts=True)
