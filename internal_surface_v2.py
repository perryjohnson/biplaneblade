"""A script to merge all the structural parts into one polygon.

Author: Perry Roth-Johnson
Last updated: October 7, 2013

"""


import lib.blade as bl
reload(bl)
from shapely.geometry import Polygon
# from shapely.validation import explain_validity
import matplotlib.pyplot as plt
from operator import attrgetter

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
for station in m.list_of_stations:
    p = station.merge_all_parts()
    good_loops = []
    for n, interior in enumerate(p.interiors):
        a = Polygon(interior).area
        if a > 10e-06:
            good_loops.append(interior)
            # c = interior.centroid.wkt
    print "Station #{0} has {1} good interior loops.".format(station.station_num, len(good_loops))
    # list_of_cx = []
    for loop in good_loops:
        cx = loop.centroid.x
        # list_of_cx.append(cx)
        print "  unsorted x-coord of centroid is {0}".format(cx)
    # sort the loops by the x-coordinate of their centroids, smallest to largest
    good_loops.sort(key=attrgetter('centroid.x'))
    for loop in good_loops:
        cx = loop.centroid.x
        print "  SORTED x-coord of centroid is {0}".format(cx)

# m.plot_selected_cross_sections(plot_edges=False, plot_parts=True)
