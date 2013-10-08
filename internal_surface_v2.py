"""A script to merge all the structural parts into one polygon.

Author: Perry Roth-Johnson
Last updated: October 7, 2013

"""


import lib.blade as bl
reload(bl)
from shapely.geometry import Polygon
import matplotlib.pyplot as plt

# --- sandia blade ------------------------------------------------------------
m = bl.MonoplaneBlade('Sandia blade SNL100-00', 'sandia_blade')
m.copy_all_airfoil_coords()

# pre-process the airfoil coordinates
for station in m.list_of_stations[:-1]:
    station.airfoil.read_coords()
    station.airfoil.scale_and_translate_coords()
    station.airfoil.split_at_LE_and_TE()
    station.airfoil.make_polygon()
    station.find_part_edges()
    # station.find_all_part_cs_coords()
    station.find_all_part_polygons()

# # plot station 14
# stn14 = m.list_of_stations[13]
# m.plot_parts(stn14)

m.plot_selected_cross_sections(plot_edges=False, plot_parts=True)
