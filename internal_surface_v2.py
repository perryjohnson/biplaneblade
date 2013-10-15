"""A script to merge all the structural parts into one polygon at each station.

Author: Perry Roth-Johnson
Last updated: October 11, 2013

"""


import lib.blade as bl
reload(bl)
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

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
        print "length of _list_of_layers =", len(station.structure._list_of_layers)
        pass
    # station.write_all_part_polygons()

# m.plot_selected_cross_sections(plot_edges=False, plot_parts=True,
    # selected_stations=[7,11,14,16,18,19,23,26,30,31,32,33])

# make a stacked bar plot of all the percent areas for each structural part
m.plot_percent_areas()

# make a stacked bar plot of all the percent masses for each structural part
m.plot_percent_masses()
