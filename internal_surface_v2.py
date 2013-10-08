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
for station in m.list_of_stations:
    station.airfoil.read_coords()
    station.airfoil.scale_and_translate_coords()
    station.airfoil.split_at_LE_and_TE()
    station.airfoil.make_polygon()
    station.find_part_edges()
    # station.find_all_part_cs_coords()
    station.find_all_part_polygons()

# plot station 14
stn14 = m.list_of_stations[13]
p = m.plot_merged_parts(stn14)
i1 = p.interiors[0]
i2 = p.interiors[1]
i3 = p.interiors[2]
i4 = p.interiors[3]

int_surfs = []
for interior in p.interiors:
    op = Polygon(interior)
    ip = op.buffer(-0.01)  # thickness of internal surface
    int_surfs.append(op.difference(ip))
stn14.structure.internal_surface_1.polygon = int_surfs[0]
stn14.structure.internal_surface_2.polygon = int_surfs[1]
stn14.structure.internal_surface_3.polygon = int_surfs[2]
stn14.structure.internal_surface_4.polygon = int_surfs[3]
m.plot_parts(stn14)