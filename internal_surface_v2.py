"""A script to merge all the structural parts into one polygon at each station.

Author: Perry Roth-Johnson
Last updated: October 11, 2013

"""


import lib.blade as bl
reload(bl)
import pandas as pd
import os

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

list_of_dicts = []
for station in m.list_of_stations:
    d = station.structure.calculate_all_percent_areas()
    list_of_dicts.append(d)

percent_areas = pd.DataFrame(list_of_dicts, index=range(1,m.number_of_stations+1))
percent_areas_path = os.path.join(m.blade_path, 'percent_areas.csv')
percent_areas.to_csv(percent_areas_path, index_label='blade station', cols=[
    'external surface (gelcoat)',
    'external surface (triax)',
    'root buildup',
    'spar cap (lower)',
    'spar cap (upper)',
    'aft panel 1 (lower)',
    'aft panel 1 (upper)',
    'aft panel 2 (lower)',
    'aft panel 2 (upper)',
    'LE panel',
    'shear web 1 (left biax)',
    'shear web 1 (foam)',
    'shear web 1 (right biax)',
    'shear web 2 (left biax)',
    'shear web 2 (foam)',
    'shear web 2 (right biax)',
    'shear web 3 (left biax)',
    'shear web 3 (foam)',
    'shear web 3 (right biax)',
    'TE reinforcement (uniax)',
    'TE reinforcement (foam)',
    'internal surface 1 (triax)',
    'internal surface 1 (resin)',
    'internal surface 2 (triax)',
    'internal surface 2 (resin)',
    'internal surface 3 (triax)',
    'internal surface 3 (resin)',
    'internal surface 4 (triax)',
    'internal surface 4 (resin)'])


# m.plot_selected_cross_sections(plot_edges=False, plot_parts=True,
    # selected_stations=[7,11,14,16,18,19,23,26,30,31,32,33])

# stn31 = m.list_of_stations[30]
# stn32 = m.list_of_stations[31]
# stn33 = m.list_of_stations[32]
# stn31.plot_parts()
# stn32.plot_parts()
# stn33.plot_parts()
