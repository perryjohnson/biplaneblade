"""A script to merge all the structural parts into one polygon at each station.

Author: Perry Roth-Johnson
Last updated: October 10, 2013

"""


import lib.blade as bl
reload(bl)

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
    # station.write_all_part_polygons()

# for station in m.list_of_stations:
#     a = station.cross_section_area()
#     print "Station #{0}, area = {1} m^2".format(station.station_num, a)
#     d = station.structure.which_parts_exist()
#     for keys, values in sorted(d.items()):
#         print "  {0} : {1}".format(keys, values)

# m.plot_selected_cross_sections(plot_edges=False, plot_parts=True)
# stn33 = m.list_of_stations[32]
# stn33.merge_all_parts(plot_flag=True, merge_internal_surface=False)
# stn33.merge_all_parts(plot_flag=True, merge_internal_surface=True)
