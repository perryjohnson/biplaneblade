import os
import matplotlib.pyplot as plt
import numpy as np
import lib.blade as bl
reload(bl)

# load the sandia blade
m = bl.MonoplaneBlade('Sandia blade SNL100-00', 'sandia_blade',
    rotate_airfoil_coords=False)

# pre-process the airfoil coordinates
for station in m.list_of_stations:
    station.airfoil.create_polygon()
    station.structure.create_all_layers()
    if station.structure.LE_panel.exists():
        station.structure.LE_panel.get_and_save_edges()
    if station.structure.spar_cap.exists():
        station.structure.spar_cap.get_and_save_edges('lower')
        station.structure.spar_cap.get_and_save_edges('upper')
    if station.structure.aft_panel_1.exists():
        station.structure.aft_panel_1.get_and_save_edges('lower')
        station.structure.aft_panel_1.get_and_save_edges('upper')
    if station.structure.aft_panel_2.exists():
        station.structure.aft_panel_2.get_and_save_edges('lower')
        station.structure.aft_panel_2.get_and_save_edges('upper')
    if station.structure.shear_web_1.exists():
        station.structure.shear_web_1.get_and_save_edges('left biax')
        station.structure.shear_web_1.get_and_save_edges('foam')
        station.structure.shear_web_1.get_and_save_edges('right biax')
    if station.structure.shear_web_2.exists():
        station.structure.shear_web_2.get_and_save_edges('left biax')
        station.structure.shear_web_2.get_and_save_edges('foam')
        station.structure.shear_web_2.get_and_save_edges('right biax')
    if station.structure.shear_web_3.exists():
        station.structure.shear_web_3.get_and_save_edges('left biax')
        station.structure.shear_web_3.get_and_save_edges('foam')
        station.structure.shear_web_3.get_and_save_edges('right biax')
    if station.structure.TE_reinforcement.exists():
        station.structure.TE_reinforcement.get_and_save_edges('uniax')
        try:
            station.structure.TE_reinforcement.get_and_save_edges('foam')
        except Warning:  # foam layer doesn't exist
            pass
    if station.structure.root_buildup.exists():
        station.structure.root_buildup.get_and_save_edges('lower left')
        station.structure.root_buildup.get_and_save_edges('lower right')
        station.structure.root_buildup.get_and_save_edges('upper right')
        station.structure.root_buildup.get_and_save_edges('upper left')
    # if station.structure.external_surface.exists():
    #     station.structure.external_surface.get_and_save_edges('gelcoat, lower left')
    #     station.structure.external_surface.get_and_save_edges('gelcoat, lower right')
    #     station.structure.external_surface.get_and_save_edges('gelcoat, upper right')
    #     station.structure.external_surface.get_and_save_edges('gelcoat, upper left')
    #     station.structure.external_surface.get_and_save_edges('triax, lower left')
    #     station.structure.external_surface.get_and_save_edges('triax, lower right')
    #     station.structure.external_surface.get_and_save_edges('triax, upper right')
    #     station.structure.external_surface.get_and_save_edges('triax, upper left')
    # station.plot_polygon_edges()

stn16 = m.list_of_stations[15]
# stn16.plot_polygon_edges()
# stn16.plot_parts()

f = open(os.path.join(stn16.station_path, 'curves.tg'), 'w')
# write the curves of the LE panel into a file
# left curve
curve_num = 1
f.write('curd {0} lp3\n'.format(curve_num))
for coord_pair in stn16.structure.LE_panel.layer[0].left:
    f.write('{0: .8f}  {1: .8f}  0.0\n'.format(coord_pair[0], coord_pair[1]))
f.write(';;\n\n')
# bottom curve
curve_num += 1
f.write('curd {0} lp3\n'.format(curve_num))
for coord_pair in stn16.structure.LE_panel.layer[0].bottom:
    f.write('{0: .8f}  {1: .8f}  0.0\n'.format(coord_pair[0], coord_pair[1]))
f.write(';;\n\n')
# right curve
curve_num += 1
f.write('curd {0} lp3\n'.format(curve_num))
for coord_pair in stn16.structure.LE_panel.layer[0].right:
    f.write('{0: .8f}  {1: .8f}  0.0\n'.format(coord_pair[0], coord_pair[1]))
f.write(';;\n\n')
# top curve
curve_num += 1
f.write('curd {0} lp3\n'.format(curve_num))
for coord_pair in stn16.structure.LE_panel.layer[0].top:
    f.write('{0: .8f}  {1: .8f}  0.0\n'.format(coord_pair[0], coord_pair[1]))
f.write(';;\n\n')
# curve_num += 1

f.close()

# m.plot_selected_cross_sections()