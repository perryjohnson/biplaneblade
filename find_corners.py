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
        station.structure.LE_panel.layer['foam'].get_and_save_edges()
    if station.structure.spar_cap.exists():
        station.structure.spar_cap.layer['lower'].get_and_save_edges()
        station.structure.spar_cap.layer['upper'].get_and_save_edges()
    if station.structure.aft_panel_1.exists():
        station.structure.aft_panel_1.layer['lower'].get_and_save_edges()
        station.structure.aft_panel_1.layer['upper'].get_and_save_edges()
    if station.structure.aft_panel_2.exists():
        station.structure.aft_panel_2.layer['lower'].get_and_save_edges()
        station.structure.aft_panel_2.layer['upper'].get_and_save_edges()
    if station.structure.shear_web_1.exists():
        station.structure.shear_web_1.layer['biax, left'].get_and_save_edges()
        station.structure.shear_web_1.layer['foam'].get_and_save_edges()
        station.structure.shear_web_1.layer['biax, right'].get_and_save_edges()
    if station.structure.shear_web_2.exists():
        station.structure.shear_web_2.layer['biax, left'].get_and_save_edges()
        station.structure.shear_web_2.layer['foam'].get_and_save_edges()
        station.structure.shear_web_2.layer['biax, right'].get_and_save_edges()
    if station.structure.shear_web_3.exists():
        station.structure.shear_web_3.layer['biax, left'].get_and_save_edges()
        station.structure.shear_web_3.layer['foam'].get_and_save_edges()
        station.structure.shear_web_3.layer['biax, right'].get_and_save_edges()
    if station.structure.TE_reinforcement.exists():
        station.structure.TE_reinforcement.layer['uniax'].get_and_save_edges()
        try:
            station.structure.TE_reinforcement.layer['foam'].get_and_save_edges()
        except KeyError:  # foam layer doesn't exist
            pass
    if station.structure.root_buildup.exists():
        station.structure.root_buildup.layer['triax, lower left'].get_and_save_edges()
        station.structure.root_buildup.layer['triax, lower right'].get_and_save_edges()
        station.structure.root_buildup.layer['triax, upper right'].get_and_save_edges()
        station.structure.root_buildup.layer['triax, upper left'].get_and_save_edges()
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

stn6 = m.list_of_stations[5]
stn6.plot_polygon_edges()
stn6.plot_parts()

stn16 = m.list_of_stations[15]
stn16.plot_polygon_edges()
stn16.plot_parts()

f = open(os.path.join(stn16.station_path, 'curves.tg'), 'w')
# write the curves of the LE panel into a file
# left curve
curve_num = 1
f.write('curd {0} lp3\n'.format(curve_num))
for coord_pair in stn16.structure.LE_panel.layer['foam'].left:
    f.write('{0: .8f}  {1: .8f}  0.0\n'.format(coord_pair[0], coord_pair[1]))
f.write(';;\n\n')
# bottom curve
curve_num += 1
f.write('curd {0} lp3\n'.format(curve_num))
for coord_pair in stn16.structure.LE_panel.layer['foam'].bottom:
    f.write('{0: .8f}  {1: .8f}  0.0\n'.format(coord_pair[0], coord_pair[1]))
f.write(';;\n\n')
# right curve
curve_num += 1
f.write('curd {0} lp3\n'.format(curve_num))
for coord_pair in stn16.structure.LE_panel.layer['foam'].right:
    f.write('{0: .8f}  {1: .8f}  0.0\n'.format(coord_pair[0], coord_pair[1]))
f.write(';;\n\n')
# top curve
curve_num += 1
f.write('curd {0} lp3\n'.format(curve_num))
for coord_pair in stn16.structure.LE_panel.layer['foam'].top:
    f.write('{0: .8f}  {1: .8f}  0.0\n'.format(coord_pair[0], coord_pair[1]))
f.write(';;\n\n')
# curve_num += 1

f.close()

# m.plot_selected_cross_sections()