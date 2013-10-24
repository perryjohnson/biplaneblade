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
    # station.plot_polygon_edges()

stn7 = m.list_of_stations[6]
stn7.plot_polygon_edges()

# m.plot_selected_cross_sections()