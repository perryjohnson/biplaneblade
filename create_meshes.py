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
    station.structure.save_all_layer_edges()
    # station.plot_layer_edges()
    # station.structure.write_all_layer_edges()
    station.structure.write_truegrid_inputfile(interrupt_flag=True)

# next steps for tomorrow:
# (1) DONE: plot meshes for all stations!   : )
# (2) fix mesh errors in TE reinforcement
#     ... may need to split LBTR curves further near sharp corners...
# (3) these stations crashed when running 'mesh.tg'
#     2, 3, 4, 6, 7, 10, 11, 12, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
#     26, 27, 28, 29, 30, 31, 32, 33, 34

m.plot_selected_cross_sections()
# m.plot_selected_cross_sections(selected_stations=range(22,34))
# stn16 = m.list_of_stations[15]
# stn16.plot_parts()