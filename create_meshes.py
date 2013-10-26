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
    station.structure.write_truegrid_inputfile()

# next steps for tomorrow:
# (1) DONE: plot meshes for all stations!   : )
# (2) fix mesh errors in TE reinforcement
#     ... may need to split LBTR curves further near sharp corners...

# m.plot_selected_cross_sections()