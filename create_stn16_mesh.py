import os
import matplotlib.pyplot as plt
import numpy as np
import lib.blade as bl
reload(bl)

# load the sandia blade
m = bl.MonoplaneBlade('Sandia blade SNL100-00', 'sandia_blade',
    rotate_airfoil_coords=False)

# pre-process the airfoil coordinates
stn16 = m.list_of_stations[15]
stn16.airfoil.create_polygon()
stn16.structure.create_all_layers()
stn16.structure.save_all_layer_edges()
stn16.structure.create_all_alternate_layers()
stn16.structure.save_all_alternate_layer_edges()
# stn16.plot_layer_edges()
stn16.structure.write_truegrid_inputfile(interrupt_flag=True)
stn16.structure.write_all_part_polygons()

stn16.plot_parts(alternate_layers=True)
