"""Write initial TrueGrid files for Sandia blade stations.

Usage
-----
start an IPython (qt)console with the pylab flag:
$ ipython qtconsole --pylab
or
$ ipython --pylab
Then, from the prompt, run this script:
|> %run create_meshes
Once you are finished looking at the meshes, you can clean up extra files:
|> %run clean
(See the 'clean.py' script in this directory for details.)

Author: Perry Roth-Johnson
Last updated: March 6, 2014

"""

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
    station.structure.create_all_alternate_layers()
    station.structure.save_all_alternate_layer_edges()
    # station.plot_layer_edges()
    station.structure.write_truegrid_inputfile(interrupt_flag=True)
    station.structure.write_all_part_polygons()

m.plot_selected_cross_sections(alternate_layers=True)
# m.list_of_stations[15].plot_parts(alternate_layers=True)
