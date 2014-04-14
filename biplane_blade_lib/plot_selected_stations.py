"""Plot structural parts for selected stations in the Sandia blade.

Usage
-----
start an IPython (qt)console with the pylab flag:
$ ipython qtconsole --pylab
or
$ ipython --pylab
Then, from the prompt, run this script:
|> %run plot_selected_stations
Once you are finished looking at the meshes, you can clean up extra files:
|> %run clean
(See the 'clean.py' script in this directory for details.)

Author: Perry Roth-Johnson
Last updated: March 10, 2014

"""

import matplotlib.pyplot as plt
import numpy as np
import lib.blade as bl
reload(bl)

# load the biplane blade
b1 = bl.BiplaneBlade(
    'biplane blade, flapwise symmetric, no stagger, rj/R=0.452, g/c=1.25',
    'biplane_blade',
    rotate_airfoil_coords=False)

# pre-process the airfoil coordinates
for station in b1.list_of_stations:
    station.airfoil.create_polygon()
    station.structure.create_all_layers()
    # station.structure.save_all_layer_edges()
    # station.structure.create_all_alternate_layers()
    # station.structure.save_all_alternate_layer_edges()
    # station.structure.write_truegrid_inputfile(interrupt_flag=True)
    # station.structure.write_all_part_polygons()

b1.plot_selected_cross_sections()

# for station in m.list_of_stations:
#     station.plot_parts(alternate_layers=False)
