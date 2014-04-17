"""Plot structural parts for selected stations in the Sandia blade.

Usage
-----
start an IPython (qt)console with the pylab flag:
$ ipython qtconsole --pylab
or
$ ipython --pylab
Then, from the prompt, run this script:
|> %run sandia_blade_lib/plot_selected_stations
Once you are finished looking at the meshes, you can clean up extra files:
|> %run clean
(See the 'clean.py' script in this directory for details.)

Author: Perry Roth-Johnson
Last updated: April 14, 2014

"""

import matplotlib.pyplot as plt
import numpy as np
import lib.blade as bl
reload(bl)

# load the sandia blade
m = bl.MonoplaneBlade('Sandia blade SNL100-00', 'sandia_blade')

# pre-process the airfoil coordinates
for station in m.list_of_stations:
    station.airfoil.create_polygon()
    station.structure.create_all_layers()
    station.structure.save_all_layer_edges()
    # station.structure.create_all_alternate_layers()
    # station.structure.save_all_alternate_layer_edges()
    # station.structure.write_truegrid_inputfile(interrupt_flag=True)
    station.structure.write_all_part_polygons()

m.plot_selected_cross_sections()

# for station in m.list_of_stations:
#     station.plot_parts(alternate_layers=False)
