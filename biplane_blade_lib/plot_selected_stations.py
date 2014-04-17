"""Plot structural parts for selected stations in the biplane blade.

Usage
-----
start an IPython (qt)console with the pylab flag:
$ ipython qtconsole --pylab
or
$ ipython --pylab
Then, from the prompt, run this script:
|> %run biplane_blade_lib/plot_selected_stations

Author: Perry Roth-Johnson
Last updated: April 17, 2014

"""

import matplotlib.pyplot as plt
import numpy as np
import lib.blade as bl
reload(bl)

# load the biplane blade
b1 = bl.BiplaneBlade(
    'biplane blade, flapwise symmetric, no stagger, rj/R=0.452, g/c=1.25',
    'biplane_blade')

# pre-process the airfoil coordinates
for station in b1.list_of_stations:
    station.airfoil.create_polygon()
    station.structure.create_all_layers()

b1.plot_selected_cross_sections()
