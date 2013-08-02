"""
Plot all the airfoils in the Sandia blade with mayavi.

Usage:
Run an IPython prompt with the pylab option (ipython qtconsole --pylab).

Author: Perry Roth-Johnson
Date: August 2, 2013

"""


import numpy as np
import scripts.blade as bl
from mayavi import mlab


# make a new figure
mlab.figure(1, size=(800, 800))
# mlab.clf()

# import the blade and its airfoil coordinates
b = bl.Blade('Sandia blade SNL100-00', 'sandia_blade')
b.copy_all_airfoil_coords()
for station in b.list_of_stations:
    station.read_airfoil_coords()

    # assemble the airfoil coordinates for mlab
    x = station.airfoil.coords['x']
    y = station.airfoil.coords['y']
    l = len(x)
    z = np.ones((l,))*station.coords.x1  # spanwise coordinate
    s = np.ones((l,))  # arbitrary scalar value, which has no meaning here (normally this is used for plotting how a scalar field that varies in 3D space)

    # make connectivity information between points to draw lines
    # taken from plotting_many_lines.py
    connections = []
    connections.append(np.vstack([np.arange(0,l-1.5),np.arange(1,l-0.5)]).T)

    # plot the airfoil on the screen
    src = mlab.pipeline.scalar_scatter(x,y,z,s)
    src.mlab_source.dataset.lines=connections
    lines = mlab.pipeline.stripper(src)
    mlab.pipeline.surface(lines, line_width=1)

# show the final plot
mlab.show()
