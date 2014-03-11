"""Write initial TrueGrid files for one Sandia blade station.

Usage
-----
start an IPython (qt)console with the pylab flag:
$ ipython qtconsole --pylab
or
$ ipython --pylab
Then, from the prompt, run this script:
|> %run create_one_mesh
Once you are finished looking at the meshes, you can clean up extra files:
|> %run clean
(See the 'clean.py' script in this directory for details.)

Author: Perry Roth-Johnson
Last updated: March 10, 2014

"""

import os
import matplotlib.pyplot as plt
import numpy as np
import lib.blade as bl
reload(bl)
from shapely.geometry import Polygon
from descartes import PolygonPatch
import lib.layer as l

# load the Sandia blade
m = bl.MonoplaneBlade('Sandia blade SNL100-00', 'sandia_blade',
    rotate_airfoil_coords=False)

# pick a station number
station_num = 1

# pre-process the station dimensions
station = m.list_of_stations[station_num-1]
station.airfoil.create_polygon()
station.structure.create_all_layers()
station.structure.save_all_layer_edges()
# station.structure.write_truegrid_inputfile(interrupt_flag=True)
station.structure.write_all_part_polygons()

# plot the parts
station.plot_parts(alternate_layers=False)

def cut_polygon(original, bounding):
    """Cut the original layer polygon with the bounding polygon."""
    return original.intersection(bounding)

def plot_polygon(p, face_color, edge_color='None'):
    """Plot a polygon on the current axes."""
    # get the current axes, so we can add polygons to the plot
    ax = plt.gca()
    patch = PolygonPatch(p, fc=face_color, ec=edge_color, alpha=0.8)
    ax.add_patch(patch)

# upper right -----------------------------------------------------------------
label = 'upper right'

# create the bounding polygon
points = [
    (0.0, 0.0),
    (3.0, 0.0),
    (3.0, 3.0),
    (0.0, 3.0)
    ]
bounding_polygon = Polygon(points)
plot_polygon(bounding_polygon, 'None', '#000000')

# cut the new layer polygons
st = station.structure

l = st.root_buildup.layer['triax']
p_new = cut_polygon(l.polygon, bounding_polygon)
plot_polygon(p_new, l.face_color)
l.parent_part.add_new_layer('triax, '+label, p_new)

l = st.external_surface.layer['triax']
p_new = cut_polygon(l.polygon, bounding_polygon)
plot_polygon(p_new, l.face_color)
l.parent_part.add_new_layer('triax, '+label, p_new, 'triax')

l = st.external_surface.layer['gelcoat']
p_new = cut_polygon(l.polygon, bounding_polygon)
plot_polygon(p_new, l.face_color)
l.parent_part.add_new_layer('gelcoat, '+label, p_new, 'gelcoat')

l = st.internal_surface_1.layer['triax']
p_new = cut_polygon(l.polygon, bounding_polygon)
plot_polygon(p_new, l.face_color)
l.parent_part.add_new_layer('triax, '+label, p_new, 'triax')

l = st.internal_surface_1.layer['resin']
p_new = cut_polygon(l.polygon, bounding_polygon)
plot_polygon(p_new, l.face_color)
l.parent_part.add_new_layer('resin, '+label, p_new, 'resin')



# show the plot
plt.show()
