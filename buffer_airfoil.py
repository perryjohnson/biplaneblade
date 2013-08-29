"""Example code that generates offset curves for an airfoil.

Main source:
https://github.com/Toblerity/Shapely/blob/master/docs/code/buffer.py

Additional sources:
https://pypi.python.org/pypi/descartes/
http://toblerity.org/shapely/manual.html#object.buffer
http://stackoverflow.com/questions/1109536/an-algorithm-for-inflating-deflating-offsetting-buffering-polygons

Author: Toblerity, Perry Roth-Johnson
Last modified: August 29, 2013

"""


import numpy as np
from matplotlib import pyplot
from shapely.geometry import LineString, Polygon
from shapely.geometry.polygon import LinearRing
from descartes import PolygonPatch
# the descartes module translates shapely objects into matplotlib objects
import lib.airfoil as af


BLUE = '#6699cc'
GRAY = '#999999'

def nparray_to_polygon(nparray):
    l = []
    for point in nparray:
        x = float(point['x'])
        y = float(point['y'])
        l.append((x,y))
    p = Polygon(l)
    return p

# import airfoil coordinates
a = af.MonoplaneAirfoil(name="NACA 64-618", filename="NACA_64-618.txt", chord=1.0, pitch_axis=0.375, twist=0.0)
a.path = "sandia_blade/airfoils/NACA_64-618.txt"
a.read_coords()
fig, ax = pyplot.subplots(figsize=(10,8))
ax.set_xlim([-0.2, 1.2])
ax.set_ylim([-0.4, 0.4])
ax.set_aspect('equal')

# convert airfoil coordinates to a polygon
p = nparray_to_polygon(a.coords)
airfoil_patch = PolygonPatch(p, fc=GRAY, ec=GRAY, alpha=0.5, zorder=1)
ax.add_patch(airfoil_patch)

# erode the airfoil polygon by the desired laminate thickness
eroded = p.buffer(-0.01)
polygon = eroded.__geo_interface__
patch2b = PolygonPatch(polygon, fc=BLUE, ec=BLUE, alpha=0.5, zorder=2)
ax.add_patch(patch2b)

# use clipping to cut out a structural part (e.g. spar cap) from `eroded`
# ...

pyplot.show()