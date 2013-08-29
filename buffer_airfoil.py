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

# set up matplotlib figure
fig, ax = pyplot.subplots(figsize=(10,8))
ax.set_xlim([-0.2, 1.2])
ax.set_ylim([-0.4, 0.4])
ax.set_aspect('equal')

# convert airfoil coordinates to a polygon
p1 = nparray_to_polygon(a.coords)
p1_patch = PolygonPatch(p1, fc=GRAY, ec=GRAY, alpha=0.5, zorder=1)
# ax.add_patch(p1_patch)

# erode the airfoil polygon by the desired laminate thickness
# eroded = p1.buffer(-0.01)
# p2 = eroded.__geo_interface__  # you can use __geo_interface__ to get the points back out, I think, if you call object.__geo_interface__ from the IPython command line, you get a dictionary of coordinates
p2 = p1.buffer(-0.01)
p2_patch = PolygonPatch(p2, fc=BLUE, ec=BLUE, alpha=0.5, zorder=2)
# ax.add_patch(p2_patch)

# cut a hole out of the middle of the airfoil polygon
# ref: http://toblerity.org/shapely/manual.html#object.difference
p3 = p1.difference(p2)
p3_patch = PolygonPatch(p3, fc=(1,0,0), ec=(1,0,0), alpha=0.7, zorder=3)
ax.add_patch(p3_patch)

# use clipping to cut out a structural part (e.g. spar cap) from `eroded`
left_box = Polygon([(-0.1,-0.2),(0.2,-0.2),(0.2,0.2),(-0.1,0.2)])
right_box = Polygon([(0.6,-0.2),(1.1,-0.2),(1.1,0.2),(0.6,0.2)])
left_box_patch = PolygonPatch(left_box, fc=(0,1,0), ec=(0,1,0), alpha=0.3, zorder=4)
right_box_patch = PolygonPatch(right_box, fc=(0,0,1), ec=(0,0,1), alpha=0.3, zorder=4)
ax.add_patch(left_box_patch)
ax.add_patch(right_box_patch)
p4 = p3.difference(left_box).difference(right_box)

# p4 is a `MultiPolygon`: it contains two polygons (upper and lower spar caps)
# extract each polygon and convert them to separate patches
# (otherwise, PolygonPatch will throw an error)
sc1 = p4.geoms[0]
sc2 = p4.geoms[1]
sc1_patch = PolygonPatch(sc1, fc=(1,1,0), ec=(1,1,0), alpha=0.5, zorder=5)
sc2_patch = PolygonPatch(sc2, fc=(0,1,1), ec=(0,1,1), alpha=0.5, zorder=5)
# now we can plot each polygon in matplotlib
ax.add_patch(sc1_patch)  # lower spar cap
ax.add_patch(sc2_patch)  # upper spar cap

# print the coordinates of the spar caps
print "lower spar cap coordinates"
print "--------------------------"
print sc1.__geo_interface__
print ""
print "upper spar cap coordinates"
print "--------------------------"
print sc2.__geo_interface__

pyplot.show()