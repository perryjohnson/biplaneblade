"""Example code that generates offset curves for the internal surface.

Main source:
https://github.com/Toblerity/Shapely/blob/master/docs/code/buffer.py

Additional sources:
https://pypi.python.org/pypi/descartes/
http://toblerity.org/shapely/manual.html#object.buffer
http://stackoverflow.com/questions/1109536/an-algorithm-for-inflating-deflating-offsetting-buffering-polygons

Author: Toblerity, Perry Roth-Johnson
Last modified: October 7, 2013

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
a = af.MonoplaneAirfoil(name="Cylinder", filename="Cylinder.txt", chord=1.0, pitch_axis=0.375, twist=0.0)
a.path = "sandia_blade/airfoils/Cylinder.txt"
a.read_coords()

# set up matplotlib figure
fig, ax = pyplot.subplots(figsize=(10,8))
ax.set_xlim([-0.2, 1.2])
ax.set_ylim([-0.6, 0.6])
ax.set_aspect('equal')

# convert airfoil coordinates to a polygon
p1 = nparray_to_polygon(a.coords)
# p1_patch = PolygonPatch(p1, fc=GRAY, ec=GRAY, alpha=0.5, zorder=1)
# ax.add_patch(p1_patch)

# erode the airfoil polygon by the desired laminate thickness
p2 = p1.buffer(-0.01)
# p2_patch = PolygonPatch(p2, fc=BLUE, ec=BLUE, alpha=0.5, zorder=2)
# ax.add_patch(p2_patch)

p3 = p2.buffer(-0.01)
# p3_patch = PolygonPatch(p3, fc=GRAY, ec=GRAY, alpha=0.5, zorder=3)
# ax.add_patch(p3_patch)

# cut out the external surface
ext_surface = p1.difference(p2)
ext_surface_patch = PolygonPatch(ext_surface, fc=GRAY, ec=GRAY, alpha=0.5, zorder=1)
ax.add_patch(ext_surface_patch)

# cut out the root buildup
root_buildup = p2.difference(p3)
root_buildup_patch = PolygonPatch(root_buildup, fc=BLUE, ec=BLUE, alpha=0.5, zorder=1)
ax.add_patch(root_buildup_patch)

# create spar caps
p4 = p3.buffer(-0.01)  # inner profile
bounding_box = Polygon([(0.2,-0.55),(0.8,-0.55),(0.8,0.55),(0.2,0.55)])
# bounding_box_patch = PolygonPatch(bounding_box, fc=(0,1,0), ec=(0,1,0), alpha=0.3, zorder=4)
# ax.add_patch(bounding_box_patch)
p5 = p3.difference(p4)  # cut out inner profile from outer profile of spar caps
spar_caps = p5.intersection(bounding_box)

# spar_caps is a `MultiPolygon`: it contains two polygons (upper and lower spar caps)
# extract each polygon and convert them to separate patches
# (otherwise, PolygonPatch will throw an error)
sc1 = spar_caps.geoms[0]
sc2 = spar_caps.geoms[1]
sc1_patch = PolygonPatch(sc1, fc=(1,0,0), ec=(1,0,0), alpha=0.8, zorder=5)
sc2_patch = PolygonPatch(sc2, fc=GRAY, ec=GRAY, alpha=0.8, zorder=5)
# now we can plot each polygon in matplotlib
ax.add_patch(sc1_patch)  # lower spar cap
ax.add_patch(sc2_patch)  # upper spar cap

# create TE reinforcement
bb2 = Polygon([(0.5,0.2),(0.5,-0.2),(1.1,-0.2),(1.1,0.2)])
TE_reinf = p5.intersection(bb2)
TE_reinf_patch = PolygonPatch(TE_reinf, fc=(1,0,1), ec=(1,0,1), alpha=0.8, zorder=5)
ax.add_patch(TE_reinf_patch)

# create internal surface (triax)
op_is_triax = p3.difference(sc1).difference(sc2).difference(TE_reinf)
ip_is_triax = op_is_triax.buffer(-0.01)
int_surface_triax = op_is_triax.difference(ip_is_triax)
int_surface_triax_patch = PolygonPatch(int_surface_triax, fc=(0.8,0.5,0.5), ec=(0.8,0.5,0.5), alpha=0.5, zorder=6)
ax.add_patch(int_surface_triax_patch)

# create internal surface (resin)
op_is_resin = ip_is_triax
ip_is_resin = op_is_resin.buffer(-0.01)
int_surface_resin = op_is_resin.difference(ip_is_resin)
int_surface_resin_patch = PolygonPatch(int_surface_resin, fc=(0.8,0.1,0.5), ec=(0.8,0.1,0.5), alpha=0.5, zorder=6)
ax.add_patch(int_surface_resin_patch)

# print the coordinates of the internal surface (triax and resin layers)
print "internal surface (triax) coordinates"
print "------------------------------------"
print int_surface_triax.__geo_interface__
print ""
print "internal surface (resin) coordinates"
print "------------------------------------"
print int_surface_resin.__geo_interface__

pyplot.show()