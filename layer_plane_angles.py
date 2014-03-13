"""Determine the layer plane angle of all the elements in a grid.

Author: Perry Roth-Johnson
Last modified: March 13, 2014

References:
http://stackoverflow.com/questions/3365171/calculating-the-angle-between-two-lines-without-having-to-calculate-the-slope/3366569#3366569
http://stackoverflow.com/questions/19295725/angle-less-than-180-between-two-segments-lines

"""

import numpy as np
import matplotlib.pyplot as plt
import lib.grid as gr
reload(gr)
import lib.abaqus_utils2 as au
reload(au)
import lib.vabs_utils as vu
from shapely.geometry import Polygon, LineString
from descartes import PolygonPatch


plt.close('all')
# create a figure
ax = plt.gcf().gca()

left_elemsets = [
    'rbtrill',
    'rbtriul',
    'estrill',
    'estriul',
    'esgelll',
    'esgelul',
    'istrill',
    'istriul',
    'isresll',
    'isresul'
    ]
right_elemsets = [
    'rbtrilr',
    'rbtriur',
    'estrilr',
    'estriur',
    'esgellr',
    'esgelur',
    'istrilr',
    'istriur',
    'isreslr',
    'isresur'
    ]

# import the initial grid object
g = au.AbaqusGrid('sandia_blade/mesh_stn01.abq', debug_flag=True)
# update the grid object with all the layer plane angles
for elem in g.list_of_elements:    
    if elem.element_set in left_elemsets:
        elem.calculate_layer_plane_angle(outer_edge_node_nums=[1,4],
            inner_edge_node_nums=[2,3])
    elif elem.element_set in right_elemsets:
        elem.calculate_layer_plane_angle(outer_edge_node_nums=[3,2],
            inner_edge_node_nums=[4,1])
    else:
        raise Warning("Element #{0} has no element set!".format(elem.elem_num))
# plot a small selection of elements to check the results
for elem in g.list_of_elements[::25]:
    elem.plot()
    print elem.elem_num, elem.element_set, elem.theta1
# show the plot
plt.xlim([-3,3])
plt.ylim([-3,3])
ax.set_aspect('equal')
plt.show()
# write the updated grid object to a VABS input file
f = vu.VabsInputFile(
    vabs_filename='sandia_blade/mesh_stn01.vabs',
    grid=g,
    debug_flag=True)
