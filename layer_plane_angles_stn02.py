"""Determine the layer plane angle of all the elements in a grid.

Author: Perry Roth-Johnson
Last modified: March 17, 2014

References:
http://stackoverflow.com/questions/3365171/calculating-the-angle-between-two-lines-without-having-to-calculate-the-slope/3366569#3366569
http://stackoverflow.com/questions/19295725/angle-less-than-180-between-two-segments-lines

"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import lib.grid as gr
reload(gr)
import lib.abaqus_utils2 as au
reload(au)
import lib.vabs_utils as vu
reload(vu)
from shapely.geometry import Polygon, LineString
from descartes import PolygonPatch


# -----------------------------------------------
# update these parameters!
station_num = 2
# -----------------------------------------------

stn_str = 'stn{0:02d}'.format(station_num)
plt.close('all')
# create a figure
plt.figure(num='Station #{0:02d}'.format(station_num))
ax = plt.gcf().gca()

list_of_unflipped_elementsets = [
    'sclower',  # correct?
    'rbtrile',  # correct: outer_edge_node_nums=[1,4], inner_edge_node_nums=[2,3]
    'rbtriscl',  # correct: outer_edge_node_nums=[1,4], inner_edge_node_nums=[2,3]
    'esgelle',  # correct: outer_edge_node_nums=[1,4], inner_edge_node_nums=[2,3]
    'esgelscl',  # correct?
    'estrile',  # correct: outer_edge_node_nums=[1,4], inner_edge_node_nums=[2,3]
    'estriscl',  # correct?
    'isresle',  # correct: outer_edge_node_nums=[1,4], inner_edge_node_nums=[2,3]
    'isresscl',  # correct?
    'istrile',  # correct: outer_edge_node_nums=[1,4], inner_edge_node_nums=[2,3]
    'istriscl'  # correct?
    ]
list_of_flipped_elementsets = [
    'scupper',  # incorrect
    'teuniax',  # incorrect
    'rbtrite',  # incorrect
    'rbtrilr',  # incorrect
    'rbtriur',  # incorrect?
    'rbtriscu',  # incorrect
    'esgelte',  # incorrect?
    'esgellr',  # incorrect?
    'esgelur',  # incorrect?
    'esgelscu',  # incorrect
    'estrite',  # incorrect?
    'estrilr',  # incorrect?
    'estriur',  # incorrect
    'estriscu',  # incorrect?
    'isreste',  # incorrect?
    'isreslr',  # incorrect?
    'isresur',  # incorrect?
    'isresscu',  # incorrect?
    'istrite',  # incorrect?
    'istrilr',  # incorrect?
    'istriur',  # incorrect?
    'istriscu'  # incorrect
    ]

# import the initial grid object
fmt_grid = 'sandia_blade/' + stn_str + '/mesh_' + stn_str + '.abq'
g = au.AbaqusGrid(fmt_grid, debug_flag=True)
# update the grid object with all the layer plane angles
for elem in g.list_of_elements:    
    if elem.element_set in list_of_unflipped_elementsets:
        elem.calculate_layer_plane_angle(outer_edge_node_nums=[1,4],
            inner_edge_node_nums=[2,3])
    elif elem.element_set in list_of_flipped_elementsets:
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
# -----------------------------------------------------------------------------
# read layers.csv to determine the number of layers
layer_file = pd.read_csv('sandia_blade/layers.csv', index_col=0)
number_of_layers = len(layer_file)
# write the updated grid object to a VABS input file
fmt_vabs = 'sandia_blade/' + stn_str + '/mesh_' + stn_str + '.vabs'
f = vu.VabsInputFile(
    vabs_filename=fmt_vabs,
    grid=g,
    material_filename='sandia_blade/materials.csv',
    layer_filename='sandia_blade/layers.csv',
    debug_flag=True,
    flags={
        'format'           : 1,
        'Timoshenko'       : 1,
        'recover'          : 0,
        'thermal'          : 0,
        'curve'            : 0,
        'oblique'          : 0,
        'trapeze'          : 0,
        'Vlasov'           : 0
    })
