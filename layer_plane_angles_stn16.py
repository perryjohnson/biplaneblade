"""Determine the layer plane angle of all the elements in a grid.

Author: Perry Roth-Johnson
Last modified: March 25, 2014

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
station_num = 16
# -----------------------------------------------

stn_str = 'stn{0:02d}'.format(station_num)
plt.close('all')
# create a figure
plt.figure(num='Station #{0:02d}'.format(station_num))
ax = plt.gcf().gca()

# element sets on the leading edge
# outer_edge_node_nums=[1,4], inner_edge_node_nums=[2,3]
list_of_LE_elementsets = [
    'lepanel', # yes
    'sw1biaxl', # yes
    'sw1foam', # yes
    'sw1biaxr', # yes
    'esleg', # yes
    'eslet', # yes
    'is1lt',
    'is1lr',
    'is1rt',
    'is1rr',
    'is2lt',
    'is2lr'
    ]

# element sets on the trailing edge
# outer_edge_node_nums=[3,2], inner_edge_node_nums=[4,1]
list_of_TE_elementsets = [
    'sw2biaxl',
    'sw2foam',
    'sw2biaxr',
    'sw3biaxl',
    'sw3foam',
    'sw3biaxr',
    'is2rt',
    'is2rr',
    'is3rt',
    'is3rr',
    'estet',
    'esteg',
    'is3lt',
    'is3lr',
    'is4lt',
    'is4lr',
    ]

# element sets on the lower surface
# outer_edge_node_nums=[2,1], inner_edge_node_nums=[3,4]
list_of_lower_elementsets = [
    'sclower',
    'ap1lower',
    'ap2lower',
    'eslowerg',
    'eslowert',
    'is2bt',
    'is2br',
    'is3bt',
    'is3br',
    'is4bt',
    'is4br',
    'teunixll',
    'teunixlm',
    'teunixlr',
    'tefoamll',
    'tefoamlm',
]

# element sets on the upper surface
# outer_edge_node_nums=[4,3], inner_edge_node_nums=[1,2]
list_of_upper_elementsets = [
    'esuppert',
    'esupperg',
    'scupper',
    'ap1upper',
    'ap2upper',
    'is2tt',
    'is2tr',
    'is3tt',
    'is3tr',
    'is4tt',
    'teunixul',
    'teunixum',
    'teunixur',
    'tefoamul',
    'tefoamum',
    'is4tr'
]

# import the initial grid object
fmt_grid = 'sandia_blade/' + stn_str + '/mesh_' + stn_str + '.abq'
g = au.AbaqusGrid(fmt_grid, debug_flag=True)
# update the grid object with all the layer plane angles
for elem in g.list_of_elements:
    if elem.element_set in list_of_LE_elementsets:
        elem.calculate_layer_plane_angle(outer_edge_node_nums=[1,4],
            inner_edge_node_nums=[2,3])
    elif elem.element_set in list_of_TE_elementsets:
        elem.calculate_layer_plane_angle(outer_edge_node_nums=[3,2],
            inner_edge_node_nums=[4,1])
    elif elem.element_set in list_of_lower_elementsets:
        elem.calculate_layer_plane_angle(outer_edge_node_nums=[2,1], 
            inner_edge_node_nums=[3,4])
    elif elem.element_set in list_of_upper_elementsets:
        elem.calculate_layer_plane_angle(outer_edge_node_nums=[4,3], 
            inner_edge_node_nums=[1,2])
    else:
        raise Warning("Element #{0} has no element set!".format(elem.elem_num))
# plot a small selection of elements to check the results
for elem in g.list_of_elements[::25]:
# for elem in g.list_of_elements[3650:3710]:
    elem.plot(label_nodes=False)
    print elem.elem_num, elem.element_set, elem.theta1
# show the plot
plt.xlim([-3,5.5])
plt.ylim([-3,3])
ax.set_aspect('equal')
print ' ------------------------'
print '  LEGEND'
print '    magenta : inner edge'
print '    blue    : outer edge'
print ' ------------------------'

# plt.figure(num='Station #{0:02d}, theta1 vs. elem_num'.format(
#     station_num))
# enum=np.arange(g.number_of_elements)+1
# theta=np.zeros(g.number_of_elements)
# elemset=[]
# for i,elem in enumerate(g.list_of_elements):
#     theta[i] = elem.theta1
#     elemset.append(elem.element_set)
# plt.plot(enum,theta)
# plt.xlabel('element number [#]')
# plt.ylabel('theta1 [deg]')
# plt.grid('on')

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
