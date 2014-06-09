"""Determine the layer plane angle of all the elements in a grid.

Author: Perry Roth-Johnson
Last modified: April 9, 2014

Usage:
1. Look through the mesh_stn11.abq file and find all the element set names.
   (Find all the lines that start with "*ELSET".)
2. Enter each of the element set names in one of the four lists below:
     (1) list_of_LE_elementsets
     (2) list_of_TE_elementsets
     (3) list_of_lower_elementsets
     (4) list_of_upper_elementsets
3. Run this script. Visually inspect the plot to make sure each of the element
   sets are in the correct list. (The blue edge should be facing outward,
   relative to the airfoil centerpoint. The magenta edge should be facing 
   inward.) If you find an element that is oriented incorrectly, note the
   element number, and look up it's element set name from the printout in the
   IPython terminal. Then, move that element set name to a different list in
   this script.
4. Repeat step 3 until your visual inspection suggests that all the edges (and
   layer plane angles) are being assigned correctly.

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
import lib.blade as bl
from shapely.geometry import Polygon, LineString
from descartes import PolygonPatch


# -----------------------------------------------
# update these parameters!
station_num = 34
skip_num = 25   # plot every 'skip_num' elements (larger values plot faster)
ES_triax_upper_tri_elem_num = [251,237,223,510,478,446]  # num of tri elems in upper external surface (triax)
ES_triax_lower_tri_elem_num = [209,195,181,387,355,323]  # num of tri elems in lower external surface (triax)
IS_triax_upper_tri_elem_num = [764,754]  # num of tri elems in upper internal surface (triax)
IS_triax_lower_tri_elem_num = [755,745]  # num of tri elems in lower internal surface (triax)
IS_resin_upper_tri_elem_num = [698,688]  # num of tri elems in upper internal surface (resin) at LE
IS_resin_lower_tri_elem_num = [693,683]  # num of tri elems in lower internal surface (resin) at LE
# -----------------------------------------------

stn_str = 'stn{0:02d}'.format(station_num)
plt.close('all')

# load the Sandia blade
m = bl.MonoplaneBlade('Sandia blade SNL100-00', 'sandia_blade')

# pre-process the station dimensions
station = m.list_of_stations[station_num-1]
station.airfoil.create_polygon()
station.structure.create_all_layers()
station.structure.save_all_layer_edges()
station.structure.write_all_part_polygons()

# plot the parts
station.plot_parts(alternate_layers=False)

# create a figure
# plt.figure(num='Station #{0:02d}'.format(station_num))
ax = plt.gcf().gca()

# element sets on the leading edge
# outer_edge_node_nums=[1,4], inner_edge_node_nums=[2,3]
list_of_LE_elementsets = [
    ]

# element sets on the trailing edge
# outer_edge_node_nums=[3,2], inner_edge_node_nums=[4,1]
list_of_TE_elementsets = [
    ]

# element sets on the lower surface
# outer_edge_node_nums=[2,1], inner_edge_node_nums=[3,4]
list_of_lower_elementsets = [
    'esgellel',
    'esgelal1',
    'esgelal2',
    'esgelfl1',
    'esglower',
    'estrilel',
    'estrial1',
    'estrial2',
    'estrifl1',
    'estlower',
    'isresal1',
    'isresfl1',
    'isrlower',
    'istlower',
    'istrial1',
    'istrifl1'
]

# element sets on the upper surface
# outer_edge_node_nums=[4,3], inner_edge_node_nums=[1,2]
list_of_upper_elementsets = [
    'esgelleu',
    'esgelau1',
    'esgelau2',
    'esgelfu1',
    'esgupper',
    'estrileu',
    'estriau1',
    'estriau2',
    'estrifu1',
    'estupper',
    'isresau1',
    'isresfu1',
    'isrupper',
    'istupper',
    'istriau1',
    'istrifu1'
]

# element sets of triangular elements on the lower surface
# outer_edge_node_nums=[2,1]
list_of_tri_lower_elementsets = [
    'estril_tri',
    'istril_tri',
    'isresl_tri'
]

# element sets of triangular elements on the upper surface
# outer_edge_node_nums=[3,2]
list_of_tri_upper_elementsets = [
    'estriu_tri',
    'istriu_tri',
    'isresu_tri'
]

# import the initial grid object
fmt_grid = 'sandia_blade/' + stn_str + '/mesh_' + stn_str + '.abq'
g = au.AbaqusGrid(fmt_grid, debug_flag=True, soft_warning=False,
    auto_parse=True)

# manually assign triangular elements into new element sets
for num in ES_triax_upper_tri_elem_num:
    g.list_of_elements[num-1].element_set = 'estriu_tri'
for num in ES_triax_lower_tri_elem_num:
    g.list_of_elements[num-1].element_set = 'estril_tri'
for num in IS_triax_upper_tri_elem_num:
    g.list_of_elements[num-1].element_set = 'istriu_tri'
for num in IS_triax_lower_tri_elem_num:
    g.list_of_elements[num-1].element_set = 'istril_tri'
for num in IS_resin_upper_tri_elem_num:
    g.list_of_elements[num-1].element_set = 'isresu_tri'
for num in IS_resin_lower_tri_elem_num:
    g.list_of_elements[num-1].element_set = 'isresl_tri'

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
    elif elem.element_set in list_of_tri_lower_elementsets:
        elem.calculate_layer_plane_angle(outer_edge_node_nums=[2,1])
    elif elem.element_set in list_of_tri_upper_elementsets:
        elem.calculate_layer_plane_angle(outer_edge_node_nums=[3,2])
    else:
        raise Warning("Element #{0} has no element set!".format(elem.elem_num))
# manually correct the theta1 values of some elements at the LE
g.list_of_elements[223-1].theta1 = g.list_of_elements[224-1].theta1
g.list_of_elements[237-1].theta1 = g.list_of_elements[238-1].theta1
g.list_of_elements[251-1].theta1 = g.list_of_elements[252-1].theta1

g.list_of_elements[209-1].theta1 = g.list_of_elements[210-1].theta1
g.list_of_elements[195-1].theta1 = g.list_of_elements[196-1].theta1
g.list_of_elements[181-1].theta1 = g.list_of_elements[182-1].theta1

# plot a small selection of elements to check the results
for elem in g.list_of_elements[::skip_num]:
    elem.plot(label_nodes=False)
    print elem.elem_num, elem.element_set, elem.theta1

for num in ES_triax_upper_tri_elem_num:
    g.list_of_elements[num-1].plot()
    g.list_of_elements[num-2].plot()
for num in ES_triax_lower_tri_elem_num:
    g.list_of_elements[num-1].plot()
    g.list_of_elements[num-2].plot()
for num in IS_triax_upper_tri_elem_num:
    g.list_of_elements[num-1].plot()
    g.list_of_elements[num-2].plot()
for num in IS_triax_lower_tri_elem_num:
    g.list_of_elements[num-1].plot()
    g.list_of_elements[num-2].plot()
for num in IS_resin_upper_tri_elem_num:
    g.list_of_elements[num-1].plot()
    g.list_of_elements[num-2].plot()
for num in IS_resin_lower_tri_elem_num:
    g.list_of_elements[num-1].plot()
    g.list_of_elements[num-2].plot()

# show the plot
plt.xlim([-0.05,0.1])
plt.ylim([-0.01,0.02])
ax.set_aspect('equal')
print ' ------------------------'
print '  LEGEND'
print '    magenta : inner edge'
print '    blue    : outer edge'
print ' ------------------------'

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
