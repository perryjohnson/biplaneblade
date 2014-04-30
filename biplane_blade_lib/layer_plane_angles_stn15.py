"""Determine the layer plane angle of all the elements in a grid.

Author: Perry Roth-Johnson
Last modified: April 30, 2014

Usage:
1. Look through the mesh_stnXX.abq file and find all the element set names.
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
station_num = 15
skip_num = 25   # plot every 'skip_num' elements (larger values plot faster)
IS4_resin_u2_tri_elem_num = 2896  # num of tri elem in int surf 3 resin upper 2
IS4_resin_l2_tri_elem_num = 2883  # num of tri elem in int surf 3 resin lower 2
TE_reinf_foam_u3_tri_elem_num = 419  # num of tri elem in TE reinf foam upper 3
TE_reinf_foam_l3_tri_elem_num = 395  # num of tri elem in TE reinf foam lower 3
# -----------------------------------------------

stn_str = 'stn{0:02d}'.format(station_num)
plt.close('all')

# load the biplane blade
b1 = bl.BiplaneBlade(
    'biplane blade, flapwise symmetric, no stagger, rj/R=0.452, g/c=1.25',
    'biplane_blade')

# pre-process the station dimensions
station = b1.list_of_stations[station_num-1]
st = station.structure
af = station.airfoil
af.create_polygon()
st.create_all_layers()
st.save_all_layer_edges()
st.write_all_part_polygons()
x3_off = af.lower_chord * af.gap_to_chord_ratio * af.gap_fraction

# plot the parts
station.plot_parts_offset(airfoil_to_plot='lower', x3_offset=x3_off)



# create a figure
ax = plt.gcf().gca()

# element sets on the leading edge
# outer_edge_node_nums=[1,4], inner_edge_node_nums=[2,3]
list_of_LE_elementsets = [
    'lepanel',
    'rbtrile',
    'esgelle',
    'estrile',
    'is1rsle',
    'is1trile',
    'sw1biaxl',
    'sw1foam',
    'sw1biaxr',
    'is1rlsw1',
    'is1tlsw1',
    'is2rrsw1',
    'is2trsw1'
    ]

# element sets on the trailing edge
# outer_edge_node_nums=[3,2], inner_edge_node_nums=[4,1]
list_of_TE_elementsets = [
    'is3reste',
    'is3trite',
    'sw2biaxl',
    'sw2foam',
    'sw2biaxr',
    'sw3biaxl',
    'sw3foam',
    'sw3biaxr',
    'is2rlsw2',
    'is2tlsw2',
    'is3rrsw2',
    'is3trsw2',
    'is3rlsw3',
    'is3tlsw3',
    'is4rrsw3',
    'is4trsw3'
    ]

# element sets on the lower surface
# outer_edge_node_nums=[2,1], inner_edge_node_nums=[3,4]
list_of_lower_elementsets = [
    'tefoaml1',
    'tefoaml2',
    'tefoaml3',
    'teunil1',
    'teunil2',
    'teunil3',
    'teunil4',
    'ap1lower',
    'ap2lower',
    'esgllap1',
    'estrlap1',
    'rbtrlap1',
    'is3rlap1',
    'is3tlap1',
    'esgllap2',
    'estrlap2',
    'is4rlap2',
    'is4tlap2',
    'sclower',
    'rbtriscl',
    'rbtrtel1',
    'rbtrtel2',
    'rbtrtel3',
    'rbtrtel4',
    'estrtel1',
    'estrtel2',
    'estrtel3',
    'estrtel4',
    'esgltel1',
    'esgltel2',
    'esgltel3',
    'esgltel4',
    'esgelscl',
    'estriscl',
    'is2rsscl',
    'is2trscl',
    'is4rtel1',
    'is4rtel2',
    'is4ttel1',
    'is4ttel2',
    'rbtrbsw1',
    'esglbsw1',
    'estrbsw1',
    'esglbsw2',
    'estrbsw2',
    'rbtrbsw2',
    'estrbsw3',
    'esglbsw3'
]

# element sets on the upper surface
# outer_edge_node_nums=[4,3], inner_edge_node_nums=[1,2]
list_of_upper_elementsets = [
    'tefoamu1',
    'tefoamu2',
    'tefoamu3',
    'teuniu1',
    'teuniu2',
    'teuniu3',
    'teuniu4',
    'ap1upper',
    'ap2upper',
    'esgluap1',
    'estruap1',
    'rbtruap1',
    'is3ruap1',
    'is3tuap1',
    'esgluap2',
    'estruap2',
    'is4ruap2',
    'is4tuap2',
    'is4rteu1',
    'is4rteu2',
    'is4tteu1',
    'is4tteu2',
    'esglasw1',
    'estrasw1',
    'rbtrasw1',
    'rbtrteu1',
    'rbtrteu2',
    'rbtrteu3',
    'rbtrteu4',
    'estrteu1',
    'estrteu2',
    'estrteu3',
    'estrteu4',
    'esglteu1',
    'esglteu2',
    'esglteu3',
    'esglteu4',
    'esglasw2',
    'estrasw2',
    'rbtrasw2',
    'esglasw3',
    'estrasw3',
    'scupper',
    'rbtriscu',
    'estriscu',
    'is2rsscu',
    'esgelscu',
    'is2trscu'
]

# element sets of triangular elements on the lower surface
# outer_edge_node_nums=[2,1]
list_of_tri_lower_elementsets = [
    'is4rtel2_tri',
    'tefoaml3_tri'
]

# element sets of triangular elements on the upper surface
# outer_edge_node_nums=[3,2]
list_of_tri_upper_elementsets = [
    'is4rteu2_tri',
    'tefoamu3_tri'
]


# import the initial grid object
fmt_grid = 'biplane_blade/' + stn_str + '/mesh_' + stn_str + '.abq'
g = au.AbaqusGrid(fmt_grid, debug_flag=True)

# manually assign two triangular elements into new element sets
g.list_of_elements[IS4_resin_u2_tri_elem_num-1].element_set = 'is4rteu2_tri'
g.list_of_elements[IS4_resin_l2_tri_elem_num-1].element_set = 'is4rtel2_tri'
g.list_of_elements[TE_reinf_foam_u3_tri_elem_num-1].element_set = 'tefoamu3_tri'
g.list_of_elements[TE_reinf_foam_l3_tri_elem_num-1].element_set = 'tefoaml3_tri'

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
# plot a small selection of elements to check the results
for elem in g.list_of_elements[::skip_num]:
    elem.plot(label_nodes=False)
    print elem.elem_num, elem.element_set, elem.theta1

g.list_of_elements[TE_reinf_foam_u3_tri_elem_num-1].plot()
g.list_of_elements[TE_reinf_foam_u3_tri_elem_num-2].plot()
g.list_of_elements[TE_reinf_foam_l3_tri_elem_num-1].plot()
g.list_of_elements[TE_reinf_foam_l3_tri_elem_num-2].plot()
g.list_of_elements[IS4_resin_u2_tri_elem_num-1].plot()
g.list_of_elements[IS4_resin_u2_tri_elem_num-2].plot()
g.list_of_elements[IS4_resin_l2_tri_elem_num-1].plot()
g.list_of_elements[IS4_resin_l2_tri_elem_num-2].plot()

# show the plot
plt.xlim([-3,4.5])
plt.ylim([-3,3])
ax.set_aspect('equal')
print ' ------------------------'
print '  LEGEND'
print '    magenta : inner edge'
print '    blue    : outer edge'
print ' ------------------------'

plt.show()
# -----------------------------------------------------------------------------
# read layers.csv to determine the number of layers
layer_file = pd.read_csv('biplane_blade/layers.csv', index_col=0)
number_of_layers = len(layer_file)
# write the updated grid object to a VABS input file
fmt_vabs = 'biplane_blade/' + stn_str + '/mesh_' + stn_str + '.vabs'
f = vu.VabsInputFile(
    vabs_filename=fmt_vabs,
    grid=g,
    material_filename='biplane_blade/materials.csv',
    layer_filename='biplane_blade/layers.csv',
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
