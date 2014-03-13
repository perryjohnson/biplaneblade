"""Determine the layer ply angle of an element in a grid.

Author: Perry Roth-Johnson
Last modified: March 12, 2014

References:
http://stackoverflow.com/questions/3365171/calculating-the-angle-between-two-lines-without-having-to-calculate-the-slope/3366569#3366569
http://stackoverflow.com/questions/19295725/angle-less-than-180-between-two-segments-lines

"""

import numpy as np
import matplotlib.pyplot as plt
import lib.grid as gr
import lib.vabs_utils as vu
from shapely.geometry import Polygon, LineString
from descartes import PolygonPatch


plt.close('all')
# create a figure
plt.figure()
ax = plt.gca()

f = vu.VabsInputFile(
    vabs_filename='sandia_blade/mesh_stn01.vabs',
    grid_filename='sandia_blade/mesh_stn01.abq',
    debug_flag=True)

def plot_line(ax, ob, color='b'):
    x, y = ob.xy
    ax.plot(x, y, color=color, alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)

# list_of_elements = f.grid.element_array[[1,250,500,750,1000,1250,1500,1750,2000,-1]]
list_of_elements = f.grid.element_array[::25]
# save a single element
for elem in list_of_elements:
    n1 = gr.Node(elem['node1'], f.grid.node_array[elem['node1']-1]['x2'], f.grid.node_array[elem['node1']-1]['x3'])
    n2 = gr.Node(elem['node2'], f.grid.node_array[elem['node2']-1]['x2'], f.grid.node_array[elem['node2']-1]['x3'])
    n3 = gr.Node(elem['node3'], f.grid.node_array[elem['node3']-1]['x2'], f.grid.node_array[elem['node3']-1]['x3'])
    n4 = gr.Node(elem['node4'], f.grid.node_array[elem['node4']-1]['x2'], f.grid.node_array[elem['node4']-1]['x3'])
    e = gr.Element(elem['elem_no'], n1, n2, n3, n4)
    print 'element #{0}'.format(e.elem_no)
    # plot the element
    p = Polygon([n1.xy, n2.xy, n3.xy, n4.xy])
    patch = PolygonPatch(p, fc='r', ec=None, alpha=0.5)
    ax.add_patch(patch)
    for i,node in enumerate([n1,n2,n3,n4]):
        # print the nodes in this element
        print node.node_no, node.xy
        # label each node on the plot
        fmt_n = '({0}) {1}'.format(i+1, node.node_no)
        ax.text(node.x2,node.x3,fmt_n)
    # plot the centroid
    (cx,cy) = p.centroid.coords.xy
    cx=cx[0]
    cy=cy[0]
    ax.scatter(cx, cy, s=50)
    if cx < 0:
        # plot the edges used to calculate the angle
        outer_edge = LineString([n1.xy, n4.xy])
        plot_line(ax, outer_edge, color='b')
        inner_edge = LineString([n2.xy, n3.xy])
        plot_line(ax, inner_edge, color='m')
        # calculate the angle of the outer edge
        outer_angle = np.arctan2(n4.x3-n1.x3, n4.x2-n1.x2)
        # calculate the angle of the inner edge
        inner_angle = np.arctan2(n3.x3-n2.x3, n3.x2-n2.x2)
    else:
        outer_edge = LineString([n3.xy, n2.xy])
        plot_line(ax, outer_edge, color='b')
        inner_edge = LineString([n4.xy, n1.xy])
        plot_line(ax, inner_edge, color='m')
        # calculate the angle of the outer edge
        outer_angle = np.arctan2(n2.x3-n3.x3, n2.x2-n3.x2)
        # calculate the angle of the inner edge
        inner_angle = np.arctan2(n1.x3-n4.x3, n1.x2-n4.x2)
    # calculate the layer ply angle by averaging the outer and inner angles,
    #   then convert to degrees and print
    layer_ply_angle = np.degrees((outer_angle + inner_angle)/2.0)
    if layer_ply_angle < 0.0:
        layer_ply_angle = layer_ply_angle + 360.0
    # save the layer ply angle as an attribute of this element
    e.theta1 = layer_ply_angle
    print 'layer ply angle =', e.theta1, 'degrees'
    print ''
    # plot the element number, and on the next line, the layer ply angle
    fmt1 = '{0}\n'.format(elem['elem_no'])
    fmt2 = r'{0:3.0f}$^\circ$'.format(layer_ply_angle)
    fmt=fmt1+fmt2
    ax.text(cx, cy, fmt)
# show the plot
plt.xlim([-3,3])
plt.ylim([-3,3])
ax.set_aspect('equal')
plt.show()
