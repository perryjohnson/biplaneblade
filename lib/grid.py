"""Create entities for a 2D unstructured grid.

Author: Perry Roth-Johnson
Last modified: March 13, 2014

"""


import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, LineString
from descartes import PolygonPatch


class Node:
    number_of_nodes = 0
    def __init__(self, node_num, x2, x3):
        Node.number_of_nodes += 1
        self.node_num = int(node_num)
        self.x2 = float(x2)
        self.x3 = float(x3)
        self.coords = (self.x2,self.x3)
        self.parent_element = None
        self.is_corner_node = None

    def __str__(self):
        return "Node #{0}: ({1:10.8f}, {2:10.8f})".format(self.node_num,
            self.x2, self.x3)


class _Element:
    number_of_elements = 0
    def __init__(self, elem_num):
        _Element.number_of_elements += 1
        self.elem_num = int(elem_num)
        self.element_set = None
        self.theta1 = None
        self.layer_num = None

    def _plot_edge(self, ax, ob, color='b'):
        x, y = ob.xy
        ax.plot(x, y, color=color, alpha=0.7, linewidth=3,
            solid_capstyle='round', zorder=2)

    def plot(self, equal_aspect_ratio=True, plot_centroid=True,
        label_nodes=True, label_element=True, plot_outer_inner_edges=True):
        """Plot this element."""
        p = Polygon([self.node1.coords, self.node2.coords, self.node3.coords,
            self.node4.coords])
        (cx,cy) = p.centroid.coords.xy
        cx=cx[0]
        cy=cy[0]
        patch = PolygonPatch(p, fc='r', ec=None, alpha=0.5)
        ax = plt.gcf().gca()
        ax.add_patch(patch)
        if equal_aspect_ratio:
            ax.set_aspect('equal')
        if plot_centroid:
            ax.scatter(cx, cy, s=50)
        if label_nodes:
            for i,node in enumerate(self.nodes[0:4]):
                ax.scatter(node.x2,node.x3,s=30)
                fmt_n = '({0}) {1}'.format(i+1, node.node_num)
                ax.text(node.x2,node.x3,fmt_n)
        if label_element:
            if self.theta1 is None:
                # only plot the element number
                ax.text(cx, cy, str(self.elem_num))
            else:
                # plot both the element number and the layer plane angle
                fmt1 = '{0}\n'.format(self.elem_num)
                fmt2 = r'{0:3.0f}$^\circ$'.format(self.theta1)
                fmt=fmt1+fmt2
                ax.text(cx, cy, fmt)
        if plot_outer_inner_edges:
            # create LineString objects along edges
            outer_edge = LineString([self._outer_edge_node0.coords,
                self._outer_edge_node1.coords])
            inner_edge = LineString([self._inner_edge_node0.coords,
                self._inner_edge_node1.coords])
            # plot the edges
            ax = plt.gcf().gca()
            self._plot_edge(ax, outer_edge, color='b')
            self._plot_edge(ax, inner_edge, color='m')

    def calculate_layer_plane_angle(self, outer_edge_node_nums=[1,4],
        inner_edge_node_nums=[2,3], plot_edges=True, plot_angle=True):
        
        # get Node objects on edges and save as non-public attributes
        self._outer_edge_node0 = self.nodes[outer_edge_node_nums[0]-1]
        self._outer_edge_node1 = self.nodes[outer_edge_node_nums[1]-1]
        self._inner_edge_node0 = self.nodes[inner_edge_node_nums[0]-1]
        self._inner_edge_node1 = self.nodes[inner_edge_node_nums[1]-1]
        # calculate the angles of the edges
        outer_angle = np.arctan2(
            self._outer_edge_node1.x3 - self._outer_edge_node0.x3,
            self._outer_edge_node1.x2 - self._outer_edge_node0.x2)
        inner_angle = np.arctan2(
            self._inner_edge_node1.x3 - self._inner_edge_node0.x3,
            self._inner_edge_node1.x2 - self._inner_edge_node0.x2)
        # calc the layer plane angle by averaging the outer and inner angles
        self.theta1 = np.degrees((outer_angle + inner_angle)/2.0)
        if self.theta1 < 0.0:
            self.theta1 += 360.0

class LinearElement(_Element):
    def __init__(self, elem_num, node1, node2, node3, node4):
        _Element.__init__(self, elem_num)
        self.node1 = node1
        self.node2 = node2
        self.node3 = node3
        self.node4 = node4
        self.nodes = (self.node1,self.node2,self.node3,self.node4)
        for node in self.nodes:
            node.parent_element = self

    def __str__(self):
        return "Element #{0}, Nodes({1}, {2}, {3}, {4})".format(self.elem_num,
            self.node1.node_num, self.node2.node_num, self.node3.node_num,
            self.node4.node_num)


class QuadraticElement(_Element):
    def __init__(self, elem_num, node1, node2, node3, node4, node5, node6,
        node7, node8):
        _Element.__init__(self, elem_num)
        self.node1 = node1
        self.node2 = node2
        self.node3 = node3
        self.node4 = node4
        self.node5 = node5
        self.node6 = node6
        self.node7 = node7
        self.node8 = node8
        self.nodes = (self.node1,self.node2,self.node3,self.node4,self.node5,
            self.node6,self.node7,self.node8)
        for node in self.nodes:
            node.parent_element = self

    def __str__(self):
        return "Element #{0}, Nodes({1}, {2}, {3}, {4}, {5}, {6}, {7}, {8})".format(
            self.elem_num, self.node1.node_num, self.node2.node_num,
            self.node3.node_num, self.node4.node_num, self.node5.node_num,
            self.node6.node_num, self.node7.node_num, self.node8.node_num)
