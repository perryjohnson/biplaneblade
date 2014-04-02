"""Create entities for a 2D unstructured grid.

Author: Perry Roth-Johnson
Last modified: April 2, 2014

"""


import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, LineString
from shapely.geometry.polygon import orient
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
    def __init__(self, elem_num, layer_num):
        _Element.number_of_elements += 1
        self.elem_num = int(elem_num)
        self.element_set = None
        self.theta1 = None
        self.layer_num = layer_num

    def swap_nodes(self, nodeA, nodeB):
        temp = nodeA
        nodeA = nodeB
        nodeB = temp
        return (nodeA, nodeB)

    def _plot_edge(self, ax, ob, color='b'):
        x, y = ob.xy
        ax.plot(x, y, color=color, alpha=0.7, linewidth=3,
            solid_capstyle='round', zorder=2)

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

class QuadrilateralLinearElement(_Element):
    """A quadrilateral linear element with four nodes.

The VABS quadrilateral element node numbering scheme is shown below:
    4-------3
    |       |
    |       |
    |       |
    1-------2

    """
    def __init__(self, elem_num, node1, node2, node3, node4, layer_num):
        _Element.__init__(self, elem_num, layer_num)
        self.node1 = node1
        self.node2 = node2
        self.node3 = node3
        self.node4 = node4
        # assign node_num=0 for nodes that are not present
        self.node5 = Node(0, 0.0, 0.0)
        self.node6 = Node(0, 0.0, 0.0)
        self.node7 = Node(0, 0.0, 0.0)
        self.node8 = Node(0, 0.0, 0.0)
        self.node9 = Node(0, 0.0, 0.0)
        self.nodes = (self.node1,self.node2,self.node3,self.node4)
        for node in self.nodes:
            node.parent_element = self
        # Calculate the middle coordinates of this element
        self.x2_middle = (self.node1.x2 + self.node2.x2
                          + self.node3.x2 + self.node4.x2)/4.0
        self.x3_middle = (self.node1.x3 + self.node2.x3
                          + self.node3.x3 + self.node4.x3)/4.0
        if not self.is_ccw():      # Try to fix a bad (CW) element
            self.reorder_nodes()   # Hopefully this makes a good (CCW) element
            print "  Trying to fix element #{0}. Reorienting...".format(
                self.elem_num)
        # Check if the element actually has CCW orientation
        if not self.is_ccw():
            fmt = "Element #{:d} is bad! Its nodes are not oriented CCW."
            raise Warning(fmt.format(self.elem_num))

    def __str__(self):
        return """Element #{0} -----
  Nodes({1}, {2}, {3}, {4})
  Layer #{9}""".format(self.elem_num,
            self.node1.node_num, self.node2.node_num, self.node3.node_num,
            self.node4.node_num, self.layer_num)

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

    def _angles(self, print_flag=False):
        """
Returns a dictionary of angles for each node in this element.

Each dictionary value represents the angle between two vectors: (1) a vector
drawn from the middle of the element to a node, and (2) a vector along the
horizontal (positive x2-axis).

Optionally, if print_flag=True, print the dictionary to the screen with nice
formatting.

        """
        from math import atan2, pi
        angle_dict={}
        # find the angle to node1
        y = self.node1.x3 - self.x3_middle
        x = self.node1.x2 - self.x2_middle
        angle_dict['node1'] = atan2(y,x) * (180.0/pi)
        if angle_dict['node1'] < 0.0:
            angle_dict['node1'] += 360.0
        # find the angle to node2
        y = self.node2.x3 - self.x3_middle
        x = self.node2.x2 - self.x2_middle
        angle_dict['node2'] = atan2(y,x) * (180.0/pi)
        if angle_dict['node2'] < 0.0:
            angle_dict['node2'] += + 360.0
        # find the angle to node3
        y = self.node3.x3 - self.x3_middle
        x = self.node3.x2 - self.x2_middle
        angle_dict['node3'] = atan2(y,x) * (180.0/pi)
        if angle_dict['node3'] < 0.0:
            angle_dict['node3'] += + 360.0
        # find the angle to node4
        y = self.node4.x3 - self.x3_middle
        x = self.node4.x2 - self.x2_middle
        angle_dict['node4'] = atan2(y,x) * (180.0/pi)
        if angle_dict['node4'] < 0.0:
            angle_dict['node4'] += + 360.0
        if print_flag:
            print "node1_angle = {:6.2f} degrees".format(angle_dict['node1'])
            print "node2_angle = {:6.2f} degrees".format(angle_dict['node2'])
            print "node3_angle = {:6.2f} degrees".format(angle_dict['node3'])
            print "node4_angle = {:6.2f} degrees".format(angle_dict['node4'])
        return angle_dict

    def is_ccw(self, print_flag=False):
        """
Returns a boolean True if the nodes of this element are arranged counter-
clockwise; False otherwise.

Usage:
from vabs_objects import *
a = Node(0,0)
b = Node(1,0)
c = Node(1,1)
d = Node(0,1)
matl = IsotropicMaterial(name='Foam', rho=0.20E+03, E=0.256E+09, nu=0.3)
L1 = Layer(matl, 45.3)
elem1 = LinearQuadElement(L1, 270, a, b, c, d)
elem1.is_ccw()

        """
        angle_dict = self._angles()         # Update the angles to each node.
        angleDiff = []
        angleDiff.append(angle_dict['node2'] - angle_dict['node1'])
        angleDiff.append(angle_dict['node3'] - angle_dict['node2'])
        angleDiff.append(angle_dict['node4'] - angle_dict['node3'])
        angleDiff.append(angle_dict['node1'] - angle_dict['node4'])
        if print_flag:
            print "Angle differences:", angleDiff
        negCounter = 0
        for i in range(len(angleDiff)):
            if angleDiff[i] < 0.0:
                negCounter += 1
        if negCounter > 1:
            flag = False
            if print_flag:
                print "this element has the WRONG orientation"
        else:
            flag = True
        return flag

    def reorder_nodes(self, print_flag=False):
        """
Reorients the nodes from a clockwise orientation to a counter-clockwise
orientation, or vice-versa.

Note: if the nodes are not oriented clockwise, this method will not fix the
nodes to have counter-clockwise orientation!

Usage:
from vabs_objects import *
a = Node(0,0)
b = Node(1,0)
c = Node(1,1)
d = Node(0,1)
matl = IsotropicMaterial(name='Foam', rho=0.20E+03, E=0.256E+09, nu=0.3)
L1 = Layer(matl, 45.3)
elem1 = LinearQuadElement(L1, 270, a, d, c, b)
elem1.is_ccw(print_flag=True)
elem1.reorder_nodes(print_flag=True)
elem1.is_ccw(print_flag=True)
elem1.reorder_nodes(print_flag=True)
elem1.is_ccw(print_flag=True)

        """
        if print_flag:
            print "BEFORE:"
            print "  node2: #{:d},    node4: #{:d}".format(self.node2.node_num,
                                                           self.node4.node_num)
        # Reorder the nodes in a counter-clockwise fashion.
        (self.node2,self.node4) = self.swap_nodes(self.node2,self.node4)
        if print_flag:
            print "AFTER:"
            print "  node2: #{:d},    node4: #{:d}".format(self.node2.node_num,
                                                           self.node4.node_num)


class QuadrilateralQuadraticElement(_Element):
    """A quadrilateral quadratic element with eight nodes.

The VABS quadrilateral element node numbering scheme is shown below:
    4---7---3
    |       |
    8       6
    |       |
    1---5---2
The central node (9) is optional.

    """
    def __init__(self, elem_num, node1, node2, node3, node4, node5, node6,
        node7, node8, layer_num, autocorrect=True):
        _Element.__init__(self, elem_num, layer_num)
        # create a polygon representation
        p = Polygon([
            node1.coords, 
            node5.coords, 
            node2.coords, 
            node6.coords, 
            node3.coords,
            node7.coords,
            node4.coords,
            node8.coords
            ])
        if autocorrect:
            p_ccw = orient(p)
            self.polygon = p_ccw
            # find each node that matches each polygon point
            node_list = []
            for point in p_ccw.exterior.coords[:-1]:
                if point == node1.coords:
                    node_list.append(node1)
                elif point == node2.coords:
                    node_list.append(node2)
                elif point == node3.coords:
                    node_list.append(node3)
                elif point == node4.coords:
                    node_list.append(node4)
                elif point == node5.coords:
                    node_list.append(node5)
                elif point == node6.coords:
                    node_list.append(node6)
                elif point == node7.coords:
                    node_list.append(node7)
                elif point == node8.coords:
                    node_list.append(node8)
                else:
                    raise Warning("No point match found!")
            # assign the nodes in a proper CCW-orientation
            self.node1 = node_list[0]
            self.node5 = node_list[1]
            self.node2 = node_list[2]
            self.node6 = node_list[3]
            self.node3 = node_list[4]
            self.node7 = node_list[5]
            self.node4 = node_list[6]
            self.node8 = node_list[7]
        else:
            self.polygon = p
            self.node1 = node1
            self.node2 = node2
            self.node3 = node3
            self.node4 = node4
            self.node5 = node5
            self.node6 = node6
            self.node7 = node7
            self.node8 = node8
        # assign node_num=0 for nodes that are not present
        self.node9 = Node(0, 0.0, 0.0)
        self.nodes = (self.node1,self.node2,self.node3,self.node4,self.node5,
            self.node6,self.node7,self.node8)
        for node in self.nodes:
            node.parent_element = self

    def __str__(self):
        return """Element #{0} -----
  Nodes({1}, {2}, {3}, {4}, {5}, {6}, {7}, {8})
  Layer #{9}
  element set: {10}
  centroid: ({11:5.3f}, {12:5.3f})
  theta1 = {13} degrees""".format(
            self.elem_num, 
            self.node1.node_num, self.node2.node_num, self.node3.node_num,
              self.node4.node_num, self.node5.node_num, self.node6.node_num, 
              self.node7.node_num, self.node8.node_num,
            self.layer_num,
            self.element_set,
            self.polygon.centroid.x, self.polygon.centroid.y,
            self.theta1)

    def plot(self, equal_aspect_ratio=True, plot_centroid=True,
        label_nodes=True, label_element=True, plot_outer_inner_edges=True):
        """Plot this element."""
        (cx,cy) = self.polygon.centroid.coords.xy
        cx=cx[0]
        cy=cy[0]
        patch = PolygonPatch(self.polygon, fc='r', ec=None, alpha=0.5)
        ax = plt.gcf().gca()
        ax.add_patch(patch)
        if equal_aspect_ratio:
            ax.set_aspect('equal')
        if plot_centroid:
            ax.scatter(cx, cy, s=50, c='m')
        if label_nodes:
            for i,node in enumerate(self.nodes):
                ax.scatter(node.x2,node.x3,s=30)
                fmt_n = '({0}) {1}'.format(i+1, node.node_num)
                ax.text(node.x2,node.x3,fmt_n)
        if label_element:
            if self.theta1 is None:
                # only plot the element number
                ax.text(cx, cy, str(self.elem_num))
            else:
                # plot both the element number and the layer plane angle
                fmt1 = '({0}) '.format(self.elem_num)
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


class TriangularLinearElement(_Element):
    """A triangular lienar element with three nodes.

The VABS triangular element node numbering scheme is shown below:
        3
       / \
      /   \  
     /     \
    1-------2

    """
    def __init__(self, elem_num, node1, node2, node3, layer_num):
        _Element.__init__(self, elem_num, layer_num)
        # create a polygon representation
        p = Polygon([
            node1.coords, 
            node2.coords, 
            node3.coords
            ])
        self.polygon = p
        self.node1 = node1
        self.node2 = node2
        self.node3 = node3
        # assign node_num=0 for nodes that are not present
        self.node4 = Node(0, 0.0, 0.0)
        self.node5 = Node(0, 0.0, 0.0)
        self.node6 = Node(0, 0.0, 0.0)
        self.node7 = Node(0, 0.0, 0.0)
        self.node8 = Node(0, 0.0, 0.0)
        self.node9 = Node(0, 0.0, 0.0)
        self.nodes = (self.node1,self.node2,self.node3)
        for node in self.nodes:
            node.parent_element = self
        
    def __str__(self):
        return """Element #{0} -----
  Nodes({1}, {2}, {3})
  Layer #{4}
  element set: {5}
  centroid: ({6:5.3f}, {7:5.3f})
  theta1 = {8} degrees""".format(
            self.elem_num, 
            self.node1.node_num, self.node2.node_num, self.node3.node_num,
            self.layer_num,
            self.element_set,
            self.polygon.centroid.x, self.polygon.centroid.y,
            self.theta1)

    # def plot(self):
    # !!! FILL IN !!!


class TriangularQuadraticElement(_Element):
    """A triangular quadratic element with six nodes.

The VABS triangular element node numbering scheme is shown below:
        3
       / \
      7   6  
     /     \
    1---5---2

    """
    def __init__(self, elem_num, node1, node2, node3, node5, node6, node7,
        layer_num, autocorrect=True):
        _Element.__init__(self, elem_num, layer_num)
        # create a polygon representation
        p = Polygon([
            node1.coords, 
            node5.coords, 
            node2.coords, 
            node6.coords, 
            node3.coords,
            node7.coords
            ])
        if autocorrect:
            p_ccw = orient(p)
            self.polygon = p_ccw
            # find each node that matches each polygon point
            node_list = []
            for point in p_ccw.exterior.coords[:-1]:
                if point == node1.coords:
                    node_list.append(node1)
                elif point == node2.coords:
                    node_list.append(node2)
                elif point == node3.coords:
                    node_list.append(node3)
                elif point == node5.coords:
                    node_list.append(node5)
                elif point == node6.coords:
                    node_list.append(node6)
                elif point == node7.coords:
                    node_list.append(node7)
                else:
                    raise Warning("No point match found!")
            # assign the nodes in a proper CCW-orientation
            self.node1 = node_list[0]
            self.node5 = node_list[1]
            self.node2 = node_list[2]
            self.node6 = node_list[3]
            self.node3 = node_list[4]
            self.node7 = node_list[5]
        else:
            self.polygon = p
            self.node1 = node1
            self.node2 = node2
            self.node3 = node3
            self.node5 = node5
            self.node6 = node6
            self.node7 = node7
        # assign node_num=0 for nodes that are not present
        self.node4 = Node(0, 0.0, 0.0)
        self.node8 = Node(0, 0.0, 0.0)
        self.node9 = Node(0, 0.0, 0.0)
        self.nodes = (self.node1,self.node2,self.node3,self.node5,self.node6,
            self.node7)
        for node in self.nodes:
            node.parent_element = self
        
    def __str__(self):
        return """Element #{0} -----
  Nodes({1}, {2}, {3}, {4}, {5}, {6})
  Layer #{7}
  element set: {8}
  centroid: ({9:5.3f}, {10:5.3f})
  theta1 = {11} degrees""".format(
            self.elem_num, 
            self.node1.node_num, self.node2.node_num, self.node3.node_num,
              self.node5.node_num, self.node6.node_num, self.node7.node_num,
            self.layer_num,
            self.element_set,
            self.polygon.centroid.x, self.polygon.centroid.y,
            self.theta1)

    def plot(self, equal_aspect_ratio=True, plot_centroid=True,
        label_nodes=True, label_element=True, plot_outer_edge=True):
        """Plot this element."""
        (cx,cy) = self.polygon.centroid.coords.xy
        cx=cx[0]
        cy=cy[0]
        patch = PolygonPatch(self.polygon, fc='r', ec=None, alpha=0.5)
        ax = plt.gcf().gca()
        ax.add_patch(patch)
        if equal_aspect_ratio:
            ax.set_aspect('equal')
        if plot_centroid:
            ax.scatter(cx, cy, s=50, c='m')
        if label_nodes:
            for i,node in enumerate(self.nodes):
                ax.scatter(node.x2,node.x3,s=30)
                if i < 3:
                    fmt_n = '({0}) {1}'.format(i+1, node.node_num)
                else:
                    fmt_n = '({0}) {1}'.format(i+2, node.node_num)
                ax.text(node.x2,node.x3,fmt_n)
        if label_element:
            if self.theta1 is None:
                # only plot the element number
                ax.text(cx, cy, str(self.elem_num))
            else:
                # plot both the element number and the layer plane angle
                fmt1 = '({0}) '.format(self.elem_num)
                fmt2 = r'{0:3.0f}$^\circ$'.format(self.theta1)
                fmt=fmt1+fmt2
                ax.text(cx, cy, fmt)
        if plot_outer_edge:
            # create LineString objects along edge
            outer_edge = LineString([self._outer_edge_node0.coords,
                self._outer_edge_node1.coords])
            # plot the edge
            ax = plt.gcf().gca()
            self._plot_edge(ax, outer_edge, color='b')

    def calculate_layer_plane_angle(self, outer_edge_node_nums=[2,3],
        inner_edge_node_nums=None, plot_edges=True, plot_angle=True):
        if inner_edge_node_nums is not None:
            raise Warning("Element #{0} does not have an inner edge defined!".format(elem.elem_num))
        # get Node objects on edges and save as non-public attributes
        self._outer_edge_node0 = self.nodes[outer_edge_node_nums[0]-1]
        self._outer_edge_node1 = self.nodes[outer_edge_node_nums[1]-1]
        # calculate the angles of the edges
        outer_angle = np.arctan2(
            self._outer_edge_node1.x3 - self._outer_edge_node0.x3,
            self._outer_edge_node1.x2 - self._outer_edge_node0.x2)
        # calc the layer plane angle by taking the outer angle
        self.theta1 = np.degrees(outer_angle)
        if self.theta1 < 0.0:
            self.theta1 += 360.0
