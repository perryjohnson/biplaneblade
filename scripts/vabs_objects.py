"""
Classes and modules that organize the data for a VABS input file.

Author: Perry Roth-Johnson
Last updated: October 24, 2012

"""


import numpy as np


class Node:
    """
Make a node (vertex) in the cross-sectional mesh.

Class ("global") attributes:
number_of_nodes - The total number of all nodes in the cross-sectional mesh.

Object ("local") attributes:
node_no - An integer that represents this unique node.
x2 - A float for the x2-coordinate of this node. (The horizontal/chordwise
    coordinate.)
x3 - A float for the x3-coordinate of this node. (The vertical/flapwise
    coordinate.)
element_list - A list of Element instances. This node belongs to each of the
    elements in element_list.

Usage:
from vabs_objects import *
a = Node(x2=1.0, x3=-0.35)
Node.number_of_nodes
a.node_no
a.x2
a.x3
print a
b = Node(node_no=15, x2=10, x3=-0.4)
Node.number_of_nodes
print b

    """
    number_of_nodes = 0
    node_list = []
    def __init__(self, x2, x3, node_no=None):
        Node.number_of_nodes += 1
        if node_no == None:
            self.node_no = Node.number_of_nodes
        else:
            self.node_no = node_no
        self.x2 = x2
        self.x3 = x3
        self.element_list = []
        Node.node_list.append(self)
    
    def __str__(self):
        L = []
        L.append('INSPECTING NODE #{:d}\n'.format(self.node_no))
        L.append('  coordinates (x2, x3) = ({:g}, {:g})\n'.format(self.x2, self.x3))
        elem_list = [elem.elem_no for elem in self.element_list]
        L.append('  belongs to elements: ' + str(elem_list))
        return "".join(L)


def assign_coordinates_to_nodes(data):
    """
Assign an array of coordinates to Node instances.

Input:
data - A comma-separated record array with three columns: 'node_no', 'x2', and
    'x3'.

Usage:
from vabs_objects import *
from matplotlib.mlab import csv2rec
f = open('node.txt', 'w')
datastring = '1,-0.836,-0.254,0.E+00\n2,-0.836,-0.25338,0.E+00\n3,-0.836,-0.25275,0.E+00\n4,-0.836,-0.25213,0.E+00\n5,-0.836,-0.2515,0.E+00\n6,-0.836,-0.25087,0.E+00\n7,-0.836,-0.25025,0.E+00\n8,-0.836,-0.24962,0.E+00\n9,-0.836,-0.249,0.E+00\n10,-0.83581,-0.254,0.E+00'
f.writelines(datastring)
f.close()
mydata = csv2rec('node.txt', names=['node_no', 'x2', 'x3', 'x1'])
assign_coordinates_to_nodes(mydata)
from shutil import os
os.remove('node.txt')
for x in Node.node_list:
    print x


    """
    for i in range(len(data)):
        Node(x2=data['x2'][i], x3=data['x3'][i])


class EmptyNode(Node):
    """
Make an empty node that can be used as a placeholder.

This is useful for elements that contain less than 9 nodes. An EmptyNode can be
assigned to the nodes not directly used by these elements. For example, a
linear quadrilateral element only uses nodes 1-4. Hence, we should assign an
EmptyNode to nodes 5-9 in this element. This is constructed to follow the
conventions detailed in the VABS Manual for Users, page 11.

Object ("local") attributes:
node_no - An integer that represents this unique node. For an EmptyNode,
    node_no = 0 always.

Usage:
from vabs_objects import *
not_a_node = EmptyNode()
not_a_node.node_no
print not_a_node

    """
    def __init__(self):
        self.node_no = 0
    
    def __str__(self):
        return 'EmptyNode:  node_no = {:d}'.format(self.node_no)


class Material:
    """
Make a material with all its material properties.

Typically, a user should not directly create instances of this class. Instead,
users should create an instance of one of the three subclasses:
IsotropicMaterial, OrthotropicMaterial, or GeneralAnisotropicMaterial. This
material (and its associated material properties) can then be assigned to
elements in the cross-sectional mesh.

Class ("global") attributes:
number_of_materials - The total number of all materials in the cross-section.

Object ("local") attributes:
material_no - An integer that represents this unique material.
name - A string for the written name of the material.
orth_flag - An integer flag (0, 1, or 2) to indicate if the material is
    isotropic (0), orthotropic (1), or general anisotropic (2)
rho - A float for the density of the material in units of [kg/m^2].

Usage:
Typically, a user should not directly create instances of this class. Instead,
users should create an instance of one of the three subclasses:
IsotropicMaterial, OrthotropicMaterial, or GeneralAnisotropicMaterial.

    """
    number_of_materials = 0
    def __init__(self, name, orth_flag, rho):
        Material.number_of_materials += 1
        self.material_no = Material.number_of_materials
        self.name = name
        self.orth_flag = int(orth_flag)
        self.rho = float(rho)


class IsotropicMaterial(Material):
    """
Make an isotropic material with three material properties (rho, E, nu).

This material (and its associated material properties) can then be 
assigned to elements in the cross-sectional mesh.

Object ("local") attributes:
material_no - An integer that represents this unique material.
name - A string for the written name of the material.
orth_flag - An integer flag (set to 0) that indicates this material is
    isotropic.
rho - A float for the density of the material in units of [kg/m^2].
E - A float for the Young's modulus of the material in units of [Pa].
nu - A float for the Poisson's ratio of the material (unitless).

Usage:
from vabs_objects import *
b = IsotropicMaterial(name='Foam', rho=0.20E+03, E=0.256E+09, nu=0.3)
print b

    """
    def __init__(self, name, rho, E, nu):
        Material.__init__(self, name, 0, rho)           # Force orth_flag = 0
        self.E = float(E)
        self.nu = float(nu)

    def __str__(self):
        L = []
        L.append('INSPECTING MATERIAL #{:d} ({:s})\n'.format(self.material_no,
                                                             self.name))
        L.append('  orth_flag = {:d} (isotropic)\n'.format(self.orth_flag))
        L.append('  rho = {:g} [kg/m^2]\n'.format(self.rho))
        L.append('  E = {:g} [Pa]\n'.format(self.E))
        L.append('  nu = {:g} [-]'.format(self.nu))
        return "".join(L)


class OrthotropicMaterial(Material):
    """
Make an orthotropic material with ten material properties (rho, E1, E2, E3,
    G12, G13, G23, nu12, nu13, nu23).

This material (and its associated material properties) can then be 
assigned to elements in the cross-sectional mesh.

Object ("local") attributes:
material_no - An integer that represents this unique material.
name - A string for the written name of the material.
orth_flag - An integer flag (set to 1) that indicates this material is
    orthotropic.
rho - A float for the density of the material.
E1 - A float for the Young's modulus of the material in the x1-direction in 
    units of [Pa].
E2 - A float for the Young's modulus of the material in the x2-direction in 
    units of [Pa].
E3 - A float for the Young's modulus of the material in the x3-direction in 
    units of [Pa].
G12 - A float for the shear modulus in the x1x2-plane in units of [Pa].
G13 - A float for the shear modulus in the x1x3-plane in units of [Pa].
G23 - A float for the shear modulus in the x2x3-plane in units of [Pa].
nu12 - A float for the Poisson's ratio of the material in the x1x2-plane
    (unitless).
nu13 - A float for the Poisson's ratio of the material in the x1x3-plane
    (unitless).
nu23 - A float for the Poisson's ratio of the material in the x2x3-plane
    (unitless).

Usage:
from vabs_objects import *
uniaxial_GFRP = OrthotropicMaterial(
    name='E-LT-5500/EP-3', rho=1920,
    E1=41.8E+09, E2=14.0E+09, E3=14.0E+09,
    G12=2.63E+09, G13=2.63E+09, G23=2.63E+09,
    nu12=0.28, nu13=0.28, nu23=0.28)
biaxial_GFRP = OrthotropicMaterial(
    name='Saertex/EP-3', rho=1780,
    E1=13.6E+09, E2=13.3E+09, E3=13.3E+09,
    G12=11.8E+09, G13=11.8E+09, G23=11.8E+09,
    nu12=0.51, nu13=0.51, nu23=0.51)
triaxial_GFRP = OrthotropicMaterial(
    name='SNL Triax', rho=1850,
    E1=27.7E+09, E2=13.65E+09, E3=13.65E+09,
    G12=7.20E+09, G13=7.20E+09, G23=7.20E+09,
    nu12=0.39, nu13=0.39, nu23=0.39)
print uniaxial_GFRP
print biaxial_GFRP
print triaxial_GFRP

    """
    def __init__(self, name, rho, E1, E2, E3, G12, G13, G23, nu12, nu13, nu23):
        Material.__init__(self, name, 1, rho)           # Force orth_flag = 1
        self.E1, self.E2, self.E3 = float(E1), float(E2), float(E3)
        self.G12, self.G13, self.G23 = float(G12), float(G13), float(G23)
        self.nu12, self.nu13, self.nu23 = float(nu12), float(nu13), float(nu23)
    
    def __str__(self):
        L = []
        L.append('INSPECTING MATERIAL #{:d} ({:s})\n'.format(self.material_no,
                                                             self.name))
        L.append('  orth_flag = {:d} (isotropic)\n'.format(self.orth_flag))
        L.append('  rho = {:g} [kg/m^2]\n'.format(self.rho))
        L.append('  E1 = {:g} [Pa]\n'.format(self.E1))
        L.append('  E2 = {:g} [Pa]\n'.format(self.E2))
        L.append('  E3 = {:g} [Pa]\n'.format(self.E3))
        L.append('  G12 = {:g} [Pa]\n'.format(self.G12))
        L.append('  G13 = {:g} [Pa]\n'.format(self.G13))
        L.append('  G23 = {:g} [Pa]\n'.format(self.G23))
        L.append('  nu12 = {:g} [-]\n'.format(self.nu12))
        L.append('  nu13 = {:g} [-]\n'.format(self.nu13))
        L.append('  nu23 = {:g} [-]'.format(self.nu23))
        return "".join(L)


class Layer:
    """
Make a layer, which is defined as a unique combination of material type and
layup angle.

Class ("global") attributes:
number_of_layers - The total number of all layers in the cross-section.

Object ("local") attributes:
layer_no - An integer that represents this unique layer.
material - An instance of the IsotropicMaterial or OrthotropicMaterial class.
theta3 - A float for the layup angle (in degrees). Acceptable values are
    -90 <= theta3 <= 90 degrees. The ply coordinate system (y1,y2,y3) is
    rotated about y3 in the right-hand sense by the amount -90 <= theta3 <= 90
    degrees to form the material system (e1,e2,e3). The range of theta3 is the
    same as is commonly defined in the field of composite materials. For more
    details, see the VABS Manual for Users, pages 8-9, figures 4-5.

Usage:
from vabs_objects import *
b = IsotropicMaterial(name='Foam', rho=0.20E+03, E=0.256E+09, nu=0.3)
L1 = Layer(b, 45.3)
print L1

    """
    number_of_layers = 0
    def __init__(self, material, theta3):
        Layer.number_of_layers += 1
        self.layer_no = Layer.number_of_layers
        self.material = material
        self.theta3 = float(theta3)

    def __str__(self):
        L = []
        L.append('INSPECTING LAYER #{:d} ----------\n'.format(self.layer_no))
        L.append('  theta3 = {:g} [deg]\n'.format(self.theta3))
        L.append(self.material.__str__() + '\n')
        L.append('-'*30)
        return "".join(L)


class Element:
    """
Make an element in the cross-sectional mesh.

Typically, a user should not directly create instances of this class. Instead,
users should create an instance of one of the four subclasses:
LinearQuadElement, QuadraticQuadElement, LinearTriElement, or
QuadraticTriElement.

Elements are made up of nodes, a layer, and a layer plane angle. There are four
types of elements: linear quadrilateral, quadratic quadrilateral, linear
triangular, and quadratic triangular. Each of these element types are defined
as subclasses of the Element baseclass.

Class ("global") attributes:
number_of_elements - The total number of all elements in the cross-section.

Object ("local") attributes:
elem_no - An integer that represents this unique element.
layer - An instance of the Layer class.
theta1 - A float for the layer plane angle (in degrees) of this element.
    Acceptable values are 0 <= theta1 <= 360 degrees. The ply coordinate system
    (y1,y2,y3) is formed by rotating the global coordinate system (x1,x2,x3) in
    the right-hand sense about x1 by the amount 0 <= theta1 <= 360 degrees. For
    example, consider a box-beam cross-section.

                            theta1 = 0
                        ---------------------
                        |    (top wall)     |
            theta1 = 90 |                   | theta1 = 270
            (left wall) |                   | (right wall)
                        |   (bottom wall)   |
                        ---------------------
                            theta1 = 180

    For the top wall, theta1=0; the left wall, theta1=90; the bottom wall,
    theta1=180; the right wall, theta1=270. For more details, see the VABS
    Manual for Users, pages 8-9, figures 4-5.

Object methods:
swap_nodes() - Swaps the assignment of two nodes in this element.

Usage:
Typically, a user should not directly create instances of this class. Instead,
users should create an instance of one of the four subclasses:
LinearQuadElement, QuadraticQuadElement, LinearTriElement, or
QuadraticTriElement.

    """
    number_of_elements = 0
    def __init__(self, layer, theta1):
        Element.number_of_elements += 1
        self.elem_no = Element.number_of_elements
        self.layer = layer
        self.theta1 = float(theta1)

    def swap_nodes(self, nodeA, nodeB):
        temp = nodeA
        nodeA = nodeB
        nodeB = temp
        return (nodeA, nodeB)


class LinearQuadElement(Element):
    """
Make a linear quadrilateral element with four nodes.

The VABS quadrilateral element node numbering scheme is shown below:
    4-------3
    |       |
    |       |
    |       |
    1-------2

Object ("local") attributes:
elem_no - An integer that represents this unique element.
layer - An instance of the Layer class.
theta1 - A float for the layer plane angle (in degrees) of this element.
    Acceptable values are 0 <= theta1 <= 360 degrees. The ply coordinate system
    (y1,y2,y3) is formed by rotating the global coordinate system (x1,x2,x3) in
    the right-hand sense about x1 by the amount 0 <= theta1 <= 360 degrees. For
    example, consider a box-beam cross-section. For the upper wall, theta1=0;
    the left wall, theta1=90; the lower wall, theta1=180; the right wall,
    theta1=270. For more details, see the VABS Manual for Users, pages 8-9,
    figures 4-5.
node1 - An instance of the Node class, which represents the unique node located
    at the bottom left corner of this element.
node2 - An instance of the Node class, which represents the unique node located
    at the bottom right corner of this element.
node3 - An instance of the Node class, which represents the unique node located
    at the top right corner of this element.
node4 - An instance of the Node class, which represents the unique node located
    at the top left corner of this element.
node5 - This node is not present in a LinearQuadElement. It is set to
    EmptyNode.
node6 - This node is not present in a LinearQuadElement. It is set to
    EmptyNode.
node7 - This node is not present in a LinearQuadElement. It is set to
    EmptyNode.
node8 - This node is not present in a LinearQuadElement. It is set to
    EmptyNode.
node9 - This node is not present in a LinearQuadElement. It is set to
    EmptyNode.
x2_middle - A float for the x2-coordinate of the middle of this element.
x3_middle - A float for the x3-coordinate of the middle of this element.

Object methods:
middle() - Returns the 'middle' x2 and x3 coordinates of this element.
_angles() - Returns a dictionary of angles for each node in this element.
is_ccw() - Returns a boolean True if the nodes of this element are arranged
    counter-clockwise; False otherwise.
reorder_nodes() - Reorders the nodes of this element in a counter-clockwise
    fashion from a clockwise fashion, or vice-versa.

Usage:
from vabs_objects import *
a = Node(0,0)
b = Node(1,0)
c = Node(1,1)
d = Node(0,1)
e = Node(2,0)
f = Node(2,1)
matl = IsotropicMaterial(name='Foam', rho=0.20E+03, E=0.256E+09, nu=0.3)
L1 = Layer(matl, 45.3)
elem1 = LinearQuadElement(L1, 270, a, b, c, d)
elem2 = LinearQuadElement(L1, 270, b, e, f, c)
print elem1
print elem2

from vabs_objects import *
a = Node(0,0)
b = Node(1,0)
c = Node(1,1)
d = Node(0,1)
matl = IsotropicMaterial(name='Foam', rho=0.20E+03, E=0.256E+09, nu=0.3)
L1 = Layer(matl, 45.3)
elem1 = LinearQuadElement(L1, 270, a, d, c, b)
elem2 = LinearQuadElement(L1, 270, a, c, b, d)

    """
    def __init__(self, layer, theta1, node1, node2, node3, node4):
        Element.__init__(self, layer, theta1)
        (self.node1, self.node2, self.node3, self.node4) = (
            node1, node2, node3, node4)
        # Force nodes 5-9 to be zero (e.g. self.node5.node_no = 0).
        not_a_node = EmptyNode()
        (self.node5, self.node6, self.node7, self.node8, self.node9) = (
            not_a_node, not_a_node, not_a_node, not_a_node, not_a_node)
        # Calculate the middle coordinates of this element
        self.x2_middle = (self.node1.x2 + self.node2.x2
                          + self.node3.x2 + self.node4.x2)/4.0
        self.x3_middle = (self.node1.x3 + self.node2.x3
                          + self.node3.x3 + self.node4.x3)/4.0
        if self.is_ccw() == False:  # Try to fix a bad (CW) element
            self.reorder_nodes()         # Hopefully this makes a good (CCW) element
            print "Trying to fix this element. Reorienting..."
        # Check if the element actually has CCW orientation
        if self.is_ccw() == False:
            print "    WARNING: Element #{:d} is bad!\n\t\tIts nodes are not oriented CCW.".format(self.elem_no)
        # else:
        #     print "    Yay, element fixed!"
        # Tell nodes 1-4 which element they belong to
        self.node1.element_list.append(self)
        self.node2.element_list.append(self)
        self.node3.element_list.append(self)
        self.node4.element_list.append(self)

    def _inspect_empty_node(self, nodeN):
        a = []
        a.append('  This node is not present in element #{:d} ({:s}):\n'.format(
            self.elem_no, LinearQuadElement.__name__))
        a.append('    ' + nodeN.__str__() + '\n')
        return "".join(a)
    
    def __str__(self):
        L = []
        L.append('INSPECTING ELEMENT #{:d} '.format(self.elem_no)
            + '-'*43 + '\n')
        L.append(self.layer.__str__() + '\n')
        L.append('NODE1********\n')
        L.append(self.node1.__str__() + '\n')
        L.append('NODE2********\n')
        L.append(self.node2.__str__() + '\n')
        L.append('NODE3********\n')
        L.append(self.node3.__str__() + '\n')
        L.append('NODE4********\n')
        L.append(self.node4.__str__() + '\n')
        L.append('NODE5********\n')
        L.append(LinearQuadElement._inspect_empty_node(self, self.node5))
        L.append('NODE6********\n')
        L.append(LinearQuadElement._inspect_empty_node(self, self.node6))
        L.append('NODE7********\n')
        L.append(LinearQuadElement._inspect_empty_node(self, self.node7))
        L.append('NODE8********\n')
        L.append(LinearQuadElement._inspect_empty_node(self, self.node8))
        L.append('NODE9********\n')
        L.append(LinearQuadElement._inspect_empty_node(self, self.node9))
        L.append('-'*65)
        return "".join(L)
    
    def middle(self):
        """
Returns the 'middle' x2 and x3 coordinates of this element.

Usage:
from vabs_objects import *
a = Node(0,0)
b = Node(1,0)
c = Node(1,1)
d = Node(0,1)
matl = IsotropicMaterial(name='Foam', rho=0.20E+03, E=0.256E+09, nu=0.3)
L1 = Layer(matl, 45.3)
elem1 = LinearQuadElement(L1, 270, a, b, c, d)
(x2m,x3m) = elem1.middle()

        """
        return (self.x2_middle, self.x3_middle)

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
            print "  node2: #{:d},    node4: #{:d}".format(self.node2.node_no,
                                                           self.node4.node_no)
        # Reorder the nodes in a counter-clockwise fashion.
        (self.node2,self.node4) = self.swap_nodes(self.node2,self.node4)
        if print_flag:
            print "AFTER:"
            print "  node2: #{:d},    node4: #{:d}".format(self.node2.node_no,
                                                           self.node4.node_no)


class QuadraticQuadElement(Element):
    """
Make a linear quadrilateral element with eight (or nine) nodes.

The VABS quadrilateral element node numbering scheme is shown below:
    4---7---3
    |       |
    8  (9)  6
    |       |
    1---5---2
The central node (9) is optional.

Object ("local") attributes:
elem_no - An integer that represents this unique element.
layer - An instance of the Layer class.
theta1 - A float for the layer plane angle (in degrees) of this element.
    Acceptable values are 0 <= theta1 <= 360 degrees. The ply coordinate system
    (y1,y2,y3) is formed by rotating the global coordinate system (x1,x2,x3) in
    the right-hand sense about x1 by the amount 0 <= theta1 <= 360 degrees. For
    example, consider a box-beam cross-section. For the upper wall, theta1=0;
    the left wall, theta1=90; the lower wall, theta1=180; the right wall,
    theta1=270. For more details, see the VABS Manual for Users, pages 8-9,
    figures 4-5.
node1 - An instance of the Node class, which represents the unique node located
    at the bottom left corner of this element.
node2 - An instance of the Node class, which represents the unique node located
    at the bottom right corner of this element.
node3 - An instance of the Node class, which represents the unique node located
    at the top right corner of this element.
node4 - An instance of the Node class, which represents the unique node located
    at the top left corner of this element.
node5 - An instance of the Node class, which represents the unique node located
    at the midpoint of the bottom side of this element.
node6 - An instance of the Node class, which represents the unique node located
    at the midpoint of the right side of this element.
node7 - An instance of the Node class, which represents the unique node located
    at the midpoint of the top side of this element.
node8 - An instance of the Node class, which represents the unique node located
    at the midpoint of the left side of this element.
node9 - This node is optional. If the user assigns node9, this is an instance
    of the Node class, which represents the unique node located at the center
    of this element. Otherwise, node9 is set to EmptyNode.
x2_middle - A float for the x2-coordinate of the middle of this element.
x3_middle - A float for the x3-coordinate of the middle of this element.

Object methods:
middle() - Returns the 'middle' x2 and x3 coordinates of this element.
_angles() - Returns a dictionary of angles for each node in this element.
is_ccw() - Returns a boolean True if the nodes of this element are arranged
    counter-clockwise; False otherwise.
reorder_nodes() - Reorients the nodes from a clockwise orientation to a
    counter-clockwise orientation, or vice-versa.

Usage:
from vabs_objects import *
a = Node(0,0)
b = Node(1,0)
c = Node(1,1)
d = Node(0,1)
e = Node(2,0)
f = Node(2,1)
g = Node(0,2)
h = Node(1,2)
i = Node(2,2)
matl = IsotropicMaterial(name='Foam', rho=0.20E+03, E=0.256E+09, nu=0.3)
L1 = Layer(matl, 45.3)
elem1 = QuadraticQuadElement(L1, 270, a, e, i, g, b, f, h, d, c)
print elem1
elem1.middle(print_flag=True)
j = Node(10,0)
k = Node(11,0)
m = Node(10,1)
n = Node(12,0)
o = Node(12,1)
p = Node(10,2)
q = Node(11,2)
r = Node(12,2)
elem2 = QuadraticQuadElement(L1, 270, j, n, r, p, k, o, q, m)
print elem2
elem2.middle(print_flag=True)

    """
    not_a_node = EmptyNode()
    def __init__(self, layer, theta1, 
                 node1, node2, node3, node4, 
                 node5, node6, node7, node8, node9=not_a_node):
        Element.__init__(self, layer, theta1)
        # Assign corner nodes.
        (self.node1, self.node2, self.node3, self.node4) = (
            node1, node2, node3, node4)
        # Assign mid-side nodes.
        (self.node5, self.node6, self.node7, self.node8) = (
            node5, node6, node7, node8)
        # Assign central node (optional).
        self.node9 = node9
        # Calculate the middle coordinates of this element
        self.x2_middle = (self.node1.x2 + self.node2.x2
                          + self.node3.x2 + self.node4.x2)/4.0
        self.x3_middle = (self.node1.x3 + self.node2.x3
                          + self.node3.x3 + self.node4.x3)/4.0
        if self.is_ccw() == False:  # Try to fix a bad (CW) element
            self.reorder_nodes()         # Hopefully this makes a good (CCW) element
            print "Trying to fix this element. Reorienting..."
        # Check if the element actually has CCW orientation
        if self.is_ccw() == False:
            print "    WARNING: Element #{:d} is bad!\n\t\tIts nodes are not oriented CCW.".format(self.elem_no)
        # else:
        #     print "    Yay, element fixed!"
        # Tell nodes 1-8 which element they belong to
        self.node1.element_list.append(self)
        self.node2.element_list.append(self)
        self.node3.element_list.append(self)
        self.node4.element_list.append(self)
        self.node5.element_list.append(self)
        self.node6.element_list.append(self)
        self.node7.element_list.append(self)
        self.node8.element_list.append(self)
        if self.node9 != QuadraticQuadElement.not_a_node:
            self.node9.element_list.append(self)
    
    def _inspect_empty_node(self, nodeN):
        a = []
        a.append('  This node is not present in element #{:d} ({:s}):\n'.format(
            self.elem_no, QuadraticQuadElement.__name__))
        a.append('    ' + nodeN.__str__() + '\n')
        return "".join(a)
    
    def __str__(self):
        L = []
        L.append('INSPECTING ELEMENT #{:d} '.format(self.elem_no)
            + '-'*43 + '\n')
        L.append(self.layer.__str__() + '\n')
        L.append('NODE1********\n')
        L.append(self.node1.__str__() + '\n')
        L.append('NODE2********\n')
        L.append(self.node2.__str__() + '\n')
        L.append('NODE3********\n')
        L.append(self.node3.__str__() + '\n')
        L.append('NODE4********\n')
        L.append(self.node4.__str__() + '\n')
        L.append('NODE5********\n')
        L.append(self.node5.__str__() + '\n')
        L.append('NODE6********\n')
        L.append(self.node6.__str__() + '\n')
        L.append('NODE7********\n')
        L.append(self.node7.__str__() + '\n')
        L.append('NODE8********\n')
        L.append(self.node8.__str__() + '\n')
        L.append('NODE9********\n')
        if self.node9 != QuadraticQuadElement.not_a_node:
            L.append(self.node9.__str__() + '\n')
        else:
            L.append(QuadraticQuadElement._inspect_empty_node(self, self.node9))
        L.append('-'*65)
        return "".join(L)

    def middle(self):
        """
Returns the 'middle' x2 and x3 coordinates of this element.

        """
        return (self.x2_middle, self.x3_middle)

    def _angles(self, print_flag=False):
        """
Returns a dictionary of angles for each node in this element.

Each dictionary value represents the angle between two vectors: (1) a vector
drawn from the middle of the element to a node, and (2) a vector along the
horizontal (positive x2-axis).

Note: This method does not deal with the central node, node9.

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
        # find the angle to node5
        y = self.node5.x3 - self.x3_middle
        x = self.node5.x2 - self.x2_middle
        angle_dict['node5'] = atan2(y,x) * (180.0/pi)
        if angle_dict['node5'] < 0.0:
            angle_dict['node5'] += + 360.0
        # find the angle to node6
        y = self.node6.x3 - self.x3_middle
        x = self.node6.x2 - self.x2_middle
        angle_dict['node6'] = atan2(y,x) * (180.0/pi)
        if angle_dict['node6'] < 0.0:
            angle_dict['node6'] += + 360.0
        # find the angle to node7
        y = self.node7.x3 - self.x3_middle
        x = self.node7.x2 - self.x2_middle
        angle_dict['node7'] = atan2(y,x) * (180.0/pi)
        if angle_dict['node7'] < 0.0:
            angle_dict['node7'] += + 360.0
        # find the angle to node8
        y = self.node8.x3 - self.x3_middle
        x = self.node8.x2 - self.x2_middle
        angle_dict['node8'] = atan2(y,x) * (180.0/pi)
        if angle_dict['node8'] < 0.0:
            angle_dict['node8'] += + 360.0
        if print_flag:
            print "node1_angle = {:6.2f} degrees".format(angle_dict['node1'])
            print "node5_angle = {:6.2f} degrees".format(angle_dict['node5'])
            print "node2_angle = {:6.2f} degrees".format(angle_dict['node2'])
            print "node6_angle = {:6.2f} degrees".format(angle_dict['node6'])
            print "node3_angle = {:6.2f} degrees".format(angle_dict['node3'])
            print "node7_angle = {:6.2f} degrees".format(angle_dict['node7'])
            print "node4_angle = {:6.2f} degrees".format(angle_dict['node4'])
            print "node8_angle = {:6.2f} degrees".format(angle_dict['node8'])
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
e = Node(2,0)
f = Node(2,1)
g = Node(0,2)
h = Node(1,2)
i = Node(2,2)
matl = IsotropicMaterial(name='Foam', rho=0.20E+03, E=0.256E+09, nu=0.3)
L1 = Layer(matl, 45.3)
elem1 = QuadraticQuadElement(L1, 270, a, e, i, g, b, f, h, d, c)
elem1.is_ccw()

        """
        angle_dict = self._angles()         # Update the angles to each node.
        angleDiff = []
        angleDiff.append(angle_dict['node5'] - angle_dict['node1'])
        angleDiff.append(angle_dict['node2'] - angle_dict['node5'])
        angleDiff.append(angle_dict['node6'] - angle_dict['node2'])
        angleDiff.append(angle_dict['node3'] - angle_dict['node6'])
        angleDiff.append(angle_dict['node7'] - angle_dict['node3'])
        angleDiff.append(angle_dict['node4'] - angle_dict['node7'])
        angleDiff.append(angle_dict['node8'] - angle_dict['node4'])
        angleDiff.append(angle_dict['node1'] - angle_dict['node8'])
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
e = Node(2,0)
f = Node(2,1)
g = Node(0,2)
h = Node(1,2)
i = Node(2,2)
matl = IsotropicMaterial(name='Foam', rho=0.20E+03, E=0.256E+09, nu=0.3)
L1 = Layer(matl, 45.3)
elem1 = QuadraticQuadElement(L1, 270, a, g, i, e, d, h, f, b, c)
elem1.is_ccw(print_flag=True)
elem1.reorder_nodes(print_flag=True)
elem1.is_ccw(print_flag=True)
elem1.reorder_nodes(print_flag=True)
elem1.is_ccw(print_flag=True)

        """
        if print_flag:
            print "BEFORE:"
            print "  node2: #{:d},    node4: #{:d}".format(self.node2.node_no,
                                                           self.node4.node_no)
            print "  node8: #{:d},    node5: #{:d}".format(self.node8.node_no,
                                                           self.node5.node_no)
            print "  node6: #{:d},    node7: #{:d}".format(self.node6.node_no,
                                                           self.node7.node_no)
        # Reorder the nodes in a counter-clockwise fashion.
        (self.node2,self.node4) = self.swap_nodes(self.node2,self.node4)
        (self.node8,self.node5) = self.swap_nodes(self.node8,self.node5)
        (self.node6,self.node7) = self.swap_nodes(self.node6,self.node7)
        if print_flag:
            print "AFTER:"
            print "  node2: #{:d},    node4: #{:d}".format(self.node2.node_no,
                                                           self.node4.node_no)
            print "  node8: #{:d},    node5: #{:d}".format(self.node8.node_no,
                                                           self.node5.node_no)
            print "  node6: #{:d},    node7: #{:d}".format(self.node6.node_no,
                                                           self.node7.node_no)