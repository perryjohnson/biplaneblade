"""Create entities for a 2D unstructured grid.

Author: Perry Roth-Johnson
Last modified: March 13, 2014

"""

class Node:
    def __init__(self, node_no, x2, x3):
        self.node_no = node_no
        self.x2 = x2
        self.x3 = x3
        self.xy = (self.x2,self.x3)
        self.parent_element = None

class Element:
    def __init__(self, elem_no, node1, node2, node3, node4):
        self.elem_no = elem_no
        self.node1 = node1
        self.node2 = node2
        self.node3 = node3
        self.node4 = node4
        self.theta1 = None
        self.node1.parent_element = self
        self.node2.parent_element = self
        self.node3.parent_element = self
        self.node4.parent_element = self
