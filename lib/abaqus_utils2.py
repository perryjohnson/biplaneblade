"""A module to parse data from an ABAQUS-formatted 2D cross-section grid file.

Authors: Perry Roth-Johnson, Phil Chiu
Last updated: March 27, 2014

"""


import re
import os
import numpy as np
import grid as gr
reload(gr)
from operator import attrgetter


class AbaqusGrid:
    """The AbaqusGrid class contains methods for parsing an ABAQUS-formatted
    2D grid file (cross-section grid).

    Usage:
    import lib.abaqus_utils2 as au
    g = au.AbaqusGrid('cs_abq.txt')
    g.number_of_nodes
    g.list_of_nodes
    g.list_of_nodes[0].node_num
    g.list_of_nodes[0].x2
    g.list_of_nodes[0].x3
    g.number_of_elements
    g.list_of_elements

    Initialization:
    AbaqusGrid(filename, debug_flag=False)
      filename - A string for the full path of the ABAQUS-formatted grid file.
      debug_flag - Optional boolean to print intermediate results to the screen.

    Public attributes:
    filename - A string for the full path of the ABAQUS-formatted grid file.
    list_of_nodes - A list of gr.Node objects, with attributes:
        node_num: A unique integer that labels this node.
        x2: A float for the x2-coordinate of this node.
        x3: A float for the x3-coordinate of this node.
    list_of_elements - A list of gr._Element objects that stores element
        connectivity.
        If the element is linear (4-noded), the list is populated with
            gr.QuadrilateralLinearElement objects, with attributes:
            layer_no: An integer for the unique layer associated with this element.
            elem_num: An integer that represents this unique element.
                                    4-------3
                                    |       |   The VABS linear quadrilateral
                                    |       |   element node numbering scheme.
                                    |       |
                                    1-------2
            node1: An integer for the unique node located at the bottom left corner
                of this element.
            node2: An integer for the unique node located at the bottom right
                corner of this element.
            node3: An integer for the unique node located at the top right corner
                of this element.
            node4: An integer for the unique node located at the top left corner of
                this element.
        If the element is quadratic (8-noded), the list is populated with
            gr.QuadrilateralQuadraticElement objects, with attributes:
            layer_no: An integer for the unique layer associated with this element.
            elem_num: An integer that represents this unique element.
                                    4---7---3
                                    |       |   The VABS quadratic quadrilateral
                                    8       6   element node numbering scheme.
                                    |       |
                                    1---5---2
            node1: An integer for the unique node located at the bottom left corner
                of this element.
            node2: An integer for the unique node located at the bottom right
                corner of this element.
            node3: An integer for the unique node located at the top right corner
                of this element.
            node4: An integer for the unique node located at the top left corner of
                this element.
            node5: An integer for the unique node located at the midpoint of the
                bottom side of this element.
            node6: An integer for the unique node located at the midpoint of the
                right side of this element.
            node7: An integer for the unique node located at the midpoint of the
                top side of this element.
            node8: An integer for the unique node located at the midpoint of the
                left side of this element.
    list_of_element_sets - A structured array that stores element orientation angles in
        2 columns:
        theta1: The orientation (in degrees) of this element. Sometimes this is
            also referred to as the "layer plane angle."
                    Consider the box cross-section below:

                                theta1 = 0
                            ---------------------
                            |    (top wall)     |
                theta1 = 90 |                   | theta1 = 270
                (left wall) |                   | (right wall)
                            |   (bottom wall)   |
                            ---------------------
                                theta1 = 180

        elem_num: An integer that represents a unique element.
    number_of_nodes - An integer for the number of nodes in the grid.
    number_of_elements - An integer for the number of elements in the grid.

    """
    def __init__(self, filename, debug_flag=False, soft_warning=False,
        auto_parse=True):
        self.filename = filename
        # attributes for self._read_file()
        self._abq_file = None
        # attributes for self._define_patterns()
        self._node_pattern = None
        self._node_header_pattern = None
        self._quadrilateral_linear_element_pattern = None
        self._quadrilateral_quadratic_element_pattern = None
        self._triangular_quadratic_element_pattern = None
        self._element_header_pattern = None
        self._elementset_header_pattern = None
        # attributes for self._find_block_starts()
        self._node_block_start = None
        self._element_block_start = None
        self._elementset_block_start = None
        # attributes for self._parse_nodes()
        self.number_of_nodes = 0
        self.list_of_nodes = []
        # attributes for self._parse_elements()
        self.number_of_elements = None
        self.list_of_elements = []
        if auto_parse:
            # parse the ABAQUS output file into grid objects
            self._parse_abaqus(debug_flag=debug_flag,
                soft_warning=soft_warning)

    def _parse_abaqus(self, debug_flag=False, soft_warning=False):
        """Parses the ABAQUS output file and saves it grid objects.

        This non-public method is automatically run when a new AbaqusGrid
        instance is created.

        """
        if debug_flag:
            print 'ABAQUS file: ' + self.filename
        self._read_file()
        if debug_flag:
            print 'STATUS: parsing the ABAQUS file...'
        self._define_patterns()
        self._find_block_starts()
        self._parse_nodes()
        self._parse_elements(debug_flag=debug_flag)
        # Sort list_of_elements by element number.
        #   This MUST happen before calling self._parse_elementsets()
        self.list_of_elements.sort(key=attrgetter('elem_num'))
        self._parse_elementsets(debug_flag=debug_flag, 
            soft_warning=soft_warning)
        if debug_flag:
            print 'list_of_nodes[0] =', self.list_of_nodes[0]
            print 'list_of_elements[0] =', self.list_of_elements[0]
            print 'number of nodes: ' + str(self.number_of_nodes)
            print 'number of elements: ' + str(self.number_of_elements)    

    def _read_file(self):
        """Saves the ABAQUS file as a list of strings (as an attribute).

        Each element in the list represents one line in the file.

        Saves:
        self._abq_file - A list of strings.

        """
        f = open(self.filename, 'r')
        self._abq_file = f.readlines()
        f.close()

    def _define_patterns(self):
        """Define regular expressions for nodes and elements in an ABAQUS file.

        Saves:
        self._node_pattern
        self._node_header_pattern
        self._quadrilateral_linear_element_pattern
        self._quadrilateral_quadratic_element_pattern
        self._triangular_quadratic_element_pattern
        self._element_header_pattern
        self._elementset_header_pattern

        """
        # node pattern --------------------------------------------------------
        self._node_pattern = re.compile(
            r'[0-9]+(,-*[0-9]+\.[0-9]*E*[+-]*[0-9]*){3}')
        # this regex pattern explained:
        # -----------------------------
        # [0-9]+ : node number
        # (,-*[0-9]+\.[0-9]*E*[+-]*[0-9]*){3} : x-coord, y-coord, z-coord
        # note: coords may be in decimal or scientific notation
        # 
        # node header pattern -------------------------------------------------
        self._node_header_pattern = re.compile(r'\*NODE.+')
        #
        # quadrilateral linear element connectivity pattern -------------------
        self._quadrilateral_linear_element_pattern = re.compile(
            r'^[0-9]+(,[0-9]+){4}$')
        # this regex pattern explained:
        # -----------------------------
        # ^ : beginning of line
        # [0-9]+ : element number
        # (,[0-9]+){4} : node1-node4
        # $ : end of line
        #
        # quadrilateral quadratic element connectivity pattern ----------------
        self._quadrilateral_quadratic_element_pattern = re.compile(
            r'^[0-9]+(,[0-9]+){8}$')
        # this regex pattern explained:
        # -----------------------------
        # [0-9]+ : element number
        # (,[0-9]+){8} : node1-node8
        #
        # triangular linear element connectivity pattern ----------------------
        self._triangular_linear_element_pattern = re.compile(
            r'^[0-9]+(,[0-9]+){3}$')
        # this regex pattern explained:
        # -----------------------------
        # ^ : beginning of line
        # [0-9]+ : element number
        # (,[0-9]+){3} : node1-node3
        # $ : end of line
        #
        # triangular quadratic element connectivity pattern -------------------
        self._triangular_quadratic_element_pattern = re.compile(
            r'^[0-9]+(,[0-9]+){6}$')
        # this regex pattern explained:
        # -----------------------------
        # [0-9]+ : element number
        # (,[0-9]+){6} : node1-node6
        #
        # element header pattern ----------------------------------------------
        self._element_header_pattern = re.compile(r'\*ELEMENT.+')
        #
        # element set header pattern ------------------------------------------
        self._elementset_header_pattern = re.compile(r'\*ELSET,ELSET=.+')
        
    def _find_block_starts(self):
        """Returns the indices (line number) for the start of the node, element
        connectivity, and element set blocks in self._abq_file.

        Saves:
        self._node_block_start
        self._element_block_start
        self._elementset_block_start

        """
        node_headers = []
        element_headers = []
        element_set_headers = []
        for i, line in enumerate(self._abq_file):
            node_header_match = self._node_header_pattern.match(line)
            element_header_match = self._element_header_pattern.match(line)
            elementset_header_match = self._elementset_header_pattern.match(line)
            if node_header_match:
                node_headers.append(i)
            elif element_header_match:
                element_headers.append(i)
            elif elementset_header_match:
                element_set_headers.append(i)
        if len(node_headers) > 0:
            self._node_block_start = node_headers[0]
        else:
            self._node_block_start = 0
            print "WARNING: node block start not found!"
        if len(element_headers) > 0:
            self._element_block_start = element_headers[0]
        else:
            self._element_block_start = 0
            print "WARNING: element block start not found!"
        if len(element_set_headers) > 0:
            self._elementset_block_start = element_set_headers[0]
        else:
            self._elementset_block_start = 0
            print "WARNING: elementset block start not found!"

    def _parse_nodes(self):
        """Save the nodes in a list of gr.Node objects.

        Saves:
        self.number_of_nodes
        self.list_of_nodes

        """
        for line in self._abq_file[self._node_block_start:self._element_block_start]:
            node_match = self._node_pattern.match(line)
            if node_match: # if we find a node
                # save the first 3 entries; drop x1 (last entry)
                (node_num, x2, x3) = line.strip().split(',')[:-1]
                n = gr.Node(node_num, x2, x3)
                self.list_of_nodes.append(n)
        self.number_of_nodes = len(self.list_of_nodes)

    def _parse_elements(self, debug_flag=False):
        """Saves the elements in a list of gr.Element objects.

        This function supports 4 element types from TrueGrid:
        (1) quadrilateral linear (4-noded) elements
        (2) quadrilateral quadratic (8-noded) elements
        (3) triangular linear (3-noded) elements
        (4) triangular quadratic (6-noded) elements

        Saves:
        self.number_of_elements
        self.list_of_elements

        """
        new_element_header_found = False
        for i, line in enumerate(self._abq_file[self._element_block_start:self._elementset_block_start]):
            quadrilateral_linear_element_match = self._quadrilateral_linear_element_pattern.match(line)
            quadrilateral_quadratic_element_match = self._quadrilateral_quadratic_element_pattern.match(line)
            triangular_linear_element_match = self._triangular_linear_element_pattern.match(line)
            triangular_quadratic_element_match = self._triangular_quadratic_element_pattern.match(line)
            element_header_match = self._element_header_pattern.match(line)
            if element_header_match:
                new_element_header_found = True
                # Extract the layer number from the element header line.
                layer_num = line.strip().split("=")[-1].split("M")[-1]
                if debug_flag:
                    print 'element header found at line ' + str(
                            i+self._element_block_start+1)
                    print 'layer #' + str(layer_num)
            elif quadrilateral_quadratic_element_match:
                if debug_flag and new_element_header_found:
                    print 'first element: #' + line.strip().split(',')[0]
                    new_element_header_found = False
                (elem_num, node1_num, node2_num, node3_num, node4_num,
                    node5_num, node6_num, node7_num,
                    node8_num) = line.strip().split(',')
                e = gr.QuadrilateralQuadraticElement(
                    elem_num = int(elem_num),
                    node1 = self.list_of_nodes[int(node1_num)-1],
                    node2 = self.list_of_nodes[int(node2_num)-1],
                    node3 = self.list_of_nodes[int(node3_num)-1],
                    node4 = self.list_of_nodes[int(node4_num)-1],
                    node5 = self.list_of_nodes[int(node5_num)-1],
                    node6 = self.list_of_nodes[int(node6_num)-1],
                    node7 = self.list_of_nodes[int(node7_num)-1],
                    node8 = self.list_of_nodes[int(node8_num)-1],
                    layer_num = int(layer_num))
                self.list_of_elements.append(e)
            elif quadrilateral_linear_element_match:
                if debug_flag and new_element_header_found:
                    print 'first element: #' + line.strip().split(',')[0]
                    new_element_header_found = False
                (elem_num, node1_num, node2_num, node3_num,
                    node4_num) = line.strip().split(',')
                e = gr.QuadrilateralLinearElement(
                    elem_num = int(elem_num),
                    node1 = self.list_of_nodes[int(node1_num)-1],
                    node2 = self.list_of_nodes[int(node2_num)-1],
                    node3 = self.list_of_nodes[int(node3_num)-1],
                    node4 = self.list_of_nodes[int(node4_num)-1],
                    layer_num = int(layer_num))
                self.list_of_elements.append(e)
            elif triangular_quadratic_element_match:
                if debug_flag and new_element_header_found:
                    print 'first element: #' + line.strip().split(',')[0]
                    new_element_header_found = False
                (elem_num, node1_num, node2_num, node3_num, node5_num,
                    node6_num, node7_num) = line.strip().split(',')
                e = gr.TriangularQuadraticElement(
                    elem_num = int(elem_num),
                    node1 = self.list_of_nodes[int(node1_num)-1],
                    node2 = self.list_of_nodes[int(node2_num)-1],
                    node3 = self.list_of_nodes[int(node3_num)-1],
                    node5 = self.list_of_nodes[int(node5_num)-1],
                    node6 = self.list_of_nodes[int(node6_num)-1],
                    node7 = self.list_of_nodes[int(node7_num)-1],
                    layer_num = int(layer_num))
                self.list_of_elements.append(e)
            elif triangular_linear_element_match:
                if debug_flag and new_element_header_found:
                    print 'first element: #' + line.strip().split(',')[0]
                    new_element_header_found = False
                (elem_num, node1_num, node2_num,
                    node3_num) = line.strip().split(',')
                e = gr.TriangularLinearElement(
                    elem_num = int(elem_num),
                    node1 = self.list_of_nodes[int(node1_num)-1],
                    node2 = self.list_of_nodes[int(node2_num)-1],
                    node3 = self.list_of_nodes[int(node3_num)-1],
                    layer_num = int(layer_num))
                self.list_of_elements.append(e)
        self.number_of_elements = len(self.list_of_elements)

    def _parse_elementsets(self, debug_flag=False, soft_warning=False):
        """Save all the element sets as attributes of their Elements.

        """
        for line in self._abq_file[self._elementset_block_start:]:
            elementset_header_match = self._elementset_header_pattern.match(line)
            if elementset_header_match:
                # Extract the elementset name
                elementset_name = line.strip().split('=')[-1]
                if debug_flag:
                    print 'element set: ' + line.strip().split('=')[-1]
            else:
                # Save the elementset name to each element in the block
                element_nums = line.strip().strip(',').split(',')
                for elem_num in element_nums:
                    # make sure the list of elements have been sorted
                    #   before assigning element sets to elements
                    if int(elem_num) != self.list_of_elements[int(elem_num)-1].elem_num:
                        if not soft_warning:
                            raise Warning("The element set '{0}' may be assigned to the wrong element (#{1}), instead of to the correct element (#{2}). In <grid>._parse_abaqus(), run:\n-->  <grid>.list_of_elements.sort(key=attrgetter('elem_num'))\nbefore calling:\n-->  <grid>._parse_elementsets(debug_flag=debug_flag)".format(elementset_name, self.list_of_elements[int(elem_num)-1].elem_num, int(elem_num)))
                        else:
                            print "The element set '{0}' may be assigned to the wrong element (#{1}), instead of to the correct element (#{2}). In <grid>._parse_abaqus(), run:\n-->  <grid>.list_of_elements.sort(key=attrgetter('elem_num'))\nbefore calling:\n-->  <grid>._parse_elementsets(debug_flag=debug_flag)".format(elementset_name, self.list_of_elements[int(elem_num)-1].elem_num, int(elem_num))
                    self.list_of_elements[int(elem_num)-1].element_set = elementset_name
