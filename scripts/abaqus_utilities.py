"""
A module to parse data from an ABAQUS-formatted 2D cross-section grid file.

Authors: Phil Chiu
         Perry Roth-Johnson
Last updated: November 2, 2012

"""

import re                        # import the regular expression module
import os                        # import the operating system module
import numpy as np               # import the the numpy module; rename it as np

class AbaqusGrid:
    """
The AbaqusGrid class contains methods for parsing an ABAQUS-formatted
2D grid file (cross-section grid).

Usage:
from abaqus_utilities import *
g = AbaqusGrid('../run_01/input/spar_station_04_linear_abq.txt')
g.node_array
g.node_array['node_no']
g.node_array['x2']
g.node_array['x3']
g.element_array
g.elementset_array
g.number_of_nodes
g.number_of_elements

Initialization:
AbaqusGrid(filename, debug_flag=False)
  filename - A string for the full path of the ABAQUS-formatted grid file.
  debug_flag - Optional boolean to print intermediate results to the screen.

Public attributes:
filename - A string for the full path of the ABAQUS-formatted grid file.
node_array - A structured array that stores nodes in 3 columns:
    node_no: A unique integer that labels this node.
    x2: A float for the x2-coordinate of this node.
    x3: A float for the x3-coordinate of this node.
element_array - A structured array that stores element connectivity.
    If the element is linear (4-noded), 6 columns are stored:
        layer_no: An integer for the unique layer associated with this element.
        elem_no: An integer that represents this unique element.
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
    If the element is quadratic (8-noded), 10 columns are stored:
        layer_no: An integer for the unique layer associated with this element.
        elem_no: An integer that represents this unique element.
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
        The VABS quadrilateral element node numbering scheme is shown below:
elementset_array - A structured array that stores element orientation angles in
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

    elem_no: An integer that represents a unique element.
number_of_nodes - An integer for the number of nodes in the grid.
number_of_elements - An integer for the number of elements in the grid.

    """


    def __init__(self, filename, debug_flag=False):
        self.filename = filename
        self._parse_abaqus(debug_flag=debug_flag)


    def _read_file(self):
        """
Saves the ABAQUS file as a list of strings (as an attribute).

Each element in the list represents one line in the file.

Saves:
self._abq_file - A list of strings.

        """
        f = open(self.filename, 'r')
        self._abq_file = f.readlines()
        f.close()


    def _define_patterns(self):
        """
Define regular expressions to search for nodes and elements in the ABAQUS file.

Saves:
self._node_pattern
self._element_pattern
self._node_header_pattern
self._element_header_pattern
self._elementset_header_pattern

        """
        self._node_pattern = re.compile(
            r'[0-9]+(,-*[0-9]+\.[0-9]*E*[+-]*[0-9]*){3}')
        #  this regex pattern explained:
        #  -----------------------------
        #  [0-9]+                               :  node number
        #  (,-*[0-9]+\.[0-9]*E*[+-]*[0-9]*){3}  :  x-coord, y-coord, z-coord
        #  note: coords may be in decimal or scientific notation

        self._node_header_pattern = re.compile(r'\*NODE.+')

        # linear element connectivity pattern:
        self._linear_element_pattern = re.compile(r'^[0-9]+(,[0-9]+){4}$')
        #  this regex pattern explained:
        #  -----------------------------
        #  ^               :  beginning of line
        #  [0-9]+          :  element number
        #  (,[0-9]+){4}    :  node1-node4
        #  $               :  end of line

        # quadratic element connectivity pattern:
        self._quadratic_element_pattern = re.compile(r'^[0-9]+(,[0-9]+){8}$')
        #  this regex pattern explained:
        #  -----------------------------
        #  [0-9]+          :  element number
        #  (,[0-9]+){8}    :  node1-node8

        self._element_header_pattern = re.compile(r'\*MATERIAL,NAME=M[0-9]+')

        # element set pattern:
        # esetPat = re.compile(r'([0-9]+,){5,16}')

        self._elementset_header_pattern = re.compile(r'\*ELSET,ELSET=.+')
        

    def _find_block_starts(self):
        """
Returns the indices (line number) for the start of the node, element
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
        self._node_block_start = node_headers[0]
        self._element_block_start = element_headers[0]
        self._elementset_block_start = element_set_headers[0]


    def _parse_nodes(self):
        """
Save the nodes in a structured array (as an attribute).

Saves:
self.number_of_nodes
self.node_array

        """
        list_of_node_tuples = []
        for line in self._abq_file[self._node_block_start:self._element_block_start]:
            node_match = self._node_pattern.match(line)
            if node_match:  # if we find a node
                # save the first 3 entries; drop x1 (last entry)
                (node_no, x2, x3) = line.strip().split(',')[:-1]
                list_of_node_tuples.append((int(node_no), float(x2), float(x3)))
        self.number_of_nodes = len(list_of_node_tuples)
        # initialize a structured array
        self.node_array = np.zeros((self.number_of_nodes,), 
            dtype=[('node_no', 'i4'), ('x2', 'f8'), ('x3', 'f8')])
        # save the nodes in the structured array
        self.node_array[:] = list_of_node_tuples


    def _parse_elements(self, debug_flag=False):
        """
Saves the elements in a structured array (as an attribute).

Note: This function only supports linear (4-noded) elements and quadratic
(8-noded) elements from TrueGrid.

Saves:
self.number_of_elements
self.element_array

        """
        list_of_element_tuples = []
        new_element_header_found = False
        for i, line in enumerate(
                self._abq_file[self._element_block_start:self._elementset_block_start]):
            linear_element_match = self._linear_element_pattern.match(line)
            quadratic_element_match = self._quadratic_element_pattern.match(line)
            element_header_match = self._element_header_pattern.match(line)
            if element_header_match:
                new_element_header_found = True
                # Extract the layer number from the element header line.
                layer_no = line.strip().split("=")[-1].split("M")[-1]
                if debug_flag:
                    print 'element header found at line ' + str(
                            i+self._element_block_start+1)
                    print 'layer #' + str(layer_no)
            elif quadratic_element_match:
                if debug_flag and new_element_header_found:
                    print 'first element: #' + line.strip().split(',')[0]
                    new_element_header_found = False
                (elem_no, node1, node2, node3, node4,
                        node5, node6, node7, node8) = line.strip().split(',')
                list_of_element_tuples.append(
                        (int(layer_no), int(elem_no), int(node1), int(node2),
                        int(node3), int(node4), int(node5), int(node6),
                        int(node7), int(node8)))
            elif linear_element_match:
                if debug_flag and new_element_header_found:
                    print 'first element: #' + line.strip().split(',')[0]
                    new_element_header_found = False
                (elem_no, node1, node2, node3, node4) = line.strip().split(',')
                list_of_element_tuples.append(
                        (int(layer_no), int(elem_no), int(node1), int(node2),
                        int(node3), int(node4)))
        self.number_of_elements = len(list_of_element_tuples)
        if len(list_of_element_tuples[0]) == 10:
            is_quadratic = True
        else:
            is_quadratic = False
        # initialize a structured array
        if is_quadratic:
            self.element_array = np.zeros(
                    (self.number_of_elements,), dtype=[('layer_no', 'i4'),
                    ('elem_no', 'i4'), ('node1', 'i4'), ('node2', 'i4'),
                    ('node3', 'i4'), ('node4', 'i4'), ('node5', 'i4'),
                    ('node6', 'i4'), ('node7', 'i4'), ('node8', 'i4')])
        else:
            self.element_array = np.zeros(
                    (self.number_of_elements,), dtype=[('layer_no', 'i4'),
                    ('elem_no', 'i4'), ('node1', 'i4'), ('node2', 'i4'),
                    ('node3', 'i4'), ('node4', 'i4')])
        # save the nodes in the structured array
        self.element_array[:] = list_of_element_tuples


    def _parse_elementsets(self, debug_flag=False):
        """
Save all the elements and their orientations (theta1) in a structured array
(as an attribute).

Saves:
self.elementset_array

FUTURE WORK: theta1_dict should be pulled out and replaced with an input file
supplied by the user...maybe add theta1_dict as a keyword argument?

KNOWN BUG: The elementsets in the ABAQUS file may be missing some elements!
seems to be more likely to happen when quadratic elements are generated with
TrueGrid. When linear elements are generated, everything is usually okay. A
quick and dirty fix was implemented in the past: search for elements in the
mesh that have not yet been assigned a theta1 value. Based on the middle
coordinates of that element, assign a theta1 value.

        """
        list_of_elementset_tuples = []
        # Define the relationship between elementset names and theta1 values
        # e.g., elementset 'scb' (bottom spar cap) has theta1=180 degrees
        theta1_dict = { 'swlbiaxl': 90,      # left wall
                        'swlfoam':  90,      # left wall
                        'swlbiaxr': 90,      # left wall
                        'swrbiaxl': 270,     # right wall
                        'swrfoam':  270,     # right wall
                        'swrbiaxr': 270,     # right wall
                        'sct':      0,       # top wall
                        'rbt':      0,       # top wall
                        'scb':      180,     # bottom wall
                        'rbb':      180 }    # bottom wall

        for line in self._abq_file[self._elementset_block_start:]:
            elementset_header_match = self._elementset_header_pattern.match(line)
            if elementset_header_match:
                # Extract the elementset name and look up its theta1 value.
                theta1 = theta1_dict[line.strip().split('=')[-1]]
                if debug_flag:
                    print 'element set: ' + line.strip().split('=')[-1]
                    print 'theta1 = ' + str(theta1) + '\n'
            else:
                elements = line.strip().strip(',').split(',')
                for elem_no in elements:
                    list_of_elementset_tuples.append((theta1,int(elem_no)))
        elementset_len = len(list_of_elementset_tuples)
        # initialize a structured array
        self.elementset_array = np.zeros((self.number_of_elements,),
            dtype=[('theta1', 'i4'), ('elem_no', 'i4')])
        # save the nodes in the structured array
        try:
            self.elementset_array[:] = list_of_elementset_tuples
        except:
            print """
<*WARNING* in abaqus_utilities.AbaqusGrid._parse_elementsets()>
    self.number_of_elements != len(list_of_elementset_tuples)
    Elements are missing from the elementset blocks in the ABAQUS file!
    ...Saving as many elements as possible in self.elementset_array
            """
            # Reinitialize the structured array with a smaller shape in order
            # to match the number of elements parsed from the elemensets.
            self.elementset_array = np.zeros((elementset_len,),
                dtype=[('theta1', 'i4'), ('elem_no', 'i4')])
            self.elementset_array[:] = list_of_elementset_tuples


    def _parse_abaqus(self, debug_flag=False):
        """
Parses the ABAQUS output file and saves the node, element, and eset structured
arrays as object attributes.

This non-public method is automatically run when a new AbaqusGrid instance is
created.

Saves:
self.number_of_nodes
self.node_array

Usage:
cd safe
import scripts.abaqus_utilities as au
g = au.AbaqusGrid('run_01/input/spar_station_04_linear_abq.txt')
g._parse_abaqus(debug_flag=True)

g.node_array
g.element_array
g.elementset_array

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
        self._parse_elementsets(debug_flag=debug_flag)
        # Sort element_array by element number.
        self.element_array.sort(order='elem_no')
        # Sort elementset_array by element number.
        self.elementset_array.sort(order='elem_no')
        if debug_flag:
            print '.node_array[0] =', self.node_array[0]
            print '.element_array[0] =', self.element_array[0]
            print '.elementset_array[0] =', self.elementset_array[0]
            print 'number of nodes: ' + str(self.number_of_nodes)
            print 'number of elements: ' + str(self.number_of_elements)