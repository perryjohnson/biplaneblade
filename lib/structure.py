"""A module for organizing structural part data for a blade station.

Note: to generate html documentation of this module, open a Windows cmd prompt:
> cd path/to/lib
> python -m pydoc -w structure
Then, open path/to/lib/structure.html in a browser.
Refs:
https://docs.python.org/2.7/library/pydoc.html
http://bytes.com/topic/python/answers/436285-how-use-pydoc

The different kinds of parts, listed from outside to inside, are:
    external surface (gelcoat, triax)
    root buildup (triax)
    LE panel (foam)
    spar caps (uniax)
    aft panel 1 (foam)
    aft panel 2 (foam)
    TE reinforcement (uniax, foam)
    shear web 1 (biax, foam, biax)
    shear web 2 (biax, foam, biax)
    shear web 3 (biax, foam, biax)
    internal surface 1 (triax, resin)
    internal surface 2 (triax, resin)
    internal surface 3 (triax, resin)
    internal surface 4 (triax, resin)

Author: Perry Roth-Johnson
Last updated: April 14, 2014

"""


import os
import numpy as np
import matplotlib.pyplot as plt
import layer as l
reload(l)
from math import isnan
from shapely.geometry import Polygon, asLineString
from shapely.ops import cascaded_union
from shapely.topology import TopologicalError
from descartes import PolygonPatch
# the descartes module translates shapely objects into matplotlib objects
from operator import attrgetter
# helps to sort lists of objects by their attributes
# ref: https://wiki.python.org/moin/HowTo/Sorting#Operator_Module_Functions


class Part:
    """Define the dimensions of a structural part."""
    def __init__(self, parent_structure, base, height):
        self.parent_structure = parent_structure
        self.base = base
        self.height = height
        self.left = None    # assigned later by <station>.find_part_edges()
        self.right = None   # assigned later by <station>.find_part_edges()
        self.layer = {}     # assigned later by <part>.create_layers()
        self.alt_layer = {} # assigned later by <part>.add_new_layer(...)
    
    def __str__(self):
        return """base:    {0} (meters)
height:  {1} (meters)""".format(self.base, self.height)
    
    def exists(self):
        """Checks if a structural part exists at this station."""
        if isnan(self.base) and isnan(self.height):
            return False
        else:
            return True
    
    def bounding_box(self, x_boundary_buffer=1.2, y_boundary_buffer=1.2,
        airfoil=None):
        """Returns a polygon of the bounding box that contains this part.

        The points of the bounding box are labeled from 1 to 4 as:

        4---3
        |   |
        1---2

        Parameters
        ----------
        x_boundary_buffer : float (default: 1.2), factor to multiply with the
            minx and maxx bound of the airfoil polygon, to stretch the bounding
            box past the left and right edges of the airfoil polygon
        y_boundary_buffer : float (default: 1.2), factor to multiply with the
            miny and maxy bound of the airfoil polygon, to stretch the bounding
            box above and below the top and bottom edges of the airfoil polygon

        """
        af = self.parent_structure.parent_station.airfoil
        if airfoil is None:
            (minx, miny, maxx, maxy) = af.polygon.bounds
        elif airfoil == 'lower':
            (minx, miny, maxx, maxy) = af.lower_polygon.bounds
        elif airfoil == 'upper':
            (minx, miny, maxx, maxy) = af.upper_polygon.bounds
        else:
            raise ValueError("Keyword `airfoil` must be None, 'lower', or 'upper'.")
        if self.left is not None and self.right is not None:
            # if the part has values for the attributes `left` and `right`
            pt1 = (self.left, miny*y_boundary_buffer)
            pt2 = (self.right, miny*y_boundary_buffer)
            pt3 = (self.right, maxy*y_boundary_buffer)
            pt4 = (self.left, maxy*y_boundary_buffer)
        else:
            # if self.left=None and self.right=None, then just use
            # the airfoil bounds to form the bounding box
            # (e.g. the root buildup doesnt have left and right edges)
            pt1 = (minx*x_boundary_buffer, miny*y_boundary_buffer)
            pt2 = (maxx*x_boundary_buffer, miny*y_boundary_buffer)
            pt3 = (maxx*x_boundary_buffer, maxy*y_boundary_buffer)
            pt4 = (minx*x_boundary_buffer, maxy*y_boundary_buffer)
        bounding_box = Polygon([pt1, pt2, pt3, pt4])
        return bounding_box
    

class ExternalSurface(Part):
    """Define triax and gelcoat dimensions of the external surface."""
    def __init__(self, parent_structure, base, height_triax, height_gelcoat):
        Part.__init__(self, parent_structure, base,
            height=(height_triax+height_gelcoat))
        self.height_triax = height_triax
        self.height_gelcoat = height_gelcoat
    
    def __str__(self):
        return """base:    {0:6.4f} (meters)
height:  {1:6.4f} (meters)
|-> height_triax:  {2:6.4f} (meters)
|-> height_gelcoat:  {3:6.4f} (meters)""".format(self.base, self.height,
    self.height_triax, self.height_gelcoat)
    
    def create_layers(self, airfoil=None):
        """Create the gelcoat and triax layers in the external surface.

        <external_surface>.layer['gelcoat'] : gelcoat layer
        <external_surface>.layer['triax'] : triax layer

        """
        st = self.parent_structure
        # create the gelcoat layer
        af = st.parent_station.airfoil
        b = st.parent_station.parent_blade
        if airfoil is None:
            op_gelcoat = af.polygon  # outer profile is the airfoil profile
        elif airfoil == 'lower':
            op_gelcoat = af.lower_polygon
        elif airfoil == 'upper':
            op_gelcoat = af.upper_polygon
        else:
            raise ValueError("Keyword `airfoil` must be None, 'lower', or 'upper'.")
        ip_gelcoat = op_gelcoat.buffer(-self.height_gelcoat)
        polygon_gelcoat = op_gelcoat.difference(ip_gelcoat)
        self.layer['gelcoat'] = l.Layer(polygon_gelcoat,
            b.dict_of_materials['gelcoat'], parent_part=self,
            name='gelcoat', face_color='#5EE54C')
        assert self.layer['gelcoat'].polygon.geom_type == 'Polygon'
        if airfoil is None:
            st._list_of_layers.append(self.layer['gelcoat'])
        elif airfoil == 'lower':
            st._list_of_lower_layers.append(self.layer['gelcoat'])
        elif airfoil == 'upper':
            st._list_of_upper_layers.append(self.layer['gelcoat'])
        # create the triax layer
        op_triax = ip_gelcoat  # outer profile is the gelcoat inner profile
        ip_triax = op_triax.buffer(-self.height_triax)
        polygon_triax = op_triax.difference(ip_triax)
        self.layer['triax'] = l.Layer(polygon_triax,
            b.dict_of_materials['triaxial GFRP'], parent_part=self,
            name='triax', face_color='#5EE54C')
        assert self.layer['triax'].polygon.geom_type == 'Polygon'
        if airfoil is None:
            st._list_of_layers.append(self.layer['triax'])
        elif airfoil == 'lower':
            st._list_of_lower_layers.append(self.layer['triax'])
        elif airfoil == 'upper':
            st._list_of_upper_layers.append(self.layer['triax'])

    def add_new_layer(self, new_name, new_polygon, material):
        """Add a new layer to the external surface part."""
        st = self.parent_structure
        b = st.parent_station.parent_blade
        if material == 'triax':
            m = b.dict_of_materials['triaxial GFRP']
        elif material == 'gelcoat':
            m = b.dict_of_materials['gelcoat']
        else:
            raise Warning("Unrecognized material requested:", material)
        self.alt_layer[new_name] = l.Layer(new_polygon, m, parent_part=self, name=new_name)
        assert self.alt_layer[new_name].polygon.geom_type == 'Polygon'

    def create_alternate_layers(self):
        """Create the alternate layers for meshing the external surface.

        The alternate layers will be used to create meshes in TrueGrid.

        If the trailing edge is smooth (sharp_TE=False), save 8 quadrilaterals:
        -----------------------------------------------------------------------
        <external_surface>.layer['gelcoat, lower left'] : gelcoat layer (lower
            left quadrant)
        <external_surface>.layer['gelcoat, lower right'] : gelcoat layer (lower
            right quadrant)
        <external_surface>.layer['gelcoat, upper right'] : gelcoat layer (upper
            right quadrant)
        <external_surface>.layer['gelcoat, upper left'] : gelcoat layer (upper
            left quadrant)
        <external_surface>.layer['triax, lower left'] : triax layer (lower left
            quadrant)
        <external_surface>.layer['triax, lower right'] : triax layer (lower
            right quadrant)
        <external_surface>.layer['triax, upper right'] : triax layer (upper
            right quadrant)
        <external_surface>.layer['triax, upper left'] : triax layer (upper left
            quadrant)

        """
        st = self.parent_structure
        b = st.parent_station.parent_blade
        sharp_TE = st.parent_station.airfoil.has_sharp_TE
        material_list = ['triax', 'gelcoat']
        if not sharp_TE:
            # cut annulus into quadrants
            for material in material_list:
                try:
                    # try to access the annulus layer
                    p = self.layer[material].polygon  # annulus
                except KeyError:
                    # if the annulus layer doesn't exist, throw a warning
                    raise Warning("<ExternalSurface>.layer[{0}] was not found.\n  Try running <ExternalSurface>.create_layers() before running <ExternalSurface>.create_alternate_layers().".format(material))
                # cut the annulus into 4 curved rectangles, w/ one in each quadrant
                bb = self.bounding_box() # get bounding boxes for each quadrant
                for (label,box) in bb.items():
                    p_quad = p.intersection(box)
                    if material == 'triax':
                        dict_key = 'triaxial GFRP'
                    else:
                        dict_key = 'gelcoat'
                    self.layer[material+', '+label] = l.Layer(p_quad,
                        b.dict_of_materials[dict_key], parent_part=self,
                        name=(material+', '+label), face_color='#5EE54C')
                    # check that the layer just created is a polygon
                    assert self.layer[material+', '+label].polygon.geom_type == 'Polygon'
                    # no need to append these polygons to <station>._list_of_layers

    def bounding_box(self, x_boundary_buffer=1.2, y_boundary_buffer=1.2):
        """Returns list of 4 polygons for bounding boxes in each quadrant.

        The bounding boxes will be used to split the root buildup into 4 curved
        polygons.

        bb[0] : lower left quadrant, (0,0) to (x_min*1.2, y_min*1.2)
        bb[1] : lower right quadrant, (0,0) to (x_max*1.2, y_min*1.2)
        bb[2] : upper right quadrant, (0,0) to (x_max*1.2, y_max*1.2)
        bb[3] : upper left quadrant, (0,0) to (x_min*1.2, y_max*1.2)

        The points of each bounding box are labeled from 1 to 4 as:

        4---3
        |   |
        1---2

        Parameters
        ----------
        x_boundary_buffer : float (default: 1.2), factor to multiply with the
            minx and maxx bound of the airfoil polygon, to stretch the bounding
            box past the left and right edges of the airfoil polygon
        y_boundary_buffer : float (default: 1.2), factor to multiply with the
            miny and maxy bound of the airfoil polygon, to stretch the bounding
            box above and below the top and bottom edges of the airfoil polygon

        """
        af = self.parent_structure.parent_station.airfoil
        bb = {}
        (minx, miny, maxx, maxy) = af.polygon.bounds
        # lower left quadrant
        pt1 = (minx*x_boundary_buffer, miny*y_boundary_buffer)
        pt2 = (0.0, miny*y_boundary_buffer)
        pt3 = (0.0, 0.0)
        pt4 = (minx*x_boundary_buffer, 0.0)
        bb['lower left'] = Polygon([pt1, pt2, pt3, pt4])
        # lower right quadrant
        pt1 = (0.0, miny*y_boundary_buffer)
        pt2 = (maxx*x_boundary_buffer, miny*y_boundary_buffer)
        pt3 = (maxx*x_boundary_buffer, 0.0)
        pt4 = (0.0, 0.0)
        bb['lower right'] = Polygon([pt1, pt2, pt3, pt4])
        # upper right quadrant
        pt1 = (0.0, 0.0)
        pt2 = (maxx*x_boundary_buffer, 0.0)
        pt3 = (maxx*x_boundary_buffer, maxy*y_boundary_buffer)
        pt4 = (0.0, maxy*y_boundary_buffer)
        bb['upper right'] = Polygon([pt1, pt2, pt3, pt4])
        # upper left quadrant
        pt1 = (minx*x_boundary_buffer, 0.0)
        pt2 = (0.0, 0.0)
        pt3 = (0.0, maxy*y_boundary_buffer)
        pt4 = (minx*x_boundary_buffer, maxy*y_boundary_buffer)
        bb['upper left'] = Polygon([pt1, pt2, pt3, pt4])
        return bb

    def get_edges(self, which_layer):
        """Returns 4 arrays of coords for each edge of the chosen layer.

        Parameters
        ----------
        which_layer : str, the desired layer, either 'lower left',
            'lower right', 'upper right', or 'upper left'

        """
        # extract the desired layer
        try:
            lyr = self.layer[which_layer]
        except KeyError:
            raise ValueError("`which_layer` must be either 'gelcoat' or 'triax'!")
            # raise ValueError("`which_layer` must be either 'gelcoat, lower left', 'gelcoat, lower right', 'gelcoat, upper right', 'gelcoat, upper left', 'triax, lower left', 'triax, lower right', 'triax, upper right', or 'triax, upper left'!")
        p = lyr.polygon  # get the polygon for this layer
        # store the polygon exterior coords as a numpy array
        a = np.array(p.exterior.coords)
        # get the x- and y-coordinates of the polygon exterior
        x = a[:,0]
        y = a[:,1]
        # find the indices where the x-coord is equal to zero
        match_x = np.nonzero(x==0.0)[0]
        # find the indices where the y-coord is equal to the right edge
        match_y = np.nonzero(y==0.0)[0]
        # group all the indices together in a sorted array
        match = np.append(match_x, match_y)
        match.sort()
        # split the polygon up at each of the corners into 4 "edges"
        edge1 = a[match[0]:match[1]+1,:]
        edge2 = a[match[1]:match[2]+1,:]
        edge3 = a[match[2]:match[3]+1,:]
        try:
            edge4 = a[match[3]:match[4]+1,:]
        except IndexError:
            edge4 = np.append(a[match[3]:,:],a[1:match[0]+1,:],axis=0)
        return (edge1, edge2, edge3, edge4)

    def get_and_save_edges(self, which_layer):
        """Identifies and saves the left, top, right, and bottom edges.

        Parameters
        ----------
        which_layer : str, the desired layer, either 'lower left',
            'lower right', 'upper right', or 'upper left'

        This method saves the LTRB edges as attributes within the layer object.

        self.layer[which_layer].left : np.array, coords for left edge
        self.layer[which_layer].top : np.array, coords for top edge
        self.layer[which_layer].right : np.array, coords for right edge
        self.layer[which_layer].bottom : np.array, coords for bottom edge

        """
        # extract the desired layer
        try:
            lyr = self.layer[which_layer]
        except KeyError:
            raise ValueError("`which_layer` must be either 'gelcoat' or 'triax'!")
            # raise ValueError("`which_layer` must be either 'gelcoat, lower left', 'gelcoat, lower right', 'gelcoat, upper right', 'gelcoat, upper left', 'triax, lower left', 'triax, lower right', 'triax, upper right', or 'triax, upper left'!")
        edges = self.get_edges(which_layer)
        # get centroids
        centroids = []
        for edge in edges:
            centroids.append(asLineString(edge).centroid)
        # determine which edges are top, bottom, left, and right
        l = range(4)  # list of indices, one for each edge
        c = np.array([[centroids[0].x, centroids[0].y],
                      [centroids[1].x, centroids[1].y],
                      [centroids[2].x, centroids[2].y],
                      [centroids[3].x, centroids[3].y]])
        cx = c[:,0]
        cy = c[:,1]
        # find centroid at x=0
        ind_x = np.nonzero(cx==0.0)[0][0]
        l.remove(ind_x)  # remove the index for the right edge
        # find centroid at y=0
        ind_y = np.nonzero(cy==0.0)[0][0]
        l.remove(ind_y)  # remove the index for the left edge
        if which_layer == 'lower left':
            lyr.right = edges[ind_x]  # right edge saved!
            lyr.left = edges[ind_y]  # left edge saved!
        elif which_layer == 'lower right':
            lyr.left = edges[ind_x]  # left edge saved!
            lyr.right = edges[ind_y]  # right edge saved!
        elif which_layer == 'upper right':
            lyr.left = edges[ind_x]  # left edge saved!
            lyr.right = edges[ind_y]  # right edge saved!
        elif which_layer == 'upper left':
            lyr.right = edges[ind_x]  # right edge saved!
            lyr.left = edges[ind_y]  # left edge saved!
        # find top and bottom edges
        if centroids[l[0]].y > centroids[l[1]].y:
            lyr.top = edges[l[0]]     # top edge saved!
            lyr.bottom = edges[l[1]]  # bottom edge saved!
        else:
            lyr.top = edges[l[1]]     # top edge saved!
            lyr.bottom = edges[l[0]]  # bottom edge saved!


class RootBuildup(Part):
    """Define triax dimensions of the root buildup."""
    def create_layers(self, airfoil=None):
        """Create the triax layer in the root buildup.

        <root_buildup>.layer['triax'] : triax layer (entire annulus)

        """
        st = self.parent_structure
        # create the triax layer
        af = st.parent_station.airfoil
        b = st.parent_station.parent_blade
        if airfoil is None:
            op = af.polygon.buffer(-st.external_surface.height)
        elif airfoil == 'lower':
            op = af.lower_polygon.buffer(-st.lower_external_surface.height)
        elif airfoil == 'upper':
            op = af.upper_polygon.buffer(-st.upper_external_surface.height)
        else:
            raise ValueError("Keyword `airfoil` must be None, 'lower', or 'upper'.")
        ip = op.buffer(-self.height)
        p = op.difference(ip)  # this polygon is like an annulus
        self.layer['triax'] = l.Layer(p, b.dict_of_materials['triaxial GFRP'],
            parent_part=self, name='triax', face_color='#BE925A')
        # check that layer['triax'] is a Polygon
        assert self.layer['triax'].polygon.geom_type == 'Polygon'
        if airfoil is None:
            st._list_of_layers.append(self.layer['triax'])
        elif airfoil == 'lower':
            st._list_of_lower_layers.append(self.layer['triax'])
        elif airfoil == 'upper':
            st._list_of_upper_layers.append(self.layer['triax'])

    def add_new_layer(self, new_name, new_polygon, material='triax'):
        """Add a new alternate layer to the root buildup part."""
        st = self.parent_structure
        b = st.parent_station.parent_blade
        if material == 'triax':
            m = b.dict_of_materials['triaxial GFRP']
        else:
            raise Warning("Unrecognized material requested:", material)
        self.alt_layer[new_name] = l.Layer(new_polygon, m, parent_part=self, name=new_name, face_color='#BE925A')
        assert self.alt_layer[new_name].polygon.geom_type == 'Polygon'

    def create_alternate_layers(self):
        """Create the alternate triax layers for meshing the root buildup.

        These alternate layers will be used to create meshes in TrueGrid.

        <root_buildup>.layer['triax, lower left'] : triax layer (lower left
            quadrant)
        <root_buildup>.layer['triax, lower right'] : triax layer (lower right
            quadrant)
        <root_buildup>.layer['triax, upper right'] : triax layer (upper right
            quadrant)
        <root_buildup>.layer['triax, upper left'] : triax layer (upper left
            quadrant)

        Note: This stores 4 curved rectangles : .layer['triax, lower left'],
            .layer['triax, lower right'], .layer['triax, upper right'],
            .layer['triax, upper left']

        """
        st = self.parent_structure
        b = st.parent_station.parent_blade
        try:
            # try to access the annulus layer
            p = self.layer['triax'].polygon  # annulus
        except KeyError:
            # if the annulus layer doesn't exist, throw a warning
            raise Warning("<RootBuildup>.layer['triax'] was not found.\n  Try running <RootBuildup>.create_layers() before running <RootBuildup>.create_alternate_layers().")
        # cut the annulus into 4 curved rectangles, w/ one in each quadrant
        bb = self.bounding_box() # get bounding boxes for each quadrant
        for (label,box) in bb.items():
            p_quad = p.intersection(box)
            self.layer['triax, '+label] = l.Layer(p_quad,
                b.dict_of_materials['triaxial GFRP'], parent_part=self,
                name=('triax, '+label), face_color='#BE925A')
            # face color is brown
            # check that the layer just created is a polygon
            assert self.layer['triax, '+label].polygon.geom_type == 'Polygon'
            # no need to append these polygons to <station>._list_of_layers

    def bounding_box(self, x_boundary_buffer=1.2, y_boundary_buffer=1.2):
        """Returns dict of 4 polygons for bounding boxes in each quadrant.

        The bounding boxes will be used to split the root buildup into 4 curved
        polygons.

        bb['lower left'] : lower left quadrant, (0,0) to (x_min*1.2, y_min*1.2)
        bb['lower right'] : lower right quadrant, (0,0) to (x_max*1.2, y_min*1.2)
        bb['upper right'] : upper right quadrant, (0,0) to (x_max*1.2, y_max*1.2)
        bb['upper left'] : upper left quadrant, (0,0) to (x_min*1.2, y_max*1.2)

        The points of each bounding box are labeled from 1 to 4 as:

        4---3
        |   |
        1---2

        Parameters
        ----------
        x_boundary_buffer : float (default: 1.2), factor to multiply with the
            minx and maxx bound of the airfoil polygon, to stretch the bounding
            box past the left and right edges of the airfoil polygon
        y_boundary_buffer : float (default: 1.2), factor to multiply with the
            miny and maxy bound of the airfoil polygon, to stretch the bounding
            box above and below the top and bottom edges of the airfoil polygon

        """
        af = self.parent_structure.parent_station.airfoil
        bb = {}
        (minx, miny, maxx, maxy) = af.polygon.bounds
        # lower left quadrant
        pt1 = (minx*x_boundary_buffer, miny*y_boundary_buffer)
        pt2 = (0.0, miny*y_boundary_buffer)
        pt3 = (0.0, 0.0)
        pt4 = (minx*x_boundary_buffer, 0.0)
        bb['lower left'] = Polygon([pt1, pt2, pt3, pt4])
        # lower right quadrant
        pt1 = (0.0, miny*y_boundary_buffer)
        pt2 = (maxx*x_boundary_buffer, miny*y_boundary_buffer)
        pt3 = (maxx*x_boundary_buffer, 0.0)
        pt4 = (0.0, 0.0)
        bb['lower right'] = Polygon([pt1, pt2, pt3, pt4])
        # upper right quadrant
        pt1 = (0.0, 0.0)
        pt2 = (maxx*x_boundary_buffer, 0.0)
        pt3 = (maxx*x_boundary_buffer, maxy*y_boundary_buffer)
        pt4 = (0.0, maxy*y_boundary_buffer)
        bb['upper right'] = Polygon([pt1, pt2, pt3, pt4])
        # upper left quadrant
        pt1 = (minx*x_boundary_buffer, 0.0)
        pt2 = (0.0, 0.0)
        pt3 = (0.0, maxy*y_boundary_buffer)
        pt4 = (minx*x_boundary_buffer, maxy*y_boundary_buffer)
        bb['upper left'] = Polygon([pt1, pt2, pt3, pt4])
        return bb


class LE_Panel(Part):
    """Define foam dimensions of the leading edge panel."""
    def create_layers(self, airfoil=None):
        """Create the foam layer in the leading edge panel.

        <LE_panel>.layer['foam'] : the only layer in the LE panel (foam)

        """
        st = self.parent_structure
        # create the foam layer
        af = st.parent_station.airfoil
        b = st.parent_station.parent_blade
        # 1. get outer profile
        if airfoil is None:
            if st.root_buildup.exists():
                op = af.polygon.buffer(-(st.external_surface.height + 
                    st.root_buildup.height))
            else:
                op = af.polygon.buffer(-st.external_surface.height)
        elif airfoil == 'lower':
            if st.lower_root_buildup.exists():
                op = af.lower_polygon.buffer(
                    -(st.lower_external_surface.height + 
                    st.lower_root_buildup.height))
            else:
                op = af.lower_polygon.buffer(-st.lower_external_surface.height)
        elif airfoil == 'upper':
            if st.upper_root_buildup.exists():
                op = af.upper_polygon.buffer(
                    -(st.upper_external_surface.height + 
                    st.upper_root_buildup.height))
            else:
                op = af.upper_polygon.buffer(-st.upper_external_surface.height)
        else:
            raise ValueError("Keyword `airfoil` must be None, 'lower', or 'upper'.")
        # 2. erode the outer profile by the part thickness
        ip = op.buffer(-self.height)
        # 3. cut out the part interior from the outer profile
        ac = op.difference(ip)
        # 4. draw a bounding box at the part edges
        if airfoil is None:
            bb = self.bounding_box()
        elif airfoil == 'lower':
            bb = self.bounding_box(y_boundary_buffer=1.0, airfoil='lower')
        elif airfoil == 'upper':
            bb = self.bounding_box(y_boundary_buffer=1.0, airfoil='upper')
        # 5. cut out the structural part
        p = ac.intersection(bb)
        self.layer['foam'] = l.Layer(p, b.dict_of_materials['foam'],
            parent_part=self, name='foam', face_color='#00A64F')
        assert self.layer['foam'].polygon.geom_type == 'Polygon'
        if airfoil is None:
            st._list_of_layers.append(self.layer['foam'])
        elif airfoil == 'lower':
            st._list_of_lower_layers.append(self.layer['foam'])
        elif airfoil == 'upper':
            st._list_of_upper_layers.append(self.layer['foam'])


class SparCap(Part):
    """Define uniax dimensions of the lower and upper spar caps."""
    def create_layers(self, airfoil=None):
        """Create the uniax layers in the lower and upper spar caps.

        <spar_cap>.layer['lower'] : lower spar cap
        <spar_cap>.layer['upper'] : upper spar cap

        """
        st = self.parent_structure
        # create the uniax layer
        af = st.parent_station.airfoil
        b = st.parent_station.parent_blade
        # 1. get outer profile
        if airfoil is None:
            if st.root_buildup.exists():
                op = af.polygon.buffer(-(st.external_surface.height + 
                    st.root_buildup.height))
            else:
                op = af.polygon.buffer(-st.external_surface.height)
        elif airfoil == 'lower':
            if st.lower_root_buildup.exists():
                op = af.lower_polygon.buffer(
                    -(st.lower_external_surface.height + 
                    st.lower_root_buildup.height))
            else:
                op = af.lower_polygon.buffer(-st.lower_external_surface.height)
        elif airfoil == 'upper':
            if st.upper_root_buildup.exists():
                op = af.upper_polygon.buffer(
                    -(st.upper_external_surface.height + 
                    st.upper_root_buildup.height))
            else:
                op = af.upper_polygon.buffer(-st.upper_external_surface.height)
        else:
            raise ValueError("Keyword `airfoil` must be None, 'lower', or 'upper'.")
        # 2. erode the outer profile by the part thickness
        ip = op.buffer(-self.height)
        # 3. cut out the part interior from the outer profile
        ac = op.difference(ip)
        # 4. draw a bounding box at the part edges
        if airfoil is None:
            bb = self.bounding_box()
        elif airfoil == 'lower':
            bb = self.bounding_box(y_boundary_buffer=1.0, airfoil='lower')
        elif airfoil == 'upper':
            bb = self.bounding_box(y_boundary_buffer=1.0, airfoil='upper')
        # 5. cut out the structural part
        p = ac.intersection(bb)
        # 6. find the lower spar cap
        if p.geoms[0].centroid.y < p.geoms[1].centroid.y:
            pl = p.geoms[0]
            pu = p.geoms[1]
        else:
            pl = p.geoms[1]
            pu = p.geoms[0]
        # 7. add the lower spar cap
        self.layer['lower'] = l.Layer(pl, b.dict_of_materials['uniaxial GFRP'],
            parent_part=self, name='lower', face_color='#00ACEF')
        assert self.layer['lower'].polygon.geom_type == 'Polygon'
        if airfoil is None:
            st._list_of_layers.append(self.layer['lower'])
        elif airfoil == 'lower':
            st._list_of_lower_layers.append(self.layer['lower'])
        elif airfoil == 'upper':
            st._list_of_upper_layers.append(self.layer['lower'])
        # 8. add the upper spar cap
        self.layer['upper'] = l.Layer(pu, b.dict_of_materials['uniaxial GFRP'],
            parent_part=self, name='upper', face_color='#00ACEF')
        assert self.layer['upper'].polygon.geom_type == 'Polygon'
        if airfoil is None:
            st._list_of_layers.append(self.layer['upper'])
        elif airfoil == 'lower':
            st._list_of_lower_layers.append(self.layer['upper'])
        elif airfoil == 'upper':
            st._list_of_upper_layers.append(self.layer['upper'])


class AftPanel(Part):
    """Define foam dimensions of the lower and upper aft panels."""
    def __init__(self, num, parent_structure, base, height):
        Part.__init__(self, parent_structure, base, height)
        self.num = num

    def create_layers(self, airfoil=None):
        """Create the foam layers in the lower and upper aft panels.

        <aft_panel>.layer['lower'] : lower aft panel
        <aft_panel>.layer['upper'] : upper aft panel

        """
        st = self.parent_structure
        # create the foam layer
        af = st.parent_station.airfoil
        b = st.parent_station.parent_blade
        # 1. get outer profile
        if airfoil is None:
            if st.root_buildup.exists():
                op = af.polygon.buffer(-(st.external_surface.height + 
                    st.root_buildup.height))
            else:
                op = af.polygon.buffer(-st.external_surface.height)
        elif airfoil == 'lower':
            if st.lower_root_buildup.exists():
                op = af.lower_polygon.buffer(
                    -(st.lower_external_surface.height + 
                    st.lower_root_buildup.height))
            else:
                op = af.lower_polygon.buffer(-st.lower_external_surface.height)
        elif airfoil == 'upper':
            if st.upper_root_buildup.exists():
                op = af.upper_polygon.buffer(
                    -(st.upper_external_surface.height + 
                    st.upper_root_buildup.height))
            else:
                op = af.upper_polygon.buffer(-st.upper_external_surface.height)
        else:
            raise ValueError("Keyword `airfoil` must be None, 'lower', or 'upper'.")
        # 2. erode the outer profile by the part thickness
        ip = op.buffer(-self.height)
        # 3. cut out the part interior from the outer profile
        ac = op.difference(ip)
        # 4. draw a bounding box at the part edges
        if airfoil is None:
            bb = self.bounding_box()
        elif airfoil == 'lower':
            bb = self.bounding_box(y_boundary_buffer=1.0, airfoil='lower')
        elif airfoil == 'upper':
            bb = self.bounding_box(y_boundary_buffer=1.0, airfoil='upper')
        # 5. cut out the structural part
        p = ac.intersection(bb)
        # 6. find the lower aft panel
        if p.geoms[0].centroid.y < p.geoms[1].centroid.y:
            pl = p.geoms[0]
            pu = p.geoms[1]
        else:
            pl = p.geoms[1]
            pu = p.geoms[0]
        # 7. add the lower aft panel
        self.layer['lower'] = l.Layer(pl, b.dict_of_materials['foam'],
            parent_part=self, name='lower', face_color='#F58612')
        assert self.layer['lower'].polygon.geom_type == 'Polygon'
        if airfoil is None:
            st._list_of_layers.append(self.layer['lower'])
        elif airfoil == 'lower':
            st._list_of_lower_layers.append(self.layer['lower'])
        elif airfoil == 'upper':
            st._list_of_upper_layers.append(self.layer['lower'])
        # 8. add the upper aft panel
        self.layer['upper'] = l.Layer(pu, b.dict_of_materials['foam'],
            parent_part=self, name='upper', face_color='#F58612')
        assert self.layer['upper'].polygon.geom_type == 'Polygon'
        if airfoil is None:
            st._list_of_layers.append(self.layer['upper'])
        elif airfoil == 'lower':
            st._list_of_lower_layers.append(self.layer['upper'])
        elif airfoil == 'upper':
            st._list_of_upper_layers.append(self.layer['upper'])


class TE_Reinforcement(Part):
    """Define uniax and foam dimensions of a trailing edge reinforcement."""
    def __init__(self, parent_structure, base, height_uniax, height_foam):
        # check for nans in height dimensions, if they exist, replace with zero
        # this prevents `float + nan = nan` from occuring
        # for example:
        #   height_uniax = 0.018, height_foam = nan
        #   height = height_uniax + height_foam
        #   height = nan
        # but, we'd rather get the result:
        #   height = height_uniax + height_foam = 0.018
        if isnan(height_foam) and isnan(height_uniax):
            # unless both are nan, let's replace the single nan with a zero
            hu = np.nan
            hf = np.nan
        else:
            if isnan(height_foam):
                hf = 0.0
            else:
                hf = height_foam
            if isnan(height_uniax):
                hu = 0.0
            else:
                hu = height_uniax
        Part.__init__(self, parent_structure, base, height=(hu+hf))
        self.height_uniax = height_uniax
        self.height_foam = height_foam
        # these 4 attributes are assigned later by self.create_alternate_layers()
        self.left_vertex = None
        self.foam_left_vertex = None
        self.uniax_left_vertex = None
        self.uniax_right_vertex = None
    
    def __str__(self):
        return """base:    {0:6.4f} (meters)
height:  {1:6.4f} (meters)
|-> height_uniax:  {2:6.4f} (meters)
|-> height_foam:   {3:6.4f} (meters)""".format(self.base, self.height,
    self.height_uniax, self.height_foam)

    def create_layers(self, airfoil=None):
        """Create the uniax and foam layers in the TE reinforcement.

        The TE reinforcement is split into one OR two regions:
        <TE_reinforcement>.layer['uniax'] : uniax layer
        <TE_reinforcement>.layer['foam'] : foam layer (optional)

        """
        st = self.parent_structure
        # create the uniax layer
        af = st.parent_station.airfoil
        b = st.parent_station.parent_blade
        # 1. get outer profile
        if airfoil is None:
            if st.root_buildup.exists():
                op_uniax = af.polygon.buffer(-(st.external_surface.height + 
                    st.root_buildup.height))
            else:
                op_uniax = af.polygon.buffer(-st.external_surface.height)
        elif airfoil == 'lower':
            if st.lower_root_buildup.exists():
                op_uniax = af.lower_polygon.buffer(
                    -(st.lower_external_surface.height + 
                    st.lower_root_buildup.height))
            else:
                op_uniax = af.lower_polygon.buffer(
                    -st.lower_external_surface.height)
        elif airfoil == 'upper':
            if st.upper_root_buildup.exists():
                op_uniax = af.upper_polygon.buffer(
                    -(st.upper_external_surface.height + 
                    st.upper_root_buildup.height))
            else:
                op_uniax = af.upper_polygon.buffer(
                    -st.upper_external_surface.height)
        else:
            raise ValueError("Keyword `airfoil` must be None, 'lower', or 'upper'.")
        # 2. erode the outer profile by the uniax thickness
        ip_uniax = op_uniax.buffer(-self.height_uniax)
        # 3. cut out the uniax layer from the outer profile
        ac_uniax = op_uniax.difference(ip_uniax)
        # 4. draw a bounding box at the TE reinforcement edges
        if airfoil is None:
            bb = self.bounding_box()
        elif airfoil == 'lower':
            bb = self.bounding_box(y_boundary_buffer=1.0, airfoil='lower')
        elif airfoil == 'upper':
            bb = self.bounding_box(y_boundary_buffer=1.0, airfoil='upper')
        # 5. cut out the uniax layer
        polygon_uniax = ac_uniax.intersection(bb)
        # 6. add the uniax layer
        self.layer['uniax'] = l.Layer(polygon_uniax,
            b.dict_of_materials['uniaxial GFRP'], parent_part=self,
            name='uniax', face_color='#F366BA')
        assert self.layer['uniax'].polygon.geom_type == 'Polygon'
        if airfoil is None:
            st._list_of_layers.append(self.layer['uniax'])
        elif airfoil == 'lower':
            st._list_of_lower_layers.append(self.layer['uniax'])
        elif airfoil == 'upper':
            st._list_of_upper_layers.append(self.layer['uniax'])
        if not isnan(self.height_foam):
            # create the foam layer
            # 1. get outer profile
            op_foam = ip_uniax  # outer profile is the uniax inner profile
            # 2. erode the outer profile by the foam thickness
            ip_foam = op_foam.buffer(-self.height_foam)
            # 3. cut out the foam layer from the outer profile
            ac_foam = op_foam.difference(ip_foam)
            # 4. cut out the foam layer with the earlier bounding box
            polygon_foam = ac_foam.intersection(bb)
            # 5. add the foam layer
            self.layer['foam'] = l.Layer(polygon_foam,
                b.dict_of_materials['foam'], parent_part=self, name='foam',
                face_color='#F366BA')
            assert self.layer['foam'].polygon.geom_type == 'Polygon'
            if airfoil is None:
                st._list_of_layers.append(self.layer['foam'])
            elif airfoil == 'lower':
                st._list_of_lower_layers.append(self.layer['foam'])
            elif airfoil == 'upper':
                st._list_of_upper_layers.append(self.layer['foam'])

    def add_new_layer(self, new_name, new_polygon, material):
        """Add a new layer to the TE reinforcement part."""
        st = self.parent_structure
        b = st.parent_station.parent_blade
        if material == 'uniax':
            m = b.dict_of_materials['uniaxial GFRP']
        elif material == 'foam':
            m = b.dict_of_materials['foam']
        else:
            raise Warning("Unrecognized material requested:", material)
        self.alt_layer[new_name] = l.Layer(new_polygon, m, parent_part=self, name=new_name, face_color='#F366BA')
        assert self.alt_layer[new_name].polygon.geom_type == 'Polygon'

    def create_alternate_layers(self):
        """Create the alternate uniax and foam layers for meshing the TE reinf.

        These alternate layers will be used to create meshes in TrueGrid.

        <TE_reinforcement>.layer['uniax, upper left'] : uniax layer (upper left
            block)
        <TE_reinforcement>.layer['uniax, upper middle'] : uniax layer (upper
            middle block)
        <TE_reinforcement>.layer['uniax, upper right'] : uniax layer (upper
            right block)
        <TE_reinforcement>.layer['uniax, lower left'] : uniax layer (lower left
            block)
        <TE_reinforcement>.layer['uniax, lower middle'] : uniax layer (lower
            middle block)
        <TE_reinforcement>.layer['uniax, lower right'] : uniax layer (lower
            right block)
        <TE_reinforcement>.layer['foam, upper left'] : foam layer (upper left
            block)
        <TE_reinforcement>.layer['foam, upper middle'] : foam layer (upper
            middle block)
        <TE_reinforcement>.layer['foam, lower left'] : foam layer (lower left
            block)
        <TE_reinforcement>.layer['foam, lower middle'] : foam layer (lower
            middle block)

        Note: This stores 8 quadrilaterals and 2 triangles:
            8 quadrilaterals : .layer['uniax, upper left'], 
                .layer['uniax, upper middle'], .layer['uniax, upper right'],
                .layer['uniax, lower left'], .layer['uniax, lower middle'],
                .layer['uniax, lower right'], .layer['foam, upper left'],
                .layer['foam, lower left']
            2 triangles : .layer['foam, upper right'],
                .layer['foam, lower right']

        """
        st = self.parent_structure
        b = st.parent_station.parent_blade
        try:
            # try to access the entire uniax layer
            u = self.layer['uniax']
        except KeyError:
            raise Warning("<TE_Reinforcement>.layer['uniax'] was not found.\n  Try running <TE_Reinforcement>.create_layers() before running <TE_Reinforcement>.create_alternate_layers().")
        try:
            # try to access the entire foam layer
            f = self.layer['foam']
        except KeyError:
            raise Warning("<TE_Reinforcement>.layer['foam'] was not found.\n  Try running <TE_Reinforcement>.create_layers() before running <TE_Reinforcement>.create_alternate_layers().")
        # get vertex on left edge of TE (self.left_vertex)
        temp = np.average(
            [self.layer['uniax'].top, self.layer['uniax'].bottom], axis=0)
        self.left_vertex = np.average(temp, axis=0)
        # get vertex on left edge of foam (self.foam_left_vertex)
        ind = np.nonzero(f.left[:,0]==f.left[:,0].max())[0][0]
        self.foam_left_vertex = f.left[ind,:]
        # get vertex on left edge of uniax (self.uniax_left_vertex)
        ind = np.nonzero(u.left[:,0]==u.left[:,0].max())[0]
        if len(ind) == 1:
            self.uniax_left_vertex = u.left[ind[0],:]
        else:
            self.uniax_left_vertex = np.average(
                [u.left[ind[0],:], u.left[ind[1],:]], axis=0)
        # get vertex on right edge of uniax (self.uniax_right_vertex)
        ind = np.nonzero(u.right[:,0]==u.right[:,0].max())[0]
        if len(ind) == 1:
            self.uniax_right_vertex = u.right[ind[0],:]
        else:
            self.uniax_right_vertex = np.average(
                [u.right[ind[0],:], u.right[ind[1],:]], axis=0)
        bb = self.alternate_bounding_boxes()
        for (label,box) in bb.items():
            # split up uniax region into alternate layers
            u_new = u.polygon.intersection(box)
            if u_new.geom_type == 'MultiPolygon':
                # if the intersection operation picked up a second polygon,
                #   keep the bigger polygon, and neglect the other polygon
                if u_new.geoms[0].area > u_new.geoms[1].area:
                    u_new = u_new.geoms[0]
                    print " Warning: At Station #{0}, throwing out smaller polygon in \n  <TE_Reinforcement>.layer['foam, {1}'] ...".format(st.parent_station.station_num, label)
                else:
                    u_new = u_new.geoms[1]
                    print " Warning: At Station #{0}, throwing out smaller polygon in \n  <TE_Reinforcement>.layer['foam, {1}'] ...".format(st.parent_station.station_num, label)
            elif u_new.geom_type == 'GeometryCollection':
                # if the intersection operation picked up a Polygon and
                #   something else (like a LineString), only keep the
                #   Polygon, and throw everything else out
                if u_new.geoms[0].geom_type == 'Polygon':
                    u_new = u_new.geoms[0]
                elif u_new.geoms[1].geom_type == 'Polygon':
                    u_new = u_new.geoms[1]
            self.layer['uniax, '+label] = l.Layer(u_new,
                b.dict_of_materials['uniaxial GFRP'], parent_part=self,
                name=('uniax, '+label))
            # check that the layer just created is a polygon
            assert self.layer['uniax, '+label].polygon.geom_type == 'Polygon'
            # no need to append these polygons to <station>._list_of_layers
            # split up foam region into alternate layers
            if not label.endswith('right'):
                # get only 'left' and 'middle' bounding boxes
                f_new = f.polygon.intersection(box)
                if f_new.geom_type == 'MultiPolygon':
                    # if the intersection operation picked up a second polygon,
                    #   keep the bigger polygon, and neglect the other polygon
                    if f_new.geoms[0].area > f_new.geoms[1].area:
                        f_new = f_new.geoms[0]
                        print " Warning: At Station #{0}, throwing out smaller polygon in \n  <TE_Reinforcement>.layer['foam, {1}'] ...".format(st.parent_station.station_num, label)
                    else:
                        f_new = f_new.geoms[1]
                        print " Warning: At Station #{0}, throwing out smaller polygon in \n  <TE_Reinforcement>.layer['foam, {1}'] ...".format(st.parent_station.station_num, label)
                elif f_new.geom_type == 'GeometryCollection':
                    # if the intersection operation picked up a Polygon and
                    #   something else (like a LineString), only keep the
                    #   Polygon, and throw everything else out
                    if f_new.geoms[0].geom_type == 'Polygon':
                        f_new = f_new.geoms[0]
                    elif f_new.geoms[1].geom_type == 'Polygon':
                        f_new = f_new.geoms[1]
                self.layer['foam, '+label] = l.Layer(f_new,
                    b.dict_of_materials['foam'], parent_part=self,
                    name=('foam, '+label))
                # check that the layer just created is a polygon
                assert self.layer['foam, '+label].polygon.geom_type == 'Polygon'
                # no need to append these polygons to <station>._list_of_layers

    def alternate_bounding_boxes(self, x_boundary_buffer=1.2,
        y_boundary_buffer=1.2):
        """Returns dict of polygons for bounding boxes for alternate layers.

        The bounding boxes will be used to split the TE reinforcement into
        several quadrilateral and triangular regions for meshing.

        This method is called by self.create_alternate_layers().

        If the foam region exists, 6 bounding boxes are returned:
        bb['upper left']
        bb['upper middle']
        bb['upper right']
        bb['lower left']
        bb['lower middle']
        bb['lower right']

        If the foam region doesn't exist, 4 bounding boxes are returned:
        bb['upper left']
        bb['upper right']
        bb['lower left']
        bb['lower right']

        The points of each bounding box are labeled from 1 to 4 as:

        4---3
        |   |
        1---2

        Parameters
        ----------
        x_boundary_buffer : float (default: 1.2), factor to multiply with the
            minx and maxx bound of the airfoil polygon, to stretch the bounding
            box past the left and right edges of the airfoil polygon
        y_boundary_buffer : float (default: 1.2), factor to multiply with the
            miny and maxy bound of the airfoil polygon, to stretch the bounding
            box above and below the top and bottom edges of the airfoil polygon
        self.left_vertex : np.array, coordinates at the midpoint of the left
            edge of the TE reinforcement. This coordinate does not sit on the
            TE reinforcement, it lies between the top and bottom "legs" of the
            TE reinforcement. It help us draw the upper left and lower left
            bounding boxes.
        self.foam_left_vertex : np.array, coordinates for the sharp corner on
            the left edge of the foam region. If the foam region doesn't exist,
            foam_left_vertex=None.
        self.uniax_left_vertex : np.array, coordinates for the sharp corner on
            the left edge of the uniax region.
        self.uniax_right_vertex : np.array, coordinates for the middle of the
            blunt TE along the right edge of the uniax region.

        """
        af = self.parent_structure.parent_station.airfoil
        bb = {}
        (minx, miny, maxx, maxy) = af.polygon.bounds
        # upper right box
        pt1 = tuple(self.uniax_left_vertex)
        pt2 = tuple(self.uniax_right_vertex)
        pt3 = (self.right, maxy*y_boundary_buffer)
        pt4 = (self.uniax_left_vertex[0], maxy*y_boundary_buffer)
        bb['upper right'] = Polygon([pt1, pt2, pt3, pt4])
        # lower right box
        pt1 = (self.uniax_left_vertex[0], miny*y_boundary_buffer)
        pt2 = (self.right, miny*y_boundary_buffer)
        pt3 = tuple(self.uniax_right_vertex)
        pt4 = tuple(self.uniax_left_vertex)
        bb['lower right'] = Polygon([pt1, pt2, pt3, pt4])
        # upper middle box
        pt1 = tuple(self.foam_left_vertex)
        pt2 = tuple(self.uniax_left_vertex)
        pt3 = (self.uniax_left_vertex[0], maxy*y_boundary_buffer)
        pt4 = (self.foam_left_vertex[0], maxy*y_boundary_buffer)
        bb['upper middle'] = Polygon([pt1, pt2, pt3, pt4])
        # lower middle box
        pt1 = (self.foam_left_vertex[0], miny*y_boundary_buffer)
        pt2 = (self.uniax_left_vertex[0], miny*y_boundary_buffer)
        pt3 = tuple(self.uniax_left_vertex)
        pt4 = tuple(self.foam_left_vertex)
        bb['lower middle'] = Polygon([pt1, pt2, pt3, pt4])
        # upper left box
        pt1 = tuple(self.left_vertex)
        pt2 = tuple(self.foam_left_vertex)
        pt3 = (self.foam_left_vertex[0], maxy*y_boundary_buffer)
        pt4 = (self.left, maxy*y_boundary_buffer)
        bb['upper left'] = Polygon([pt1, pt2, pt3, pt4])
        # lower left box
        pt1 = (self.left, miny*y_boundary_buffer)
        pt2 = (self.foam_left_vertex[0], miny*y_boundary_buffer)
        pt3 = tuple(self.foam_left_vertex)
        pt4 = tuple(self.left_vertex)
        bb['lower left'] = Polygon([pt1, pt2, pt3, pt4])
        return bb


class ShearWeb(Part):
    """Define the biax (skin) and foam (core) dimensions of a shear web.

    Parameters
    ----------
    parent_structure : Structure object, the parent structure of this shear web
    base_biax : float, base of biax material in shear web
    base_foam : float, base of foam material in shear web
    x2 : float, chordwise distance from pitch axis to edge of shear web
    height : NaN (DO NOT SPECIFY), height of shear web

    """
    def __init__(self, num, parent_structure, base_biax, base_foam, x2,
        height=np.nan):
        Part.__init__(self, parent_structure, base=(2.0*base_biax+base_foam),
            height=height)
        self.base_biax = base_biax
        self.base_foam = base_foam
        self.x2 = x2
        self.cs_coords = None   # assigned later by <station>.find_SW_cs_coords()
        self.num = num

    def __str__(self):
        return """base:    {0:6.4f} (meters)
|-> base_biax:  {1:6.4f} (meters)
|-> base_foam:  {2:6.4f} (meters)
height:  {3} (meters)
x2:      {4:6.4f} (meters)""".format(self.base, self.base_biax,
    self.base_foam, self.height, self.x2)

    def bounding_box(self, y_boundary_buffer=1.2, airfoil=None):
        """Returns 3 bounding boxes for the biax and foam regions of the SW.

        The points of each bounding box are labeled from 1 to 4 as:

        4---3
        |   |
        1---2

        """
        af = self.parent_structure.parent_station.airfoil
        if airfoil is None:
            (minx, miny, maxx, maxy) = af.polygon.bounds
        elif airfoil == 'lower':
            (minx, miny, maxx, maxy) = af.lower_polygon.bounds
        elif airfoil == 'upper':
            (minx, miny, maxx, maxy) = af.upper_polygon.bounds
        else:
            raise ValueError("Keyword `airfoil` must be None, 'lower', or 'upper'.")
        # left biax bounding box
        pt1 = (self.left, miny*y_boundary_buffer)
        pt2 = (self.left+self.base_biax, miny*y_boundary_buffer)
        pt3 = (self.left+self.base_biax, maxy*y_boundary_buffer)
        pt4 = (self.left, maxy*y_boundary_buffer)
        left_biax_bb = Polygon([pt1, pt2, pt3, pt4])
        # foam bounding box
        pt1 = (self.left+self.base_biax, miny*y_boundary_buffer)
        pt2 = (self.right-self.base_biax, miny*y_boundary_buffer)
        pt3 = (self.right-self.base_biax, maxy*y_boundary_buffer)
        pt4 = (self.left+self.base_biax, maxy*y_boundary_buffer)
        foam_bb = Polygon([pt1, pt2, pt3, pt4])
        # right biax bounding box
        pt1 = (self.right-self.base_biax, miny*y_boundary_buffer)
        pt2 = (self.right, miny*y_boundary_buffer)
        pt3 = (self.right, maxy*y_boundary_buffer)
        pt4 = (self.right-self.base_biax, maxy*y_boundary_buffer)
        right_biax_bb = Polygon([pt1, pt2, pt3, pt4])
        return (left_biax_bb, foam_bb, right_biax_bb)
    
    def create_layers(self, airfoil=None):
        """Create the biax and foam layers in this shear web.

        <shear_web>.layer['biax, left'] is the left biax layer
        <shear_web>.layer['foam'] is the foam layer
        <shear_web>.layer['biax, right'] is the right biax layer

        """
        st = self.parent_structure
        # create the foam layer
        af = st.parent_station.airfoil
        b = st.parent_station.parent_blade
        # 1. get outer profile
        if airfoil is None:
            if st.root_buildup.exists():
                op = af.polygon.buffer(-(st.external_surface.height + 
                    st.root_buildup.height))
            else:
                op = af.polygon.buffer(-st.external_surface.height)
        elif airfoil == 'lower':
            if st.lower_root_buildup.exists():
                op = af.lower_polygon.buffer(
                    -(st.lower_external_surface.height + 
                    st.lower_root_buildup.height))
            else:
                op = af.lower_polygon.buffer(-st.lower_external_surface.height)
        elif airfoil == 'upper':
            if st.upper_root_buildup.exists():
                op = af.upper_polygon.buffer(
                    -(st.upper_external_surface.height + 
                    st.upper_root_buildup.height))
            else:
                op = af.upper_polygon.buffer(-st.upper_external_surface.height)
        else:
            raise ValueError("Keyword `airfoil` must be None, 'lower', or 'upper'.")
        # 2. get bounding boxes for the biax and foam regions
        if airfoil is None:
            (bb_left_biax, bb_foam, bb_right_biax) = self.bounding_box()
        elif airfoil == 'lower':
            (bb_left_biax, bb_foam, bb_right_biax) = self.bounding_box(
                y_boundary_buffer=1.0, airfoil='lower')
        elif airfoil == 'upper':
            (bb_left_biax, bb_foam, bb_right_biax) = self.bounding_box(
                y_boundary_buffer=1.0, airfoil='upper')
        # 3. cut out the structural parts
        p_left_biax = op.intersection(bb_left_biax)
        p_foam = op.intersection(bb_foam)
        p_right_biax = op.intersection(bb_right_biax)
        # 4. add the left biax layer
        self.layer['biax, left'] = l.Layer(p_left_biax,
            b.dict_of_materials['biaxial GFRP'], parent_part=self,
            name='biax, left', face_color='#FFF100')
        assert self.layer['biax, left'].polygon.geom_type == 'Polygon'
        if airfoil is None:
            st._list_of_layers.append(self.layer['biax, left'])
        elif airfoil == 'lower':
            st._list_of_lower_layers.append(self.layer['biax, left'])
        elif airfoil == 'upper':
            st._list_of_upper_layers.append(self.layer['biax, left'])
        # 5. add the foam layer
        self.layer['foam'] = l.Layer(p_foam, b.dict_of_materials['foam'],
            parent_part=self, name='foam', face_color='#FFF100')
        assert self.layer['foam'].polygon.geom_type == 'Polygon'
        if airfoil is None:
            st._list_of_layers.append(self.layer['foam'])
        elif airfoil == 'lower':
            st._list_of_lower_layers.append(self.layer['foam'])
        elif airfoil == 'upper':
            st._list_of_upper_layers.append(self.layer['foam'])
        # 6. add the right biax layer
        self.layer['biax, right'] = l.Layer(p_right_biax,
            b.dict_of_materials['biaxial GFRP'], parent_part=self,
            name='biax, right', face_color='#FFF100')
        assert self.layer['biax, right'].polygon.geom_type == 'Polygon'
        if airfoil is None:
            st._list_of_layers.append(self.layer['biax, right'])
        elif airfoil == 'lower':
            st._list_of_lower_layers.append(self.layer['biax, right'])
        elif airfoil == 'upper':
            st._list_of_upper_layers.append(self.layer['biax, right'])


class InternalSurface(Part):
    """Define triax and resin dimensions of the internal surface."""
    def __init__(self, num, parent_structure, base, height_triax,
        height_resin):
        Part.__init__(self, parent_structure, base,
            height=(height_triax+height_resin))
        self.height_triax = height_triax
        self.height_resin = height_resin
        self.num = num
    
    def __str__(self):
        return """base:    {0:6.4f} (meters)
height:  {1:6.4f} (meters)
|-> height_triax:  {2:6.4f} (meters)
|-> height_resin:  {3:6.4f} (meters)""".format(self.base, self.height,
    self.height_triax, self.height_resin)

    def interior_loop(self, p, area_threshold=10e-06, debug_flag=False):
        """Get the interior loop for this internal surface.

        Returns a polygon object of the requested interior loop.

        Parameters
        ----------
        p : polygon, merged polygon for all parts in this station
        area_threshold : float (default: 10e-06), threshold value to decide if
            an interior loop is valid (when the loop area > area_threshold)

        """
        st = self.parent_structure
        # p = st.merge_all_polygons(plot_flag=False)
        # find interior loops that have a 'big' area (> area_threshold)
        good_loops = []
        for interior in p.interiors:
            a = Polygon(interior).area
            if a > area_threshold:
                good_loops.append(interior)
        if debug_flag:
            print "Station #{0} has {1} good interior loops.".format(
                self.station_num, len(good_loops))
        # sort loops by x-coordinate of their centroids, smallest to largest
        if len(good_loops) > 1:
            good_loops.sort(key=attrgetter('centroid.x'))
        # get the interior loop for the desired internal surface
        #   internal_surface_1 ==> self.num=1 ==> good_loops[0]
        #   internal_surface_2 ==> self.num=2 ==> good_loops[1]
        #   internal_surface_3 ==> self.num=3 ==> good_loops[2]
        #   internal_surface_4 ==> self.num=4 ==> good_loops[3]
        loop = Polygon(good_loops[self.num-1])
        return loop

    def create_layers(self, merged_polygon, airfoil=None):
        """Create the triax and resin layers in the internal surface.

        <internal_surface>.layer['triax'] : triax layer
        <internal_surface>.layer['resin'] : resin layer

        """
        st = self.parent_structure
        b = st.parent_station.parent_blade
        # triax region
        op_triax = self.interior_loop(merged_polygon)
        ip_triax = op_triax.buffer(-self.height_triax)
        polygon_triax = op_triax.difference(ip_triax)
        self.layer['triax'] = l.Layer(polygon_triax,
            b.dict_of_materials['triaxial GFRP'], parent_part=self,
            name='triax', face_color='#999999')
        assert self.layer['triax'].polygon.geom_type == 'Polygon'
        if airfoil is None:
            st._list_of_layers.append(self.layer['triax'])
        elif airfoil == 'lower':
            st._list_of_lower_layers.append(self.layer['triax'])
        elif airfoil == 'upper':
            st._list_of_upper_layers.append(self.layer['triax'])
        else:
            raise ValueError("Keyword `airfoil` must be None, 'lower', or 'upper'.")
        # resin region
        op_resin = ip_triax
        ip_resin = op_resin.buffer(-self.height_resin)
        polygon_resin = op_resin.difference(ip_resin)
        self.layer['resin'] = l.Layer(polygon_resin,
            b.dict_of_materials['resin'], parent_part=self, name='resin',
            face_color='#999999')
        assert self.layer['resin'].polygon.geom_type == 'Polygon'
        if airfoil is None:
            st._list_of_layers.append(self.layer['resin'])
        elif airfoil == 'lower':
            st._list_of_lower_layers.append(self.layer['resin'])
        elif airfoil == 'upper':
            st._list_of_upper_layers.append(self.layer['resin'])

    def add_new_layer(self, new_name, new_polygon, material):
        """Add a new alternate layer to the internal surface part."""
        st = self.parent_structure
        b = st.parent_station.parent_blade
        if material == 'triax':
            m = b.dict_of_materials['triaxial GFRP']
        elif material == 'resin':
            m = b.dict_of_materials['resin']
        else:
            raise Warning("Unrecognized material requested:", material)
        self.alt_layer[new_name] = l.Layer(new_polygon, m, parent_part=self, name=new_name)
        assert self.alt_layer[new_name].polygon.geom_type == 'Polygon'


class MonoplaneStructure:
    """Define the monoplane laminate schedule (internal dimensions)."""
    def __init__(self, h_RB, b_SC, h_SC, b_SW1_biax, b_SW1_foam, x2_SW1,
                 b_SW2_biax, b_SW2_foam, x2_SW2, b_SW3_biax, b_SW3_foam,
                 x2_SW3, b_TE_reinf, h_TE_reinf_uniax, h_TE_reinf_foam,
                 h_LE_panel, h_aft_panel_1, h_aft_panel_2, h_int_surf_1_triax,
                 h_int_surf_1_resin, h_int_surf_2_triax, h_int_surf_2_resin,
                 h_int_surf_3_triax, h_int_surf_3_resin, h_int_surf_4_triax,
                 h_int_surf_4_resin, h_ext_surf_triax, h_ext_surf_gelcoat,
                 parent_station):
        self.parent_station = parent_station
        self._list_of_layers = []
        self._dict_of_edge_nums = {}
        self.truegrid_input_filename = 'mesh_stn{0:02d}_start.tg'.format(self.parent_station.station_num)
        self.root_buildup = RootBuildup(
            parent_structure = self,
            base = np.nan,
            height = h_RB)
        self.spar_cap = SparCap(
            parent_structure = self,
            base = b_SC,
            height = h_SC)
        self.shear_web_1 = ShearWeb(
            num = 1,
            parent_structure = self,
            base_biax = b_SW1_biax,
            base_foam = b_SW1_foam,
            x2 = x2_SW1)
        self.shear_web_2 = ShearWeb(
            num = 2,
            parent_structure = self,
            base_biax = b_SW2_biax,
            base_foam = b_SW2_foam,
            x2 = x2_SW2)
        self.shear_web_3 = ShearWeb(
            num = 3,
            parent_structure = self,
            base_biax = b_SW3_biax,
            base_foam = b_SW3_foam,
            x2 = x2_SW3)
        self.TE_reinforcement = TE_Reinforcement(
            parent_structure = self,
            base = b_TE_reinf,
            height_uniax = h_TE_reinf_uniax,
            height_foam = h_TE_reinf_foam)
        self.LE_panel = LE_Panel(
            parent_structure = self,
            base = np.nan,
            height = h_LE_panel)
        self.aft_panel_1 = AftPanel(
            num = 1,
            parent_structure = self,
            base = np.nan,
            height = h_aft_panel_1)
        self.aft_panel_2 = AftPanel(
            num = 2,
            parent_structure = self,
            base = np.nan,
            height = h_aft_panel_2)
        self.internal_surface_1 = InternalSurface(
            num = 1,
            parent_structure = self,
            base = np.nan,
            height_triax = h_int_surf_1_triax,
            height_resin = h_int_surf_1_resin)
        self.internal_surface_2 = InternalSurface(
            num = 2,
            parent_structure = self,
            base = np.nan,
            height_triax = h_int_surf_2_triax,
            height_resin = h_int_surf_2_resin)
        self.internal_surface_3 = InternalSurface(
            num = 3,
            parent_structure = self,
            base = np.nan,
            height_triax = h_int_surf_3_triax,
            height_resin = h_int_surf_3_resin)
        self.internal_surface_4 = InternalSurface(
            num = 4,
            parent_structure = self,
            base = np.nan,
            height_triax = h_int_surf_4_triax,
            height_resin = h_int_surf_4_resin)
        self.external_surface = ExternalSurface(
            parent_structure = self,
            base = np.nan,
            height_triax = h_ext_surf_triax,
            height_gelcoat = h_ext_surf_gelcoat)
        self.area = None
        self.mass = None

    def __str__(self):
        """Returns a string of all the internal dimensions for this structure."""
        s = ''
        s += "--- ROOT BUILDUP ---\n"
        s += str(self.root_buildup) + '\n'
        s += "--- SPAR CAP ---\n"
        s += str(self.spar_cap) + '\n'
        s += "--- SHEAR WEB 1 ---\n"
        s += str(self.shear_web_1) + '\n'
        s += "--- SHEAR WEB 2 ---\n"
        s += str(self.shear_web_2) + '\n'
        s += "--- SHEAR WEB 3 ---\n"
        s += str(self.shear_web_3) + '\n'
        s += "---TE REINFORCEMENT ---\n"
        s += str(self.TE_reinforcement) + '\n'
        s += "--- LE PANEL ---\n"
        s += str(self.LE_panel) + '\n'
        s += "--- AFT PANEL 1 ---\n"
        s += str(self.aft_panel_1) + '\n'
        s += "--- AFT PANEL 2 ---\n"
        s += str(self.aft_panel_2) + '\n'
        s += "--- INTERNAL SURFACE 1 ---\n"
        s += str(self.internal_surface_1) + '\n'
        s += "--- INTERNAL SURFACE 2 ---\n"
        s += str(self.internal_surface_2) + '\n'
        s += "--- INTERNAL SURFACE 3 ---\n"
        s += str(self.internal_surface_3) + '\n'
        s += "--- INTERNAL SURFACE 4 ---\n"
        s += str(self.internal_surface_4) + '\n'
        s += "--- EXTERNAL SURFACE ---\n"
        s += str(self.external_surface) + '\n'
        return s

    def which_parts_exist(self):
        """Check which structural parts exist at this monoplane station.

        Returns a dictionary of booleans.

        """
        d = {'root buildup': self.root_buildup.exists(),
             'spar cap': self.spar_cap.exists(),
             'shear web 1': self.shear_web_1.exists(),
             'shear web 2': self.shear_web_2.exists(),
             'shear web 3': self.shear_web_3.exists(),
             'TE reinforcement': self.TE_reinforcement.exists(),
             'LE panel': self.LE_panel.exists(),
             'aft panel 1': self.aft_panel_1.exists(),
             'aft panel 2': self.aft_panel_2.exists(),
             'internal surface 1': self.internal_surface_1.exists(),
             'internal surface 2': self.internal_surface_2.exists(),
             'internal surface 3': self.internal_surface_3.exists(),
             'internal surface 4': self.internal_surface_4.exists(),
             'external surface': self.external_surface.exists()}
        return d

    def create_all_layers(self):
        """Create polygons for single material layers of each structural part."""
        if self.external_surface.exists():
            self.external_surface.create_layers()
        if self.root_buildup.exists():
            self.root_buildup.create_layers()
        if self.LE_panel.exists():
            self.LE_panel.create_layers()
        if self.spar_cap.exists():
            self.spar_cap.create_layers()
        if self.aft_panel_1.exists():
            self.aft_panel_1.create_layers()
        if self.aft_panel_2.exists():
            self.aft_panel_2.create_layers()
        if self.shear_web_1.exists():
            self.shear_web_1.create_layers()
        if self.shear_web_2.exists():
            self.shear_web_2.create_layers()
        if self.shear_web_3.exists():
            self.shear_web_3.create_layers()
        if self.TE_reinforcement.exists():
            self.TE_reinforcement.create_layers()
        mp = self.merge_all_polygons()
        if self.internal_surface_1.exists():
            self.internal_surface_1.create_layers(mp)
        if self.internal_surface_2.exists():
            self.internal_surface_2.create_layers(mp)
        if self.internal_surface_3.exists():
            self.internal_surface_3.create_layers(mp)
        if self.internal_surface_4.exists():
            self.internal_surface_4.create_layers(mp)

    def create_all_alternate_layers(self):
        """Create alternate layers for meshing certain structural parts."""
        if self.root_buildup.exists():
            self.root_buildup.create_alternate_layers()
        if self.TE_reinforcement.exists():
            if self.parent_station.airfoil.has_sharp_TE:
                self.TE_reinforcement.create_alternate_layers()

    def merge_all_polygons(self, plot_flag=False):
        """Merges all the layer polygons in this structure into one polygon.

        NOTE: internal surface polygons are NOT merged!

        """
        stn = self.parent_station
        if plot_flag:
            fig, ax = plt.subplots()
            ax.set_title("Station #{0}, {1}, {2}% span".format(stn.station_num,
                stn.airfoil.name, stn.coords.x1))
            ax.set_aspect('equal')
            ax.grid('on')
            ax.set_xlabel('x2 [meters]')
            ax.set_ylabel('x3 [meters]')
            patch = PolygonPatch(stn.airfoil.polygon, fc='None', ec='#999999',
                alpha=0.8)
            ax.add_patch(patch)
            (minx, miny, maxx, maxy) = stn.airfoil.polygon.bounds
            ax.set_xlim([minx*1.2,maxx*1.2])
            ax.set_ylim([miny*1.2,maxy*1.2])
        # merge everything
        list_of_polygons = []
        for layer in self._list_of_layers:
            list_of_polygons.append(layer.polygon)
        try:
            p = cascaded_union(list_of_polygons)
        except ValueError:
            # gather all the parts
            p = self.external_surface.layer['gelcoat'].polygon
            p = p.union(self.external_surface.layer['triax'].polygon)
            if self.root_buildup.exists():
                RB = self.root_buildup.layer['triax'].polygon
                p = p.union(RB)
            if self.LE_panel.exists():
                LE = self.LE_panel.layer['foam'].polygon
                p = p.union(LE)
            if self.spar_cap.exists():
                sc_l = self.spar_cap.layer['lower'].polygon
                try:
                    p = p.union(sc_l)
                except TopologicalError:
                    print " [Warning] could not merge lower spar cap in Station #{0} ... skipping!".format(self.parent_station.station_num)
                sc_u = self.spar_cap.layer['upper'].polygon
                try:
                    p = p.union(sc_u)
                except TopologicalError:
                    print " [Warning] could not merge upper spar cap in Station #{0}".format(self.parent_station.station_num)
                    try:
                        p = sc_u.union(p)
                        print " ... SUCCESSFULLY merged upper spar cap on the second try!"
                    except TopologicalError:
                        print " ... on second try, still COULD NOT merge upper spar cap!"
            if self.aft_panel_1.exists():
                aft1_u = self.aft_panel_1.layer['upper'].polygon
                aft1_l = self.aft_panel_1.layer['lower'].polygon
                p = p.union(aft1_u)
                p = p.union(aft1_l)
            if self.aft_panel_2.exists():
                aft2_u = self.aft_panel_2.layer['upper'].polygon
                aft2_l = self.aft_panel_2.layer['lower'].polygon
                p = p.union(aft2_u)
                p = p.union(aft2_l)
            if self.TE_reinforcement.exists():
                TE_uniax = self.TE_reinforcement.layer['uniax'].polygon
                try:
                    p = p.union(TE_uniax)
                except TopologicalError:
                    print " [Warning] could not merge uniax layer of TE reinforcement in Station #{0}".format(self.parent_station.station_num)
                    try:
                        p = TE_uniax.union(p)
                        print " ... SUCCESSFULLY merged uniax layer of TE reinforcement on the second try!"
                    except TopologicalError:
                        print " ... on second try, still COULD NOT merge uniax layer of TE reinforcement!"
                        try:
                            p = cascaded_union([p, TE_uniax])
                        except ValueError:
                            print " ... on third try, still COULD NOT merge uniax layer of TE reinforcement!"
                try:
                    TE_foam = self.TE_reinforcement.layer['foam'].polygon
                    try:
                        p = p.union(TE_foam)
                    except TopologicalError:
                        print " [Warning] could not merge foam layer of TE reinforcement in Station #{0} ... skipping!".format(self.parent_station.station_num)
                except ValueError:
                    print " foam layer of TE reinforcement does not exist in Station #{0}".format(self.parent_station.station_num)
            if self.shear_web_1.exists():
                sw1 = self.shear_web_1.layer['biax, left'].polygon.union(self.shear_web_1.layer['foam'].polygon)
                sw1 = sw1.union(self.shear_web_1.layer['biax, right'].polygon)
                p = p.union(sw1)
            if self.shear_web_2.exists():
                sw2 = self.shear_web_2.layer['biax, left'].polygon.union(self.shear_web_2.layer['foam'].polygon)
                sw2 = sw2.union(self.shear_web_2.layer['biax, right'].polygon)
                p = p.union(sw2)
            if self.shear_web_3.exists():
                sw3 = self.shear_web_3.layer['biax, left'].polygon.union(self.shear_web_3.layer['foam'].polygon)
                sw3 = sw3.union(self.shear_web_3.layer['biax, right'].polygon)
                p = p.union(sw3)
        if plot_flag:
            # plot the merged polygon
            patch2 = PolygonPatch(p, fc='#4000FF', ec = '#000000', alpha=0.8)
            ax.add_patch(patch2)
            plt.show()
        return p

    def calculate_area(self):
        """Add the area of all polygons in this station."""
        a = 0
        for layer in self._list_of_layers:
            a += layer.polygon.area
        self.area = a
        return a

    def calculate_mass(self):
        """Add the mass (per unit length) of all polygons in this station."""
        m = 0
        for layer in self._list_of_layers:
            m += layer.mass
        self.mass = m
        return m

    def calculate_all_percent_areas(self, print_flag=False):
        """Calculate the percent areas of all parts in this station.

        Returns a dictionary of area fractions for each structural part.

        """
        d = {}
        if print_flag:
            print " ----- STATION #{0} -----".format(self.parent_station.station_num)
        # d['blade station'] = self.parent_station.station_num
        if self.external_surface.exists():
            d['external surface (gelcoat)'] = self.external_surface.layer['gelcoat'].area_fraction()
            d['external surface (triax)'] = self.external_surface.layer['triax'].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, external surface, gelcoat".format(
                    self.external_surface.layer['gelcoat'].area_fraction())
                print "  {0:5.1%} area, external surface, triax".format(
                    self.external_surface.layer['triax'].area_fraction())
        if self.root_buildup.exists():
            d['root buildup'] = self.root_buildup.layer['triax'].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, root buildup".format(
                    self.root_buildup.layer['triax'].area_fraction())
        if self.spar_cap.exists():
            d['spar cap (lower)'] = self.spar_cap.layer['lower'].area_fraction()
            d['spar cap (upper)'] = self.spar_cap.layer['upper'].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, spar cap, lower".format(
                    self.spar_cap.layer['lower'].area_fraction())
                print "  {0:5.1%} area, spar cap, upper".format(
                    self.spar_cap.layer['upper'].area_fraction())
        if self.aft_panel_1.exists():
            d['aft panel 1 (lower)'] = self.aft_panel_1.layer['lower'].area_fraction()
            d['aft panel 1 (upper)'] = self.aft_panel_1.layer['upper'].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, aft panel 1, lower".format(
                    self.aft_panel_1.layer['lower'].area_fraction())
                print "  {0:5.1%} area, aft panel 1, upper".format(
                    self.aft_panel_1.layer['upper'].area_fraction())
        if self.aft_panel_2.exists():
            d['aft panel 2 (lower)'] = self.aft_panel_2.layer['lower'].area_fraction()
            d['aft panel 2 (upper)'] = self.aft_panel_2.layer['upper'].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, aft panel 2, lower".format(
                    self.aft_panel_2.layer['lower'].area_fraction())
                print "  {0:5.1%} area, aft panel 2, upper".format(
                    self.aft_panel_2.layer['upper'].area_fraction())
        if self.LE_panel.exists():
            d['LE panel'] = self.LE_panel.layer['foam'].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, LE panel".format(
                    self.LE_panel.layer['foam'].area_fraction())
        if self.shear_web_1.exists():
            d['shear web 1 (left biax)'] = self.shear_web_1.layer['biax, left'].area_fraction()
            d['shear web 1 (foam)'] = self.shear_web_1.layer['foam'].area_fraction()
            d['shear web 1 (right biax)'] = self.shear_web_1.layer['biax, right'].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, shear web 1, left biax".format(
                    self.shear_web_1.layer['biax, left'].area_fraction())
                print "  {0:5.1%} area, shear web 1, foam".format(
                    self.shear_web_1.layer['foam'].area_fraction())
                print "  {0:5.1%} area, shear web 1, right biax".format(
                    self.shear_web_1.layer['biax, right'].area_fraction())
        if self.shear_web_2.exists():
            d['shear web 2 (left biax)'] = self.shear_web_2.layer['biax, left'].area_fraction()
            d['shear web 2 (foam)'] = self.shear_web_2.layer['foam'].area_fraction()
            d['shear web 2 (right biax)'] = self.shear_web_2.layer['biax, right'].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, shear web 2, left biax".format(
                    self.shear_web_2.layer['biax, left'].area_fraction())
                print "  {0:5.1%} area, shear web 2, foam".format(
                    self.shear_web_2.layer['foam'].area_fraction())
                print "  {0:5.1%} area, shear web 2, right biax".format(
                    self.shear_web_2.layer['biax, right'].area_fraction())
        if self.shear_web_3.exists():
            d['shear web 3 (left biax)'] = self.shear_web_3.layer['biax, left'].area_fraction()
            d['shear web 3 (foam)'] = self.shear_web_3.layer['foam'].area_fraction()
            d['shear web 3 (right biax)'] = self.shear_web_3.layer['biax, right'].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, shear web 3, left biax".format(
                    self.shear_web_3.layer['biax, left'].area_fraction())
                print "  {0:5.1%} area, shear web 3, foam".format(
                    self.shear_web_3.layer['foam'].area_fraction())
                print "  {0:5.1%} area, shear web 3, right biax".format(
                    self.shear_web_3.layer['biax, right'].area_fraction())
        if self.TE_reinforcement.exists():
            d['TE reinforcement (uniax)'] = self.TE_reinforcement.layer['uniax'].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, TE reinforcement, uniax".format(
                    self.TE_reinforcement.layer['uniax'].area_fraction())
            try:
                d['TE reinforcement (foam)'] = self.TE_reinforcement.layer['foam'].area_fraction()
                if print_flag:
                    print "  {0:5.1%} area, TE reinforcement, foam".format(
                        self.TE_reinforcement.layer['foam'].area_fraction())
            except KeyError:
                # the foam layer doesn't exist in TE reinf at this station
                d['TE reinforcement (foam)'] = 0.0
        if self.internal_surface_1.exists():
            d['internal surface 1 (triax)'] = self.internal_surface_1.layer['triax'].area_fraction()
            d['internal surface 1 (resin)'] = self.internal_surface_1.layer['resin'].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, internal surface 1, triax".format(
                    self.internal_surface_1.layer['triax'].area_fraction())
                print "  {0:5.1%} area, internal surface 1, resin".format(
                    self.internal_surface_1.layer['resin'].area_fraction())
        if self.internal_surface_2.exists():
            d['internal surface 2 (triax)'] = self.internal_surface_2.layer['triax'].area_fraction()
            d['internal surface 2 (resin)'] = self.internal_surface_2.layer['resin'].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, internal surface 2, triax".format(
                    self.internal_surface_2.layer['triax'].area_fraction())
                print "  {0:5.1%} area, internal surface 2, resin".format(
                    self.internal_surface_2.layer['resin'].area_fraction())
        if self.internal_surface_3.exists():
            d['internal surface 3 (triax)'] = self.internal_surface_3.layer['triax'].area_fraction()
            d['internal surface 3 (resin)'] = self.internal_surface_3.layer['resin'].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, internal surface 3, triax".format(
                    self.internal_surface_3.layer['triax'].area_fraction())
                print "  {0:5.1%} area, internal surface 3, resin".format(
                    self.internal_surface_3.layer['resin'].area_fraction())
        if self.internal_surface_4.exists():
            d['internal surface 4 (triax)'] = self.internal_surface_4.layer['triax'].area_fraction()
            d['internal surface 4 (resin)'] = self.internal_surface_4.layer['resin'].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, internal surface 4, triax".format(
                    self.internal_surface_4.layer['triax'].area_fraction())
                print "  {0:5.1%} area, internal surface 4, resin".format(
                    self.internal_surface_4.layer['resin'].area_fraction())
        return d

    def calculate_all_percent_masses(self, print_flag=False):
        """Calculate the mass fractions of all parts in this station.

        Returns a dictionary of mass (per unit length) fractions for each
        structural part.

        """
        d = {}
        if print_flag:
            print " ----- STATION #{0} -----".format(self.parent_station.station_num)
        # d['blade station'] = self.parent_station.station_num
        if self.external_surface.exists():
            d['external surface (gelcoat)'] = self.external_surface.layer['gelcoat'].mass_fraction()
            d['external surface (triax)'] = self.external_surface.layer['triax'].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, external surface, gelcoat".format(
                    self.external_surface.layer['gelcoat'].mass_fraction())
                print "  {0:5.1%} mass/length, external surface, triax".format(
                    self.external_surface.layer['triax'].mass_fraction())
        if self.root_buildup.exists():
            d['root buildup'] = self.root_buildup.layer['triax'].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, root buildup".format(
                    self.root_buildup.layer['triax'].mass_fraction())
        if self.spar_cap.exists():
            d['spar cap (lower)'] = self.spar_cap.layer['lower'].mass_fraction()
            d['spar cap (upper)'] = self.spar_cap.layer['upper'].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, spar cap, lower".format(
                    self.spar_cap.layer['lower'].mass_fraction())
                print "  {0:5.1%} mass/length, spar cap, upper".format(
                    self.spar_cap.layer['upper'].mass_fraction())
        if self.aft_panel_1.exists():
            d['aft panel 1 (lower)'] = self.aft_panel_1.layer['lower'].mass_fraction()
            d['aft panel 1 (upper)'] = self.aft_panel_1.layer['upper'].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, aft panel 1, lower".format(
                    self.aft_panel_1.layer['lower'].mass_fraction())
                print "  {0:5.1%} mass/length, aft panel 1, upper".format(
                    self.aft_panel_1.layer['upper'].mass_fraction())
        if self.aft_panel_2.exists():
            d['aft panel 2 (lower)'] = self.aft_panel_2.layer['lower'].mass_fraction()
            d['aft panel 2 (upper)'] = self.aft_panel_2.layer['upper'].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, aft panel 2, lower".format(
                    self.aft_panel_2.layer['lower'].mass_fraction())
                print "  {0:5.1%} mass/length, aft panel 2, upper".format(
                    self.aft_panel_2.layer['upper'].mass_fraction())
        if self.LE_panel.exists():
            d['LE panel'] = self.LE_panel.layer['foam'].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, LE panel".format(
                    self.LE_panel.layer['foam'].mass_fraction())
        if self.shear_web_1.exists():
            d['shear web 1 (left biax)'] = self.shear_web_1.layer['biax, left'].mass_fraction()
            d['shear web 1 (foam)'] = self.shear_web_1.layer['foam'].mass_fraction()
            d['shear web 1 (right biax)'] = self.shear_web_1.layer['biax, right'].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, shear web 1, left biax".format(
                    self.shear_web_1.layer['biax, left'].mass_fraction())
                print "  {0:5.1%} mass/length, shear web 1, foam".format(
                    self.shear_web_1.layer['foam'].mass_fraction())
                print "  {0:5.1%} mass/length, shear web 1, right biax".format(
                    self.shear_web_1.layer['biax, right'].mass_fraction())
        if self.shear_web_2.exists():
            d['shear web 2 (left biax)'] = self.shear_web_2.layer['biax, left'].mass_fraction()
            d['shear web 2 (foam)'] = self.shear_web_2.layer['foam'].mass_fraction()
            d['shear web 2 (right biax)'] = self.shear_web_2.layer['biax, right'].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, shear web 2, left biax".format(
                    self.shear_web_2.layer['biax, left'].mass_fraction())
                print "  {0:5.1%} mass/length, shear web 2, foam".format(
                    self.shear_web_2.layer['foam'].mass_fraction())
                print "  {0:5.1%} mass/length, shear web 2, right biax".format(
                    self.shear_web_2.layer['biax, right'].mass_fraction())
        if self.shear_web_3.exists():
            d['shear web 3 (left biax)'] = self.shear_web_3.layer['biax, left'].mass_fraction()
            d['shear web 3 (foam)'] = self.shear_web_3.layer['foam'].mass_fraction()
            d['shear web 3 (right biax)'] = self.shear_web_3.layer['biax, right'].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, shear web 3, left biax".format(
                    self.shear_web_3.layer['biax, left'].mass_fraction())
                print "  {0:5.1%} mass/length, shear web 3, foam".format(
                    self.shear_web_3.layer['foam'].mass_fraction())
                print "  {0:5.1%} mass/length, shear web 3, right biax".format(
                    self.shear_web_3.layer['biax, right'].mass_fraction())
        if self.TE_reinforcement.exists():
            d['TE reinforcement (uniax)'] = self.TE_reinforcement.layer['uniax'].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, TE reinforcement, uniax".format(
                    self.TE_reinforcement.layer['uniax'].mass_fraction())
            try:
                d['TE reinforcement (foam)'] = self.TE_reinforcement.layer['foam'].mass_fraction()
                if print_flag:
                    print "  {0:5.1%} mass/length, TE reinforcement, foam".format(
                        self.TE_reinforcement.layer['foam'].mass_fraction())
            except KeyError:
                # the foam layer doesn't exist in TE reinf at this station
                d['TE reinforcement (foam)'] = 0.0
        if self.internal_surface_1.exists():
            d['internal surface 1 (triax)'] = self.internal_surface_1.layer['triax'].mass_fraction()
            d['internal surface 1 (resin)'] = self.internal_surface_1.layer['resin'].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, internal surface 1, triax".format(
                    self.internal_surface_1.layer['triax'].mass_fraction())
                print "  {0:5.1%} mass/length, internal surface 1, resin".format(
                    self.internal_surface_1.layer['resin'].mass_fraction())
        if self.internal_surface_2.exists():
            d['internal surface 2 (triax)'] = self.internal_surface_2.layer['triax'].mass_fraction()
            d['internal surface 2 (resin)'] = self.internal_surface_2.layer['resin'].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, internal surface 2, triax".format(
                    self.internal_surface_2.layer['triax'].mass_fraction())
                print "  {0:5.1%} mass/length, internal surface 2, resin".format(
                    self.internal_surface_2.layer['resin'].mass_fraction())
        if self.internal_surface_3.exists():
            d['internal surface 3 (triax)'] = self.internal_surface_3.layer['triax'].mass_fraction()
            d['internal surface 3 (resin)'] = self.internal_surface_3.layer['resin'].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, internal surface 3, triax".format(
                    self.internal_surface_3.layer['triax'].mass_fraction())
                print "  {0:5.1%} mass/length, internal surface 3, resin".format(
                    self.internal_surface_3.layer['resin'].mass_fraction())
        if self.internal_surface_4.exists():
            d['internal surface 4 (triax)'] = self.internal_surface_4.layer['triax'].mass_fraction()
            d['internal surface 4 (resin)'] = self.internal_surface_4.layer['resin'].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, internal surface 4, triax".format(
                    self.internal_surface_4.layer['triax'].mass_fraction())
                print "  {0:5.1%} mass/length, internal surface 4, resin".format(
                    self.internal_surface_4.layer['resin'].mass_fraction())
        return d

    def write_all_part_polygons(self):
        """Write the coordinates of all structural parts to `station_path`s."""
        if self.external_surface.exists():
            self.external_surface.layer['gelcoat'].write_polygon_edges()
            self.external_surface.layer['triax'].write_polygon_edges()
        if self.root_buildup.exists():
            self.root_buildup.layer['triax'].write_polygon_edges()
        if self.spar_cap.exists():
            self.spar_cap.layer['lower'].write_polygon_edges()
            self.spar_cap.layer['upper'].write_polygon_edges()
        if self.aft_panel_1.exists():
            self.aft_panel_1.layer['lower'].write_polygon_edges()
            self.aft_panel_1.layer['upper'].write_polygon_edges()
        if self.aft_panel_2.exists():
            self.aft_panel_2.layer['lower'].write_polygon_edges()
            self.aft_panel_2.layer['upper'].write_polygon_edges()
        if self.LE_panel.exists():
            self.LE_panel.layer['foam'].write_polygon_edges()
        if self.shear_web_1.exists():
            self.shear_web_1.layer['biax, left'].write_polygon_edges()
            self.shear_web_1.layer['foam'].write_polygon_edges()
            self.shear_web_1.layer['biax, right'].write_polygon_edges()
        if self.shear_web_2.exists():
            self.shear_web_2.layer['biax, left'].write_polygon_edges()
            self.shear_web_2.layer['foam'].write_polygon_edges()
            self.shear_web_2.layer['biax, right'].write_polygon_edges()
        if self.shear_web_3.exists():
            self.shear_web_3.layer['biax, left'].write_polygon_edges()
            self.shear_web_3.layer['foam'].write_polygon_edges()
            self.shear_web_3.layer['biax, right'].write_polygon_edges()
        if self.TE_reinforcement.exists():
            self.TE_reinforcement.layer['uniax'].write_polygon_edges()
            try:
                self.TE_reinforcement.layer['foam'].write_polygon_edges()
            except KeyError:
                # the foam layer doesn't exist
                pass
        if self.internal_surface_1.exists():
            self.internal_surface_1.layer['triax'].write_polygon_edges()
            self.internal_surface_1.layer['resin'].write_polygon_edges()
        if self.internal_surface_2.exists():
            self.internal_surface_2.layer['triax'].write_polygon_edges()
            self.internal_surface_2.layer['resin'].write_polygon_edges()
        if self.internal_surface_3.exists():
            self.internal_surface_3.layer['triax'].write_polygon_edges()
            self.internal_surface_3.layer['resin'].write_polygon_edges()
        if self.internal_surface_4.exists():
            self.internal_surface_4.layer['triax'].write_polygon_edges()
            self.internal_surface_4.layer['resin'].write_polygon_edges()

    def save_all_layer_edges(self):
        """Save all layer edges as layer attributes.

        Identifies and saves the left, top, right, and bottom (LTRB) edges for
        each layer.

        This method saves LTRB edges as attributes within each layer object.

        <structure>.<part>.<layer>.left : np.array, coords for left edge
        <structure>.<part>.<layer>.top : np.array, coords for top edge
        <structure>.<part>.<layer>.right : np.array, coords for right edge
        <structure>.<part>.<layer>.bottom : np.array, coords for bottom edge

        Note: External and internal surfaces have not yet been implemented!

        """
        if self.LE_panel.exists():
            self.LE_panel.layer['foam'].get_and_save_edges()
        if self.spar_cap.exists():
            self.spar_cap.layer['lower'].get_and_save_edges()
            self.spar_cap.layer['upper'].get_and_save_edges()
        if self.aft_panel_1.exists():
            self.aft_panel_1.layer['lower'].get_and_save_edges()
            self.aft_panel_1.layer['upper'].get_and_save_edges()
        if self.aft_panel_2.exists():
            self.aft_panel_2.layer['lower'].get_and_save_edges()
            self.aft_panel_2.layer['upper'].get_and_save_edges()
        if self.shear_web_1.exists():
            self.shear_web_1.layer['biax, left'].get_and_save_edges()
            self.shear_web_1.layer['foam'].get_and_save_edges()
            self.shear_web_1.layer['biax, right'].get_and_save_edges()
        if self.shear_web_2.exists():
            self.shear_web_2.layer['biax, left'].get_and_save_edges()
            self.shear_web_2.layer['foam'].get_and_save_edges()
            self.shear_web_2.layer['biax, right'].get_and_save_edges()
        if self.shear_web_3.exists():
            self.shear_web_3.layer['biax, left'].get_and_save_edges()
            self.shear_web_3.layer['foam'].get_and_save_edges()
            self.shear_web_3.layer['biax, right'].get_and_save_edges()
        if self.TE_reinforcement.exists():
            self.TE_reinforcement.layer['uniax'].get_and_save_edges()
            try:
                self.TE_reinforcement.layer['foam'].get_and_save_edges()
            except KeyError:  # foam layer doesn't exist
                pass

    def save_all_alternate_layer_edges(self):
        """Save all alternate layer edges as layer attributes.

        Identifies and saves the left, top, right, and bottom (LTRB) edges for
        each layer.

        This method saves LTRB edges as attributes within each layer object.

        <structure>.<part>.<layer>.left : np.array, coords for left edge
        <structure>.<part>.<layer>.top : np.array, coords for top edge
        <structure>.<part>.<layer>.right : np.array, coords for right edge
        <structure>.<part>.<layer>.bottom : np.array, coords for bottom edge

        Note: External and internal surfaces have not yet been implemented!

        """
        if self.root_buildup.exists():
            self.root_buildup.layer['triax, lower left'].get_and_save_edges()
            self.root_buildup.layer['triax, lower right'].get_and_save_edges()
            self.root_buildup.layer['triax, upper right'].get_and_save_edges()
            self.root_buildup.layer['triax, upper left'].get_and_save_edges()
        if (self.TE_reinforcement.exists() and
            self.parent_station.airfoil.has_sharp_TE):
            te = self.TE_reinforcement
            te.layer['uniax, upper left'].get_and_save_edges()
            te.layer['uniax, upper right'].get_and_save_edges()
            te.layer['uniax, lower left'].get_and_save_edges()
            te.layer['uniax, lower right'].get_and_save_edges()
            te.layer['uniax, upper middle'].get_and_save_edges()
            te.layer['uniax, lower middle'].get_and_save_edges()
            te.layer['foam, upper left'].get_and_save_edges()
            te.layer['foam, upper middle'].get_and_save_edges()
            te.layer['foam, lower left'].get_and_save_edges()
            te.layer['foam, lower middle'].get_and_save_edges()
        # if self.external_surface.exists():
        #     es = self.external_surface
        #     if not self.parent_station.airfoil.has_sharp_TE:
        #         es.layer['gelcoat, lower left'].get_and_save_edges()
        #         es.layer['gelcoat, lower right'].get_and_save_edges()
        #         es.layer['gelcoat, upper right'].get_and_save_edges()
        #         es.layer['gelcoat, upper left'].get_and_save_edges()
        #         es.layer['triax, lower left'].get_and_save_edges()
        #         es.layer['triax, lower right'].get_and_save_edges()
        #         es.layer['triax, upper right'].get_and_save_edges()
        #         es.layer['triax, upper left'].get_and_save_edges()
        # if self.internal_surface_1.exists():
        #     # do something...

    def write_truegrid_header(self, outputfile_type='abaqus'):
        """Create a TrueGrid input file and write the header.

        This file is formatted as a TrueGrid input file (*.tg).

        """
        stn = self.parent_station
        b = stn.parent_blade
        separator = "c " + "-"*40 + "\n"
        f = open(os.path.join(stn.station_path, self.truegrid_input_filename),
            'w')
        f.write("c {0} ".format(b.name) + "-"*40 + "\n")
        f.write("c Station #{0:02d}\n".format(stn.station_num))
        f.write("\n")
        f.write(separator)
        f.write("\n")
        f.write("c set the name of the mesh output file\n")
        f.write("mof stn{0:02d}/mesh_stn{0:02d}.abq\n".format(stn.station_num))
        f.write("\n")
        f.write("c set the output file type\n")
        f.write(outputfile_type + "\n")
        f.write("\n")
        f.write("c set some {0}-specific parameters\n".format(outputfile_type))
        f.write("c   assign a title to the problem\n")
        f.write("title blade station #{0}\n".format(stn.station_num))
        f.write("c   build 2nd-order quadrilateral elements (with mid-side nodes)\n")
        f.write("quadratic\n")
        f.write("\n")
        f.write("c write symbolic 'pbs' commands instead of absolute 'pb' commands\n")
        f.write("cooref symbolic\n")
        f.write("\n")
        f.write(separator)
        f.write("\n")
        f.write("c set some global parameters\n")
        f.write("para t_elem 4;      c num of elements across laminate thickness\n")
        f.write("para l_elemLE 120;  c num of elements along leading edge\n")
        f.write("para l_elemSP 30;   c num of elements along spar cap\n")
        f.write("para l_elemAP 30;   c num of elements along aft panel\n")
        f.write("para l_elemTE 80;   c num of elements along trailing edge\n")
        f.write("\n")
        f.write(separator)
        f.write("\n")
        f.close()

    def write_truegrid_footer(self, interrupt_flag=False):
        """Write the footer for the TrueGrid input file."""
        stn = self.parent_station
        f = open(os.path.join(stn.station_path, self.truegrid_input_filename),
            'a')
        if interrupt_flag:
            f.write("interrupt\n")
        f.write("c merge all the individual block meshes into a single mesh\n")
        f.write("merge\n")
        f.write("c display all 3D curves\n")
        f.write("dacd\n")
        f.write("c display numbers of defined 3D curves\n")
        f.write("labels crv\n")
        f.write("c display the mesh (wireframe)\n")
        f.write("disp\n")
        f.write("\n")
        if interrupt_flag:
            f.write("interrupt\n")
        f.write("c delete redundant nodes at part boundaries\n")
        f.write("stp 0.0001\n")
        f.write("c the minimum space between nodes at biax plies is(0.836-0.833)/8 = 0.000375\n")
        f.write("c   choose a tolerance smaller than this (0.0001)\n")
        f.write("c   the command stp deletes nodes who are near each other (w/in the tolerance)\n")
        f.write("c display the mesh (filled)\n")
        f.write("tvv\n")
        f.write("\n")
        if interrupt_flag:
            f.write("interrupt\n")
        f.write("c write the mesh to an output file\n")
        f.write("write\n")
        f.write("\n")
        f.write("c exit\n")
        f.close()

    def write_all_layer_edges(self):
        """Write the coordinates of all layer edges to `station_path`.

        This file is formatted as a TrueGrid input file (*.tg).

        Saves self._dict_of_edge_nums, a dictionary of edge names for curve ID
        numbers. For example, some entries might look like:

        _dict_of_edge_nums['TE_Reinforcement; uniax; right'] = 5
        _dict_of_edge_nums['LE_Panel; foam; left'] = 10
        _dict_of_edge_nums['ShearWeb1; biax, left; top'] = 15
        _dict_of_edge_nums['AftPanel2; foam; bottom'] = 19
        _dict_of_edge_nums['SparCap; uniax; right'] = 22

        """
        stn = self.parent_station
        start_edge_num = 1
        f = open(os.path.join(stn.station_path, self.truegrid_input_filename),
            'a')
        # Procedure for each structural part:
        # 1. write edges for a layer
        # 2. increment start_edge_num by 4 edges (left, bottom, top, right)
        # 3. append d to self._dict_of_edge_nums
        if self.root_buildup.exists():
            f.write("c root buildup " + "-"*20 + "\n")
            # remove the unnecessary 'triax' layer, which is not used
            #   for meshing
            # make a new copy of the dictionary, so we don't mutate the
            #   original 'layer' dictionary
            rb_dict = self.root_buildup.layer.copy()
            rb_dict.pop('triax')
            for (layer_name, layer_obj) in rb_dict.items():
                d = layer_obj.write_layer_edges(f, start_edge_num)
                start_edge_num += 4
                self._dict_of_edge_nums.update(d)
        if self.spar_cap.exists():
            f.write("c spar cap " + "-"*20 + "\n")
            for (layer_name, layer_obj) in self.spar_cap.layer.items():
                d = layer_obj.write_layer_edges(f, start_edge_num)
                start_edge_num += 4
                self._dict_of_edge_nums.update(d)
        if self.aft_panel_1.exists():
            f.write("c aft panel #1 " + "-"*20 + "\n")
            for (layer_name, layer_obj) in self.aft_panel_1.layer.items():
                d = layer_obj.write_layer_edges(f, start_edge_num)
                start_edge_num += 4
                self._dict_of_edge_nums.update(d)
        if self.aft_panel_2.exists():
            f.write("c aft panel #2 " + "-"*20 + "\n")
            for (layer_name, layer_obj) in self.aft_panel_2.layer.items():
                d = layer_obj.write_layer_edges(f, start_edge_num)
                start_edge_num += 4
                self._dict_of_edge_nums.update(d)
        if self.LE_panel.exists():
            f.write("c LE panel " + "-"*20 + "\n")
            for (layer_name, layer_obj) in self.LE_panel.layer.items():
                d = layer_obj.write_layer_edges(f, start_edge_num)
                start_edge_num += 4
                self._dict_of_edge_nums.update(d)
        if self.shear_web_1.exists():
            f.write("c shear web #1 " + "-"*20 + "\n")
            for (layer_name, layer_obj) in self.shear_web_1.layer.items():
                d = layer_obj.write_layer_edges(f, start_edge_num)
                start_edge_num += 4
                self._dict_of_edge_nums.update(d)
        if self.shear_web_2.exists():
            f.write("c shear web #2 " + "-"*20 + "\n")
            for (layer_name, layer_obj) in self.shear_web_2.layer.items():
                d = layer_obj.write_layer_edges(f, start_edge_num)
                start_edge_num += 4
                self._dict_of_edge_nums.update(d)
        if self.shear_web_3.exists():
            f.write("c shear web #3 " + "-"*20 + "\n")
            for (layer_name, layer_obj) in self.shear_web_3.layer.items():
                d = layer_obj.write_layer_edges(f, start_edge_num)
                start_edge_num += 4
                self._dict_of_edge_nums.update(d)
        if self.TE_reinforcement.exists():
            f.write("c TE reinforcement " + "-"*20 + "\n")
            te_dict = self.TE_reinforcement.layer.copy()
            # make a new copy of the dictionary, so we don't mutate the
            #   original 'layer' dictionary
            if self.parent_station.airfoil.has_sharp_TE:
                # remove the unnecessary 'uniax' and 'foam' layers, which are
                #   not used for meshing
                te_dict.pop('uniax')
                te_dict.pop('foam')
            for (layer_name, layer_obj) in te_dict.items():
                if layer_obj.right is None:
                    # this is a triangular region
                    d = layer_obj.write_layer_edges(f, start_edge_num,
                        triangular_region=True)
                else:
                    d = layer_obj.write_layer_edges(f, start_edge_num)
                start_edge_num += len(d)
                self._dict_of_edge_nums.update(d)
        f.close()

    def write_all_alt_layer_edges(self, alt_TE_reinforcement=False,
        soft_warning=False):
        """Write the coordinates of all layer edges to `station_path`.

        This file is formatted as a TrueGrid input file (*.tg).

        """
        stn = self.parent_station
        start_edge_num = 1
        f = open(os.path.join(stn.station_path, self.truegrid_input_filename),
            'a')
        if self.root_buildup.exists():
            f.write("c root buildup " + "-"*40 + "\n")
            sd = sorted(self.root_buildup.alt_layer.items())
            for (layer_name, layer_obj) in sd:
                f.write("c " + layer_name + " " + "-"*5 + "\n")
                layer_obj.write_alt_layer_edges2(f, start_edge_num)
                if len(layer_obj.edges) > 4:
                    fmt = "More than 4 edges found in layer 'RootBuildup; {0}'!"
                    if soft_warning:
                        print "*** Warning: " + fmt.format(layer_name)
                    else:
                        raise Warning(fmt.format(layer_name))
                elif len(layer_obj.edges) == 4:
                    pass
                elif len(layer_obj.edges) == 3:
                    fmt2 = "Only 3 edges found in layer 'RootBuildup; {0}'!"
                    print "*** Warning: " + fmt2.format(layer_name)
                else:
                    fmt3 = "Layer 'RootBuildup; {0}' does not have 3 or 4 edges!"
                    raise Warning(fmt3.format(layer_name))
                start_edge_num += 4
        if self.external_surface.exists():
            f.write("c external surface " + "-"*40 + "\n")
            sd = sorted(self.external_surface.alt_layer.items())
            for (layer_name, layer_obj) in sd:
                f.write("c " + layer_name + " " + "-"*5 + "\n")
                layer_obj.write_alt_layer_edges2(f, start_edge_num)
                if len(layer_obj.edges) > 4:
                    fmt = "More than 4 edges found in layer 'ExternalSurface; {0}'!"
                    if soft_warning:
                        print "*** Warning: " + fmt.format(layer_name)
                    else:
                        raise Warning(fmt.format(layer_name))
                elif len(layer_obj.edges) == 4:
                    pass
                elif len(layer_obj.edges) == 3:
                    fmt2 = "Only 3 edges found in layer 'ExternalSurface; {0}'!"
                    print "*** Warning: " + fmt2.format(layer_name)
                else:
                    fmt3 = "Layer 'ExternalSurface; {0}' does not have 3 or 4 edges!"
                    raise Warning(fmt3.format(layer_name))
                start_edge_num += 4
        if self.internal_surface_1.exists():
            f.write("c internal surface 1 " + "-"*40 + "\n")
            sd = sorted(self.internal_surface_1.alt_layer.items())
            for (layer_name, layer_obj) in sd:
                f.write("c " + layer_name + " " + "-"*5 + "\n")
                layer_obj.write_alt_layer_edges2(f, start_edge_num)
                if len(layer_obj.edges) > 4:
                    fmt = "More than 4 edges found in layer 'InternalSurface1; {0}'!"
                    if soft_warning:
                        print "*** Warning: " + fmt.format(layer_name)
                    else:
                        raise Warning(fmt.format(layer_name))
                elif len(layer_obj.edges) == 4:
                    pass
                elif len(layer_obj.edges) == 3:
                    fmt2 = "Only 3 edges found in layer 'InternalSurface1; {0}'!"
                    print "*** Warning: " + fmt2.format(layer_name)
                else:
                    fmt3 = "Layer 'InternalSurface1; {0}' does not have 3 or 4 edges!"
                    raise Warning(fmt3.format(layer_name))
                start_edge_num += 4
        if self.internal_surface_2.exists():
            f.write("c internal surface 2 " + "-"*40 + "\n")
            sd = sorted(self.internal_surface_2.alt_layer.items())
            for (layer_name, layer_obj) in sd:
                f.write("c " + layer_name + " " + "-"*5 + "\n")
                layer_obj.write_alt_layer_edges2(f, start_edge_num)
                if len(layer_obj.edges) > 4:
                    fmt = "More than 4 edges found in layer 'InternalSurface2; {0}'!"
                    if soft_warning:
                        print "*** Warning: " + fmt.format(layer_name)
                    else:
                        raise Warning(fmt.format(layer_name))
                elif len(layer_obj.edges) == 4:
                    pass
                elif len(layer_obj.edges) == 3:
                    fmt2 = "Only 3 edges found in layer 'InternalSurface2; {0}'!"
                    print "*** Warning: " + fmt2.format(layer_name)
                else:
                    fmt3 = "Layer 'InternalSurface2; {0}' does not have 3 or 4 edges!"
                    raise Warning(fmt3.format(layer_name))
                start_edge_num += 4
        if self.internal_surface_3.exists():
            f.write("c internal surface 3 " + "-"*40 + "\n")
            sd = sorted(self.internal_surface_3.alt_layer.items())
            for (layer_name, layer_obj) in sd:
                f.write("c " + layer_name + " " + "-"*5 + "\n")
                layer_obj.write_alt_layer_edges2(f, start_edge_num)
                if len(layer_obj.edges) > 4:
                    fmt = "More than 4 edges found in layer 'InternalSurface3; {0}'!"
                    if soft_warning:
                        print "*** Warning: " + fmt.format(layer_name)
                    else:
                        raise Warning(fmt.format(layer_name))
                elif len(layer_obj.edges) == 4:
                    pass
                elif len(layer_obj.edges) == 3:
                    fmt2 = "Only 3 edges found in layer 'InternalSurface3; {0}'!"
                    print "*** Warning: " + fmt2.format(layer_name)
                else:
                    fmt3 = "Layer 'InternalSurface3; {0}' does not have 3 or 4 edges!"
                    raise Warning(fmt3.format(layer_name))
                start_edge_num += 4
        if self.internal_surface_4.exists():
            f.write("c internal surface 4 " + "-"*40 + "\n")
            sd = sorted(self.internal_surface_4.alt_layer.items())
            for (layer_name, layer_obj) in sd:
                f.write("c " + layer_name + " " + "-"*5 + "\n")
                layer_obj.write_alt_layer_edges2(f, start_edge_num)
                if len(layer_obj.edges) > 4:
                    fmt = "More than 4 edges found in layer 'InternalSurface4; {0}'!"
                    if soft_warning:
                        print "*** Warning: " + fmt.format(layer_name)
                    else:
                        raise Warning(fmt.format(layer_name))
                elif len(layer_obj.edges) == 4:
                    pass
                elif len(layer_obj.edges) == 3:
                    fmt2 = "Only 3 edges found in layer 'InternalSurface4; {0}'!"
                    print "*** Warning: " + fmt2.format(layer_name)
                else:
                    fmt3 = "Layer 'InternalSurface4; {0}' does not have 3 or 4 edges!"
                    raise Warning(fmt3.format(layer_name))
                start_edge_num += 4
        if self.TE_reinforcement.exists() and alt_TE_reinforcement:
            f.write("c TE reinforcement " + "-"*40 + "\n")
            sd = sorted(self.TE_reinforcement.alt_layer.items())
            for (layer_name, layer_obj) in sd:
                f.write("c " + layer_name + " " + "-"*5 + "\n")
                layer_obj.write_alt_layer_edges2(f, start_edge_num)
                if len(layer_obj.edges) > 4:
                    fmt = "More than 4 edges found in layer 'TE_Reinforcement; {0}'!"
                    if soft_warning:
                        print "*** Warning: " + fmt.format(layer_name)
                    else:
                        raise Warning(fmt.format(layer_name))
                elif len(layer_obj.edges) == 4:
                    pass
                elif len(layer_obj.edges) == 3:
                    fmt2 = "Only 3 edges found in layer 'TE_Reinforcement; {0}'!"
                    print "*** Warning: " + fmt2.format(layer_name)
                else:
                    fmt3 = "Layer 'TE_Reinforcement; {0}' does not have 3 or 4 edges!"
                    raise Warning(fmt3.format(layer_name))
                start_edge_num += 4
        # if self.spar_cap.exists():
        #     f.write("c spar cap " + "-"*20 + "\n")
        #     for (layer_name, layer_obj) in self.spar_cap.layer.items():
        #         d = layer_obj.write_layer_edges(f, start_edge_num)
        #         start_edge_num += 4
        #         self._dict_of_edge_nums.update(d)
        # if self.aft_panel_1.exists():
        #     f.write("c aft panel #1 " + "-"*20 + "\n")
        #     for (layer_name, layer_obj) in self.aft_panel_1.layer.items():
        #         d = layer_obj.write_layer_edges(f, start_edge_num)
        #         start_edge_num += 4
        #         self._dict_of_edge_nums.update(d)
        # if self.aft_panel_2.exists():
        #     f.write("c aft panel #2 " + "-"*20 + "\n")
        #     for (layer_name, layer_obj) in self.aft_panel_2.layer.items():
        #         d = layer_obj.write_layer_edges(f, start_edge_num)
        #         start_edge_num += 4
        #         self._dict_of_edge_nums.update(d)
        # if self.LE_panel.exists():
        #     f.write("c LE panel " + "-"*20 + "\n")
        #     for (layer_name, layer_obj) in self.LE_panel.layer.items():
        #         d = layer_obj.write_layer_edges(f, start_edge_num)
        #         start_edge_num += 4
        #         self._dict_of_edge_nums.update(d)
        # if self.shear_web_1.exists():
        #     f.write("c shear web #1 " + "-"*20 + "\n")
        #     for (layer_name, layer_obj) in self.shear_web_1.layer.items():
        #         d = layer_obj.write_layer_edges(f, start_edge_num)
        #         start_edge_num += 4
        #         self._dict_of_edge_nums.update(d)
        # if self.shear_web_2.exists():
        #     f.write("c shear web #2 " + "-"*20 + "\n")
        #     for (layer_name, layer_obj) in self.shear_web_2.layer.items():
        #         d = layer_obj.write_layer_edges(f, start_edge_num)
        #         start_edge_num += 4
        #         self._dict_of_edge_nums.update(d)
        # if self.shear_web_3.exists():
        #     f.write("c shear web #3 " + "-"*20 + "\n")
        #     for (layer_name, layer_obj) in self.shear_web_3.layer.items():
        #         d = layer_obj.write_layer_edges(f, start_edge_num)
        #         start_edge_num += 4
        #         self._dict_of_edge_nums.update(d)
        # if self.TE_reinforcement.exists():
        #     f.write("c TE reinforcement " + "-"*20 + "\n")
        #     te_dict = self.TE_reinforcement.layer.copy()
        #     # make a new copy of the dictionary, so we don't mutate the
        #     #   original 'layer' dictionary
        #     if self.parent_station.airfoil.has_sharp_TE:
        #         # remove the unnecessary 'uniax' and 'foam' layers, which are
        #         #   not used for meshing
        #         te_dict.pop('uniax')
        #         te_dict.pop('foam')
        #     for (layer_name, layer_obj) in te_dict.items():
        #         if layer_obj.right is None:
        #             # this is a triangular region
        #             d = layer_obj.write_layer_edges(f, start_edge_num,
        #                 triangular_region=True)
        #         else:
        #             d = layer_obj.write_layer_edges(f, start_edge_num)
        #         start_edge_num += len(d)
        #         self._dict_of_edge_nums.update(d)
        f.close()
        return start_edge_num

    def write_block_mesh(self, dict_key_prefix, i_cells, j_cells):
        """Write the commands for creating ONE block mesh for ONE layer.

        This file is formatted as a Truegrid input file (*.tg).

        Usage:
        write_block_mesh(
            dict_key_prefix='RootBuildup; triax, lower left',
            i_cells=60,
            j_cells=2)

        """
        stn = self.parent_station
        d = self._dict_of_edge_nums
        f = open(os.path.join(stn.station_path, self.truegrid_input_filename),
            'a')
        f.write("c make a block mesh for: {0}\n".format(dict_key_prefix))
        # create the default block mesh
        f.write("block 1 {0}; 1 {1}; -1;\n".format(i_cells,j_cells))
        f.write("-0.1 0.1; -0.1 0.1; 0;\n")
        # attach the left edge
        f.write("c left edge (i=1, j varies)\n")
        f.write("cure 1 1 1 1 2 1 {0}\n".format(
            d['{0}; left'.format(dict_key_prefix)]))
        # attach the bottom edge
        f.write("c bottom edge (j=1, i varies)\n")
        f.write("cure 1 1 1 2 1 1 {0}\n".format(
            d['{0}; bottom'.format(dict_key_prefix)]))
        # attach the right edge
        try:
            f.write("c right edge (i=2, j varies)\n")
            f.write("cure 2 1 1 2 2 1 {0}\n".format(
                d['{0}; right'.format(dict_key_prefix)]))
        except KeyError:
            # there is no right edge (this is a triangular region)
            pass
        # attach the top edge
        f.write("c top edge (j=2, i varies)\n")
        f.write("cure 1 2 1 2 2 1 {0}\n".format(
            d['{0}; top'.format(dict_key_prefix)]))
        f.write("\n")
        f.close()

    def write_all_block_meshes(self, interrupt_flag=False):
        """Write the commands for creating all TrueGrid block meshes.

        This file is formatted as a TrueGrid input file (*.tg).

        """
        if self.root_buildup.exists():
            self.write_block_mesh(
                dict_key_prefix='RootBuildup; triax, lower left',
                i_cells=60,
                j_cells=2)
            self.write_block_mesh(
                dict_key_prefix='RootBuildup; triax, lower right',
                i_cells=60,
                j_cells=2)
            self.write_block_mesh(
                dict_key_prefix='RootBuildup; triax, upper right',
                i_cells=60,
                j_cells=2)
            self.write_block_mesh(
                dict_key_prefix='RootBuildup; triax, upper left',
                i_cells=60,
                j_cells=2)
        if self.spar_cap.exists():
            self.write_block_mesh(
                dict_key_prefix='SparCap; upper',
                i_cells=30,
                j_cells=2)
            self.write_block_mesh(
                dict_key_prefix='SparCap; lower',
                i_cells=30,
                j_cells=2)
        if self.aft_panel_1.exists():
            self.write_block_mesh(
                dict_key_prefix='AftPanel1; upper',
                i_cells=18,
                j_cells=2)
            self.write_block_mesh(
                dict_key_prefix='AftPanel1; lower',
                i_cells=18,
                j_cells=2)
        if self.aft_panel_2.exists():
            self.write_block_mesh(
                dict_key_prefix='AftPanel2; upper',
                i_cells=18,
                j_cells=2)
            self.write_block_mesh(
                dict_key_prefix='AftPanel2; lower',
                i_cells=18,
                j_cells=2)
        if self.LE_panel.exists():
            self.write_block_mesh(
                dict_key_prefix='LE_Panel; foam',
                i_cells=2,
                j_cells=60)
        if self.shear_web_1.exists():
            self.write_block_mesh(
                dict_key_prefix='ShearWeb1; biax, left',
                i_cells=2,
                j_cells=18)
            self.write_block_mesh(
                dict_key_prefix='ShearWeb1; foam',
                i_cells=2,
                j_cells=18)
            self.write_block_mesh(
                dict_key_prefix='ShearWeb1; biax, right',
                i_cells=2,
                j_cells=18)
        if self.shear_web_2.exists():
            self.write_block_mesh(
                dict_key_prefix='ShearWeb2; biax, left',
                i_cells=2,
                j_cells=18)
            self.write_block_mesh(
                dict_key_prefix='ShearWeb2; foam',
                i_cells=2,
                j_cells=18)
            self.write_block_mesh(
                dict_key_prefix='ShearWeb2; biax, right',
                i_cells=2,
                j_cells=18)
        if self.shear_web_3.exists():
            self.write_block_mesh(
                dict_key_prefix='ShearWeb3; biax, left',
                i_cells=2,
                j_cells=18)
            self.write_block_mesh(
                dict_key_prefix='ShearWeb3; foam',
                i_cells=2,
                j_cells=18)
            self.write_block_mesh(
                dict_key_prefix='ShearWeb3; biax, right',
                i_cells=2,
                j_cells=18)
        if self.TE_reinforcement.exists():
            if self.parent_station.airfoil.has_sharp_TE:
                self.write_block_mesh(
                    dict_key_prefix='TE_Reinforcement; uniax, upper left',
                    i_cells=10,
                    j_cells=2)
                self.write_block_mesh(
                    dict_key_prefix='TE_Reinforcement; uniax, upper middle',
                    i_cells=10,
                    j_cells=2)
                self.write_block_mesh(
                    dict_key_prefix='TE_Reinforcement; uniax, upper right',
                    i_cells=10,
                    j_cells=2)
                self.write_block_mesh(
                    dict_key_prefix='TE_Reinforcement; uniax, lower left',
                    i_cells=10,
                    j_cells=2)
                self.write_block_mesh(
                    dict_key_prefix='TE_Reinforcement; uniax, lower middle',
                    i_cells=10,
                    j_cells=2)
                self.write_block_mesh(
                    dict_key_prefix='TE_Reinforcement; uniax, lower right',
                    i_cells=10,
                    j_cells=2)
                self.write_block_mesh(
                    dict_key_prefix='TE_Reinforcement; foam, upper left',
                    i_cells=10,
                    j_cells=2)
                self.write_block_mesh(
                    dict_key_prefix='TE_Reinforcement; foam, upper middle',
                    i_cells=10,
                    j_cells=2)
                self.write_block_mesh(
                    dict_key_prefix='TE_Reinforcement; foam, lower left',
                    i_cells=10,
                    j_cells=2)
                self.write_block_mesh(
                    dict_key_prefix='TE_Reinforcement; foam, lower middle',
                    i_cells=10,
                    j_cells=2)
            else:
                # uniax layer
                self.write_block_mesh(
                    dict_key_prefix='TE_Reinforcement; uniax',
                    i_cells=2,
                    j_cells=30)
                # check if the foam layer exists
                d = self._dict_of_edge_nums
                if 'TE_Reinforcement; foam; left' in d.keys():
                    self.write_block_mesh(
                        dict_key_prefix='TE_Reinforcement; foam',
                        i_cells=2,
                        j_cells=30)

    def write_truegrid_inputfile(self, interrupt_flag=False,
        additional_layers=[], alt_TE_reinforcement=False, soft_warning=False):
        """Write the TrueGrid input file in `station_path`."""
        self.write_truegrid_header()
        # self.write_all_layer_edges()
        start_edge_num = self.write_all_alt_layer_edges(alt_TE_reinforcement=alt_TE_reinforcement,
            soft_warning=soft_warning)
        if len(additional_layers) > 0:
            stn = self.parent_station
            f = open(os.path.join(stn.station_path,
                self.truegrid_input_filename), 'a')
            for layer in additional_layers:
                part_name = layer.parent_part.__class__.__name__  # part name
                layer_name = layer.name  # layer name
                if part_name in ['ShearWeb', 'AftPanel', 'InternalSurface']:
                    part_num = layer.parent_part.num
                    fmt = 'c {0}{1}; {2}'.format(part_name, part_num, layer_name)
                else:
                    fmt = 'c {0}; {1}'.format(part_name, layer_name)
                fmt = fmt + '-'*40 + '\n'
                f.write(fmt.format(part_name, layer_name))
                layer.write_layer_edges(f, start_edge_num)
                start_edge_num += 4
            f.close()
        # self.write_all_block_meshes()
        self.write_truegrid_footer(interrupt_flag=interrupt_flag)
        print " Wrote TrueGrid input file for Station #{0}.".format(
            self.parent_station.station_num)


class BiplaneStructure:
    """Define the biplane laminate schedule (internal dimensions)."""
    def __init__(self, h_RB, b_SC, h_SC, b_SW1_biax, b_SW1_foam, x2_SW1,
                 b_SW2_biax, b_SW2_foam, x2_SW2, b_SW3_biax, b_SW3_foam,
                 x2_SW3, b_TE_reinf, h_TE_reinf_uniax, h_TE_reinf_foam,
                 h_LE_panel, h_aft_panel_1, h_aft_panel_2, h_int_surf_1_triax,
                 h_int_surf_1_resin, h_int_surf_2_triax, h_int_surf_2_resin,
                 h_int_surf_3_triax, h_int_surf_3_resin, h_int_surf_4_triax,
                 h_int_surf_4_resin, h_ext_surf_triax, h_ext_surf_gelcoat,
                 h_RB_u, b_SC_u, h_SC_u, b_SW1_biax_u, b_SW1_foam_u, x2_SW1_u,
                 b_SW2_biax_u, b_SW2_foam_u, x2_SW2_u, b_SW3_biax_u,
                 b_SW3_foam_u, x2_SW3_u, b_TE_reinf_u, h_TE_reinf_uniax_u,
                 h_TE_reinf_foam_u, h_LE_panel_u, h_aft_panel_1_u,
                 h_aft_panel_2_u, h_int_surf_1_triax_u, h_int_surf_1_resin_u,
                 h_int_surf_2_triax_u, h_int_surf_2_resin_u,
                 h_int_surf_3_triax_u, h_int_surf_3_resin_u,
                 h_int_surf_4_triax_u, h_int_surf_4_resin_u,
                 h_ext_surf_triax_u, h_ext_surf_gelcoat_u, parent_station):
        self.parent_station = parent_station
        self._list_of_lower_layers = []
        self._list_of_upper_layers = []
        # self._dict_of_edge_nums = {}
        self.truegrid_input_filename = 'mesh_stn{0:02d}_start.tg'.format(self.parent_station.station_num)
        self.lower_root_buildup = RootBuildup(
            parent_structure = self,
            base = np.nan,
            height = h_RB)
        self.lower_spar_cap = SparCap(
            parent_structure = self,
            base = b_SC,
            height = h_SC)
        self.lower_shear_web_1 = ShearWeb(
            num = 1,
            parent_structure = self,
            base_biax = b_SW1_biax,
            base_foam = b_SW1_foam,
            x2 = x2_SW1)
        self.lower_shear_web_2 = ShearWeb(
            num = 2,
            parent_structure = self,
            base_biax = b_SW2_biax,
            base_foam = b_SW2_foam,
            x2 = x2_SW2)
        self.lower_shear_web_3 = ShearWeb(
            num = 3,
            parent_structure = self,
            base_biax = b_SW3_biax,
            base_foam = b_SW3_foam,
            x2 = x2_SW3)
        self.lower_TE_reinforcement = TE_Reinforcement(
            parent_structure = self,
            base = b_TE_reinf,
            height_uniax = h_TE_reinf_uniax,
            height_foam = h_TE_reinf_foam)
        self.lower_LE_panel = LE_Panel(
            parent_structure = self,
            base = np.nan,
            height = h_LE_panel)
        self.lower_aft_panel_1 = AftPanel(
            num = 1,
            parent_structure = self,
            base = np.nan,
            height = h_aft_panel_1)
        self.lower_aft_panel_2 = AftPanel(
            num = 2,
            parent_structure = self,
            base = np.nan,
            height = h_aft_panel_2)
        self.lower_internal_surface_1 = InternalSurface(
            num = 1,
            parent_structure = self,
            base = np.nan,
            height_triax = h_int_surf_1_triax,
            height_resin = h_int_surf_1_resin)
        self.lower_internal_surface_2 = InternalSurface(
            num = 2,
            parent_structure = self,
            base = np.nan,
            height_triax = h_int_surf_2_triax,
            height_resin = h_int_surf_2_resin)
        self.lower_internal_surface_3 = InternalSurface(
            num = 3,
            parent_structure = self,
            base = np.nan,
            height_triax = h_int_surf_3_triax,
            height_resin = h_int_surf_3_resin)
        self.lower_internal_surface_4 = InternalSurface(
            num = 4,
            parent_structure = self,
            base = np.nan,
            height_triax = h_int_surf_4_triax,
            height_resin = h_int_surf_4_resin)
        self.lower_external_surface = ExternalSurface(
            parent_structure = self,
            base = np.nan,
            height_triax = h_ext_surf_triax,
            height_gelcoat = h_ext_surf_gelcoat)
        self.lower_area = None
        self.lower_mass = None
        self.upper_root_buildup = RootBuildup(
            parent_structure = self,
            base = np.nan,
            height = h_RB)
        self.upper_spar_cap = SparCap(
            parent_structure = self,
            base = b_SC,
            height = h_SC)
        self.upper_shear_web_1 = ShearWeb(
            num = 1,
            parent_structure = self,
            base_biax = b_SW1_biax,
            base_foam = b_SW1_foam,
            x2 = x2_SW1)
        self.upper_shear_web_2 = ShearWeb(
            num = 2,
            parent_structure = self,
            base_biax = b_SW2_biax,
            base_foam = b_SW2_foam,
            x2 = x2_SW2)
        self.upper_shear_web_3 = ShearWeb(
            num = 3,
            parent_structure = self,
            base_biax = b_SW3_biax,
            base_foam = b_SW3_foam,
            x2 = x2_SW3)
        self.upper_TE_reinforcement = TE_Reinforcement(
            parent_structure = self,
            base = b_TE_reinf,
            height_uniax = h_TE_reinf_uniax,
            height_foam = h_TE_reinf_foam)
        self.upper_LE_panel = LE_Panel(
            parent_structure = self,
            base = np.nan,
            height = h_LE_panel)
        self.upper_aft_panel_1 = AftPanel(
            num = 1,
            parent_structure = self,
            base = np.nan,
            height = h_aft_panel_1)
        self.upper_aft_panel_2 = AftPanel(
            num = 2,
            parent_structure = self,
            base = np.nan,
            height = h_aft_panel_2)
        self.upper_internal_surface_1 = InternalSurface(
            num = 1,
            parent_structure = self,
            base = np.nan,
            height_triax = h_int_surf_1_triax,
            height_resin = h_int_surf_1_resin)
        self.upper_internal_surface_2 = InternalSurface(
            num = 2,
            parent_structure = self,
            base = np.nan,
            height_triax = h_int_surf_2_triax,
            height_resin = h_int_surf_2_resin)
        self.upper_internal_surface_3 = InternalSurface(
            num = 3,
            parent_structure = self,
            base = np.nan,
            height_triax = h_int_surf_3_triax,
            height_resin = h_int_surf_3_resin)
        self.upper_internal_surface_4 = InternalSurface(
            num = 4,
            parent_structure = self,
            base = np.nan,
            height_triax = h_int_surf_4_triax,
            height_resin = h_int_surf_4_resin)
        self.upper_external_surface = ExternalSurface(
            parent_structure = self,
            base = np.nan,
            height_triax = h_ext_surf_triax,
            height_gelcoat = h_ext_surf_gelcoat)
        self.upper_area = None
        self.upper_mass = None
        self.area = None
        self.mass = None

    def __str__(self):
        """Returns a string of all the internal dimensions for this structure."""
        s = ''
        s += "****** Lower Airfoil ******\n"
        s += "  --- ROOT BUILDUP ---\n"
        s += "  " + str(self.lower_root_buildup) + '\n'
        s += "  --- SPAR CAP ---\n"
        s += "  " + str(self.lower_spar_cap) + '\n'
        s += "  --- SHEAR WEB 1 ---\n"
        s += "  " + str(self.lower_shear_web_1) + '\n'
        s += "  --- SHEAR WEB 2 ---\n"
        s += "  " + str(self.lower_shear_web_2) + '\n'
        s += "  --- SHEAR WEB 3 ---\n"
        s += "  " + str(self.lower_shear_web_3) + '\n'
        s += "  ---TE REINFORCEMENT ---\n"
        s += "  " + str(self.lower_TE_reinforcement) + '\n'
        s += "  --- LE PANEL ---\n"
        s += "  " + str(self.lower_LE_panel) + '\n'
        s += "  --- AFT PANEL 1 ---\n"
        s += "  " + str(self.lower_aft_panel_1) + '\n'
        s += "  --- AFT PANEL 2 ---\n"
        s += "  " + str(self.lower_aft_panel_2) + '\n'
        s += "  --- INTERNAL SURFACE 1 ---\n"
        s += "  " + str(self.lower_internal_surface_1) + '\n'
        s += "  --- INTERNAL SURFACE 2 ---\n"
        s += "  " + str(self.lower_internal_surface_2) + '\n'
        s += "  --- INTERNAL SURFACE 3 ---\n"
        s += "  " + str(self.lower_internal_surface_3) + '\n'
        s += "  --- INTERNAL SURFACE 4 ---\n"
        s += "  " + str(self.lower_internal_surface_4) + '\n'
        s += "  --- EXTERNAL SURFACE ---\n"
        s += "  " + str(self.lower_external_surface) + '\n'
        s += "****** Upper Airfoil ******\n"
        s += "  --- ROOT BUILDUP ---\n"
        s += "  " + str(self.upper_root_buildup) + '\n'
        s += "  --- SPAR CAP ---\n"
        s += "  " + str(self.upper_spar_cap) + '\n'
        s += "  --- SHEAR WEB 1 ---\n"
        s += "  " + str(self.upper_shear_web_1) + '\n'
        s += "  --- SHEAR WEB 2 ---\n"
        s += "  " + str(self.upper_shear_web_2) + '\n'
        s += "  --- SHEAR WEB 3 ---\n"
        s += "  " + str(self.upper_shear_web_3) + '\n'
        s += "  ---TE REINFORCEMENT ---\n"
        s += "  " + str(self.upper_TE_reinforcement) + '\n'
        s += "  --- LE PANEL ---\n"
        s += "  " + str(self.upper_LE_panel) + '\n'
        s += "  --- AFT PANEL 1 ---\n"
        s += "  " + str(self.upper_aft_panel_1) + '\n'
        s += "  --- AFT PANEL 2 ---\n"
        s += "  " + str(self.upper_aft_panel_2) + '\n'
        s += "  --- INTERNAL SURFACE 1 ---\n"
        s += "  " + str(self.upper_internal_surface_1) + '\n'
        s += "  --- INTERNAL SURFACE 2 ---\n"
        s += "  " + str(self.upper_internal_surface_2) + '\n'
        s += "  --- INTERNAL SURFACE 3 ---\n"
        s += "  " + str(self.upper_internal_surface_3) + '\n'
        s += "  --- INTERNAL SURFACE 4 ---\n"
        s += "  " + str(self.upper_internal_surface_4) + '\n'
        s += "  --- EXTERNAL SURFACE ---\n"
        s += "  " + str(self.upper_external_surface) + '\n'
        return s

    def which_parts_exist(self):
        """Check which structural parts exist at this biplane station.

        Returns a dictionary of booleans.

        """
        d = {'lower root buildup': self.lower_root_buildup.exists(),
             'lower spar cap': self.lower_spar_cap.exists(),
             'lower shear web 1': self.lower_shear_web_1.exists(),
             'lower shear web 2': self.lower_shear_web_2.exists(),
             'lower shear web 3': self.lower_shear_web_3.exists(),
             'lower TE reinforcement': self.lower_TE_reinforcement.exists(),
             'lower LE panel': self.lower_LE_panel.exists(),
             'lower aft panel 1': self.lower_aft_panel_1.exists(),
             'lower aft panel 2': self.lower_aft_panel_2.exists(),
             'lower internal surface': self.lower_internal_surface.exists(),
             'lower external surface': self.lower_external_surface.exists(),
             'upper root buildup': self.upper_root_buildup.exists(),
             'upper spar cap': self.upper_spar_cap.exists(),
             'upper shear web 1': self.upper_shear_web_1.exists(),
             'upper shear web 2': self.upper_shear_web_2.exists(),
             'upper shear web 3': self.upper_shear_web_3.exists(),
             'upper TE reinforcement': self.upper_TE_reinforcement.exists(),
             'upper LE panel': self.upper_LE_panel.exists(),
             'upper aft panel 1': self.upper_aft_panel_1.exists(),
             'upper aft panel 2': self.upper_aft_panel_2.exists(),
             'upper internal surface': self.upper_internal_surface.exists(),
             'upper external surface': self.upper_external_surface.exists()}
        return d

    def create_all_layers(self):
        # fix !!!
        """Create polygons for single material layers of each structural part."""
        # lower layers --------------------------------------------------------
        if self.lower_external_surface.exists():
            self.lower_external_surface.create_layers(airfoil='lower')
        if self.lower_root_buildup.exists():
            self.lower_root_buildup.create_layers(airfoil='lower')
        if self.lower_spar_cap.exists():
            self.lower_spar_cap.create_layers(airfoil='lower')
        if self.lower_aft_panel_1.exists():
            self.lower_aft_panel_1.create_layers(airfoil='lower')
        if self.lower_aft_panel_2.exists():
            self.lower_aft_panel_2.create_layers(airfoil='lower')
        if self.lower_LE_panel.exists():
            self.lower_LE_panel.create_layers(airfoil='lower')
        if self.lower_shear_web_1.exists():
            self.lower_shear_web_1.create_layers(airfoil='lower')
        if self.lower_shear_web_2.exists():
            self.lower_shear_web_2.create_layers(airfoil='lower')
        if self.lower_shear_web_3.exists():
            self.lower_shear_web_3.create_layers(airfoil='lower')
        if self.lower_TE_reinforcement.exists():
            self.lower_TE_reinforcement.create_layers(airfoil='lower')
        lmp = self.merge_all_polygons(airfoil='lower')
        if self.lower_internal_surface_1.exists():
            self.lower_internal_surface_1.create_layers(lmp, airfoil='lower')
        if self.lower_internal_surface_2.exists():
            self.lower_internal_surface_2.create_layers(lmp, airfoil='lower')
        if self.lower_internal_surface_3.exists():
            self.lower_internal_surface_3.create_layers(lmp, airfoil='lower')
        if self.lower_internal_surface_4.exists():
            self.lower_internal_surface_4.create_layers(lmp, airfoil='lower')
        # upper layers --------------------------------------------------------
        if self.upper_external_surface.exists():
            self.upper_external_surface.create_layers(airfoil='upper')
        if self.upper_root_buildup.exists():
            self.upper_root_buildup.create_layers(airfoil='upper')
        if self.upper_spar_cap.exists():
            self.upper_spar_cap.create_layers(airfoil='upper')
        if self.upper_aft_panel_1.exists():
            self.upper_aft_panel_1.create_layers(airfoil='upper')
        if self.upper_aft_panel_2.exists():
            self.upper_aft_panel_2.create_layers(airfoil='upper')
        if self.upper_LE_panel.exists():
            self.upper_LE_panel.create_layers(airfoil='upper')
        if self.upper_shear_web_1.exists():
            self.upper_shear_web_1.create_layers(airfoil='upper')
        if self.upper_shear_web_2.exists():
            self.upper_shear_web_2.create_layers(airfoil='upper')
        if self.upper_shear_web_3.exists():
            self.upper_shear_web_3.create_layers(airfoil='upper')
        if self.upper_TE_reinforcement.exists():
            self.upper_TE_reinforcement.create_layers(airfoil='upper')
        ump = self.merge_all_polygons(airfoil='upper')
        if self.upper_internal_surface_1.exists():
            self.upper_internal_surface_1.create_layers(ump, airfoil='upper')
        if self.upper_internal_surface_2.exists():
            self.upper_internal_surface_2.create_layers(ump, airfoil='upper')
        if self.upper_internal_surface_3.exists():
            self.upper_internal_surface_3.create_layers(ump, airfoil='upper')
        if self.upper_internal_surface_4.exists():
            self.upper_internal_surface_4.create_layers(ump, airfoil='upper')

    def merge_all_polygons(self, airfoil, plot_flag=False):
        """Merges all the layer polygons in this structure into one polygon.

        NOTE: internal surface polygons are NOT merged!

        """
        stn = self.parent_station
        if plot_flag:
            fig, ax = plt.subplots()
            ax.set_title("Station #{0}, {1}, {2}% span".format(stn.station_num,
                stn.airfoil.name, stn.coords.x1))
            ax.set_aspect('equal')
            ax.grid('on')
            ax.set_xlabel('x2 [meters]')
            ax.set_ylabel('x3 [meters]')
            patch = PolygonPatch(stn.airfoil.polygon, fc='None', ec='#999999',
                alpha=0.8)
            ax.add_patch(patch)
            (minx, miny, maxx, maxy) = stn.airfoil.polygon.bounds
            ax.set_xlim([minx*1.2,maxx*1.2])
            ax.set_ylim([miny*1.2,maxy*1.2])
        # merge everything
        list_of_polygons = []
        if airfoil is None:
            this_list = self._list_of_layers
        elif airfoil == 'lower':
            this_list = self._list_of_lower_layers
        elif airfoil == 'upper':
            this_list = self._list_of_upper_layers
        else:
            raise ValueError("Keyword `airfoil` must be None, 'lower', or 'upper'.")
        for layer in this_list:
            list_of_polygons.append(layer.polygon)
        try:
            p = cascaded_union(list_of_polygons)
        except ValueError:
            # gather all the parts
            if airfoil is None:
                p = self.external_surface.layer['gelcoat'].polygon
                p = p.union(self.external_surface.layer['triax'].polygon)
                if self.root_buildup.exists():
                    RB = self.root_buildup.layer['triax'].polygon
                    p = p.union(RB)
                if self.LE_panel.exists():
                    LE = self.LE_panel.layer['foam'].polygon
                    p = p.union(LE)
                if self.spar_cap.exists():
                    sc_l = self.spar_cap.layer['lower'].polygon
                    try:
                        p = p.union(sc_l)
                    except TopologicalError:
                        print " [Warning] could not merge lower spar cap in Station #{0} ... skipping!".format(self.parent_station.station_num)
                    sc_u = self.spar_cap.layer['upper'].polygon
                    try:
                        p = p.union(sc_u)
                    except TopologicalError:
                        print " [Warning] could not merge upper spar cap in Station #{0}".format(self.parent_station.station_num)
                        try:
                            p = sc_u.union(p)
                            print " ... SUCCESSFULLY merged upper spar cap on the second try!"
                        except TopologicalError:
                            print " ... on second try, still COULD NOT merge upper spar cap!"
                if self.aft_panel_1.exists():
                    aft1_u = self.aft_panel_1.layer['upper'].polygon
                    aft1_l = self.aft_panel_1.layer['lower'].polygon
                    p = p.union(aft1_u)
                    p = p.union(aft1_l)
                if self.aft_panel_2.exists():
                    aft2_u = self.aft_panel_2.layer['upper'].polygon
                    aft2_l = self.aft_panel_2.layer['lower'].polygon
                    p = p.union(aft2_u)
                    p = p.union(aft2_l)
                if self.TE_reinforcement.exists():
                    TE_uniax = self.TE_reinforcement.layer['uniax'].polygon
                    try:
                        p = p.union(TE_uniax)
                    except TopologicalError:
                        print " [Warning] could not merge uniax layer of TE reinforcement in Station #{0}".format(self.parent_station.station_num)
                        try:
                            p = TE_uniax.union(p)
                            print " ... SUCCESSFULLY merged uniax layer of TE reinforcement on the second try!"
                        except TopologicalError:
                            print " ... on second try, still COULD NOT merge uniax layer of TE reinforcement!"
                            try:
                                p = cascaded_union([p, TE_uniax])
                            except ValueError:
                                print " ... on third try, still COULD NOT merge uniax layer of TE reinforcement!"
                    try:
                        TE_foam = self.TE_reinforcement.layer['foam'].polygon
                        try:
                            p = p.union(TE_foam)
                        except TopologicalError:
                            print " [Warning] could not merge foam layer of TE reinforcement in Station #{0} ... skipping!".format(self.parent_station.station_num)
                    except ValueError:
                        print " foam layer of TE reinforcement does not exist in Station #{0}".format(self.parent_station.station_num)
                if self.shear_web_1.exists():
                    sw1 = self.shear_web_1.layer['biax, left'].polygon.union(self.shear_web_1.layer['foam'].polygon)
                    sw1 = sw1.union(self.shear_web_1.layer['biax, right'].polygon)
                    p = p.union(sw1)
                if self.shear_web_2.exists():
                    sw2 = self.shear_web_2.layer['biax, left'].polygon.union(self.shear_web_2.layer['foam'].polygon)
                    sw2 = sw2.union(self.shear_web_2.layer['biax, right'].polygon)
                    p = p.union(sw2)
                if self.shear_web_3.exists():
                    sw3 = self.shear_web_3.layer['biax, left'].polygon.union(self.shear_web_3.layer['foam'].polygon)
                    sw3 = sw3.union(self.shear_web_3.layer['biax, right'].polygon)
                    p = p.union(sw3)
            elif airfoil == 'lower':
                p = self.lower_external_surface.layer['gelcoat'].polygon
                p = p.union(self.lower_external_surface.layer['triax'].polygon)
                if self.lower_root_buildup.exists():
                    RB = self.lower_root_buildup.layer['triax'].polygon
                    p = p.union(RB)
                if self.lower_LE_panel.exists():
                    LE = self.lower_LE_panel.layer['foam'].polygon
                    p = p.union(LE)
                if self.lower_spar_cap.exists():
                    sc_l = self.lower_spar_cap.layer['lower'].polygon
                    try:
                        p = p.union(sc_l)
                    except TopologicalError:
                        print " [Warning] could not merge lower spar cap in Station #{0} ... skipping!".format(self.parent_station.station_num)
                    sc_u = self.lower_spar_cap.layer['upper'].polygon
                    try:
                        p = p.union(sc_u)
                    except TopologicalError:
                        print " [Warning] could not merge upper spar cap in Station #{0}".format(self.parent_station.station_num)
                        try:
                            p = sc_u.union(p)
                            print " ... SUCCESSFULLY merged upper spar cap on the second try!"
                        except TopologicalError:
                            print " ... on second try, still COULD NOT merge upper spar cap!"
                if self.lower_aft_panel_1.exists():
                    aft1_u = self.lower_aft_panel_1.layer['upper'].polygon
                    aft1_l = self.lower_aft_panel_1.layer['lower'].polygon
                    p = p.union(aft1_u)
                    p = p.union(aft1_l)
                if self.lower_aft_panel_2.exists():
                    aft2_u = self.lower_aft_panel_2.layer['upper'].polygon
                    aft2_l = self.lower_aft_panel_2.layer['lower'].polygon
                    p = p.union(aft2_u)
                    p = p.union(aft2_l)
                if self.lower_TE_reinforcement.exists():
                    TE_uniax = self.lower_TE_reinforcement.layer['uniax'].polygon
                    try:
                        p = p.union(TE_uniax)
                    except TopologicalError:
                        print " [Warning] could not merge uniax layer of TE reinforcement in Station #{0}".format(self.parent_station.station_num)
                        try:
                            p = TE_uniax.union(p)
                            print " ... SUCCESSFULLY merged uniax layer of TE reinforcement on the second try!"
                        except TopologicalError:
                            print " ... on second try, still COULD NOT merge uniax layer of TE reinforcement!"
                            try:
                                p = cascaded_union([p, TE_uniax])
                            except ValueError:
                                print " ... on third try, still COULD NOT merge uniax layer of TE reinforcement!"
                    try:
                        TE_foam = self.lower_TE_reinforcement.layer['foam'].polygon
                        try:
                            p = p.union(TE_foam)
                        except TopologicalError:
                            print " [Warning] could not merge foam layer of TE reinforcement in Station #{0} ... skipping!".format(self.parent_station.station_num)
                    except ValueError:
                        print " foam layer of TE reinforcement does not exist in Station #{0}".format(self.parent_station.station_num)
                if self.lower_shear_web_1.exists():
                    sw1 = self.lower_shear_web_1.layer['biax, left'].polygon.union(self.lower_shear_web_1.layer['foam'].polygon)
                    sw1 = sw1.union(self.lower_shear_web_1.layer['biax, right'].polygon)
                    p = p.union(sw1)
                if self.lower_shear_web_2.exists():
                    sw2 = self.lower_shear_web_2.layer['biax, left'].polygon.union(self.lower_shear_web_2.layer['foam'].polygon)
                    sw2 = sw2.union(self.lower_shear_web_2.layer['biax, right'].polygon)
                    p = p.union(sw2)
                if self.lower_shear_web_3.exists():
                    sw3 = self.lower_shear_web_3.layer['biax, left'].polygon.union(self.lower_shear_web_3.layer['foam'].polygon)
                    sw3 = sw3.union(self.lower_shear_web_3.layer['biax, right'].polygon)
                    p = p.union(sw3)
            elif airfoil == 'upper':
                p = self.upper_external_surface.layer['gelcoat'].polygon
                p = p.union(self.upper_external_surface.layer['triax'].polygon)
                if self.upper_root_buildup.exists():
                    RB = self.upper_root_buildup.layer['triax'].polygon
                    p = p.union(RB)
                if self.upper_LE_panel.exists():
                    LE = self.upper_LE_panel.layer['foam'].polygon
                    p = p.union(LE)
                if self.upper_spar_cap.exists():
                    sc_l = self.upper_spar_cap.layer['lower'].polygon
                    try:
                        p = p.union(sc_l)
                    except TopologicalError:
                        print " [Warning] could not merge lower spar cap in Station #{0} ... skipping!".format(self.parent_station.station_num)
                    sc_u = self.upper_spar_cap.layer['upper'].polygon
                    try:
                        p = p.union(sc_u)
                    except TopologicalError:
                        print " [Warning] could not merge upper spar cap in Station #{0}".format(self.parent_station.station_num)
                        try:
                            p = sc_u.union(p)
                            print " ... SUCCESSFULLY merged upper spar cap on the second try!"
                        except TopologicalError:
                            print " ... on second try, still COULD NOT merge upper spar cap!"
                if self.upper_aft_panel_1.exists():
                    aft1_u = self.upper_aft_panel_1.layer['upper'].polygon
                    aft1_l = self.upper_aft_panel_1.layer['lower'].polygon
                    p = p.union(aft1_u)
                    p = p.union(aft1_l)
                if self.upper_aft_panel_2.exists():
                    aft2_u = self.upper_aft_panel_2.layer['upper'].polygon
                    aft2_l = self.upper_aft_panel_2.layer['lower'].polygon
                    p = p.union(aft2_u)
                    p = p.union(aft2_l)
                if self.upper_TE_reinforcement.exists():
                    TE_uniax = self.upper_TE_reinforcement.layer['uniax'].polygon
                    try:
                        p = p.union(TE_uniax)
                    except TopologicalError:
                        print " [Warning] could not merge uniax layer of TE reinforcement in Station #{0}".format(self.parent_station.station_num)
                        try:
                            p = TE_uniax.union(p)
                            print " ... SUCCESSFULLY merged uniax layer of TE reinforcement on the second try!"
                        except TopologicalError:
                            print " ... on second try, still COULD NOT merge uniax layer of TE reinforcement!"
                            try:
                                p = cascaded_union([p, TE_uniax])
                            except ValueError:
                                print " ... on third try, still COULD NOT merge uniax layer of TE reinforcement!"
                    try:
                        TE_foam = self.upper_TE_reinforcement.layer['foam'].polygon
                        try:
                            p = p.union(TE_foam)
                        except TopologicalError:
                            print " [Warning] could not merge foam layer of TE reinforcement in Station #{0} ... skipping!".format(self.parent_station.station_num)
                    except ValueError:
                        print " foam layer of TE reinforcement does not exist in Station #{0}".format(self.parent_station.station_num)
                if self.upper_shear_web_1.exists():
                    sw1 = self.upper_shear_web_1.layer['biax, left'].polygon.union(self.upper_shear_web_1.layer['foam'].polygon)
                    sw1 = sw1.union(self.upper_shear_web_1.layer['biax, right'].polygon)
                    p = p.union(sw1)
                if self.upper_shear_web_2.exists():
                    sw2 = self.upper_shear_web_2.layer['biax, left'].polygon.union(self.upper_shear_web_2.layer['foam'].polygon)
                    sw2 = sw2.union(self.upper_shear_web_2.layer['biax, right'].polygon)
                    p = p.union(sw2)
                if self.upper_shear_web_3.exists():
                    sw3 = self.upper_shear_web_3.layer['biax, left'].polygon.union(self.upper_shear_web_3.layer['foam'].polygon)
                    sw3 = sw3.union(self.upper_shear_web_3.layer['biax, right'].polygon)
                    p = p.union(sw3)
        if plot_flag:
            # plot the merged polygon
            patch2 = PolygonPatch(p, fc='#4000FF', ec = '#000000', alpha=0.8)
            ax.add_patch(patch2)
            plt.show()
        return p

    def write_all_part_polygons(self):
        """Write the coordinates of all structural parts to `station_path`s."""
        # lower layers --------------------------------------------------------
        if self.lower_external_surface.exists():
            self.lower_external_surface.layer['gelcoat'].write_polygon_edges(airfoil='lower')
            self.lower_external_surface.layer['triax'].write_polygon_edges(airfoil='lower')
        if self.lower_root_buildup.exists():
            self.lower_root_buildup.layer['triax'].write_polygon_edges(airfoil='lower')
        if self.lower_spar_cap.exists():
            self.lower_spar_cap.layer['lower'].write_polygon_edges(airfoil='lower')
            self.lower_spar_cap.layer['upper'].write_polygon_edges(airfoil='lower')
        if self.lower_aft_panel_1.exists():
            self.lower_aft_panel_1.layer['lower'].write_polygon_edges(airfoil='lower')
            self.lower_aft_panel_1.layer['upper'].write_polygon_edges(airfoil='lower')
        if self.lower_aft_panel_2.exists():
            self.lower_aft_panel_2.layer['lower'].write_polygon_edges(airfoil='lower')
            self.lower_aft_panel_2.layer['upper'].write_polygon_edges(airfoil='lower')
        if self.lower_LE_panel.exists():
            self.lower_LE_panel.layer['foam'].write_polygon_edges(airfoil='lower')
        if self.lower_shear_web_1.exists():
            self.lower_shear_web_1.layer['biax, left'].write_polygon_edges(airfoil='lower')
            self.lower_shear_web_1.layer['foam'].write_polygon_edges(airfoil='lower')
            self.lower_shear_web_1.layer['biax, right'].write_polygon_edges(airfoil='lower')
        if self.lower_shear_web_2.exists():
            self.lower_shear_web_2.layer['biax, left'].write_polygon_edges(airfoil='lower')
            self.lower_shear_web_2.layer['foam'].write_polygon_edges(airfoil='lower')
            self.lower_shear_web_2.layer['biax, right'].write_polygon_edges(airfoil='lower')
        if self.lower_shear_web_3.exists():
            self.lower_shear_web_3.layer['biax, left'].write_polygon_edges(airfoil='lower')
            self.lower_shear_web_3.layer['foam'].write_polygon_edges(airfoil='lower')
            self.lower_shear_web_3.layer['biax, right'].write_polygon_edges(airfoil='lower')
        if self.lower_TE_reinforcement.exists():
            self.lower_TE_reinforcement.layer['uniax'].write_polygon_edges(airfoil='lower')
            try:
                self.lower_TE_reinforcement.layer['foam'].write_polygon_edges(airfoil='lower')
            except KeyError:
                # the foam layer doesn't exist
                pass
        if self.lower_internal_surface_1.exists():
            self.lower_internal_surface_1.layer['triax'].write_polygon_edges(airfoil='lower')
            self.lower_internal_surface_1.layer['resin'].write_polygon_edges(airfoil='lower')
        if self.lower_internal_surface_2.exists():
            self.lower_internal_surface_2.layer['triax'].write_polygon_edges(airfoil='lower')
            self.lower_internal_surface_2.layer['resin'].write_polygon_edges(airfoil='lower')
        if self.lower_internal_surface_3.exists():
            self.lower_internal_surface_3.layer['triax'].write_polygon_edges(airfoil='lower')
            self.lower_internal_surface_3.layer['resin'].write_polygon_edges(airfoil='lower')
        if self.lower_internal_surface_4.exists():
            self.lower_internal_surface_4.layer['triax'].write_polygon_edges(airfoil='lower')
            self.lower_internal_surface_4.layer['resin'].write_polygon_edges(airfoil='lower')
        # upper layers --------------------------------------------------------
        if self.upper_external_surface.exists():
            self.upper_external_surface.layer['gelcoat'].write_polygon_edges(airfoil='upper')
            self.upper_external_surface.layer['triax'].write_polygon_edges(airfoil='upper')
        if self.upper_root_buildup.exists():
            self.upper_root_buildup.layer['triax'].write_polygon_edges(airfoil='upper')
        if self.upper_spar_cap.exists():
            self.upper_spar_cap.layer['lower'].write_polygon_edges(airfoil='upper')
            self.upper_spar_cap.layer['upper'].write_polygon_edges(airfoil='upper')
        if self.upper_aft_panel_1.exists():
            self.upper_aft_panel_1.layer['lower'].write_polygon_edges(airfoil='upper')
            self.upper_aft_panel_1.layer['upper'].write_polygon_edges(airfoil='upper')
        if self.upper_aft_panel_2.exists():
            self.upper_aft_panel_2.layer['lower'].write_polygon_edges(airfoil='upper')
            self.upper_aft_panel_2.layer['upper'].write_polygon_edges(airfoil='upper')
        if self.upper_LE_panel.exists():
            self.upper_LE_panel.layer['foam'].write_polygon_edges(airfoil='upper')
        if self.upper_shear_web_1.exists():
            self.upper_shear_web_1.layer['biax, left'].write_polygon_edges(airfoil='upper')
            self.upper_shear_web_1.layer['foam'].write_polygon_edges(airfoil='upper')
            self.upper_shear_web_1.layer['biax, right'].write_polygon_edges(airfoil='upper')
        if self.upper_shear_web_2.exists():
            self.upper_shear_web_2.layer['biax, left'].write_polygon_edges(airfoil='upper')
            self.upper_shear_web_2.layer['foam'].write_polygon_edges(airfoil='upper')
            self.upper_shear_web_2.layer['biax, right'].write_polygon_edges(airfoil='upper')
        if self.upper_shear_web_3.exists():
            self.upper_shear_web_3.layer['biax, left'].write_polygon_edges(airfoil='upper')
            self.upper_shear_web_3.layer['foam'].write_polygon_edges(airfoil='upper')
            self.upper_shear_web_3.layer['biax, right'].write_polygon_edges(airfoil='upper')
        if self.upper_TE_reinforcement.exists():
            self.upper_TE_reinforcement.layer['uniax'].write_polygon_edges(airfoil='upper')
            try:
                self.upper_TE_reinforcement.layer['foam'].write_polygon_edges(airfoil='upper')
            except KeyError:
                # the foam layer doesn't exist
                pass
        if self.upper_internal_surface_1.exists():
            self.upper_internal_surface_1.layer['triax'].write_polygon_edges(airfoil='upper')
            self.upper_internal_surface_1.layer['resin'].write_polygon_edges(airfoil='upper')
        if self.upper_internal_surface_2.exists():
            self.upper_internal_surface_2.layer['triax'].write_polygon_edges(airfoil='upper')
            self.upper_internal_surface_2.layer['resin'].write_polygon_edges(airfoil='upper')
        if self.upper_internal_surface_3.exists():
            self.upper_internal_surface_3.layer['triax'].write_polygon_edges(airfoil='upper')
            self.upper_internal_surface_3.layer['resin'].write_polygon_edges(airfoil='upper')
        if self.upper_internal_surface_4.exists():
            self.upper_internal_surface_4.layer['triax'].write_polygon_edges(airfoil='upper')
            self.upper_internal_surface_4.layer['resin'].write_polygon_edges(airfoil='upper')

    def save_all_layer_edges(self):
        """Save all layer edges as layer attributes.

        Identifies and saves the left, top, right, and bottom (LTRB) edges for
        each layer.

        This method saves LTRB edges as attributes within each layer object.

        <structure>.<part>.<layer>.left : np.array, coords for left edge
        <structure>.<part>.<layer>.top : np.array, coords for top edge
        <structure>.<part>.<layer>.right : np.array, coords for right edge
        <structure>.<part>.<layer>.bottom : np.array, coords for bottom edge

        Note: External and internal surfaces have not yet been implemented!

        """
        # lower layers --------------------------------------------------------
        if self.lower_LE_panel.exists():
            self.lower_LE_panel.layer['foam'].get_and_save_edges()
        if self.lower_spar_cap.exists():
            self.lower_spar_cap.layer['lower'].get_and_save_edges()
            self.lower_spar_cap.layer['upper'].get_and_save_edges()
        if self.lower_aft_panel_1.exists():
            self.lower_aft_panel_1.layer['lower'].get_and_save_edges()
            self.lower_aft_panel_1.layer['upper'].get_and_save_edges()
        if self.lower_aft_panel_2.exists():
            self.lower_aft_panel_2.layer['lower'].get_and_save_edges()
            self.lower_aft_panel_2.layer['upper'].get_and_save_edges()
        if self.lower_shear_web_1.exists():
            self.lower_shear_web_1.layer['biax, left'].get_and_save_edges()
            self.lower_shear_web_1.layer['foam'].get_and_save_edges()
            self.lower_shear_web_1.layer['biax, right'].get_and_save_edges()
        if self.lower_shear_web_2.exists():
            self.lower_shear_web_2.layer['biax, left'].get_and_save_edges()
            self.lower_shear_web_2.layer['foam'].get_and_save_edges()
            self.lower_shear_web_2.layer['biax, right'].get_and_save_edges()
        if self.lower_shear_web_3.exists():
            self.lower_shear_web_3.layer['biax, left'].get_and_save_edges()
            self.lower_shear_web_3.layer['foam'].get_and_save_edges()
            self.lower_shear_web_3.layer['biax, right'].get_and_save_edges()
        if self.lower_TE_reinforcement.exists():
            self.lower_TE_reinforcement.layer['uniax'].get_and_save_edges()
            try:
                self.lower_TE_reinforcement.layer['foam'].get_and_save_edges()
            except KeyError:  # foam layer doesn't exist
                pass
        # upper layers --------------------------------------------------------
        if self.upper_LE_panel.exists():
            self.upper_LE_panel.layer['foam'].get_and_save_edges()
        if self.upper_spar_cap.exists():
            self.upper_spar_cap.layer['lower'].get_and_save_edges()
            self.upper_spar_cap.layer['upper'].get_and_save_edges()
        if self.upper_aft_panel_1.exists():
            self.upper_aft_panel_1.layer['lower'].get_and_save_edges()
            self.upper_aft_panel_1.layer['upper'].get_and_save_edges()
        if self.upper_aft_panel_2.exists():
            self.upper_aft_panel_2.layer['lower'].get_and_save_edges()
            self.upper_aft_panel_2.layer['upper'].get_and_save_edges()
        if self.upper_shear_web_1.exists():
            self.upper_shear_web_1.layer['biax, left'].get_and_save_edges()
            self.upper_shear_web_1.layer['foam'].get_and_save_edges()
            self.upper_shear_web_1.layer['biax, right'].get_and_save_edges()
        if self.upper_shear_web_2.exists():
            self.upper_shear_web_2.layer['biax, left'].get_and_save_edges()
            self.upper_shear_web_2.layer['foam'].get_and_save_edges()
            self.upper_shear_web_2.layer['biax, right'].get_and_save_edges()
        if self.upper_shear_web_3.exists():
            self.upper_shear_web_3.layer['biax, left'].get_and_save_edges()
            self.upper_shear_web_3.layer['foam'].get_and_save_edges()
            self.upper_shear_web_3.layer['biax, right'].get_and_save_edges()
        if self.upper_TE_reinforcement.exists():
            self.upper_TE_reinforcement.layer['uniax'].get_and_save_edges()
            try:
                self.upper_TE_reinforcement.layer['foam'].get_and_save_edges()
            except KeyError:  # foam layer doesn't exist
                pass

    def calculate_area(self):
        """Add the area of all polygons in this station."""
        a = 0
        for layer in self._list_of_lower_layers:
            a += layer.polygon.area
        for layer in self._list_of_upper_layers:
            a += layer.polygon.area
        self.area = a
        return a

    def calculate_mass(self):
        """Add the mass (per unit length) of all polygons in this station."""
        m = 0
        for layer in self._list_of_lower_layers:
            m += layer.mass
        for layer in self._list_of_upper_layers:
            m += layer.mass
        self.mass = m
        return m
        
    def calculate_all_percent_areas(self):
        """Calculate the percent areas of all parts in this station.

        Returns a dictionary of area fractions for each structural part.

        """
        # lower airfoil -------------------------------------------------------
        dl = {}  # dictionary for lower airfoils
        if self.lower_external_surface.exists():
            dl['external surface (gelcoat)'] = self.lower_external_surface.layer['gelcoat'].area_fraction()
            dl['external surface (triax)'] = self.lower_external_surface.layer['triax'].area_fraction()
        if self.lower_root_buildup.exists():
            dl['root buildup'] = self.lower_root_buildup.layer['triax'].area_fraction()
        if self.lower_spar_cap.exists():
            dl['spar cap (lower)'] = self.lower_spar_cap.layer['lower'].area_fraction()
            dl['spar cap (upper)'] = self.lower_spar_cap.layer['upper'].area_fraction()
        if self.lower_aft_panel_1.exists():
            dl['aft panel 1 (lower)'] = self.lower_aft_panel_1.layer['lower'].area_fraction()
            dl['aft panel 1 (upper)'] = self.lower_aft_panel_1.layer['upper'].area_fraction()
        if self.lower_aft_panel_2.exists():
            dl['aft panel 2 (lower)'] = self.lower_aft_panel_2.layer['lower'].area_fraction()
            dl['aft panel 2 (upper)'] = self.lower_aft_panel_2.layer['upper'].area_fraction()
        if self.lower_LE_panel.exists():
            dl['LE panel'] = self.lower_LE_panel.layer['foam'].area_fraction()
        if self.lower_shear_web_1.exists():
            dl['shear web 1 (left biax)'] = self.lower_shear_web_1.layer['biax, left'].area_fraction()
            dl['shear web 1 (foam)'] = self.lower_shear_web_1.layer['foam'].area_fraction()
            dl['shear web 1 (right biax)'] = self.lower_shear_web_1.layer['biax, right'].area_fraction()
        if self.lower_shear_web_2.exists():
            dl['shear web 2 (left biax)'] = self.lower_shear_web_2.layer['biax, left'].area_fraction()
            dl['shear web 2 (foam)'] = self.lower_shear_web_2.layer['foam'].area_fraction()
            dl['shear web 2 (right biax)'] = self.lower_shear_web_2.layer['biax, right'].area_fraction()
        if self.lower_shear_web_3.exists():
            dl['shear web 3 (left biax)'] = self.lower_shear_web_3.layer['biax, left'].area_fraction()
            dl['shear web 3 (foam)'] = self.lower_shear_web_3.layer['foam'].area_fraction()
            dl['shear web 3 (right biax)'] = self.lower_shear_web_3.layer['biax, right'].area_fraction()
        if self.lower_TE_reinforcement.exists():
            dl['TE reinforcement (uniax)'] = self.lower_TE_reinforcement.layer['uniax'].area_fraction()
            try:
                dl['TE reinforcement (foam)'] = self.lower_TE_reinforcement.layer['foam'].area_fraction()
            except KeyError:
                # the foam layer doesn't exist in TE reinf at this station
                dl['TE reinforcement (foam)'] = 0.0
        if self.lower_internal_surface_1.exists():
            dl['internal surface 1 (triax)'] = self.lower_internal_surface_1.layer['triax'].area_fraction()
            dl['internal surface 1 (resin)'] = self.lower_internal_surface_1.layer['resin'].area_fraction()
        if self.lower_internal_surface_2.exists():
            dl['internal surface 2 (triax)'] = self.lower_internal_surface_2.layer['triax'].area_fraction()
            dl['internal surface 2 (resin)'] = self.lower_internal_surface_2.layer['resin'].area_fraction()
        if self.lower_internal_surface_3.exists():
            dl['internal surface 3 (triax)'] = self.lower_internal_surface_3.layer['triax'].area_fraction()
            dl['internal surface 3 (resin)'] = self.lower_internal_surface_3.layer['resin'].area_fraction()
        if self.lower_internal_surface_4.exists():
            dl['internal surface 4 (triax)'] = self.lower_internal_surface_4.layer['triax'].area_fraction()
            dl['internal surface 4 (resin)'] = self.lower_internal_surface_4.layer['resin'].area_fraction()
        # upper airfoil -------------------------------------------------------
        du = {}  # dictionary for upper airfoils
        if self.lower_external_surface.exists():
            du['external surface (gelcoat)'] = self.lower_external_surface.layer['gelcoat'].area_fraction()
            du['external surface (triax)'] = self.lower_external_surface.layer['triax'].area_fraction()
        if self.lower_root_buildup.exists():
            du['root buildup'] = self.lower_root_buildup.layer['triax'].area_fraction()
        if self.lower_spar_cap.exists():
            du['spar cap (lower)'] = self.lower_spar_cap.layer['lower'].area_fraction()
            du['spar cap (upper)'] = self.lower_spar_cap.layer['upper'].area_fraction()
        if self.lower_aft_panel_1.exists():
            du['aft panel 1 (lower)'] = self.lower_aft_panel_1.layer['lower'].area_fraction()
            du['aft panel 1 (upper)'] = self.lower_aft_panel_1.layer['upper'].area_fraction()
        if self.lower_aft_panel_2.exists():
            du['aft panel 2 (lower)'] = self.lower_aft_panel_2.layer['lower'].area_fraction()
            du['aft panel 2 (upper)'] = self.lower_aft_panel_2.layer['upper'].area_fraction()
        if self.lower_LE_panel.exists():
            du['LE panel'] = self.lower_LE_panel.layer['foam'].area_fraction()
        if self.lower_shear_web_1.exists():
            du['shear web 1 (left biax)'] = self.lower_shear_web_1.layer['biax, left'].area_fraction()
            du['shear web 1 (foam)'] = self.lower_shear_web_1.layer['foam'].area_fraction()
            du['shear web 1 (right biax)'] = self.lower_shear_web_1.layer['biax, right'].area_fraction()
        if self.lower_shear_web_2.exists():
            du['shear web 2 (left biax)'] = self.lower_shear_web_2.layer['biax, left'].area_fraction()
            du['shear web 2 (foam)'] = self.lower_shear_web_2.layer['foam'].area_fraction()
            du['shear web 2 (right biax)'] = self.lower_shear_web_2.layer['biax, right'].area_fraction()
        if self.lower_shear_web_3.exists():
            du['shear web 3 (left biax)'] = self.lower_shear_web_3.layer['biax, left'].area_fraction()
            du['shear web 3 (foam)'] = self.lower_shear_web_3.layer['foam'].area_fraction()
            du['shear web 3 (right biax)'] = self.lower_shear_web_3.layer['biax, right'].area_fraction()
        if self.lower_TE_reinforcement.exists():
            du['TE reinforcement (uniax)'] = self.lower_TE_reinforcement.layer['uniax'].area_fraction()
            try:
                du['TE reinforcement (foam)'] = self.lower_TE_reinforcement.layer['foam'].area_fraction()
            except KeyError:
                # the foam layer doesn't exist in TE reinf at this station
                du['TE reinforcement (foam)'] = 0.0
        if self.lower_internal_surface_1.exists():
            du['internal surface 1 (triax)'] = self.lower_internal_surface_1.layer['triax'].area_fraction()
            du['internal surface 1 (resin)'] = self.lower_internal_surface_1.layer['resin'].area_fraction()
        if self.lower_internal_surface_2.exists():
            du['internal surface 2 (triax)'] = self.lower_internal_surface_2.layer['triax'].area_fraction()
            du['internal surface 2 (resin)'] = self.lower_internal_surface_2.layer['resin'].area_fraction()
        if self.lower_internal_surface_3.exists():
            du['internal surface 3 (triax)'] = self.lower_internal_surface_3.layer['triax'].area_fraction()
            du['internal surface 3 (resin)'] = self.lower_internal_surface_3.layer['resin'].area_fraction()
        if self.lower_internal_surface_4.exists():
            du['internal surface 4 (triax)'] = self.lower_internal_surface_4.layer['triax'].area_fraction()
            du['internal surface 4 (resin)'] = self.lower_internal_surface_4.layer['resin'].area_fraction()
        # add up the contributions from the upper and lower airfoils for the entire station
        d = {}
        if self.lower_external_surface.exists() and self.upper_external_surface.exists():
            d['external surface (gelcoat)'] = dl['external surface (gelcoat)'] + du['external surface (gelcoat)']
            d['external surface (triax)'] = dl['external surface (triax)'] + du['external surface (triax)']
        if self.lower_root_buildup.exists() and self.upper_root_buildup.exists():
            d['root buildup'] = dl['root buildup'] + du['root buildup']
        if self.lower_spar_cap.exists() and self.upper_spar_cap.exists():
            d['spar cap (lower)'] = dl['spar cap (lower)'] + du['spar cap (lower)']
            d['spar cap (upper)'] = dl['spar cap (upper)'] + du['spar cap (upper)']
        if self.lower_aft_panel_1.exists() and self.upper_aft_panel_1.exists():
            d['aft panel 1 (lower)'] = dl['aft panel 1 (lower)'] + du['aft panel 1 (lower)']
            d['aft panel 1 (upper)'] = dl['aft panel 1 (upper)'] + du['aft panel 1 (upper)']
        if self.lower_aft_panel_2.exists() and self.upper_aft_panel_2.exists():
            d['aft panel 2 (lower)'] = dl['aft panel 2 (lower)'] + du['aft panel 2 (lower)']
            d['aft panel 2 (upper)'] = dl['aft panel 2 (upper)'] + du['aft panel 2 (upper)']
        if self.lower_LE_panel.exists() and self.upper_LE_panel.exists():
            d['LE panel'] = dl['LE panel'] + du['LE panel']
        if self.lower_shear_web_1.exists() and self.upper_shear_web_1.exists():
            d['shear web 1 (left biax)'] = dl['shear web 1 (left biax)'] + du['shear web 1 (left biax)']
            d['shear web 1 (foam)'] = dl['shear web 1 (foam)'] + du['shear web 1 (foam)']
            d['shear web 1 (right biax)'] = dl['shear web 1 (right biax)'] + du['shear web 1 (right biax)']
        if self.lower_shear_web_2.exists() and self.upper_shear_web_2.exists():
            d['shear web 2 (left biax)'] = dl['shear web 2 (left biax)'] + du['shear web 2 (left biax)']
            d['shear web 2 (foam)'] = dl['shear web 2 (foam)'] + du['shear web 2 (foam)']
            d['shear web 2 (right biax)'] = dl['shear web 2 (right biax)'] + du['shear web 2 (right biax)']
        if self.lower_shear_web_3.exists() and self.upper_shear_web_3.exists():
            d['shear web 3 (left biax)'] = dl['shear web 3 (left biax)'] + du['shear web 3 (left biax)']
            d['shear web 3 (foam)'] = dl['shear web 3 (foam)'] + du['shear web 3 (foam)']
            d['shear web 3 (right biax)'] = dl['shear web 3 (right biax)'] + du['shear web 3 (right biax)']
        if self.lower_TE_reinforcement.exists() and self.upper_TE_reinforcement.exists():
            d['TE reinforcement (uniax)'] = dl['TE reinforcement (uniax)'] + du['TE reinforcement (uniax)']
            try:
                d['TE reinforcement (foam)'] = dl['TE reinforcement (foam)'] + du['TE reinforcement (foam)']
            except KeyError:
                # foam layer doesn't exist
                pass
        if self.lower_internal_surface_1.exists() and self.upper_internal_surface_1.exists():
            d['internal surface 1 (triax)'] = dl['internal surface 1 (triax)'] + du['internal surface 1 (triax)']
            d['internal surface 1 (resin)'] = dl['internal surface 1 (resin)'] + du['internal surface 1 (resin)']
        if self.lower_internal_surface_2.exists() and self.upper_internal_surface_2.exists():
            d['internal surface 2 (triax)'] = dl['internal surface 2 (triax)'] + du['internal surface 2 (triax)']
            d['internal surface 2 (resin)'] = dl['internal surface 2 (resin)'] + du['internal surface 2 (resin)']
        if self.lower_internal_surface_3.exists() and self.upper_internal_surface_3.exists():
            d['internal surface 3 (triax)'] = dl['internal surface 3 (triax)'] + du['internal surface 3 (triax)']
            d['internal surface 3 (resin)'] = dl['internal surface 3 (resin)'] + du['internal surface 3 (resin)']
        if self.lower_internal_surface_4.exists() and self.upper_internal_surface_4.exists():
            d['internal surface 4 (triax)'] = dl['internal surface 4 (triax)'] + du['internal surface 4 (triax)']
            d['internal surface 4 (resin)'] = dl['internal surface 4 (resin)'] + du['internal surface 4 (resin)']
        return d

    def calculate_all_percent_masses(self):
        """Calculate the percent masses of all parts in this station.

        Returns a dictionary of mass fractions for each structural part.

        """
        # lower airfoil -------------------------------------------------------
        dl = {}  # dictionary for lower airfoils
        if self.lower_external_surface.exists():
            dl['external surface (gelcoat)'] = self.lower_external_surface.layer['gelcoat'].mass_fraction()
            dl['external surface (triax)'] = self.lower_external_surface.layer['triax'].mass_fraction()
        if self.lower_root_buildup.exists():
            dl['root buildup'] = self.lower_root_buildup.layer['triax'].mass_fraction()
        if self.lower_spar_cap.exists():
            dl['spar cap (lower)'] = self.lower_spar_cap.layer['lower'].mass_fraction()
            dl['spar cap (upper)'] = self.lower_spar_cap.layer['upper'].mass_fraction()
        if self.lower_aft_panel_1.exists():
            dl['aft panel 1 (lower)'] = self.lower_aft_panel_1.layer['lower'].mass_fraction()
            dl['aft panel 1 (upper)'] = self.lower_aft_panel_1.layer['upper'].mass_fraction()
        if self.lower_aft_panel_2.exists():
            dl['aft panel 2 (lower)'] = self.lower_aft_panel_2.layer['lower'].mass_fraction()
            dl['aft panel 2 (upper)'] = self.lower_aft_panel_2.layer['upper'].mass_fraction()
        if self.lower_LE_panel.exists():
            dl['LE panel'] = self.lower_LE_panel.layer['foam'].mass_fraction()
        if self.lower_shear_web_1.exists():
            dl['shear web 1 (left biax)'] = self.lower_shear_web_1.layer['biax, left'].mass_fraction()
            dl['shear web 1 (foam)'] = self.lower_shear_web_1.layer['foam'].mass_fraction()
            dl['shear web 1 (right biax)'] = self.lower_shear_web_1.layer['biax, right'].mass_fraction()
        if self.lower_shear_web_2.exists():
            dl['shear web 2 (left biax)'] = self.lower_shear_web_2.layer['biax, left'].mass_fraction()
            dl['shear web 2 (foam)'] = self.lower_shear_web_2.layer['foam'].mass_fraction()
            dl['shear web 2 (right biax)'] = self.lower_shear_web_2.layer['biax, right'].mass_fraction()
        if self.lower_shear_web_3.exists():
            dl['shear web 3 (left biax)'] = self.lower_shear_web_3.layer['biax, left'].mass_fraction()
            dl['shear web 3 (foam)'] = self.lower_shear_web_3.layer['foam'].mass_fraction()
            dl['shear web 3 (right biax)'] = self.lower_shear_web_3.layer['biax, right'].mass_fraction()
        if self.lower_TE_reinforcement.exists():
            dl['TE reinforcement (uniax)'] = self.lower_TE_reinforcement.layer['uniax'].mass_fraction()
            try:
                dl['TE reinforcement (foam)'] = self.lower_TE_reinforcement.layer['foam'].mass_fraction()
            except KeyError:
                # the foam layer doesn't exist in TE reinf at this station
                dl['TE reinforcement (foam)'] = 0.0
        if self.lower_internal_surface_1.exists():
            dl['internal surface 1 (triax)'] = self.lower_internal_surface_1.layer['triax'].mass_fraction()
            dl['internal surface 1 (resin)'] = self.lower_internal_surface_1.layer['resin'].mass_fraction()
        if self.lower_internal_surface_2.exists():
            dl['internal surface 2 (triax)'] = self.lower_internal_surface_2.layer['triax'].mass_fraction()
            dl['internal surface 2 (resin)'] = self.lower_internal_surface_2.layer['resin'].mass_fraction()
        if self.lower_internal_surface_3.exists():
            dl['internal surface 3 (triax)'] = self.lower_internal_surface_3.layer['triax'].mass_fraction()
            dl['internal surface 3 (resin)'] = self.lower_internal_surface_3.layer['resin'].mass_fraction()
        if self.lower_internal_surface_4.exists():
            dl['internal surface 4 (triax)'] = self.lower_internal_surface_4.layer['triax'].mass_fraction()
            dl['internal surface 4 (resin)'] = self.lower_internal_surface_4.layer['resin'].mass_fraction()
        # upper airfoil -------------------------------------------------------
        du = {}  # dictionary for upper airfoils
        if self.lower_external_surface.exists():
            du['external surface (gelcoat)'] = self.lower_external_surface.layer['gelcoat'].mass_fraction()
            du['external surface (triax)'] = self.lower_external_surface.layer['triax'].mass_fraction()
        if self.lower_root_buildup.exists():
            du['root buildup'] = self.lower_root_buildup.layer['triax'].mass_fraction()
        if self.lower_spar_cap.exists():
            du['spar cap (lower)'] = self.lower_spar_cap.layer['lower'].mass_fraction()
            du['spar cap (upper)'] = self.lower_spar_cap.layer['upper'].mass_fraction()
        if self.lower_aft_panel_1.exists():
            du['aft panel 1 (lower)'] = self.lower_aft_panel_1.layer['lower'].mass_fraction()
            du['aft panel 1 (upper)'] = self.lower_aft_panel_1.layer['upper'].mass_fraction()
        if self.lower_aft_panel_2.exists():
            du['aft panel 2 (lower)'] = self.lower_aft_panel_2.layer['lower'].mass_fraction()
            du['aft panel 2 (upper)'] = self.lower_aft_panel_2.layer['upper'].mass_fraction()
        if self.lower_LE_panel.exists():
            du['LE panel'] = self.lower_LE_panel.layer['foam'].mass_fraction()
        if self.lower_shear_web_1.exists():
            du['shear web 1 (left biax)'] = self.lower_shear_web_1.layer['biax, left'].mass_fraction()
            du['shear web 1 (foam)'] = self.lower_shear_web_1.layer['foam'].mass_fraction()
            du['shear web 1 (right biax)'] = self.lower_shear_web_1.layer['biax, right'].mass_fraction()
        if self.lower_shear_web_2.exists():
            du['shear web 2 (left biax)'] = self.lower_shear_web_2.layer['biax, left'].mass_fraction()
            du['shear web 2 (foam)'] = self.lower_shear_web_2.layer['foam'].mass_fraction()
            du['shear web 2 (right biax)'] = self.lower_shear_web_2.layer['biax, right'].mass_fraction()
        if self.lower_shear_web_3.exists():
            du['shear web 3 (left biax)'] = self.lower_shear_web_3.layer['biax, left'].mass_fraction()
            du['shear web 3 (foam)'] = self.lower_shear_web_3.layer['foam'].mass_fraction()
            du['shear web 3 (right biax)'] = self.lower_shear_web_3.layer['biax, right'].mass_fraction()
        if self.lower_TE_reinforcement.exists():
            du['TE reinforcement (uniax)'] = self.lower_TE_reinforcement.layer['uniax'].mass_fraction()
            try:
                du['TE reinforcement (foam)'] = self.lower_TE_reinforcement.layer['foam'].mass_fraction()
            except KeyError:
                # the foam layer doesn't exist in TE reinf at this station
                du['TE reinforcement (foam)'] = 0.0
        if self.lower_internal_surface_1.exists():
            du['internal surface 1 (triax)'] = self.lower_internal_surface_1.layer['triax'].mass_fraction()
            du['internal surface 1 (resin)'] = self.lower_internal_surface_1.layer['resin'].mass_fraction()
        if self.lower_internal_surface_2.exists():
            du['internal surface 2 (triax)'] = self.lower_internal_surface_2.layer['triax'].mass_fraction()
            du['internal surface 2 (resin)'] = self.lower_internal_surface_2.layer['resin'].mass_fraction()
        if self.lower_internal_surface_3.exists():
            du['internal surface 3 (triax)'] = self.lower_internal_surface_3.layer['triax'].mass_fraction()
            du['internal surface 3 (resin)'] = self.lower_internal_surface_3.layer['resin'].mass_fraction()
        if self.lower_internal_surface_4.exists():
            du['internal surface 4 (triax)'] = self.lower_internal_surface_4.layer['triax'].mass_fraction()
            du['internal surface 4 (resin)'] = self.lower_internal_surface_4.layer['resin'].mass_fraction()
        # add up the contributions from the upper and lower airfoils for the entire station
        d = {}
        if self.lower_external_surface.exists() and self.upper_external_surface.exists():
            d['external surface (gelcoat)'] = dl['external surface (gelcoat)'] + du['external surface (gelcoat)']
            d['external surface (triax)'] = dl['external surface (triax)'] + du['external surface (triax)']
        if self.lower_root_buildup.exists() and self.upper_root_buildup.exists():
            d['root buildup'] = dl['root buildup'] + du['root buildup']
        if self.lower_spar_cap.exists() and self.upper_spar_cap.exists():
            d['spar cap (lower)'] = dl['spar cap (lower)'] + du['spar cap (lower)']
            d['spar cap (upper)'] = dl['spar cap (upper)'] + du['spar cap (upper)']
        if self.lower_aft_panel_1.exists() and self.upper_aft_panel_1.exists():
            d['aft panel 1 (lower)'] = dl['aft panel 1 (lower)'] + du['aft panel 1 (lower)']
            d['aft panel 1 (upper)'] = dl['aft panel 1 (upper)'] + du['aft panel 1 (upper)']
        if self.lower_aft_panel_2.exists() and self.upper_aft_panel_2.exists():
            d['aft panel 2 (lower)'] = dl['aft panel 2 (lower)'] + du['aft panel 2 (lower)']
            d['aft panel 2 (upper)'] = dl['aft panel 2 (upper)'] + du['aft panel 2 (upper)']
        if self.lower_LE_panel.exists() and self.upper_LE_panel.exists():
            d['LE panel'] = dl['LE panel'] + du['LE panel']
        if self.lower_shear_web_1.exists() and self.upper_shear_web_1.exists():
            d['shear web 1 (left biax)'] = dl['shear web 1 (left biax)'] + du['shear web 1 (left biax)']
            d['shear web 1 (foam)'] = dl['shear web 1 (foam)'] + du['shear web 1 (foam)']
            d['shear web 1 (right biax)'] = dl['shear web 1 (right biax)'] + du['shear web 1 (right biax)']
        if self.lower_shear_web_2.exists() and self.upper_shear_web_2.exists():
            d['shear web 2 (left biax)'] = dl['shear web 2 (left biax)'] + du['shear web 2 (left biax)']
            d['shear web 2 (foam)'] = dl['shear web 2 (foam)'] + du['shear web 2 (foam)']
            d['shear web 2 (right biax)'] = dl['shear web 2 (right biax)'] + du['shear web 2 (right biax)']
        if self.lower_shear_web_3.exists() and self.upper_shear_web_3.exists():
            d['shear web 3 (left biax)'] = dl['shear web 3 (left biax)'] + du['shear web 3 (left biax)']
            d['shear web 3 (foam)'] = dl['shear web 3 (foam)'] + du['shear web 3 (foam)']
            d['shear web 3 (right biax)'] = dl['shear web 3 (right biax)'] + du['shear web 3 (right biax)']
        if self.lower_TE_reinforcement.exists() and self.upper_TE_reinforcement.exists():
            d['TE reinforcement (uniax)'] = dl['TE reinforcement (uniax)'] + du['TE reinforcement (uniax)']
            try:
                d['TE reinforcement (foam)'] = dl['TE reinforcement (foam)'] + du['TE reinforcement (foam)']
            except KeyError:
                # foam layer doesn't exist
                pass
        if self.lower_internal_surface_1.exists() and self.upper_internal_surface_1.exists():
            d['internal surface 1 (triax)'] = dl['internal surface 1 (triax)'] + du['internal surface 1 (triax)']
            d['internal surface 1 (resin)'] = dl['internal surface 1 (resin)'] + du['internal surface 1 (resin)']
        if self.lower_internal_surface_2.exists() and self.upper_internal_surface_2.exists():
            d['internal surface 2 (triax)'] = dl['internal surface 2 (triax)'] + du['internal surface 2 (triax)']
            d['internal surface 2 (resin)'] = dl['internal surface 2 (resin)'] + du['internal surface 2 (resin)']
        if self.lower_internal_surface_3.exists() and self.upper_internal_surface_3.exists():
            d['internal surface 3 (triax)'] = dl['internal surface 3 (triax)'] + du['internal surface 3 (triax)']
            d['internal surface 3 (resin)'] = dl['internal surface 3 (resin)'] + du['internal surface 3 (resin)']
        if self.lower_internal_surface_4.exists() and self.upper_internal_surface_4.exists():
            d['internal surface 4 (triax)'] = dl['internal surface 4 (triax)'] + du['internal surface 4 (triax)']
            d['internal surface 4 (resin)'] = dl['internal surface 4 (resin)'] + du['internal surface 4 (resin)']
        return d
