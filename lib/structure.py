"""A module for organizing structural part data for a blade station.

Author: Perry Roth-Johnson
Last updated: October 11, 2013

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
    
    def __str__(self):
        return """base:    {0} (meters)
height:  {1} (meters)""".format(self.base, self.height)
    
    def exists(self):
        """Checks if a structural part exists at this station."""
        if isnan(self.base) and isnan(self.height):
            return False
        else:
            return True
    
    def bounding_box(self, x_boundary_buffer=1.2, y_boundary_buffer=1.2):
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
        (minx, miny, maxx, maxy) = af.polygon.bounds
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
    
    def create_layers(self, debug_flag=False):
        """Create the gelcoat and triax layers in the external surface.

        <external_surface>.layer['gelcoat'] : gelcoat layer
        <external_surface>.layer['triax'] : triax layer
        # --
        # <external_surface>.layer[2] : gelcoat layer (lower left quadrant)
        # <external_surface>.layer[3] : gelcoat layer (lower right quadrant)
        # <external_surface>.layer[4] : gelcoat layer (upper right quadrant)
        # <external_surface>.layer[5] : gelcoat layer (upper left quadrant)
        # --
        # <external_surface>.layer[6] : triax layer (lower left quadrant)
        # <external_surface>.layer[7] : triax layer (lower right quadrant)
        # <external_surface>.layer[8] : triax layer (upper right quadrant)
        # <external_surface>.layer[9] : triax layer (upper left quadrant)

        # Note: This stores 2 versions of the same external surface.
        # (1) entire annulus : .layer[0], .layer[1]
        # (2) 8 curved rectangles : .layer[2], .layer[3], .layer[4], .layer[5],
        #     .layer[6], .layer[7], .layer[8], .layer[9]

        """
        st = self.parent_structure
        if self.exists():
            # create the gelcoat layer
            af = st.parent_station.airfoil
            b = st.parent_station.parent_blade
            op_gelcoat = af.polygon  # outer profile is the airfoil profile
            ip_gelcoat = op_gelcoat.buffer(-self.height_gelcoat)
            polygon_gelcoat = op_gelcoat.difference(ip_gelcoat)
            self.layer['gelcoat'] = l.Layer(polygon_gelcoat,
                b.dict_of_materials['gelcoat'], parent_part=self)
            assert self.layer['gelcoat'].polygon.geom_type == 'Polygon'
            st._list_of_layers.append(self.layer['gelcoat'])
            # create the triax layer
            op_triax = ip_gelcoat  # outer profile is the gelcoat inner profile
            ip_triax = op_triax.buffer(-self.height_triax)
            polygon_triax = op_triax.difference(ip_triax)
            self.layer['triax'] = l.Layer(polygon_triax,
                b.dict_of_materials['triaxial GFRP'], parent_part=self)
            assert self.layer['triax'].polygon.geom_type == 'Polygon'
            st._list_of_layers.append(self.layer['triax'])
            # # now, we cut the annulus into 8 curved rectangles (4 per layer),
            # #   with one in each quadrant
            # bb = self.bounding_box() # get bounding boxes for each quadrant
            # # gelcoat layer
            # print "Station #{0} ----------".format(st.parent_station.station_num)
            # for i, box in enumerate(bb):
            #     p_quad = polygon_gelcoat.intersection(box)
            #     self.layer.append(l.Layer(p_quad,
            #         b.dict_of_materials['gelcoat'], parent_part=self))
            #     # check that layer i+2 (start counting at 2) is a polygon
            #     print "Layer #{0}, {1}".format((i+2), self.layer[i+2].polygon.geom_type)
            #     assert self.layer[i+2].polygon.geom_type == 'Polygon'
            #     # no need to append these polygons to <station>._list_of_layers
            # # triax layer
            # for i, box in enumerate(bb):
            #     p_quad = polygon_triax.intersection(box)
            #     self.layer.append(l.Layer(p_quad,
            #         b.dict_of_materials['triaxial GFRP'], parent_part=self))
            #     # check that layer i+6 (start counting at 6) is a polygon
            #     assert self.layer[i+6].polygon.geom_type == 'Polygon'
            #     # no need to append these polygons to <station>._list_of_layers
        else:
            if debug_flag:
                print " The external surface for Station #{0} does not exist!\n  No layers created.".format(st.parent_station.station_num)

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
        bb = []
        (minx, miny, maxx, maxy) = af.polygon.bounds
        # lower left quadrant
        pt1 = (minx*x_boundary_buffer, miny*y_boundary_buffer)
        pt2 = (0.0, miny*y_boundary_buffer)
        pt3 = (0.0, 0.0)
        pt4 = (minx*x_boundary_buffer, 0.0)
        bounding_box = Polygon([pt1, pt2, pt3, pt4])
        bb.append(bounding_box)
        # lower right quadrant
        pt1 = (0.0, miny*y_boundary_buffer)
        pt2 = (maxx*x_boundary_buffer, miny*y_boundary_buffer)
        pt3 = (maxx*x_boundary_buffer, 0.0)
        pt4 = (0.0, 0.0)
        bounding_box = Polygon([pt1, pt2, pt3, pt4])
        bb.append(bounding_box)
        # upper right quadrant
        pt1 = (0.0, 0.0)
        pt2 = (maxx*x_boundary_buffer, 0.0)
        pt3 = (maxx*x_boundary_buffer, maxy*y_boundary_buffer)
        pt4 = (0.0, maxy*y_boundary_buffer)
        bounding_box = Polygon([pt1, pt2, pt3, pt4])
        bb.append(bounding_box)
        # upper left quadrant
        pt1 = (minx*x_boundary_buffer, 0.0)
        pt2 = (0.0, 0.0)
        pt3 = (0.0, maxy*y_boundary_buffer)
        pt4 = (minx*x_boundary_buffer, maxy*y_boundary_buffer)
        bounding_box = Polygon([pt1, pt2, pt3, pt4])
        bb.append(bounding_box)
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
    def create_layers(self, debug_flag=False):
        """Create the triax layer in the root buildup.

        <root_buildup>.layer['triax, annulus'] : triax layer (entire annulus)
        --
        <root_buildup>.layer['triax, lower left'] : triax layer (lower left
            quadrant)
        <root_buildup>.layer['triax, lower right'] : triax layer (lower right
            quadrant)
        <root_buildup>.layer['triax, upper right'] : triax layer (upper right
            quadrant)
        <root_buildup>.layer['triax, upper left'] : triax layer (upper left
            quadrant)

        Note: This stores 2 versions of the same root buildup.
        (1) entire annulus : .layer['triax, annulus']
        (2) 4 curved rectangles : .layer['triax, lower left'],
            .layer['triax, lower right'], .layer['triax, upper right'],
            .layer['triax, upper left']

        """
        st = self.parent_structure
        if self.exists():
            # create the triax layer
            af = st.parent_station.airfoil
            b = st.parent_station.parent_blade
            op = af.polygon.buffer(-st.external_surface.height)
            ip = op.buffer(-self.height)
            p = op.difference(ip)  # this polygon is like an annulus
            self.layer['triax, annulus'] = l.Layer(p,
                b.dict_of_materials['triaxial GFRP'], parent_part=self)
            # check that layer['triax, annulus'] is a Polygon
            assert self.layer['triax, annulus'].polygon.geom_type == 'Polygon'
            st._list_of_layers.append(self.layer['triax, annulus'])
            # now, we cut the annulus into 4 curved rectangles, with one in
            #   each quadrant
            bb = self.bounding_box() # get bounding boxes for each quadrant
            for (label,box) in bb.items():
                p_quad = p.intersection(box)
                self.layer['triax, '+label] = l.Layer(p_quad,
                    b.dict_of_materials['triaxial GFRP'], parent_part=self,
                    name=('triax, '+label))
                # check that the layer just created is a polygon
                assert self.layer['triax, '+label].polygon.geom_type == 'Polygon'
                # no need to append these polygons to <station>._list_of_layers
        else:
            if debug_flag:
                print " The root buildup for Station #{0} does not exist!\n  No layers created.".format(st.parent_station.station_num)

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
    def create_layers(self, debug_flag=False):
        """Create the foam layer in the leading edge panel.

        <LE_panel>.layer['foam'] : the only layer in the LE panel (foam)

        """
        st = self.parent_structure
        if self.exists():
            # create the foam layer
            af = st.parent_station.airfoil
            b = st.parent_station.parent_blade
            # 1. get outer profile
            if st.root_buildup.exists():
                op = af.polygon.buffer(-(st.external_surface.height + 
                    st.root_buildup.height))
            else:
                op = af.polygon.buffer(-st.external_surface.height)
            # 2. erode the outer profile by the part thickness
            ip = op.buffer(-self.height)
            # 3. cut out the part interior from the outer profile
            ac = op.difference(ip)
            # 4. draw a bounding box at the part edges
            bb = self.bounding_box()
            # 5. cut out the structural part
            p = ac.intersection(bb)
            self.layer['foam'] = l.Layer(p, b.dict_of_materials['foam'],
                parent_part=self)
            assert self.layer['foam'].polygon.geom_type == 'Polygon'
            st._list_of_layers.append(self.layer['foam'])
        else:
            if debug_flag:
                print " The LE panel for Station #{0} does not exist!\n  No layers created.".format(st.parent_station.station_num)


class SparCap(Part):
    """Define uniax dimensions of the lower and upper spar caps."""
    def create_layers(self, debug_flag=False):
        """Create the uniax layers in the lower and upper spar caps.

        <spar_cap>.layer['lower'] : lower spar cap
        <spar_cap>.layer['upper'] : upper spar cap

        """
        st = self.parent_structure
        if self.exists():
            # create the uniax layer
            af = st.parent_station.airfoil
            b = st.parent_station.parent_blade
            # 1. get outer profile
            if st.root_buildup.exists():
                op = af.polygon.buffer(-(st.external_surface.height + 
                    st.root_buildup.height))
            else:
                op = af.polygon.buffer(-st.external_surface.height)
            # 2. erode the outer profile by the part thickness
            ip = op.buffer(-self.height)
            # 3. cut out the part interior from the outer profile
            ac = op.difference(ip)
            # 4. draw a bounding box at the part edges
            bb = self.bounding_box()
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
                parent_part=self)
            assert self.layer['lower'].polygon.geom_type == 'Polygon'
            st._list_of_layers.append(self.layer['lower'])
            # 8. add the upper spar cap
            self.layer['upper'] = l.Layer(pu, b.dict_of_materials['uniaxial GFRP'],
                parent_part=self)
            assert self.layer['upper'].polygon.geom_type == 'Polygon'
            st._list_of_layers.append(self.layer['upper'])
        else:
            if debug_flag:
                print " The spar caps for Station #{0} do not exist!\n  No layers created.".format(st.parent_station.station_num)


class AftPanel(Part):
    """Define foam dimensions of the lower and upper aft panels."""
    def create_layers(self, debug_flag=False):
        """Create the foam layers in the lower and upper aft panels.

        <aft_panel>.layer['lower'] : lower aft panel
        <aft_panel>.layer['upper'] : upper aft panel

        """
        st = self.parent_structure
        if self.exists():
            # create the foam layer
            af = st.parent_station.airfoil
            b = st.parent_station.parent_blade
            # 1. get outer profile
            if st.root_buildup.exists():
                op = af.polygon.buffer(-(st.external_surface.height + 
                    st.root_buildup.height))
            else:
                op = af.polygon.buffer(-st.external_surface.height)
            # 2. erode the outer profile by the part thickness
            ip = op.buffer(-self.height)
            # 3. cut out the part interior from the outer profile
            ac = op.difference(ip)
            # 4. draw a bounding box at the part edges
            bb = self.bounding_box()
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
                parent_part=self)
            assert self.layer['lower'].polygon.geom_type == 'Polygon'
            st._list_of_layers.append(self.layer['lower'])
            # 8. add the upper aft panel
            self.layer['upper'] = l.Layer(pu, b.dict_of_materials['foam'],
                parent_part=self)
            assert self.layer['upper'].polygon.geom_type == 'Polygon'
            st._list_of_layers.append(self.layer['upper'])
        else:
            if debug_flag:
                print " The aft panels for Station #{0} do not exist!\n  No layers created.".format(st.parent_station.station_num)


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
    
    def __str__(self):
        return """base:    {0:6.4f} (meters)
height:  {1:6.4f} (meters)
|-> height_uniax:  {2:6.4f} (meters)
|-> height_foam:   {3:6.4f} (meters)""".format(self.base, self.height,
    self.height_uniax, self.height_foam)

    def create_layers(self, debug_flag=False):
        """Create the uniax and foam layers in the TE reinforcement.

        The TE reinforcement is split into one OR two regions:
        <TE_reinforcement>.layer['uniax'] : uniax layer
        <TE_reinforcement>.layer['foam'] : foam layer (optional)

        """
        st = self.parent_structure
        if self.exists():
            # create the uniax layer
            af = st.parent_station.airfoil
            b = st.parent_station.parent_blade
            # 1. get outer profile
            if st.root_buildup.exists():
                op_uniax = af.polygon.buffer(-(st.external_surface.height + 
                    st.root_buildup.height))
            else:
                op_uniax = af.polygon.buffer(-st.external_surface.height)
            # 2. erode the outer profile by the uniax thickness
            ip_uniax = op_uniax.buffer(-self.height_uniax)
            # 3. cut out the uniax layer from the outer profile
            ac_uniax = op_uniax.difference(ip_uniax)
            # 4. draw a bounding box at the TE reinforcement edges
            bb = self.bounding_box()
            # 5. cut out the uniax layer
            polygon_uniax = ac_uniax.intersection(bb)
            # 6. add the uniax layer
            self.layer['uniax'] = l.Layer(polygon_uniax,
                b.dict_of_materials['uniaxial GFRP'], parent_part=self)
            assert self.layer['uniax'].polygon.geom_type == 'Polygon'
            st._list_of_layers.append(self.layer['uniax'])
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
                    b.dict_of_materials['foam'], parent_part=self)
                assert self.layer['foam'].polygon.geom_type == 'Polygon'
                st._list_of_layers.append(self.layer['foam'])
        else:
            if debug_flag:
                print " The TE reinforcement for Station #{0} does not exist!\n  No layers created.".format(st.parent_station.station_num)


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
    def __init__(self, parent_structure, base_biax, base_foam, x2,
        height=np.nan):
        Part.__init__(self, parent_structure, base=(2.0*base_biax+base_foam),
            height=height)
        self.base_biax = base_biax
        self.base_foam = base_foam
        self.x2 = x2
        self.cs_coords = None   # assigned later by <station>.find_SW_cs_coords()

    def __str__(self):
        return """base:    {0:6.4f} (meters)
|-> base_biax:  {1:6.4f} (meters)
|-> base_foam:  {2:6.4f} (meters)
height:  {3} (meters)
x2:      {4:6.4f} (meters)""".format(self.base, self.base_biax,
    self.base_foam, self.height, self.x2)

    def bounding_box(self, y_boundary_buffer=1.2):
        """Returns 3 bounding boxes for the biax and foam regions of the SW.

        The points of each bounding box are labeled from 1 to 4 as:

        4---3
        |   |
        1---2

        """
        af = self.parent_structure.parent_station.airfoil
        (minx, miny, maxx, maxy) = af.polygon.bounds
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
    
    def create_layers(self, debug_flag=False):
        """Create the biax and foam layers in this shear web.

        <shear_web>.layer['biax, left'] is the left biax layer
        <shear_web>.layer['foam'] is the foam layer
        <shear_web>.layer['biax, right'] is the right biax layer

        """
        st = self.parent_structure
        if self.exists():
            # create the foam layer
            af = st.parent_station.airfoil
            b = st.parent_station.parent_blade
            # 1. get outer profile
            if st.root_buildup.exists():
                op = af.polygon.buffer(-(st.external_surface.height + 
                    st.root_buildup.height))
            else:
                op = af.polygon.buffer(-st.external_surface.height)
            # 2. get bounding boxes for the biax and foam regions
            (bb_left_biax, bb_foam, bb_right_biax) = self.bounding_box()
            # 3. cut out the structural parts
            p_left_biax = op.intersection(bb_left_biax)
            p_foam = op.intersection(bb_foam)
            p_right_biax = op.intersection(bb_right_biax)
            # 4. add the left biax layer
            self.layer['biax, left'] = l.Layer(p_left_biax,
                b.dict_of_materials['biaxial GFRP'], parent_part=self,
                name='biax, left')
            assert self.layer['biax, left'].polygon.geom_type == 'Polygon'
            st._list_of_layers.append(self.layer['biax, left'])
            # 5. add the foam layer
            self.layer['foam'] = l.Layer(p_foam, b.dict_of_materials['foam'],
                parent_part=self, name='foam')
            assert self.layer['foam'].polygon.geom_type == 'Polygon'
            st._list_of_layers.append(self.layer['foam'])
            # 6. add the right biax layer
            self.layer['biax, right'] = l.Layer(p_right_biax,
                b.dict_of_materials['biaxial GFRP'], parent_part=self,
                name='biax, right')
            assert self.layer['biax, right'].polygon.geom_type == 'Polygon'
            st._list_of_layers.append(self.layer['biax, right'])
        else:
            if debug_flag:
                print " The shear web for Station #{0} does not exist!\n  No layers created.".format(st.parent_station.station_num)


class InternalSurface(Part):
    """Define triax and resin dimensions of the internal surface."""
    def __init__(self, parent_structure, base, height_triax, height_resin,
        internal_surface_num):
        Part.__init__(self, parent_structure, base,
            height=(height_triax+height_resin))
        self.height_triax = height_triax
        self.height_resin = height_resin
        self.internal_surface_num = internal_surface_num
    
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
        #   internal_surface_1 ==> internal_surface_num=1 ==> good_loops[0]
        #   internal_surface_2 ==> internal_surface_num=2 ==> good_loops[1]
        #   internal_surface_3 ==> internal_surface_num=3 ==> good_loops[2]
        #   internal_surface_4 ==> internal_surface_num=4 ==> good_loops[3]
        loop = Polygon(good_loops[self.internal_surface_num-1])
        return loop

    def create_layers(self, merged_polygon, debug_flag=False):
        """Create the triax and resin layers in the internal surface.

        <internal_surface>.layer['triax'] : triax layer
        <internal_surface>.layer['resin'] : resin layer

        """
        st = self.parent_structure
        if self.exists():
            b = st.parent_station.parent_blade
            # triax region
            op_triax = self.interior_loop(merged_polygon)
            ip_triax = op_triax.buffer(-self.height_triax)
            polygon_triax = op_triax.difference(ip_triax)
            self.layer['triax'] = l.Layer(polygon_triax,
                b.dict_of_materials['triaxial GFRP'], parent_part=self)
            assert self.layer['triax'].polygon.geom_type == 'Polygon'
            st._list_of_layers.append(self.layer['triax'])
            # resin region
            op_resin = ip_triax
            ip_resin = op_resin.buffer(-self.height_resin)
            polygon_resin = op_resin.difference(ip_resin)
            self.layer['resin'] = l.Layer(polygon_resin,
                b.dict_of_materials['resin'], parent_part=self)
            assert self.layer['resin'].polygon.geom_type == 'Polygon'
            st._list_of_layers.append(self.layer['resin'])
        else:
            if debug_flag:
                print " The internal surface for Station #{0} does not exist!\n  No layers created.".format(st.parent_station.station_num)


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
        self.root_buildup = RootBuildup(
            parent_structure = self,
            base = np.nan,
            height = h_RB)
        self.spar_cap = SparCap(
            parent_structure = self,
            base = b_SC,
            height = h_SC)
        self.shear_web_1 = ShearWeb(
            parent_structure = self,
            base_biax = b_SW1_biax,
            base_foam = b_SW1_foam,
            x2 = x2_SW1)
        self.shear_web_2 = ShearWeb(
            parent_structure = self,
            base_biax = b_SW2_biax,
            base_foam = b_SW2_foam,
            x2 = x2_SW2)
        self.shear_web_3 = ShearWeb(
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
            parent_structure = self,
            base = np.nan,
            height = h_aft_panel_1)
        self.aft_panel_2 = AftPanel(
            parent_structure = self,
            base = np.nan,
            height = h_aft_panel_2)
        self.internal_surface_1 = InternalSurface(
            parent_structure = self,
            base = np.nan,
            height_triax = h_int_surf_1_triax,
            height_resin = h_int_surf_1_resin,
            internal_surface_num = 1)
        self.internal_surface_2 = InternalSurface(
            parent_structure = self,
            base = np.nan,
            height_triax = h_int_surf_2_triax,
            height_resin = h_int_surf_2_resin,
            internal_surface_num = 2)
        self.internal_surface_3 = InternalSurface(
            parent_structure = self,
            base = np.nan,
            height_triax = h_int_surf_3_triax,
            height_resin = h_int_surf_3_resin,
            internal_surface_num = 3)
        self.internal_surface_4 = InternalSurface(
            parent_structure = self,
            base = np.nan,
            height_triax = h_int_surf_4_triax,
            height_resin = h_int_surf_4_resin,
            internal_surface_num = 4)
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
        self.external_surface.create_layers()
        self.root_buildup.create_layers()
        self.LE_panel.create_layers()
        self.spar_cap.create_layers()
        self.aft_panel_1.create_layers()
        self.aft_panel_2.create_layers()
        self.shear_web_1.create_layers()
        self.shear_web_2.create_layers()
        self.shear_web_3.create_layers()
        self.TE_reinforcement.create_layers()
        mp = self.merge_all_polygons()
        self.internal_surface_1.create_layers(mp)
        self.internal_surface_2.create_layers(mp)
        self.internal_surface_3.create_layers(mp)
        self.internal_surface_4.create_layers(mp)

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
                RB = self.root_buildup.layer['triax, annulus'].polygon
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
            d['external surface (gelcoat)'] = self.external_surface.layer[0].area_fraction()
            d['external surface (triax)'] = self.external_surface.layer[1].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, external surface, gelcoat".format(self.external_surface.layer[0].area_fraction())
                print "  {0:5.1%} area, external surface, triax".format(self.external_surface.layer[1].area_fraction())
        if self.root_buildup.exists():
            d['root buildup'] = self.root_buildup.layer[0].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, root buildup".format(self.root_buildup.layer[0].area_fraction())
        if self.spar_cap.exists():
            d['spar cap (lower)'] = self.spar_cap.layer[0].area_fraction()
            d['spar cap (upper)'] = self.spar_cap.layer[1].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, spar cap, lower".format(self.spar_cap.layer[0].area_fraction())
                print "  {0:5.1%} area, spar cap, upper".format(self.spar_cap.layer[1].area_fraction())
        if self.aft_panel_1.exists():
            d['aft panel 1 (lower)'] = self.aft_panel_1.layer[0].area_fraction()
            d['aft panel 1 (upper)'] = self.aft_panel_1.layer[1].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, aft panel 1, lower".format(self.aft_panel_1.layer[0].area_fraction())
                print "  {0:5.1%} area, aft panel 1, upper".format(self.aft_panel_1.layer[1].area_fraction())
        if self.aft_panel_2.exists():
            d['aft panel 2 (lower)'] = self.aft_panel_2.layer[0].area_fraction()
            d['aft panel 2 (upper)'] = self.aft_panel_2.layer[1].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, aft panel 2, lower".format(self.aft_panel_2.layer[0].area_fraction())
                print "  {0:5.1%} area, aft panel 2, upper".format(self.aft_panel_2.layer[1].area_fraction())
        if self.LE_panel.exists():
            d['LE panel'] = self.LE_panel.layer[0].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, LE panel".format(self.LE_panel.layer[0].area_fraction())
        if self.shear_web_1.exists():
            d['shear web 1 (left biax)'] = self.shear_web_1.layer[0].area_fraction()
            d['shear web 1 (foam)'] = self.shear_web_1.layer[1].area_fraction()
            d['shear web 1 (right biax)'] = self.shear_web_1.layer[2].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, shear web 1, left biax".format(self.shear_web_1.layer[0].area_fraction())
                print "  {0:5.1%} area, shear web 1, foam".format(self.shear_web_1.layer[1].area_fraction())
                print "  {0:5.1%} area, shear web 1, right biax".format(self.shear_web_1.layer[2].area_fraction())
        if self.shear_web_2.exists():
            d['shear web 2 (left biax)'] = self.shear_web_2.layer[0].area_fraction()
            d['shear web 2 (foam)'] = self.shear_web_2.layer[1].area_fraction()
            d['shear web 2 (right biax)'] = self.shear_web_2.layer[2].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, shear web 2, left biax".format(self.shear_web_2.layer[0].area_fraction())
                print "  {0:5.1%} area, shear web 2, foam".format(self.shear_web_2.layer[1].area_fraction())
                print "  {0:5.1%} area, shear web 2, right biax".format(self.shear_web_2.layer[2].area_fraction())
        if self.shear_web_3.exists():
            d['shear web 3 (left biax)'] = self.shear_web_3.layer[0].area_fraction()
            d['shear web 3 (foam)'] = self.shear_web_3.layer[1].area_fraction()
            d['shear web 3 (right biax)'] = self.shear_web_3.layer[2].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, shear web 3, left biax".format(self.shear_web_3.layer[0].area_fraction())
                print "  {0:5.1%} area, shear web 3, foam".format(self.shear_web_3.layer[1].area_fraction())
                print "  {0:5.1%} area, shear web 3, right biax".format(self.shear_web_3.layer[2].area_fraction())
        if self.TE_reinforcement.exists():
            d['TE reinforcement (uniax)'] = self.TE_reinforcement.layer[0].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, TE reinforcement, uniax".format(self.TE_reinforcement.layer[0].area_fraction())
            try:
                d['TE reinforcement (foam)'] = self.TE_reinforcement.layer[1].area_fraction()
                if print_flag:
                    print "  {0:5.1%} area, TE reinforcement, foam".format(self.TE_reinforcement.layer[1].area_fraction())
            except IndexError:
                pass
        if self.internal_surface_1.exists():
            d['internal surface 1 (triax)'] = self.internal_surface_1.layer[0].area_fraction()
            d['internal surface 1 (resin)'] = self.internal_surface_1.layer[1].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, internal surface 1, triax".format(self.internal_surface_1.layer[0].area_fraction())
                print "  {0:5.1%} area, internal surface 1, resin".format(self.internal_surface_1.layer[1].area_fraction())
        if self.internal_surface_2.exists():
            d['internal surface 2 (triax)'] = self.internal_surface_2.layer[0].area_fraction()
            d['internal surface 2 (resin)'] = self.internal_surface_2.layer[1].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, internal surface 2, triax".format(self.internal_surface_2.layer[0].area_fraction())
                print "  {0:5.1%} area, internal surface 2, resin".format(self.internal_surface_2.layer[1].area_fraction())
        if self.internal_surface_3.exists():
            d['internal surface 3 (triax)'] = self.internal_surface_3.layer[0].area_fraction()
            d['internal surface 3 (resin)'] = self.internal_surface_3.layer[1].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, internal surface 3, triax".format(self.internal_surface_3.layer[0].area_fraction())
                print "  {0:5.1%} area, internal surface 3, resin".format(self.internal_surface_3.layer[1].area_fraction())
        if self.internal_surface_4.exists():
            d['internal surface 4 (triax)'] = self.internal_surface_4.layer[0].area_fraction()
            d['internal surface 4 (resin)'] = self.internal_surface_4.layer[1].area_fraction()
            if print_flag:
                print "  {0:5.1%} area, internal surface 4, triax".format(self.internal_surface_4.layer[0].area_fraction())
                print "  {0:5.1%} area, internal surface 4, resin".format(self.internal_surface_4.layer[1].area_fraction())
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
            d['external surface (gelcoat)'] = self.external_surface.layer[0].mass_fraction()
            d['external surface (triax)'] = self.external_surface.layer[1].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, external surface, gelcoat".format(self.external_surface.layer[0].mass_fraction())
                print "  {0:5.1%} mass/length, external surface, triax".format(self.external_surface.layer[1].mass_fraction())
        if self.root_buildup.exists():
            d['root buildup'] = self.root_buildup.layer[0].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, root buildup".format(self.root_buildup.layer[0].mass_fraction())
        if self.spar_cap.exists():
            d['spar cap (lower)'] = self.spar_cap.layer[0].mass_fraction()
            d['spar cap (upper)'] = self.spar_cap.layer[1].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, spar cap, lower".format(self.spar_cap.layer[0].mass_fraction())
                print "  {0:5.1%} mass/length, spar cap, upper".format(self.spar_cap.layer[1].mass_fraction())
        if self.aft_panel_1.exists():
            d['aft panel 1 (lower)'] = self.aft_panel_1.layer[0].mass_fraction()
            d['aft panel 1 (upper)'] = self.aft_panel_1.layer[1].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, aft panel 1, lower".format(self.aft_panel_1.layer[0].mass_fraction())
                print "  {0:5.1%} mass/length, aft panel 1, upper".format(self.aft_panel_1.layer[1].mass_fraction())
        if self.aft_panel_2.exists():
            d['aft panel 2 (lower)'] = self.aft_panel_2.layer[0].mass_fraction()
            d['aft panel 2 (upper)'] = self.aft_panel_2.layer[1].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, aft panel 2, lower".format(self.aft_panel_2.layer[0].mass_fraction())
                print "  {0:5.1%} mass/length, aft panel 2, upper".format(self.aft_panel_2.layer[1].mass_fraction())
        if self.LE_panel.exists():
            d['LE panel'] = self.LE_panel.layer[0].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, LE panel".format(self.LE_panel.layer[0].mass_fraction())
        if self.shear_web_1.exists():
            d['shear web 1 (left biax)'] = self.shear_web_1.layer[0].mass_fraction()
            d['shear web 1 (foam)'] = self.shear_web_1.layer[1].mass_fraction()
            d['shear web 1 (right biax)'] = self.shear_web_1.layer[2].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, shear web 1, left biax".format(self.shear_web_1.layer[0].mass_fraction())
                print "  {0:5.1%} mass/length, shear web 1, foam".format(self.shear_web_1.layer[1].mass_fraction())
                print "  {0:5.1%} mass/length, shear web 1, right biax".format(self.shear_web_1.layer[2].mass_fraction())
        if self.shear_web_2.exists():
            d['shear web 2 (left biax)'] = self.shear_web_2.layer[0].mass_fraction()
            d['shear web 2 (foam)'] = self.shear_web_2.layer[1].mass_fraction()
            d['shear web 2 (right biax)'] = self.shear_web_2.layer[2].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, shear web 2, left biax".format(self.shear_web_2.layer[0].mass_fraction())
                print "  {0:5.1%} mass/length, shear web 2, foam".format(self.shear_web_2.layer[1].mass_fraction())
                print "  {0:5.1%} mass/length, shear web 2, right biax".format(self.shear_web_2.layer[2].mass_fraction())
        if self.shear_web_3.exists():
            d['shear web 3 (left biax)'] = self.shear_web_3.layer[0].mass_fraction()
            d['shear web 3 (foam)'] = self.shear_web_3.layer[1].mass_fraction()
            d['shear web 3 (right biax)'] = self.shear_web_3.layer[2].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, shear web 3, left biax".format(self.shear_web_3.layer[0].mass_fraction())
                print "  {0:5.1%} mass/length, shear web 3, foam".format(self.shear_web_3.layer[1].mass_fraction())
                print "  {0:5.1%} mass/length, shear web 3, right biax".format(self.shear_web_3.layer[2].mass_fraction())
        if self.TE_reinforcement.exists():
            d['TE reinforcement (uniax)'] = self.TE_reinforcement.layer[0].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, TE reinforcement, uniax".format(self.TE_reinforcement.layer[0].mass_fraction())
            try:
                d['TE reinforcement (foam)'] = self.TE_reinforcement.layer[1].mass_fraction()
                if print_flag:
                    print "  {0:5.1%} mass/length, TE reinforcement, foam".format(self.TE_reinforcement.layer[1].mass_fraction())
            except IndexError:
                pass
        if self.internal_surface_1.exists():
            d['internal surface 1 (triax)'] = self.internal_surface_1.layer[0].mass_fraction()
            d['internal surface 1 (resin)'] = self.internal_surface_1.layer[1].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, internal surface 1, triax".format(self.internal_surface_1.layer[0].mass_fraction())
                print "  {0:5.1%} mass/length, internal surface 1, resin".format(self.internal_surface_1.layer[1].mass_fraction())
        if self.internal_surface_2.exists():
            d['internal surface 2 (triax)'] = self.internal_surface_2.layer[0].mass_fraction()
            d['internal surface 2 (resin)'] = self.internal_surface_2.layer[1].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, internal surface 2, triax".format(self.internal_surface_2.layer[0].mass_fraction())
                print "  {0:5.1%} mass/length, internal surface 2, resin".format(self.internal_surface_2.layer[1].mass_fraction())
        if self.internal_surface_3.exists():
            d['internal surface 3 (triax)'] = self.internal_surface_3.layer[0].mass_fraction()
            d['internal surface 3 (resin)'] = self.internal_surface_3.layer[1].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, internal surface 3, triax".format(self.internal_surface_3.layer[0].mass_fraction())
                print "  {0:5.1%} mass/length, internal surface 3, resin".format(self.internal_surface_3.layer[1].mass_fraction())
        if self.internal_surface_4.exists():
            d['internal surface 4 (triax)'] = self.internal_surface_4.layer[0].mass_fraction()
            d['internal surface 4 (resin)'] = self.internal_surface_4.layer[1].mass_fraction()
            if print_flag:
                print "  {0:5.1%} mass/length, internal surface 4, triax".format(self.internal_surface_4.layer[0].mass_fraction())
                print "  {0:5.1%} mass/length, internal surface 4, resin".format(self.internal_surface_4.layer[1].mass_fraction())
        return d

    def write_all_part_polygons(self):
        """Write the coordinates of all structural parts to `station_path`s."""
        stn = self.parent_station
        if self.external_surface.exists():
            f = open(os.path.join(stn.station_path,'external_surface.txt'), 'w')
            f.write("# gelcoat region\n")
            f.write("# --------------\n")
            f.write(str(self.external_surface.layer[0].polygon.__geo_interface__))
            f.write("\n\n")
            f.write("# triax region\n")
            f.write("# ------------\n")
            f.write(str(self.external_surface.layer[1].polygon.__geo_interface__))
            f.close()
        if self.root_buildup.exists():
            f = open(os.path.join(stn.station_path,'root_buildup.txt'), 'w')
            f.write(str(self.root_buildup.layer[0].polygon.__geo_interface__))
            f.close()
        if self.spar_cap.exists():
            f = open(os.path.join(stn.station_path,'spar_cap.txt'), 'w')
            f.write("# lower spar cap\n")
            f.write("# --------------\n")
            f.write(str(self.spar_cap.layer[0].polygon.__geo_interface__))
            f.write("\n\n")
            f.write("# upper spar cap\n")
            f.write("# --------------\n")
            f.write(str(self.spar_cap.layer[1].polygon.__geo_interface__))
            f.close()
        if self.aft_panel_1.exists():
            f = open(os.path.join(stn.station_path,'aft_panel_1.txt'), 'w')
            f.write("# lower aft panel #1\n")
            f.write("# ------------------\n")
            f.write(str(self.aft_panel_1.layer[0].polygon.__geo_interface__))
            f.write("\n\n")
            f.write("# upper aft panel #1\n")
            f.write("# ------------------\n")
            f.write(str(self.aft_panel_1.layer[1].polygon.__geo_interface__))
            f.close()
        if self.aft_panel_2.exists():
            f = open(os.path.join(stn.station_path,'aft_panel_2.txt'), 'w')
            f.write("# lower aft panel #2\n")
            f.write("# ------------------\n")
            f.write(str(self.aft_panel_2.layer[0].polygon.__geo_interface__))
            f.write("\n\n")
            f.write("# upper aft panel #2\n")
            f.write("# ------------------\n")
            f.write(str(self.aft_panel_2.layer[1].polygon.__geo_interface__))
            f.close()
        if self.LE_panel.exists():
            f = open(os.path.join(stn.station_path,'LE_panel.txt'), 'w')
            f.write(str(self.LE_panel.layer[0].polygon.__geo_interface__))
            f.close()
        if self.shear_web_1.exists():
            f = open(os.path.join(stn.station_path,'shear_web_1.txt'), 'w')
            f.write("# left biax region\n")
            f.write("# ----------------\n")
            f.write(str(self.shear_web_1.layer[0].polygon.__geo_interface__))
            f.write("\n\n")
            f.write("# foam region\n")
            f.write("# -----------\n")
            f.write(str(self.shear_web_1.layer[1].polygon.__geo_interface__))
            f.write("\n\n")
            f.write("# right biax region\n")
            f.write("# -----------------\n")
            f.write(str(self.shear_web_1.layer[2].polygon.__geo_interface__))
            f.close()
        if self.shear_web_2.exists():
            f = open(os.path.join(stn.station_path,'shear_web_2.txt'), 'w')
            f.write("# left biax region\n")
            f.write("# ----------------\n")
            f.write(str(self.shear_web_2.layer[0].polygon.__geo_interface__))
            f.write("\n\n")
            f.write("# foam region\n")
            f.write("# -----------\n")
            f.write(str(self.shear_web_2.layer[1].polygon.__geo_interface__))
            f.write("\n\n")
            f.write("# right biax region\n")
            f.write("# -----------------\n")
            f.write(str(self.shear_web_2.layer[2].polygon.__geo_interface__))
            f.close()
        if self.shear_web_3.exists():
            f = open(os.path.join(stn.station_path,'shear_web_3.txt'), 'w')
            f.write("# left biax region\n")
            f.write("# ----------------\n")
            f.write(str(self.shear_web_3.layer[0].polygon.__geo_interface__))
            f.write("\n\n")
            f.write("# foam region\n")
            f.write("# -----------\n")
            f.write(str(self.shear_web_3.layer[1].polygon.__geo_interface__))
            f.write("\n\n")
            f.write("# right biax region\n")
            f.write("# -----------------\n")
            f.write(str(self.shear_web_3.layer[2].polygon.__geo_interface__))
            f.close()
        if self.TE_reinforcement.exists():
            f = open(os.path.join(stn.station_path,'TE_reinforcement.txt'), 'w')
            f.write("# uniax region\n")
            f.write("# ------------\n")
            f.write(str(self.TE_reinforcement.layer[0].polygon.__geo_interface__))
            f.write("\n\n")
            f.write("# foam region\n")
            f.write("# -----------\n")
            try:
                f.write(str(self.TE_reinforcement.layer[1].polygon.__geo_interface__))
            except IndexError:
                f.write("# ...the foam region doesn't exist!")
            f.close()
        if self.internal_surface_1.exists():
            f = open(os.path.join(stn.station_path,'internal_surface_1.txt'), 'w')
            f.write("# triax region\n")
            f.write("--------------\n")
            f.write(str(self.internal_surface_1.layer[0].polygon.__geo_interface__))
            f.write("\n\n")
            f.write("# resin region\n")
            f.write("--------------\n")
            f.write(str(self.internal_surface_1.layer[1].polygon.__geo_interface__))
            f.close()
        if self.internal_surface_2.exists():
            f = open(os.path.join(stn.station_path,'internal_surface_2.txt'), 'w')
            f.write("# triax region\n")
            f.write("--------------\n")
            f.write(str(self.internal_surface_2.layer[0].polygon.__geo_interface__))
            f.write("\n\n")
            f.write("# resin region\n")
            f.write("--------------\n")
            f.write(str(self.internal_surface_2.layer[1].polygon.__geo_interface__))
            f.close()
        if self.internal_surface_3.exists():
            f = open(os.path.join(stn.station_path,'internal_surface_3.txt'), 'w')
            f.write("# triax region\n")
            f.write("--------------\n")
            f.write(str(self.internal_surface_3.layer[0].polygon.__geo_interface__))
            f.write("\n\n")
            f.write("# resin region\n")
            f.write("--------------\n")
            f.write(str(self.internal_surface_3.layer[1].polygon.__geo_interface__))
            f.close()
        if self.internal_surface_4.exists():
            f = open(os.path.join(stn.station_path,'internal_surface_4.txt'), 'w')
            f.write("# triax region\n")
            f.write("--------------\n")
            f.write(str(self.internal_surface_4.layer[0].polygon.__geo_interface__))
            f.write("\n\n")
            f.write("# resin region\n")
            f.write("--------------\n")
            f.write(str(self.internal_surface_4.layer[1].polygon.__geo_interface__))
            f.close()


class BiplaneStructure:
    """Define the biplane laminate schedule (internal dimensions)."""
    def __init__(self, h_RB, b_SC, h_SC, b_SW1_biax, b_SW1_foam, x2_SW1,
                 b_SW2_biax, b_SW2_foam, x2_SW2, b_SW3_biax, b_SW3_foam,
                 x2_SW3, b_TE_reinf, h_TE_reinf_uniax, h_TE_reinf_foam,
                 h_LE_panel, h_aft_panel_1, h_aft_panel_2, h_int_surf_triax,
                 h_int_surf_resin, h_ext_surf_triax, h_ext_surf_gelcoat,
                 h_RB_u, b_SC_u, h_SC_u, b_SW1_biax_u, b_SW1_foam_u, x2_SW1_u,
                 b_SW2_biax_u, b_SW2_foam_u, x2_SW2_u, b_SW3_biax_u,
                 b_SW3_foam_u, x2_SW3_u, b_TE_reinf_u, h_TE_reinf_uniax_u,
                 h_TE_reinf_foam_u, h_LE_panel_u, h_aft_panel_1_u,
                 h_aft_panel_2_u, h_int_surf_triax_u, h_int_surf_resin_u,
                 h_ext_surf_triax_u, h_ext_surf_gelcoat_u):
        self.lower_root_buildup = Part(np.nan, h_RB)
        self.lower_spar_cap = Part(b_SC, h_SC)
        self.lower_shear_web_1 = ShearWeb(b_SW1_biax, b_SW1_foam, x2_SW1)
        self.lower_shear_web_2 = ShearWeb(b_SW2_biax, b_SW2_foam, x2_SW2)
        self.lower_shear_web_3 = ShearWeb(b_SW3_biax, b_SW3_foam, x2_SW3)
        self.lower_TE_reinforcement = TE_Reinforcement(b_TE_reinf, h_TE_reinf_uniax, 
                                                 h_TE_reinf_foam)
        self.lower_LE_panel = Part(np.nan, h_LE_panel)
        self.lower_aft_panel_1 = Part(np.nan, h_aft_panel_1)
        self.lower_aft_panel_2 = Part(np.nan, h_aft_panel_2)
        self.lower_internal_surface = InternalSurface(np.nan, h_int_surf_triax, h_int_surf_resin)
        self.lower_external_surface = ExternalSurface(np.nan, h_ext_surf_triax, h_ext_surf_gelcoat)
        self.upper_root_buildup = Part(np.nan, h_RB_u)
        self.upper_spar_cap = Part(b_SC_u, h_SC_u)
        self.upper_shear_web_1 = ShearWeb(b_SW1_biax_u, b_SW1_foam_u, x2_SW1_u)
        self.upper_shear_web_2 = ShearWeb(b_SW2_biax_u, b_SW2_foam_u, x2_SW2_u)
        self.upper_shear_web_3 = ShearWeb(b_SW3_biax_u, b_SW3_foam_u, x2_SW3_u)
        self.upper_TE_reinforcement = TE_Reinforcement(b_TE_reinf_u,
            h_TE_reinf_uniax_u, h_TE_reinf_foam_u)
        self.upper_LE_panel = Part(np.nan, h_LE_panel_u)
        self.upper_aft_panel_1 = Part(np.nan, h_aft_panel_1_u)
        self.upper_aft_panel_2 = Part(np.nan, h_aft_panel_2_u)
        self.upper_internal_surface = InternalSurface(np.nan,
            h_int_surf_triax_u, h_int_surf_resin_u)
        self.upper_external_surface = ExternalSurface(np.nan,
            h_ext_surf_triax_u, h_ext_surf_gelcoat_u)

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
        s += "  --- INTERNAL SURFACE ---\n"
        s += "  " + str(self.lower_internal_surface) + '\n'
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
        s += "  --- INTERNAL SURFACE ---\n"
        s += "  " + str(self.upper_internal_surface) + '\n'
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