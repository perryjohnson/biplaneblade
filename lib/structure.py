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


import numpy as np
import matplotlib.pyplot as plt
import layer as l
from math import isnan
from shapely.geometry import Polygon
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
        self.layer = []     # assigned later by <part>.create_layers()
    
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
        # self.layer = [None, None]
    
    def __str__(self):
        return """base:    {0:6.4f} (meters)
height:  {1:6.4f} (meters)
|-> height_triax:  {2:6.4f} (meters)
|-> height_gelcoat:  {3:6.4f} (meters)""".format(self.base, self.height,
    self.height_triax, self.height_gelcoat)
    
    def create_layers(self, debug_flag=False):
        """Create the gelcoat and triax layers in the external surface."""
        st = self.parent_structure
        if self.exists():
            # create the gelcoat layer
            af = st.parent_station.airfoil
            b = st.parent_station.parent_blade
            op_gelcoat = af.polygon  # outer profile is the airfoil profile
            ip_gelcoat = op_gelcoat.buffer(-self.height_gelcoat)
            polygon_gelcoat = op_gelcoat.difference(ip_gelcoat)
            self.layer.append(l.Layer(polygon_gelcoat,
                b.dict_of_materials['gelcoat'], parent_part=self))
            assert self.layer[0].polygon.geom_type == 'Polygon'
            st.list_of_polygons.append(self.layer[0].polygon)
            # create the triax layer
            op_triax = ip_gelcoat  # outer profile is the gelcoat inner profile
            ip_triax = op_triax.buffer(-self.height_triax)
            polygon_triax = op_triax.difference(ip_triax)
            self.layer.append(l.Layer(polygon_triax,
                b.dict_of_materials['triaxial GFRP'], parent_part=self))
            assert self.layer[1].polygon.geom_type == 'Polygon'
            st.list_of_polygons.append(self.layer[1].polygon)
        else:
            if debug_flag:
                print " The external surface for Station #{0} does not exist!\n  No layers created.".format(st.parent_station.station_num)


class RootBuildup(Part):
    """Define triax dimensions of the root buildup."""
    def create_layers(self, debug_flag=False):
        """Create the triax layer in the root buildup."""
        st = self.parent_structure
        if self.exists():
            # create the triax layer
            af = st.parent_station.airfoil
            b = st.parent_station.parent_blade
            op = af.polygon.buffer(-st.external_surface.height)
            ip = op.buffer(-self.height)
            p = op.difference(ip)
            self.layer.append(l.Layer(p, b.dict_of_materials['triaxial GFRP'],
                parent_part=self))
            assert self.layer[0].polygon.geom_type == 'Polygon'
            st.list_of_polygons.append(self.layer[0].polygon)
        else:
            if debug_flag:
                print " The root buildup for Station #{0} does not exist!\n  No layers created.".format(st.parent_station.station_num)


class LE_Panel(Part):
    """Define foam dimensions of the leading edge panel."""
    def create_layers(self, debug_flag=False):
        """Create the foam layer in the leading edge panel."""
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
            self.layer.append(l.Layer(p, b.dict_of_materials['foam'],
                parent_part=self))
            assert self.layer[0].polygon.geom_type == 'Polygon'
            st.list_of_polygons.append(self.layer[0].polygon)
        else:
            if debug_flag:
                print " The LE panel for Station #{0} does not exist!\n  No layers created.".format(st.parent_station.station_num)


class SparCap(Part):
    """Define uniax dimensions of the lower and upper spar caps."""
    def create_layers(self, debug_flag=False):
        """Create the uniax layers in the lower and upper spar caps.

        <spar_cap>.layer[0] is the lower spar cap
        <spar_cap>.layer[1] is the upper spar cap

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
            self.layer.append(l.Layer(pl, b.dict_of_materials['uniaxial GFRP'],
                parent_part=self))
            assert self.layer[0].polygon.geom_type == 'Polygon'
            st.list_of_polygons.append(self.layer[0].polygon)
            # 8. add the upper spar cap
            self.layer.append(l.Layer(pu, b.dict_of_materials['uniaxial GFRP'],
                parent_part=self))
            assert self.layer[1].polygon.geom_type == 'Polygon'
            st.list_of_polygons.append(self.layer[1].polygon)
        else:
            if debug_flag:
                print " The spar caps for Station #{0} do not exist!\n  No layers created.".format(st.parent_station.station_num)


class AftPanel(Part):
    """Define foam dimensions of the lower and upper aft panels."""
    def create_layers(self, debug_flag=False):
        """Create the foam layers in the lower and upper aft panels.

        <aft_panel>.layer[0] is the lower aft panel
        <aft_panel>.layer[1] is the upper aft panel

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
            self.layer.append(l.Layer(pl, b.dict_of_materials['foam'],
                parent_part=self))
            assert self.layer[0].polygon.geom_type == 'Polygon'
            st.list_of_polygons.append(self.layer[0].polygon)
            # 8. add the upper aft panel
            self.layer.append(l.Layer(pu, b.dict_of_materials['foam'],
                parent_part=self))
            assert self.layer[1].polygon.geom_type == 'Polygon'
            st.list_of_polygons.append(self.layer[1].polygon)
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
        <TE_reinforcement>.layer[0] is made of uniax
        <TE_reinforcement>.layer[1] is made of foam ... (optional layer)

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
            self.layer.append(l.Layer(polygon_uniax,
                b.dict_of_materials['uniaxial GFRP'], parent_part=self))
            assert self.layer[0].polygon.geom_type == 'Polygon'
            st.list_of_polygons.append(self.layer[0].polygon)
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
                self.layer.append(l.Layer(polygon_foam,
                    b.dict_of_materials['foam'], parent_part=self))
                assert self.layer[1].polygon.geom_type == 'Polygon'
                st.list_of_polygons.append(self.layer[1].polygon)
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

        <shear_web>.layer[0] is the left biax layer
        <shear_web>.layer[1] is the foam layer
        <shear_web>.layer[2] is the right biax layer

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
            self.layer.append(l.Layer(p_left_biax,
                b.dict_of_materials['biaxial GFRP'], parent_part=self))
            assert self.layer[0].polygon.geom_type == 'Polygon'
            st.list_of_polygons.append(self.layer[0].polygon)
            # 5. add the foam layer
            self.layer.append(l.Layer(p_foam, b.dict_of_materials['foam'],
                parent_part=self))
            assert self.layer[1].polygon.geom_type == 'Polygon'
            st.list_of_polygons.append(self.layer[1].polygon)
            # 6. add the right biax layer
            self.layer.append(l.Layer(p_right_biax,
                b.dict_of_materials['biaxial GFRP'], parent_part=self))
            assert self.layer[2].polygon.geom_type == 'Polygon'
            st.list_of_polygons.append(self.layer[2].polygon)
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
        """Create the triax and resin layers in the internal surface."""
        st = self.parent_structure
        if self.exists():
            b = st.parent_station.parent_blade
            # triax region
            op_triax = self.interior_loop(merged_polygon)
            ip_triax = op_triax.buffer(-self.height_triax)
            polygon_triax = op_triax.difference(ip_triax)
            self.layer.append(l.Layer(polygon_triax,
                b.dict_of_materials['triaxial GFRP'], parent_part=self))
            assert self.layer[0].polygon.geom_type == 'Polygon'
            st.list_of_polygons.append(self.layer[0].polygon)
            # resin region
            op_resin = ip_triax
            ip_resin = op_resin.buffer(-self.height_resin)
            polygon_resin = op_resin.difference(ip_resin)
            self.layer.append(l.Layer(polygon_resin,
                b.dict_of_materials['resin'], parent_part=self))
            assert self.layer[1].polygon.geom_type == 'Polygon'
            st.list_of_polygons.append(self.layer[1].polygon)
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
        self.list_of_polygons = []
        self.area = None

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
        try:
            p = cascaded_union(self.list_of_polygons)
        except ValueError:
            # gather all the parts
            p = self.external_surface.layer[0].polygon
            p = p.union(self.external_surface.layer[1].polygon)
            if self.root_buildup.exists():
                RB = self.root_buildup.layer[0].polygon
                p = p.union(RB)
            if self.LE_panel.exists():
                LE = self.LE_panel.layer[0].polygon
                p = p.union(LE)
            if self.spar_cap.exists():
                sc_l = self.spar_cap.layer[0].polygon
                try:
                    p = p.union(sc_l)
                except TopologicalError:
                    print " [Warning] could not merge lower spar cap in Station #{0} ... skipping!".format(self.parent_station.station_num)
                sc_u = self.spar_cap.layer[1].polygon
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
                aft1_u = self.aft_panel_1.layer[0].polygon
                aft1_l = self.aft_panel_1.layer[1].polygon
                p = p.union(aft1_u)
                p = p.union(aft1_l)
            if self.aft_panel_2.exists():
                aft2_u = self.aft_panel_2.layer[0].polygon
                aft2_l = self.aft_panel_2.layer[1].polygon
                p = p.union(aft2_u)
                p = p.union(aft2_l)
            if self.TE_reinforcement.exists():
                TE_uniax = self.TE_reinforcement.layer[0].polygon
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
                    TE_foam = self.TE_reinforcement.layer[1].polygon
                    try:
                        p = p.union(TE_foam)
                    except TopologicalError:
                        print " [Warning] could not merge foam layer of TE reinforcement in Station #{0} ... skipping!".format(self.parent_station.station_num)
                except ValueError:
                    print " foam layer of TE reinforcement does not exist in Station #{0}".format(self.parent_station.station_num)
            if self.shear_web_1.exists():
                sw1 = self.shear_web_1.layer[0].polygon.union(self.shear_web_1.layer[1].polygon)
                sw1 = sw1.union(self.shear_web_1.layer[2].polygon)
                p = p.union(sw1)
            if self.shear_web_2.exists():
                sw2 = self.shear_web_2.layer[0].polygon.union(self.shear_web_2.layer[1].polygon)
                sw2 = sw2.union(self.shear_web_2.layer[2].polygon)
                p = p.union(sw2)
            if self.shear_web_3.exists():
                sw3 = self.shear_web_3.layer[0].polygon.union(self.shear_web_3.layer[1].polygon)
                sw3 = sw3.union(self.shear_web_3.layer[2].polygon)
                p = p.union(sw3)
            # if merge_internal_surface:
            #     # also merge the internal surface with the other structural parts
            #     if self.internal_surface_1.exists():
            #         IS1 = self.internal_surface_1.polygon_triax
            #         IS1 = IS1.union(self.internal_surface_1.polygon_resin)
            #         list_of_parts.append(IS1)
            #     if self.internal_surface_2.exists():
            #         IS2 = self.internal_surface_2.polygon_triax
            #         IS2 = IS2.union(self.internal_surface_2.polygon_resin)
            #         list_of_parts.append(IS2)
            #     if self.internal_surface_3.exists():
            #         IS3 = self.internal_surface_3.polygon_triax
            #         IS3 = IS3.union(self.internal_surface_3.polygon_resin)
            #         list_of_parts.append(IS3)
            #     if self.internal_surface_4.exists():
            #         IS4 = self.internal_surface_4.polygon_triax
            #         IS4 = IS4.union(self.internal_surface_4.polygon_resin)
            #         list_of_parts.append(IS4)
            # try merging everything again
            # try:
            #     p = cascaded_union(list_of_parts)
            # except ValueError:
            #     raise Warning("Could not merge all polygons in Station #{0}".format(stn.station_num))
        if plot_flag:
            # plot the merged polygon
            patch2 = PolygonPatch(p, fc='#4000FF', ec = '#000000', alpha=0.8)
            ax.add_patch(patch2)
            plt.show()
        return p

    def calculate_area(self):
        """Add the area of all polygons in this station, save it to self.area."""
        a = 0
        for p in self.list_of_polygons:
            a += p.area
        self.area = a
        return a

    def calculate_all_percent_areas(self):
        """Calculate the percent areas of all parts in this station."""
        print " ----- STATION #{0} -----".format(self.parent_station.station_num)
        if self.external_surface.exists():
            print "  {0:5.1%} area, external surface, gelcoat".format(self.external_surface.layer[0].area_fraction())
            print "  {0:5.1%} area, external surface, triax".format(self.external_surface.layer[1].area_fraction())
        if self.root_buildup.exists():
            print "  {0:5.1%} area, root buildup".format(self.root_buildup.layer[0].area_fraction())
        if self.spar_cap.exists():
            print "  {0:5.1%} area, spar cap, lower".format(self.spar_cap.layer[0].area_fraction())
            print "  {0:5.1%} area, spar cap, upper".format(self.spar_cap.layer[1].area_fraction())
        if self.aft_panel_1.exists():
            print "  {0:5.1%} area, aft panel 1, lower".format(self.aft_panel_1.layer[0].area_fraction())
            print "  {0:5.1%} area, aft panel 1, upper".format(self.aft_panel_1.layer[1].area_fraction())
        if self.aft_panel_2.exists():
            print "  {0:5.1%} area, aft panel 2, lower".format(self.aft_panel_2.layer[0].area_fraction())
            print "  {0:5.1%} area, aft panel 2, upper".format(self.aft_panel_2.layer[1].area_fraction())
        if self.LE_panel.exists():
            print "  {0:5.1%} area, LE panel".format(self.LE_panel.layer[0].area_fraction())
        if self.shear_web_1.exists():
            print "  {0:5.1%} area, shear web 1, left biax".format(self.shear_web_1.layer[0].area_fraction())
            print "  {0:5.1%} area, shear web 1, foam".format(self.shear_web_1.layer[1].area_fraction())
            print "  {0:5.1%} area, shear web 1, right biax".format(self.shear_web_1.layer[2].area_fraction())
        if self.shear_web_2.exists():
            print "  {0:5.1%} area, shear web 2, left biax".format(self.shear_web_2.layer[0].area_fraction())
            print "  {0:5.1%} area, shear web 2, foam".format(self.shear_web_2.layer[1].area_fraction())
            print "  {0:5.1%} area, shear web 2, right biax".format(self.shear_web_2.layer[2].area_fraction())
        if self.shear_web_3.exists():
            print "  {0:5.1%} area, shear web 3, left biax".format(self.shear_web_3.layer[0].area_fraction())
            print "  {0:5.1%} area, shear web 3, foam".format(self.shear_web_3.layer[1].area_fraction())
            print "  {0:5.1%} area, shear web 3, right biax".format(self.shear_web_3.layer[2].area_fraction())
        if self.TE_reinforcement.exists():
            print "  {0:5.1%} area, TE reinforcement, uniax".format(self.TE_reinforcement.layer[0].area_fraction())
            try:
                print "  {0:5.1%} area, TE reinforcement, foam".format(self.TE_reinforcement.layer[1].area_fraction())
            except IndexError:
                pass
        if self.internal_surface_1.exists():
            print "  {0:5.1%} area, internal surface 1, triax".format(self.internal_surface_1.layer[0].area_fraction())
            print "  {0:5.1%} area, internal surface 1, resin".format(self.internal_surface_1.layer[1].area_fraction())
        if self.internal_surface_2.exists():
            print "  {0:5.1%} area, internal surface 2, triax".format(self.internal_surface_2.layer[0].area_fraction())
            print "  {0:5.1%} area, internal surface 2, resin".format(self.internal_surface_2.layer[1].area_fraction())
        if self.internal_surface_3.exists():
            print "  {0:5.1%} area, internal surface 3, triax".format(self.internal_surface_3.layer[0].area_fraction())
            print "  {0:5.1%} area, internal surface 3, resin".format(self.internal_surface_3.layer[1].area_fraction())
        if self.internal_surface_4.exists():
            print "  {0:5.1%} area, internal surface 4, triax".format(self.internal_surface_4.layer[0].area_fraction())
            print "  {0:5.1%} area, internal surface 4, resin".format(self.internal_surface_4.layer[1].area_fraction())


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