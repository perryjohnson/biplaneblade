"""A module for organizing geometrical data for a blade station.

Author: Perry Roth-Johnson
Last updated: August 6, 2013

"""


import os
import numpy as np
import matplotlib.pyplot as plt
import transformation as tf
reload(tf)
from coordinates import *
from airfoil import *
from structure import *
from shapely.geometry import Polygon
from descartes import PolygonPatch
# the descartes module translates shapely objects into matplotlib objects


class _Station:
    """Define a station for a wind turbine blade.

    The _Station base class is not intended for use.
    Use MonoplaneStation or BiplaneStation instead.

    This class also contains methods to split the airfoil curve into separate
    segments for each structural part: spar caps, shear webs, etc.

    Usage
    -----
    import pandas as pd
    import station as stn
    df = pd.read_csv('Sandia_blade.csv', index_col=0)
    s5 = stn._Station(df.ix[5])  # import station 5

    """
    logfile_name = 'station.log'
    number_of_stations = 0
    def __init__(self, stn_series, blade_path):
        """Create a new blade station.

        Parameters
        ---------
        stn_series : pandas.Series, properties for this station
        blade_path: string, the local target directory for storing blade data

        Attributes
        ----------
        .station_num : int, the blade station number
        .station_path : str, local directory for storing this station's data
        .logf : file handle, log file for all Station operations
        .coords
            .x1 : float, spanwise coordinate (meters)
            .x2 : float, edgewise coordinate (meters)
            .x3 : float, flapwise coordinate (meters)
        .airfoil
            .name : str, the airfoil name
            .filename : str, the airfoil filename
            .path : str, the airfoil path
            .pitch_axis : float, the chord fraction distance between the
                leading edge and the pitch axis (unitless)
            .chord : float, the chord length (meters)
            .twist : float, the twist about the x1 axis (degrees)
            .coords : numpy array, the airfoil coordinates (scaled by the
                .chord and .pitch_axis dimensions)
                [note: created by .read_coords()]
            .LE_index : int, the index of .coords for the leading edge of this
                airfoil
            .suction : numpy array, the (scaled) airfoil coordinates of the
                suction surface
            .pressure : numpy array, the (scaled) airfoil coordinates of the
                pressure surface
            .polygon : shapely.Polygon, a closed polygon representation of the
                entire airfoil surface
        .structure
            .root_buildup
                .base=np.nan
                .height : float, the root buildup height (meters)
            .spar_cap
                .base : float, the spar cap base (meters)
                .height : float, the spar cap height (meters)
                .left : float, the left edge -- chordwise coord (meters)
                .right : float, the right edge -- chordwise coord (meters)
                .lower_coords : numpy array, coordinates for the lower spar cap
                .upper_coords : numpy array, coordinates for the upper spar cap
            .shear_web_1
                .base : float, the shear web #1 total base (meters)
                .base_biax : float, the shear web #1 base for biax (meters)
                .base_foam : float, the shear web #1 base for foam (meters)
                .x2 : float, dist from pitch axis to edge of SW #1 (meters)
                .height=np.nan
                .left : float, the left edge -- chordwise coord (meters)
                .right : float, the right edge -- chordwise coord (meters)
                .cs_coords : numpy array, the 4 coordinates for the corners of
                    the cross-section of the shear web at this station, ordered
                    as [lower left, lower right, upper right, upper left]
            .shear_web_2
                .base : float, the shear web #2 total base (meters)
                .base_biax : float, the shear web #2 base for biax (meters)
                .base_foam : float, the shear web #2 base for foam (meters)
                .x2 : float, dist from pitch axis to edge of SW #2 (meters)
                .height=np.nan
                .left : float, the left edge -- chordwise coord (meters)
                .right : float, the right edge -- chordwise coord (meters)
                .cs_coords : numpy array, the 4 coordinates for the corners of
                    the cross-section of the shear web at this station, ordered
                    as [lower left, lower right, upper right, upper left]
            .shear_web_3
                .base : float, the shear web #3 total base (meters)
                .base_biax : float, the shear web #3 base for biax (meters)
                .base_foam : float, the shear web #3 base for foam (meters)
                .x2 : float, dist from pitch axis to edge of SW #3 (meters)
                .height=np.nan
                .left : float, the left edge -- chordwise coord (meters)
                .right : float, the right edge -- chordwise coord (meters)
                .cs_coords : numpy array, the 4 coordinates for the corners of
                    the cross-section of the shear web at this station, ordered
                    as [lower left, lower right, upper right, upper left]
            .TE_reinforcement
                .base : float, the trailing edge reinforcement base (meters)
                .height_uniax : float, the TE reinf height for uniax (meters)
                .height_foam : float, the TE reinf height for foam (meters)
                .height : float, the TE reinf total height (meters)
                .left : float, the left edge -- chordwise coord (meters)
                .right : float, the right edge -- chordwise coord (meters)
            .LE_panel
                .base=np.nan
                .height : the leading edge panel height (meters)
                .left : float, the left edge -- chordwise coord (meters)
                .right : float, the right edge -- chordwise coord (meters)
            .aft_panel
                .base=np.nan
                .height : the aft panel height (meters)
                .left : float, the left edge -- chordwise coord (meters)
                .right : float, the right edge -- chordwise coord (meters)
            .internal_surface
                .base=np.nan
                .height_triax : float, the internal surface height for triax (meters)
                .height_resin : float, the internal surface height for resin (meters)
                .height : float, the internal surface total height (meters)
            .external_surface
                .base=np.nan
                .height_triax : float, the external surface height for triax (meters)
                .height_gelcoat : float, the external surface height for gelcoat (meters)
                .height : float, the external surface total height (meters)

        Methods
        -------
        .structure
            .<Part>
                .exists() : bool, check if a Part exists at this station
        .airfoil
            .read_coords() : Read the airfoil coordinates and scale wrt the
                airfoil dims. Create a new attribute, <Station>.airfoil.coords,
                a numpy array of airfoil coordinates.
            .scale_coords()
            .rotate_coords()
            .split_at_LE_and_TE() : Split the airfoil curve into suction and
                pressure segments. Create new attributes: 
                <Station>.airfoil.LE_index, <Station>.airfoil.suction, and 
                <Station>.airfoil.pressure
            .plot_coords() : plot the airfoil coordinates of this station

        Usage
        -----    
        _Station(b._df.ix[5], 'sandia_blade')
        # this creates station #5 of the Sandia blade
        # _df is a pandas DataFrame containing properties of all blade stations
        # _df.ix[5] gets the Series object for station #5 from DataFrame _df
        # Note: _Stations are usually not created directly. New _Stations are
        # usually created by the _Blade class.

        """
        _Station.number_of_stations += 1
        self.station_num = _Station.number_of_stations
        self.station_path = os.path.join(blade_path, 'stn{0:02d}'.format(self.station_num))
        try:
            os.mkdir(self.station_path)
        except WindowsError:
            print "[WindowsError] The station path '{0}' already exists!".format(os.path.split(self.station_path)[-1])
        self.logf = open(_Station.logfile_name, "a")
        self.logf.write("............(Created blade station #{0})............\n".format(self.station_num))
        print " Created blade station #{0}".format(self.station_num)
        self.coords = Coordinates(stn_series['x1'], 
                                  stn_series['x2'], 
                                  stn_series['x3'])
        self.logf.write("****** COORDINATES ******\n")
        self.logf.write(str(self.coords) + '\n')
        self.logf.flush()
        self.logf.close()

    def __del__(self):
        _Station.number_of_stations = _Station.number_of_stations - 1
        print " Station deleted, and now _Station.number_of_stations = {0}".format(_Station.number_of_stations)

    def create_plot(self, legend_flag=False):
        """Create a plot for this station.

        Returns handles to the figure and its axes: (fig, axes)

        Several settings are applied ---------
        Title : Station #[num], [airfoil name], [num]% span
        Aspect ratio : equal
        Grid : on
        x-label : x2 [meters]
        y-label : x3 [meters]

        """
        af = self.airfoil
        fig, axes = plt.subplots()
        axes.set_title("Station #{0}, {1}, {2}% span".format(self.station_num, af.name, self.coords.x1))
        axes.set_aspect('equal')
        axes.grid('on')
        axes.set_xlabel('x2 [meters]')
        axes.set_ylabel('x3 [meters]')
        if legend_flag:
            axes.legend(loc='center')
        return (fig, axes)

    def show_plot(self):
        """Show the plot."""
        plt.show()

    def save_plot(self, fig):
        """Save the plot in the station path as a PNG file: stnXX.png"""
        fname = os.path.join(self.station_path, 'stn{0:02d}.png'.format(self.station_num))
        fig.savefig(fname)

    # note: use <Part>.exists() method to decide whether or not to split the airfoil curve for a particular Part
    # check these parts:
    # leading edge panel
    # shear webs (1,2,3)
    # spar caps
    # aft panels
    # TE reinforcement

    # note: keep implementing methods from airfoil_utils.py into this Station class!!!!


class MonoplaneStation(_Station):
    """Define a monoplane station for a wind turbine blade."""
    def __init__(self, stn_series, blade_path):
        """Create a new biplane station for a biplane blade."""
        _Station.__init__(self, stn_series, blade_path)
        self.type = 'monoplane'
        self.airfoil = MonoplaneAirfoil(
            name=stn_series['airfoil'],
            filename=stn_series['airfoil']+'.txt',
            chord=stn_series['chord'],
            pitch_axis=stn_series['pitch axis'],
            twist=stn_series['twist'])
        self.logf = open(_Station.logfile_name, "a")
        self.logf.write("****** AIRFOIL AND CHORD PROPERTIES ******\n")
        self.logf.write(str(self.airfoil) + '\n')
        self.structure = MonoplaneStructure(
            h_RB=stn_series['root buildup height'],
            b_SC=stn_series['spar cap base'],
            h_SC=stn_series['spar cap height'],
            b_SW1_biax=stn_series['shear web 1 base biax'],
            b_SW1_foam=stn_series['shear web 1 base foam'],
            x2_SW1=stn_series['shear web 1 x2'],
            b_SW2_biax=stn_series['shear web 2 base biax'],
            b_SW2_foam=stn_series['shear web 2 base foam'],
            x2_SW2=stn_series['shear web 2 x2'],
            b_SW3_biax=stn_series['shear web 3 base biax'],
            b_SW3_foam=stn_series['shear web 3 base foam'],
            x2_SW3=stn_series['shear web 3 x2'],
            b_TE_reinf=stn_series['TE reinf base'],
            h_TE_reinf_uniax=stn_series['TE reinf height uniax'],
            h_TE_reinf_foam=stn_series['TE reinf height foam'],
            h_LE_panel=stn_series['LE panel height'],
            h_aft_panel_1=stn_series['aft panel 1 height'],
            h_aft_panel_2=stn_series['aft panel 2 height'],
            h_int_surf_triax=stn_series['internal surface height triax'],
            h_int_surf_resin=stn_series['internal surface height resin'],
            h_ext_surf_triax=stn_series['external surface height triax'],
            h_ext_surf_gelcoat=stn_series['external surface height gelcoat'])
        self.logf.write("****** LAMINATE SCHEDULE ******\n")
        self.logf.write(str(self.structure) + '\n')
        self.logf.flush()
        self.logf.close()

    def find_all_part_cs_coords(self):
        """Find the corners of the cross-sections for each structural part.

        Saves cross-section coordinates (in meters) as the '.cs_coords' 
        attribute (a numpy array) within each Part instance (OOP style).

        NOTE: only shear webs have been implemented so far!

        """
        st = self.structure
        af = self.airfoil
        # if st.spar_cap.exists():
        #     st.spar_cap.left = -st.spar_cap.base/2.0
        #     st.spar_cap.right = st.spar_cap.base/2.0
        # if st.TE_reinforcement.exists():
        #     st.TE_reinforcement.left = -af.pitch_axis*af.chord+af.chord-st.TE_reinforcement.base
        #     st.TE_reinforcement.right = -af.pitch_axis*af.chord+af.chord
        if st.shear_web_1.exists():
            ((x1,y1),(x4,y4)) = af.find_part_edge_coords(st.shear_web_1.left)
            ((x2,y2),(x3,y3)) = af.find_part_edge_coords(st.shear_web_1.right)
            st.shear_web_1.cs_coords = np.array([[x1,y1],  # 1 (lower left)
                                                 [x2,y2],  # 2 (lower right)
                                                 [x3,y3],  # 3 (upper right)
                                                 [x4,y4]]) # 4 (upper left)
        if st.shear_web_2.exists():
            ((x1,y1),(x4,y4)) = af.find_part_edge_coords(st.shear_web_2.left)
            ((x2,y2),(x3,y3)) = af.find_part_edge_coords(st.shear_web_2.right)
            st.shear_web_2.cs_coords = np.array([[x1,y1],  # 1 (lower left)
                                                 [x2,y2],  # 2 (lower right)
                                                 [x3,y3],  # 3 (upper right)
                                                 [x4,y4]]) # 4 (upper left)
        if st.shear_web_3.exists():
            ((x1,y1),(x4,y4)) = af.find_part_edge_coords(st.shear_web_3.left)
            ((x2,y2),(x3,y3)) = af.find_part_edge_coords(st.shear_web_3.right)
            st.shear_web_3.cs_coords = np.array([[x1,y1],  # 1 (lower left)
                                                 [x2,y2],  # 2 (lower right)
                                                 [x3,y3],  # 3 (upper right)
                                                 [x4,y4]]) # 4 (upper left)
        # if st.LE_panel.exists():
        #     st.LE_panel.left = -af.pitch_axis*af.chord
        #     if st.shear_web_1.exists():
        #         st.LE_panel.right = st.shear_web_1.left
        #     elif st.spar_cap.exists():
        #         st.LE_panel.right = st.spar_cap.left
        #     else:
        #         st.LE_panel.right = np.nan
        #         raise Warning("'LE panel, right' is undefined for station #{0}".format(self.station_num))
        # if st.aft_panel.exists():
        #     if st.shear_web_2.exists():
        #         st.aft_panel.left = st.shear_web_2.right
        #     elif st.spar_cap.exists():
        #         st.aft_panel.left = st.spar_cap.right
        #     else:
        #         st.aft_panel.left = np.nan
        #         raise Warning("'aft panel, left' is undefined for station #{0}".format(self.station_num))
        #     if st.TE_reinforcement.exists():
        #         st.aft_panel.right = st.TE_reinforcement.left
        #     else:
        #         st.aft_panel.right = np.nan
        #         raise Warning("'aft panel, right' is undefined for station #{0}".format(self.station_num))

    def find_part_edges(self):
        """Find the edges of each structural part in this monoplane station.

        Saves coordinates (in meters) as '.left' and '.right' attributes
        (floats) within each Part instance (OOP style).

        """
        st = self.structure
        af = self.airfoil
        if st.spar_cap.exists():
            st.spar_cap.left = -st.spar_cap.base/2.0
            st.spar_cap.right = st.spar_cap.base/2.0
        if st.TE_reinforcement.exists():
            st.TE_reinforcement.left = -af.pitch_axis*af.chord+af.chord-st.TE_reinforcement.base
            st.TE_reinforcement.right = -af.pitch_axis*af.chord+af.chord
        if st.shear_web_1.exists():
            st.shear_web_1.right = st.shear_web_1.x2
            st.shear_web_1.left = st.shear_web_1.x2-st.shear_web_1.base
        if st.shear_web_2.exists():
            st.shear_web_2.left = st.shear_web_2.x2
            st.shear_web_2.right = st.shear_web_2.x2+st.shear_web_2.base
        if st.shear_web_3.exists():
            st.shear_web_3.left = st.shear_web_3.x2
            st.shear_web_3.right = st.shear_web_3.x2+st.shear_web_3.base
        if st.LE_panel.exists():
            st.LE_panel.left = -af.pitch_axis*af.chord
            if st.shear_web_1.exists():
                st.LE_panel.right = st.shear_web_1.left
            elif st.spar_cap.exists():
                st.LE_panel.right = st.spar_cap.left
            else:
                st.LE_panel.right = np.nan
                raise Warning("'LE panel, right' is undefined for station #{0}".format(self.station_num))
        if st.aft_panel_1.exists():
            if st.shear_web_2.exists():
                st.aft_panel_1.left = st.shear_web_2.right
            elif st.spar_cap.exists():
                st.aft_panel_1.left = st.spar_cap.right
            else:
                st.aft_panel_1.left = np.nan
                raise Warning("'aft panel 1, left' is undefined for station #{0}".format(self.station_num))
            if st.shear_web_3.exists():
                st.aft_panel_1.right = st.shear_web_3.left
            elif st.TE_reinforcement.exists():
                st.aft_panel_1.right = st.TE_reinforcement.left
            else:
                st.aft_panel_1.right = np.nan
                raise Warning("'aft panel 1, right' is undefined for station #{0}".format(self.station_num))
        if st.aft_panel_2.exists():
            if st.shear_web_3.exists():
                st.aft_panel_2.left = st.shear_web_3.right
            else:
                st.aft_panel_2.left = np.nan
                raise Warning("'aft panel 2, left' is undefined for station #{0}".format(self.station_num))
            if st.TE_reinforcement.exists():
                st.aft_panel_2.right = st.TE_reinforcement.left
            else:
                st.aft_panel_2.right = np.nan
                raise Warning("'aft panel 2, right' is undefined for station #{0}".format(self.station_num))

    def plot_part_edges(self, axes):
        """Plot color block for each structural part region.

        Each color block spans the plot from top to bottom.

        Uses coordinates saved as attributes within each Part instance
        (OOP style) by <Station>.find_part_edges().

        Must run <Station>.find_part_edges() first.

        KNOWN BUG: this doesn't work after rotating the airfoil coordinates.
        (This feature will not be implemented.)

        """
        st = self.structure
        try:
            if st.spar_cap.exists():
                axes.axvspan(st.spar_cap.left, st.spar_cap.right, facecolor='cyan', edgecolor='cyan', alpha=0.7)
            if st.TE_reinforcement.exists():
                axes.axvspan(st.TE_reinforcement.left, st.TE_reinforcement.right, facecolor='pink', edgecolor='pink', alpha=0.7)
            if st.LE_panel.exists():
                axes.axvspan(st.LE_panel.left, st.LE_panel.right, facecolor='magenta', edgecolor='magenta', alpha=0.7)
            if st.aft_panel_1.exists():
                axes.axvspan(st.aft_panel_1.left, st.aft_panel_1.right, facecolor='orange', edgecolor='orange', alpha=0.7)
            if st.aft_panel_2.exists():
                axes.axvspan(st.aft_panel_2.left, st.aft_panel_2.right, facecolor='orange', edgecolor='orange', alpha=0.7)
            if st.shear_web_1.exists():
                axes.axvspan(st.shear_web_1.left, st.shear_web_1.right, facecolor='green', edgecolor='green')
            if st.shear_web_2.exists():
                axes.axvspan(st.shear_web_2.left, st.shear_web_2.right, facecolor='green', edgecolor='green')
            if st.shear_web_3.exists():
                axes.axvspan(st.shear_web_3.left, st.shear_web_3.right, facecolor='green', edgecolor='green')
        except AttributeError:
            raise AttributeError("Part edges (.left and .right) have not been defined yet!\n  Try running <Station>.find_part_edges() first.")

    def plot_polygon(self, polygon, axes, face_color=(1,0,0),
        edge_color=(1,0,0), alpha=0.5):
        """Plot a polygon in a matplotlib figure."""
        patch = PolygonPatch(polygon, fc=face_color, ec=edge_color, alpha=alpha)
        axes.add_patch(patch)

    def erode_spar_cap_thickness(self):
        """Returns a polygon of the airfoil profile, eroded by the spar cap height"""
        t = self.structure.spar_cap.height
        p2 = self.airfoil.polygon.buffer(-t)
        return p2

    def cut_out_spar_cap_interior(self, interior_polygon):
        """Returns a polygon of the airfoil profile with the interior boundary of the spar cap cut out"""
        return self.airfoil.polygon.difference(interior_polygon)

    def spar_cap_bounding_box(self, y_boundary_buffer=1.2):
        """Returns a polygon of the bounding box that contains the spar caps.

        The points of the bounding box are labeled from 1 to 4 as:

        4---3
        |   |
        1---2

        """
        (minx, miny, maxx, maxy) = self.airfoil.polygon.bounds
        pt1 = (self.structure.spar_cap.left, miny*y_boundary_buffer)
        pt2 = (self.structure.spar_cap.right, miny*y_boundary_buffer)
        pt3 = (self.structure.spar_cap.right, maxy*y_boundary_buffer)
        pt4 = (self.structure.spar_cap.left, maxy*y_boundary_buffer)
        bounding_box = Polygon([pt1, pt2, pt3, pt4])
        return bounding_box

    def cut_out_spar_caps(self, airfoil_cutout, bounding_box):
        """Returns two polygons: the upper and lower spar caps.

        Spar caps are obtained by finding the intersection of two polygons:
        `airfoil_cutout` and `bounding_box`.

        """
        p4 = airfoil_cutout.intersection(bounding_box)
        # p4 is a `MultiPolygon` that contains two polygons
        # (p4 contains the upper and lower spar caps)
        # extract each polygon and convert them to separate patches
        # (otherwise, PolygonPatch will throw an error)
        scL = p4.geoms[0]  # lower spar cap
        scU = p4.geoms[1]  # upper spar cap
        return scL, scU

    def extract_and_plot_spar_caps(self, axes):
        """Extract the spar caps from the blade definition.

        Saves the spar cap polygon coordinates as attributes:
        <station>.structure.spar_cap
            .lower_coords : numpy array, lower spar cap coordinates
            .upper_coords : numpy array, upper spar cap coordinates

        """
        ip = self.erode_spar_cap_thickness()
        # self.plot_polygon(ip, axes, face_color='#6699cc', edge_color='#6699cc')
        ac = self.cut_out_spar_cap_interior(ip)
        # self.plot_polygon(ac, axes, alpha=0.7)
        bb = self.spar_cap_bounding_box()
        # self.plot_polygon(bb, axes, face_color=(0,1,0), edge_color=(0,1,0), alpha=0.3)
        (scL,scU) = self.cut_out_spar_caps(ac, bb)
        self.plot_polygon(scL, axes, face_color=(1,1,0), edge_color=(1,1,0),
            alpha=0.7)
        self.plot_polygon(scU, axes, face_color=(0,1,1), edge_color=(0,1,1),
            alpha=0.7)
        self.structure.spar_cap.lower_coords = np.array(
            scL.__geo_interface__['coordinates'][0])
        self.structure.spar_cap.upper_coords = np.array(
            scU.__geo_interface__['coordinates'][0])


class BiplaneStation(_Station):
    """Define a biplane station for a biplane wind turbine blade."""
    def __init__(self, stn_series, blade_path):
        """Create a new biplane station for a biplane blade."""
        _Station.__init__(self, stn_series, blade_path)
        self.type = 'biplane'
        self.airfoil = BiplaneAirfoil(
            name=stn_series['airfoil']+'_biplane',
            name_L=stn_series['airfoil'],
            filename_L=stn_series['airfoil']+'.txt',
            chord_L=stn_series['chord'],
            SW_ref_pt_L=stn_series['lower SW ref pt fraction'],
            name_U=stn_series['airfoil upper'],
            filename_U=stn_series['airfoil upper']+'.txt',
            chord_U=stn_series['chord'],
            SW_ref_pt_U=stn_series['upper SW ref pt fraction'],
            pitch_axis=stn_series['pitch axis'],
            twist=stn_series['twist'],
            gap_to_chord_ratio=stn_series['gap-to-chord ratio'],
            gap_fraction=stn_series['gap fraction'],
            stagger_to_chord_ratio=stn_series['stagger-to-chord ratio'])
        self.logf = open(_Station.logfile_name, "a")
        self.logf.write("****** AIRFOIL AND CHORD PROPERTIES ******\n")
        self.logf.write(str(self.airfoil) + '\n')
        self.structure = BiplaneStructure(
            h_RB=stn_series['root buildup height'],
            b_SC=stn_series['spar cap base'],
            h_SC=stn_series['spar cap height'],
            b_SW1_biax=stn_series['shear web 1 base biax'],
            b_SW1_foam=stn_series['shear web 1 base foam'],
            x2_SW1=stn_series['shear web 1 x2'],
            b_SW2_biax=stn_series['shear web 2 base biax'],
            b_SW2_foam=stn_series['shear web 2 base foam'],
            x2_SW2=stn_series['shear web 2 x2'],
            b_SW3_biax=stn_series['shear web 3 base biax'],
            b_SW3_foam=stn_series['shear web 3 base foam'],
            x2_SW3=stn_series['shear web 3 x2'],
            b_TE_reinf=stn_series['TE reinf base'],
            h_TE_reinf_uniax=stn_series['TE reinf height uniax'],
            h_TE_reinf_foam=stn_series['TE reinf height foam'],
            h_LE_panel=stn_series['LE panel height'],
            h_aft_panel_1=stn_series['aft panel 1 height'],
            h_aft_panel_2=stn_series['aft panel 2 height'],
            h_int_surf_triax=stn_series['internal surface height triax'],
            h_int_surf_resin=stn_series['internal surface height resin'],
            h_ext_surf_triax=stn_series['external surface height triax'],
            h_ext_surf_gelcoat=stn_series['external surface height gelcoat'],
            h_RB_u=stn_series['root buildup height upper'],
            b_SC_u=stn_series['spar cap base upper'],
            h_SC_u=stn_series['spar cap height upper'],
            b_SW1_biax_u=stn_series['shear web 1 base biax upper'],
            b_SW1_foam_u=stn_series['shear web 1 base foam upper'],
            x2_SW1_u=stn_series['shear web 1 x2 upper'],
            b_SW2_biax_u=stn_series['shear web 2 base biax upper'],
            b_SW2_foam_u=stn_series['shear web 2 base foam upper'],
            x2_SW2_u=stn_series['shear web 2 x2 upper'],
            b_SW3_biax_u=stn_series['shear web 3 base biax upper'],
            b_SW3_foam_u=stn_series['shear web 3 base foam upper'],
            x2_SW3_u=stn_series['shear web 3 x2 upper'],
            b_TE_reinf_u=stn_series['TE reinf base upper'],
            h_TE_reinf_uniax_u=stn_series['TE reinf height uniax upper'],
            h_TE_reinf_foam_u=stn_series['TE reinf height foam upper'],
            h_LE_panel_u=stn_series['LE panel height upper'],
            h_aft_panel_1_u=stn_series['aft panel 1 height upper'],
            h_aft_panel_2_u=stn_series['aft panel 2 height upper'],
            h_int_surf_triax_u=stn_series['internal surface height triax upper'],
            h_int_surf_resin_u=stn_series['internal surface height resin upper'],
            h_ext_surf_triax_u=stn_series['external surface height triax upper'],
            h_ext_surf_gelcoat_u=stn_series['external surface height gelcoat upper'])
        self.logf.write("****** LAMINATE SCHEDULE ******\n")
        self.logf.write(str(self.structure) + '\n')
        self.logf.flush()
        self.logf.close()

    def find_part_edges(self):
        """Find the edges of each structural part in this biplane station.

        Saves coordinates (in meters) as '.left' and '.right' attributes
        (floats) within each Part instance (OOP style).

        """
        st = self.structure
        af = self.airfoil
        # upper airfoil
        upper_refpt = -(af.pitch_axis*af.total_chord)+(af.upper_SW_ref_pt_fraction*af.upper_chord)
        if st.upper_spar_cap.exists():
            st.upper_spar_cap.left = upper_refpt - st.upper_spar_cap.base/2.0
            st.upper_spar_cap.right = upper_refpt + st.upper_spar_cap.base/2.0
        if st.upper_TE_reinforcement.exists():
            st.upper_TE_reinforcement.left = upper_refpt - af.upper_SW_ref_pt_fraction*af.upper_chord+af.upper_chord-st.upper_TE_reinforcement.base
            st.upper_TE_reinforcement.right = upper_refpt - af.upper_SW_ref_pt_fraction*af.upper_chord+af.upper_chord
        if st.upper_shear_web_1.exists():
            st.upper_shear_web_1.right = upper_refpt + st.upper_shear_web_1.x2
            st.upper_shear_web_1.left = upper_refpt + st.upper_shear_web_1.x2-st.upper_shear_web_1.base
        if st.upper_shear_web_2.exists():
            st.upper_shear_web_2.left = upper_refpt + st.upper_shear_web_2.x2
            st.upper_shear_web_2.right = upper_refpt + st.upper_shear_web_2.x2+st.upper_shear_web_2.base
        if st.upper_shear_web_3.exists():
            st.upper_shear_web_3.left = upper_refpt + st.upper_shear_web_3.x2
            st.upper_shear_web_3.right = upper_refpt + st.upper_shear_web_3.x2+st.upper_shear_web_3.base
        if st.upper_LE_panel.exists():
            st.upper_LE_panel.left = upper_refpt - af.upper_SW_ref_pt_fraction*af.upper_chord
            if st.upper_shear_web_1.exists():
                st.upper_LE_panel.right = st.upper_shear_web_1.left
            elif st.upper_spar_cap.exists():
                st.upper_LE_panel.right = st.upper_spar_cap.left
            else:
                st.upper_LE_panel.right = np.nan
                raise Warning("'LE panel, right' is undefined for station #{0}".format(self.station_num))
        if st.upper_aft_panel_1.exists():
            if st.upper_shear_web_2.exists():
                st.upper_aft_panel_1.left = st.upper_shear_web_2.right
            elif st.upper_spar_cap.exists():
                st.upper_aft_panel_1.left = st.upper_spar_cap.right
            else:
                st.upper_aft_panel_1.left = np.nan
                raise Warning("'aft panel 1, left' is undefined for station #{0}".format(self.station_num))
            if st.upper_shear_web_3.exists():
                st.upper_aft_panel_1.right = st.upper_shear_web_3.left
            elif st.upper_TE_reinforcement.exists():
                st.upper_aft_panel_1.right = st.upper_TE_reinforcement.left
            else:
                st.upper_aft_panel_1.right = np.nan
                raise Warning("'aft panel 1, right' is undefined for station #{0}".format(self.station_num))
        if st.upper_aft_panel_2.exists():
            if st.upper_shear_web_3.exists():
                st.upper_aft_panel_2.left = st.upper_shear_web_3.right
            else:
                st.upper_aft_panel_2.left = np.nan
                raise Warning("'aft panel 2, left' is undefined for station #{0}".format(self.station_num))
            if st.upper_TE_reinforcement.exists():
                st.upper_aft_panel_2.right = st.upper_TE_reinforcement.left
            else:
                st.upper_aft_panel_2.right = np.nan
                raise Warning("'aft panel 2, right' is undefined for station #{0}".format(self.station_num))
        # lower airfoil
        lower_refpt = -(af.pitch_axis*af.total_chord)+(af.stagger)+(af.lower_SW_ref_pt_fraction*af.lower_chord)
        if st.lower_spar_cap.exists():
            st.lower_spar_cap.left = lower_refpt - st.lower_spar_cap.base/2.0
            st.lower_spar_cap.right = lower_refpt + st.lower_spar_cap.base/2.0
        if st.lower_TE_reinforcement.exists():
            st.lower_TE_reinforcement.left = lower_refpt - af.lower_SW_ref_pt_fraction*af.lower_chord+af.lower_chord-st.lower_TE_reinforcement.base
            st.lower_TE_reinforcement.right = lower_refpt - af.lower_SW_ref_pt_fraction*af.lower_chord+af.lower_chord
        if st.lower_shear_web_1.exists():
            st.lower_shear_web_1.right = lower_refpt + st.lower_shear_web_1.x2
            st.lower_shear_web_1.left = lower_refpt + st.lower_shear_web_1.x2-st.lower_shear_web_1.base
        if st.lower_shear_web_2.exists():
            st.lower_shear_web_2.left = lower_refpt + st.lower_shear_web_2.x2
            st.lower_shear_web_2.right = lower_refpt + st.lower_shear_web_2.x2+st.lower_shear_web_2.base
        if st.lower_shear_web_3.exists():
            st.lower_shear_web_3.left = lower_refpt + st.lower_shear_web_3.x2
            st.lower_shear_web_3.right = lower_refpt + st.lower_shear_web_3.x2+st.lower_shear_web_3.base
        if st.lower_LE_panel.exists():
            st.lower_LE_panel.left = lower_refpt - af.lower_SW_ref_pt_fraction*af.lower_chord
            if st.lower_shear_web_1.exists():
                st.lower_LE_panel.right = st.lower_shear_web_1.left
            elif st.lower_spar_cap.exists():
                st.lower_LE_panel.right = st.lower_spar_cap.left
            else:
                st.lower_LE_panel.right = np.nan
                raise Warning("'LE panel, right' is undefined for station #{0}".format(self.station_num))
        if st.lower_aft_panel_1.exists():
            if st.lower_shear_web_2.exists():
                st.lower_aft_panel_1.left = st.lower_shear_web_2.right
            elif st.lower_spar_cap.exists():
                st.lower_aft_panel_1.left = st.lower_spar_cap.right
            else:
                st.lower_aft_panel_1.left = np.nan
                raise Warning("'aft panel 1, left' is undefined for station #{0}".format(self.station_num))
            if st.lower_shear_web_3.exists():
                st.lower_aft_panel_1.right = st.lower_shear_web_3.left
            elif st.lower_TE_reinforcement.exists():
                st.lower_aft_panel_1.right = st.lower_TE_reinforcement.left
            else:
                st.lower_aft_panel_1.right = np.nan
                raise Warning("'aft panel 1, right' is undefined for station #{0}".format(self.station_num))
        if st.lower_aft_panel_2.exists():
            if st.lower_shear_web_3.exists():
                st.lower_aft_panel_2.left = st.lower_shear_web_3.right
            else:
                st.lower_aft_panel_2.left = np.nan
                raise Warning("'aft panel 2, left' is undefined for station #{0}".format(self.station_num))
            if st.lower_TE_reinforcement.exists():
                st.lower_aft_panel_2.right = st.lower_TE_reinforcement.left
            else:
                st.lower_aft_panel_2.right = np.nan
                raise Warning("'aft panel 2, right' is undefined for station #{0}".format(self.station_num))

    def plot_part_edges(self, axes):
        """Plot color block for each structural part region.

        Each color block spans the plot from top to bottom.

        Uses coordinates saved as attributes within each Part instance
        (OOP style) by <Station>.find_part_edges().

        Must run <Station>.find_part_edges() first.

        KNOWN BUG: this doesn't work after rotating the airfoil coordinates.
        (This feature will not be implemented.)

        """
        st = self.structure
        # upper airfoil
        try:
            if st.upper_spar_cap.exists():
                axes.axvspan(st.upper_spar_cap.left, st.upper_spar_cap.right, ymin=0.5, facecolor='cyan', edgecolor='cyan', alpha=0.7)
            if st.upper_TE_reinforcement.exists():
                axes.axvspan(st.upper_TE_reinforcement.left, st.upper_TE_reinforcement.right, ymin=0.5, facecolor='pink', edgecolor='pink', alpha=0.7)
            if st.upper_LE_panel.exists():
                axes.axvspan(st.upper_LE_panel.left, st.upper_LE_panel.right, ymin=0.5, facecolor='magenta', edgecolor='magenta', alpha=0.7)
            if st.upper_aft_panel_1.exists():
                axes.axvspan(st.upper_aft_panel_1.left, st.upper_aft_panel_1.right, ymin=0.5, facecolor='orange', edgecolor='orange', alpha=0.7)
            if st.upper_aft_panel_2.exists():
                axes.axvspan(st.upper_aft_panel_2.left, st.upper_aft_panel_2.right, ymin=0.5, facecolor='orange', edgecolor='orange', alpha=0.7)
            if st.upper_shear_web_1.exists():
                axes.axvspan(st.upper_shear_web_1.left, st.upper_shear_web_1.right, ymin=0.5, facecolor='green', edgecolor='green')
            if st.upper_shear_web_2.exists():
                axes.axvspan(st.upper_shear_web_2.left, st.upper_shear_web_2.right, ymin=0.5, facecolor='green', edgecolor='green')
            if st.upper_shear_web_3.exists():
                axes.axvspan(st.upper_shear_web_3.left, st.upper_shear_web_3.right, ymin=0.5, facecolor='green', edgecolor='green')
        except AttributeError:
            raise AttributeError("Part edges (.left and .right) have not been defined yet!\n  Try running <Station>.find_part_edges() first.")
        # lower airfoil
        try:
            if st.lower_spar_cap.exists():
                axes.axvspan(st.lower_spar_cap.left, st.lower_spar_cap.right, ymax=0.5, facecolor='cyan', edgecolor='cyan', alpha=0.7)
            if st.lower_TE_reinforcement.exists():
                axes.axvspan(st.lower_TE_reinforcement.left, st.lower_TE_reinforcement.right, ymax=0.5, facecolor='pink', edgecolor='pink', alpha=0.7)
            if st.lower_LE_panel.exists():
                axes.axvspan(st.lower_LE_panel.left, st.lower_LE_panel.right, ymax=0.5, facecolor='magenta', edgecolor='magenta', alpha=0.7)
            if st.lower_aft_panel_1.exists():
                axes.axvspan(st.lower_aft_panel_1.left, st.lower_aft_panel_1.right, ymax=0.5, facecolor='orange', edgecolor='orange', alpha=0.7)
            if st.lower_aft_panel_2.exists():
                axes.axvspan(st.lower_aft_panel_2.left, st.lower_aft_panel_2.right, ymax=0.5, facecolor='orange', edgecolor='orange', alpha=0.7)
            if st.lower_shear_web_1.exists():
                axes.axvspan(st.lower_shear_web_1.left, st.lower_shear_web_1.right, ymax=0.5, facecolor='green', edgecolor='green')
            if st.lower_shear_web_2.exists():
                axes.axvspan(st.lower_shear_web_2.left, st.lower_shear_web_2.right, ymax=0.5, facecolor='green', edgecolor='green')
            if st.lower_shear_web_3.exists():
                axes.axvspan(st.lower_shear_web_3.left, st.lower_shear_web_3.right, ymax=0.5, facecolor='green', edgecolor='green')
        except AttributeError:
            raise AttributeError("Part edges (.left and .right) have not been defined yet!\n  Try running <Station>.find_part_edges() first.")

    def find_all_part_cs_coords(self):
        """Find the corners of the cross-sections for each structural part.

        Saves cross-section coordinates (in meters) as the '.cs_coords' 
        attribute (a numpy array) within each Part instance (OOP style).

        NOTE: only shear webs have been implemented so far!

        """
        st = self.structure
        af = self.airfoil
        # if st.spar_cap.exists():
        #     st.spar_cap.left = -st.spar_cap.base/2.0
        #     st.spar_cap.right = st.spar_cap.base/2.0
        # if st.TE_reinforcement.exists():
        #     st.TE_reinforcement.left = -af.pitch_axis*af.chord+af.chord-st.TE_reinforcement.base
        #     st.TE_reinforcement.right = -af.pitch_axis*af.chord+af.chord
        # lower airfoil
        if st.lower_shear_web_1.exists():
            ((x1,y1),(x4,y4)) = af.find_part_edge_coords(st.lower_shear_web_1.left, airfoil='lower')
            ((x2,y2),(x3,y3)) = af.find_part_edge_coords(st.lower_shear_web_1.right, airfoil='lower')
            st.lower_shear_web_1.cs_coords = np.array([[x1,y1],  # 1 (lower left)
                                                 [x2,y2],  # 2 (lower right)
                                                 [x3,y3],  # 3 (upper right)
                                                 [x4,y4]]) # 4 (upper left)
        if st.lower_shear_web_2.exists():
            ((x1,y1),(x4,y4)) = af.find_part_edge_coords(st.lower_shear_web_2.left, airfoil='lower')
            ((x2,y2),(x3,y3)) = af.find_part_edge_coords(st.lower_shear_web_2.right, airfoil='lower')
            st.lower_shear_web_2.cs_coords = np.array([[x1,y1],  # 1 (lower left)
                                                 [x2,y2],  # 2 (lower right)
                                                 [x3,y3],  # 3 (upper right)
                                                 [x4,y4]]) # 4 (upper left)
        if st.lower_shear_web_3.exists():
            ((x1,y1),(x4,y4)) = af.find_part_edge_coords(st.lower_shear_web_3.left, airfoil='lower')
            ((x2,y2),(x3,y3)) = af.find_part_edge_coords(st.lower_shear_web_3.right, airfoil='lower')
            st.lower_shear_web_3.cs_coords = np.array([[x1,y1],  # 1 (lower left)
                                                 [x2,y2],  # 2 (lower right)
                                                 [x3,y3],  # 3 (upper right)
                                                 [x4,y4]]) # 4 (upper left)
        # upper airfoil
        if st.upper_shear_web_1.exists():
            ((x1,y1),(x4,y4)) = af.find_part_edge_coords(st.upper_shear_web_1.left, airfoil='upper')
            ((x2,y2),(x3,y3)) = af.find_part_edge_coords(st.upper_shear_web_1.right, airfoil='upper')
            st.upper_shear_web_1.cs_coords = np.array([[x1,y1],  # 1 (lower left)
                                                 [x2,y2],  # 2 (lower right)
                                                 [x3,y3],  # 3 (upper right)
                                                 [x4,y4]]) # 4 (upper left)
        if st.upper_shear_web_2.exists():
            ((x1,y1),(x4,y4)) = af.find_part_edge_coords(st.upper_shear_web_2.left, airfoil='upper')
            ((x2,y2),(x3,y3)) = af.find_part_edge_coords(st.upper_shear_web_2.right, airfoil='upper')
            st.upper_shear_web_2.cs_coords = np.array([[x1,y1],  # 1 (lower left)
                                                 [x2,y2],  # 2 (lower right)
                                                 [x3,y3],  # 3 (upper right)
                                                 [x4,y4]]) # 4 (upper left)
        if st.upper_shear_web_3.exists():
            ((x1,y1),(x4,y4)) = af.find_part_edge_coords(st.upper_shear_web_3.left, airfoil='upper')
            ((x2,y2),(x3,y3)) = af.find_part_edge_coords(st.upper_shear_web_3.right, airfoil='upper')
            st.upper_shear_web_3.cs_coords = np.array([[x1,y1],  # 1 (lower left)
                                                 [x2,y2],  # 2 (lower right)
                                                 [x3,y3],  # 3 (upper right)
                                                 [x4,y4]]) # 4 (upper left)
        # if st.LE_panel.exists():
        #     st.LE_panel.left = -af.pitch_axis*af.chord
        #     if st.shear_web_1.exists():
        #         st.LE_panel.right = st.shear_web_1.left
        #     elif st.spar_cap.exists():
        #         st.LE_panel.right = st.spar_cap.left
        #     else:
        #         st.LE_panel.right = np.nan
        #         raise Warning("'LE panel, right' is undefined for station #{0}".format(self.station_num))
        # if st.aft_panel.exists():
        #     if st.shear_web_2.exists():
        #         st.aft_panel.left = st.shear_web_2.right
        #     elif st.spar_cap.exists():
        #         st.aft_panel.left = st.spar_cap.right
        #     else:
        #         st.aft_panel.left = np.nan
        #         raise Warning("'aft panel, left' is undefined for station #{0}".format(self.station_num))
        #     if st.TE_reinforcement.exists():
        #         st.aft_panel.right = st.TE_reinforcement.left
        #     else:
        #         st.aft_panel.right = np.nan
        #         raise Warning("'aft panel, right' is undefined for station #{0}".format(self.station_num))
