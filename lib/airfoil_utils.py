"""A module that splits an airfoil curve into separate segments for each
structural component: spar caps, shear webs, aft panels, etc.

*****
Note: the `airfoil_utils` module will be deprecated. All methods from 
`airfoil_utils` will be moved into the `Station` class in `station.py`
*****

Author: Perry Roth-Johnson
Last updated: July 24, 2013

"""


import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as ipl


dtype = [('x', 'f8'), ('y', 'f8')]


class SegmentedAirfoil:
    """Splits airfoil coordinates into segments for each structural component.

    Usage:
    import blade as bl
    import airfoil_utils as af
    b = bl.Blade('Sandia blade SNL100-00', 'sandia_blade')
    b.copy_all_airfoil_coords()
    stn26 = b.list_of_stations[25]
    a26 = af.Airfoil(stn26)
    a26.plot_all_segments()

    """
    def __init__(self, station):
        # set some global parameters
        self.name = station.airfoil.name
        # global parameters ---------------------------------------------------
        self.set_airfoil_dims(station.airfoil)
        self.set_structure_dims(station.structure)
        # read in the coordinates ---------------------------------------------
        self.read_airfoil_coords(station.airfoil.path)

    def set_airfoil_dims(self, airfoil):
        """Set the airfoil dimensions (chord, pitch axis)."""
        self.chord = airfoil.chord
        self.pitch_axis = airfoil.pitch_axis

    def set_structure_dims(self, structure):
        """Set the structural dimensions (e.g. spar cap height)."""
        self.hSC = structure.spar_cap.height
        self.bSC = structure.spar_cap.base
        self.bSW = structure.shear_web_1.base
        self.bSWbiax = structure.shear_web_1.base_biax
        self.bSWfoam = structure.shear_web_1.base_foam
        # **TODO** implement all three shear webs! ****************************
        # self.bSWbiax = structure.shear_web_2.base_biax
        # self.bSWfoam = structure.shear_web_2.base_foam
        # self.bSWbiax = structure.shear_web_3.base_biax
        # self.bSWfoam = structure.shear_web_3.base_foam
        self.bTEreinf = structure.TE_reinforcement.base
        self.hTEuniax = structure.TE_reinforcement.height_uniax
        self.hTEfoam = structure.TE_reinforcement.height_foam
        self.hLEpanel = structure.LE_panel.height
        self.hAFTpanel = structure.aft_panel.height

    def read_airfoil_coords(self, airfoil_path):
        """Read the airfoil coordinates and scale wrt the airfoil dims."""
        self.coords = np.loadtxt(airfoil_path, dtype=dtype)
        # translate the airfoil horizontally, so the pitch axis is at the origin
        self.coords['x'] = self.coords['x'] - self.pitch_axis
        # scale the airfoil to the full chord length
        self.coords['x'] = self.coords['x'] * self.chord
        self.coords['y'] = self.coords['y'] * self.chord
        try:
            # save the upper and lower airfoil surfaces
            #   index of coordinate at leading edge
            self.LE_index = np.nonzero(self.coords['y']==0.0)[0][-1]
            self.upper = self.coords[self.LE_index:]
            self.lower = self.coords[:self.LE_index+1]
        except:
            print " [Warning] Couldn't find LE index!"
            # write this to the log! ****************************************

    def plot_thickness_vs_chord(self):
        """Plots the thickness of the airfoil versus the chord length."""
        plt.figure()
        plt.axes().set_aspect('equal')
        t = []
        for i in range(len(self.lower)):
            t.append(self.upper[::-1]['y'][i] - self.lower['y'][i])
        plt.plot(self.lower['x'],t,'b+-')
        plt.show()
        return t

    def save_region_edges(self):
        """Save the x-coordinates of the region edges for each structural component."""
        # LE panel
        self.LE_panel_left = -self.pitch_axis*self.chord
        self.LE_panel_right = -self.bSC/2.0-self.bSW
        # spar caps
        self.SC_left = -self.bSC/2.0
        self.SC_right = self.bSC/2.0
        # shear webs
        self.fwd_SW_left = -self.bSC/2.0-self.bSW
        self.fwd_SW_foam_left = self.fwd_SW_left+self.bSWbiax
        self.fwd_SW_foam_right = self.fwd_SW_foam_left+self.bSWfoam
        self.fwd_SW_right = -self.bSC/2.0
        self.rear_SW_left = self.bSC/2.0
        self.rear_SW_foam_left = self.rear_SW_left+self.bSWbiax
        self.rear_SW_foam_right = self.rear_SW_foam_left+self.bSWfoam
        self.rear_SW_right = self.bSC/2.0+self.bSW
        # aft panels
        self.aft_panel_left = self.bSC/2.0+self.bSW
        self.aft_panel_right = -self.pitch_axis*self.chord+self.chord-self.bTEreinf
        # TE reinforcement
        self.TE_reinf_left = -self.pitch_axis*self.chord+self.chord-self.bTEreinf
        self.TE_reinf_right = -self.pitch_axis*self.chord+self.chord

    def plot_region_edges(self):
        """Plots color blocks to mark off the regions for each structural
        component. Each color block spans the plot from top to bottom.

        """
        # LE panel
        plt.axvspan(self.LE_panel_left, self.LE_panel_right, facecolor='magenta', edgecolor='magenta')
        # spar caps
        plt.axvspan(self.SC_left, self.SC_right, facecolor='cyan', edgecolor='cyan')
        # shear webs
        #   biax
        plt.axvspan(self.fwd_SW_left, self.fwd_SW_foam_left, facecolor='green', edgecolor='green')
        plt.axvspan(self.fwd_SW_foam_right, self.fwd_SW_right, facecolor='green', edgecolor='green')
        plt.axvspan(self.rear_SW_left, self.rear_SW_foam_left, facecolor='green', edgecolor='green')
        plt.axvspan(self.rear_SW_foam_right, self.rear_SW_right, facecolor='green', edgecolor='green')
        #   foam
        plt.axvspan(self.fwd_SW_foam_left, self.fwd_SW_foam_right, facecolor='yellow', edgecolor='yellow')
        plt.axvspan(self.rear_SW_foam_left, self.rear_SW_foam_right, facecolor='yellow', edgecolor='yellow')
        # aft panels
        plt.axvspan(self.aft_panel_left, self.aft_panel_right, facecolor='orange', edgecolor='orange')
        # TE reinforcement
        plt.axvspan(self.TE_reinf_left, self.TE_reinf_right, facecolor='pink', edgecolor='pink')

    def find_edge_coords(self, x_edge):
        """Find the airfoil coordinates at the edges of each structural component.

        Returns two coordinate pairs as tuples, one coordinate pair for the lower
        surface (x_edge,y_edge_lower), and another for the upper surface of the airfoil
        (x_edge,y_edge_upper).

        """
        # lower airfoil surface
        index_right = np.nonzero(self.lower['x']>x_edge)[0][-1]
        index_left = index_right + 1
        f = ipl.interp1d(self.lower[index_right:index_left+1][::-1]['x'],
                         self.lower[index_right:index_left+1][::-1]['y'])
        y_edge_lower = float(f(x_edge))
        # plt.plot(x_edge,y_edge_lower,'ro')
        temp = np.append(self.lower[:index_left],
                         np.array((x_edge,y_edge_lower),
                                  dtype=dtype))
        self.lower = np.append(temp, self.lower[index_left:])
        # ---------------------
        # upper airfoil surface
        index_right = np.nonzero(self.upper['x']>x_edge)[0][0]
        index_left = index_right - 1
        f = ipl.interp1d(self.upper[index_left:index_right+1]['x'],
                         self.upper[index_left:index_right+1]['y'])
        y_edge_upper = float(f(x_edge))
        # plt.plot(x_edge,y_edge_upper,'gs')
        temp = np.append(self.upper[:index_right],
                         np.array((x_edge,y_edge_upper),
                                  dtype=dtype))
        self.upper = np.append(temp, self.upper[index_right:])
        return ((x_edge,y_edge_lower),(x_edge,y_edge_upper))

    def save_edge_coords(self):
        """Use the method find_edge_coords() to save edge coordinates for all structural
        components as attributes in the Airfoil class.

        """
        self.aft_panel_right_coords = self.find_edge_coords(self.aft_panel_right)
        self.aft_panel_left_coords = self.find_edge_coords(self.aft_panel_left)
        self.SC_right_coords = self.find_edge_coords(self.SC_right)
        self.SC_left_coords = self.find_edge_coords(self.SC_left)
        self.fwd_SW_left_coords = self.find_edge_coords(self.fwd_SW_left)
        self.fwd_SW_foam_left_coords = self.find_edge_coords(self.fwd_SW_foam_left)
        self.fwd_SW_foam_right_coords = self.find_edge_coords(self.fwd_SW_foam_right)
        # fwd_SW_right_coords are the same as SC_left_coords
        self.fwd_SW_right_coords = self.SC_left_coords
        # rear_SW_left_coords are the same as SC_right_coords
        self.rear_SW_left_coords = self.SC_right_coords
        self.rear_SW_foam_left_coords = self.find_edge_coords(self.rear_SW_foam_left)
        self.rear_SW_foam_right_coords = self.find_edge_coords(self.rear_SW_foam_right)
        # rear_SW_right_coords are the same as aft_panel_left_coords
        self.rear_SW_right_coords = self.aft_panel_left_coords
        # LE_panel_right_coords are the same as fwd_SW_left_coords
        self.LE_panel_right_coords = self.fwd_SW_left_coords
        # TE_reinf_left_coords are the same as aft_panel_right_coords
        self.TE_reinf_left_coords = self.aft_panel_right_coords

    def extract_segment_along_airfoil_profile(self, left_coords, right_coords):
        """Extract a single segment along the airfoil profile, given the left and right
        edge coordinates for that segment.

        Returns two numpy arrays of coordinates, one for the lower segment and another
        for the upper segment.

        """
        # find the indices of the segment edges
        left_index_lower = np.nonzero(self.lower==np.array(left_coords[0],dtype=dtype))[0][0]
        right_index_lower = np.nonzero(self.lower==np.array(right_coords[0],dtype=dtype))[0][0]
        left_index_upper = np.nonzero(self.upper==np.array(left_coords[1],dtype=dtype))[0][0]
        right_index_upper = np.nonzero(self.upper==np.array(right_coords[1],dtype=dtype))[0][0]
        # save the segments
        lower_segment = self.lower[right_index_lower:left_index_lower+1][::-1]
        upper_segment = self.upper[left_index_upper:right_index_upper+1]
        return (lower_segment, upper_segment)

    def extract_LE_segment(self):
        """Extract the segment along the leading edge panel.

        Returns a numpy array of coordinates for the leading edge segment.

        """
        # find the indices
        lower_index = np.nonzero(self.lower==np.array(self.LE_panel_right_coords[0],dtype=dtype))[0][0]
        upper_index = np.nonzero(self.upper==np.array(self.LE_panel_right_coords[1],dtype=dtype))[0][0]
        # save the segment
        lower_segment = self.lower[lower_index:]
        upper_segment = self.upper[:upper_index+1]
        LE_segment = np.append(lower_segment, upper_segment[1:])
        return LE_segment

    def extract_fwd_SW_segments(self):
        """Extract the segments along the biax and foam parts of the forward shear web."""
        # forward biax
        (self.fwd_SW_fwdbiax_lower_segment, self.fwd_SW_fwdbiax_upper_segment) = self.extract_segment_along_airfoil_profile(self.fwd_SW_left_coords, self.fwd_SW_foam_left_coords)
        # foam core
        (self.fwd_SW_foam_lower_segment, self.fwd_SW_foam_upper_segment) = self.extract_segment_along_airfoil_profile(self.fwd_SW_foam_left_coords, self.fwd_SW_foam_right_coords)
        # rear biax
        (self.fwd_SW_rearbiax_lower_segment, self.fwd_SW_rearbiax_upper_segment) = self.extract_segment_along_airfoil_profile(self.fwd_SW_foam_right_coords, self.fwd_SW_right_coords)

    def extract_rear_SW_segments(self):
        """Extract the segments along the biax and foam parts of the rear shear web."""
        # forward biax
        (self.rear_SW_fwdbiax_lower_segment, self.rear_SW_fwdbiax_upper_segment) = self.extract_segment_along_airfoil_profile(self.rear_SW_left_coords, self.rear_SW_foam_left_coords)
        # foam core
        (self.rear_SW_foam_lower_segment, self.rear_SW_foam_upper_segment) = self.extract_segment_along_airfoil_profile(self.rear_SW_foam_left_coords, self.rear_SW_foam_right_coords)
        # rear biax
        (self.rear_SW_rearbiax_lower_segment, self.rear_SW_rearbiax_upper_segment) = self.extract_segment_along_airfoil_profile(self.rear_SW_foam_right_coords, self.rear_SW_right_coords)

    def extract_TE_segments(self):
        """Extract the segments along the trailing edge reinforcement.

        Returns five numpy arrays of coordinates:
        (1) lower main segment: The main segment for the trailing edge reinforcement,
            along the lower airfoil surface.
        (2) upper main segment: The main segment for the trailing edge reinforcement,
            along the upper airfoil surface.
        (3) lower tip segment: An additional segment near the tip of the trailing edge
            reinforcement, along the lower airfoil surface.
        (4) upper sharp segment: An additional segment near the tip of the trailing
            edge reinforcement, along the upper airfoil surface.
        (5) inner surface segment: An additional segment near the tip of the trailing
            edge reinforcement, along the interface between the upper and lower
            laminates.

        """
        # find the normal vectors on the upper and lower parts of the TE_segment
        laminate_thickness = self.hTEfoam+self.hTEuniax
        # upper sharp trailing edge segment ----------------------------------
        x0 = self.upper['x'][-2]
        y0 = self.upper['y'][-2]
        x1 = self.upper['x'][-1]
        y1 = self.upper['y'][-1]
        tang_upper = np.array([x1-x0, y1-y0])
        tang_upper = tang_upper/np.linalg.norm(tang_upper)
        tang_upper_x = tang_upper[0]
        tang_upper_y = tang_upper[1]
        norm_upper = -np.array([-tang_upper_y, tang_upper_x])
        normpt = np.array([x0, y0]) + norm_upper*laminate_thickness
        x3 = normpt[0]
        y3 = normpt[1]
        # plt.plot([x0, x3], [y0, y3], 'go-')
        x5 = x3 + tang_upper_x
        y5 = y3 + tang_upper_y
        # plt.plot([x3, x5], [y3, y5], 'go:')
        # lower sharp trailing edge segment ----------------------------------
        x2 = self.lower['x'][1]
        y2 = self.lower['y'][1]
        x7 = self.lower['x'][0]
        y7 = self.lower['y'][0]
        tang_lower = np.array([x7-x2, y7-y2])
        tang_lower = tang_lower/np.linalg.norm(tang_lower)
        tang_lower_x = tang_lower[0]
        tang_lower_y = tang_lower[1]
        norm_lower = np.array([-tang_lower_y, tang_lower_x])
        normpt = np.array([x2, y2]) + norm_lower*laminate_thickness
        x4 = normpt[0]
        y4 = normpt[1]
        # plt.plot([x2, normpt[0]], [y2, normpt[1]], 'go-')
        x6 = x4 + tang_lower_x
        y6 = y4 + tang_lower_y
        # plt.plot([x4, x6], [y4, y6], 'go:')
        # find the intersection of the TE laminate thicknesses ---------------
        m64 = (y6-y4)/(x6-x4)
        m53 = (y5-y3)/(x5-x3)
        x_int = (m64*x4 - y4 - m53*x3 + y3)/(m64 - m53)
        y_int = m64*(x_int - x4) + y4
        # plt.plot([x_int], [y_int], 'yp')
        # plot normals above and below the intersection point ----------------
        normpt = np.array([x_int, y_int]) - norm_upper*laminate_thickness
        x9 = normpt[0]
        y9 = normpt[1]
        # plt.plot([x_int, x9], [y_int, y9], 'c^:')
        normpt = np.array([x_int, y_int]) - norm_lower*laminate_thickness
        x10 = normpt[0]
        y10 = normpt[1]
        # plt.plot([x_int, x10], [y_int, y10], 'c^:')
        # find the midpoint of the TE thickness ------------------------------
        x8 = x1
        y8 = abs(y1-y7)/2.0 + y7
        # plot the inner surface of the TE reinforcement ---------------------
        # plt.plot([x3,x_int,x8], [y3,y_int,y8], 'mp-')
        # plt.plot([x4,x_int,x8], [y4,y_int,y8], 'mp-')
        # --------------------------------------------------------------------
        # append new coords from normals above and below the intersection pt
        temp = np.append(self.lower[0], np.array((x10,y10),dtype=dtype))
        self.lower = np.append(temp, self.lower[1:])
        temp = np.append(np.array((x9,y9),dtype=dtype), self.upper[-1])
        self.upper = np.append(self.upper[:-1], temp)
        # find the indices at the left boundary of the TE reinforcement
        lower_index = np.nonzero(self.lower==np.array(self.TE_reinf_left_coords[0],dtype=dtype))[0][0]
        upper_index = np.nonzero(self.upper==np.array(self.TE_reinf_left_coords[1],dtype=dtype))[0][0]
        # save the main segments
        lower_main_segment = self.lower[1:lower_index+1]
        upper_main_segment = self.upper[upper_index:-1]
        # save the additional segments near the sharp trailing edge
        lower_sharp_segment = self.lower[:2]
        upper_sharp_segment = self.upper[-2:]
        inner_surf_segment = np.array([(x_int,y_int),(x8,y8)], dtype=dtype)
        return (lower_main_segment, upper_main_segment, lower_sharp_segment, upper_sharp_segment, inner_surf_segment)

    def extract_all_segments_along_airfoil_profile(self):
        """Extract all segments along the airfoil profile.

        Segments for each structural component are extracted in this order:
        (1) spar caps
        (2) aft panels
        (3) forward shear web
        (4) rear shear web
        (5) leading edge panel
        (6) trailing edge reinforcement

        """
        (self.SC_lower_segment, self.SC_upper_segment) = self.extract_segment_along_airfoil_profile(self.SC_left_coords, self.SC_right_coords)
        (self.aft_panel_lower_segment, self.aft_panel_upper_segment) = self.extract_segment_along_airfoil_profile(self.aft_panel_left_coords, self.aft_panel_right_coords)
        (self.fwd_SW_lower_segment, self.fwd_SW_upper_segment) = self.extract_segment_along_airfoil_profile(self.fwd_SW_left_coords, self.fwd_SW_right_coords)
        (self.rear_SW_lower_segment, self.rear_SW_upper_segment) = self.extract_segment_along_airfoil_profile(self.rear_SW_left_coords, self.rear_SW_right_coords)
        # self.extract_fwd_SW_segments()
        # self.extract_rear_SW_segments()
        self.LE_segment = self.extract_LE_segment()
        (self.TE_lower_main_segment, self.TE_upper_main_segment, self.TE_lower_sharp_segment, self.TE_upper_sharp_segment, self.TE_inner_surf_segment) = self.extract_TE_segments()

    def plot_all_segments(self, plot_title=''):
        """Plot all the segments along the airfoil profile with different symbols."""
        plt.figure()
        plt.title(plot_title)
        plt.axes().set_aspect('equal')
        self.save_region_edges()
        self.plot_region_edges()
        # self.save_edge_coords()
        # self.extract_all_segments_along_airfoil_profile()
        # plt.plot(self.LE_segment['x'],self.LE_segment['y'],'ko-')
        # plt.plot(self.fwd_SW_upper_segment['x'],self.fwd_SW_upper_segment['y'],'bs-')
        # plt.plot(self.fwd_SW_lower_segment['x'],self.fwd_SW_lower_segment['y'],'rs-')
        # plt.plot(self.SC_upper_segment['x'],self.SC_upper_segment['y'],'m^-')
        # plt.plot(self.SC_lower_segment['x'],self.SC_lower_segment['y'],'g^-')
        # plt.plot(self.rear_SW_upper_segment['x'],self.rear_SW_upper_segment['y'],'bo-')
        # plt.plot(self.rear_SW_lower_segment['x'],self.rear_SW_lower_segment['y'],'ro-')
        # plt.plot(self.aft_panel_upper_segment['x'],self.aft_panel_upper_segment['y'],'ms-')
        # plt.plot(self.aft_panel_lower_segment['x'],self.aft_panel_lower_segment['y'],'gs-')
        # plt.plot(self.TE_upper_main_segment['x'],self.TE_upper_main_segment['y'],'b^-')
        # plt.plot(self.TE_lower_main_segment['x'],self.TE_lower_main_segment['y'],'r^-')
        # plt.plot(self.TE_upper_sharp_segment['x'],self.TE_upper_sharp_segment['y'],'mo-')
        # plt.plot(self.TE_lower_sharp_segment['x'],self.TE_lower_sharp_segment['y'],'go-')
        # plt.plot(self.TE_inner_surf_segment['x'],self.TE_inner_surf_segment['y'],'c^--')
        plt.show()
        