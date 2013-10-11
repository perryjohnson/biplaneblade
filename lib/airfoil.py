"""A module for organizing airfoil data for a blade station.

Author: Perry Roth-Johnson
Last updated: August 7, 2013

"""


import numpy as np
import transformation as tf
import scipy.interpolate as ipl
from shapely.geometry import Polygon


class _Airfoil:
    """Define an airfoil (external dimensions).

    The _Airfoil base class is not intended for use.
    Use MonoplaneAirfoil or BiplaneAirfoil instead.

    """
    def __init__(self, name, pitch_axis, twist):
        self.name = name
        self.pitch_axis = pitch_axis  # units [-]  (chord fraction)
        self.twist = twist            # units [deg]


class MonoplaneAirfoil(_Airfoil):
    """Define a monoplane airfoil (external dimensions)."""
    def __init__(self, name, filename, chord, pitch_axis, twist):
        _Airfoil.__init__(self, name, pitch_axis, twist)
        self.filename = filename
        self.path = None        # assigned later by Blade.copy_airfoil_coords()
        self.chord = chord      # units [m]
        self.coords = None      # assigned later by read_coords()
        self.LE_index = None    # assigned later by split_at_LE_and_TE()
        self.suction = None     # assigned later by split_at_LE_and_TE()
        self.pressure = None    # assigned later by split_at_LE_and_TE()
        self.polygon = None     # assigned later by create_polygon()

    def __str__(self):
        return """Monoplane Airfoil ---
Airfoil:     {0}
Filename:    {1}
Chord:       {2:6.4f} (meters)
Pitch axis:  {3:6.4f} (chord fraction)
Twist:       {4:6.4f} (degrees)""".format(self.name, self.filename,
            self.chord, self.pitch_axis, self.twist)

    def read_coords(self, comment_char='#'):
        """Read the airfoil coordinates into memory from a file in station_path.

        Creates a new attribute for this airfoil: <Airfoil>.coords,
        which is a numpy array of airfoil coordinates.

        Note
        ----
        If the trailing edge of the airfoil is a thin feature, it must have a
        finite thickness. Check the airfoil coordinates files to make sure this
        is true.

        For "airfoils" with thick trailing edges (transition, ellipse, or
        cylinder), it is okay to have a trailing edge with zero thickness.

        """
        try:
            self.coords = np.loadtxt(self.path, dtype=[('x', 'f8'), ('y', 'f8')], 
                comments=comment_char)
        except IOError:
            raise IOError("Airfoil file does not exist yet!\n  Run <Blade>.copy_all_airfoil_coords() first.")

    def scale_coords(self, scale_factor):
        """Scale the airfoil coordinates by a scale factor.

        Must run <Airfoil>.read_coords() first.

        Parameters
        ----------
        scale_factor : float, amount to scale the airfoil coords by

        Usage
        -----
        # scale all airfoils in blade 'm' to the specified chord length
        for station in m.list_of_stations:
            station.airfoil.scale_coords(scale_factor=station.airfoil.chord)

        """
        self.coords['x'] = self.coords['x'] * scale_factor
        self.coords['y'] = self.coords['y'] * scale_factor

    def translate_coords_chordwise(self, x, direc):
        """Translate the airfoil along the chordwise direction.

        Parameters
        ----------
        x : float, absolute amount to translate the airfoil coords by
        direc : str (either 'fwd' or 'aft'), desired direction of translation

        Usage
        -----
        # translate all airfoils in blade 'm', so pitch axis is at origin
        for station in m.list_of_stations:
            station.airfoil.translate_coords_chordwise(
                x=station.airfoil.pitch_axis*station.airfoil.chord,
                direc='fwd')

        """
        if direc != 'fwd' and direc != 'aft':
            raise ValueError("keyword 'direc' must be either 'fwd' or 'aft'.")
        if direc == 'fwd':
            self.coords['x'] = self.coords['x'] - abs(x)
        elif direc == 'aft':
            self.coords['x'] = self.coords['x'] + abs(x)

    def scale_and_translate_coords(self):
        """Scale and translate the monoplane airfoil coords.

        First, scale up the airfoil by the chord length.
        Second, shift the airfoil forward by (pitch axis fraction)*(chord)

        """
        self.scale_coords(scale_factor=self.chord)
        self.translate_coords_chordwise(x=self.pitch_axis*self.chord,
            direc='fwd')

    def rotate_coords(self):
        """Rotate the airfoil coordinates wrt the local twist angle.

        Must run <Airfoil>.read_coords() and <Airfoil>.scale_coords() first.

        """
        for point in self.coords:
            (x,y) = point
            (point['x'], point['y']) = tf.rotate_coord_pair(x, y, self.twist)

    def split_at_LE_and_TE(self):
        """Split the monoplane airfoil curve into suction and pressure segments.

        Must be run after <Airfoil>.scale_and_translate_coords().

        """
        try:
            temp_list = np.nonzero(self.coords['y']==0.0)[0]
        except AttributeError:
            raise AttributeError("{0} coordinates for station #{1} haven't been read!\n  You need to first read in the coordinates with <Station>.airfoil.read_coords().".format(self.name, self.station_num))
        else:
            # drop zeros from the list (which correspond to the TE, not the LE)
            temp_list = temp_list[np.nonzero(temp_list)[0]]
            # grab the first item in the list
            # (the last item will also correspond to the TE, not the LE)
            try:
                self.LE_index = temp_list[0]
            except IndexError:
                raise IndexError("Leading edge index was not found!")
            self.suction = self.coords[self.LE_index:]
            self.pressure = self.coords[:self.LE_index+1]

    def create_polygon(self):
        """Convert the numpy array of coordinates into a polygon object.

        This allows us to use offset and clipping methods in the Shapely module

        """
        l = []
        for point in self.coords:
            x = float(point['x'])
            y = float(point['y'])
            l.append((x,y))
        self.polygon = Polygon(l)

    def plot_coords(self, axes, split_flag=False):
        """Plot the monoplane airfoil coordinates of this station."""
        if split_flag:
            try:
                axes.plot(self.suction['x'], self.suction['y'], 'bo-', 
                    label='suction surface')
                axes.plot(self.pressure['x'], self.pressure['y'], 'rs-', 
                    label='pressure surface')
            except AttributeError:
                raise AttributeError("Suction and pressure surface {0} coordinates\n  haven't been read!\n  You need to first run <Airfoil>.split_at_LE_and_TE().".format(self.name))
        else:
            try:
                axes.plot(self.coords['x'], self.coords['y'])
            except AttributeError:
                raise AttributeError("{0} coordinates haven't been read!\n  You need to first read in the coordinates with <Station>.airfoil.read_coords().".format(self.name))

    def find_part_edge_coords(self, x_edge):
        """Find the airfoil coordinates at the edges of each structural part.

        Returns two coordinate pairs as tuples, one coordinate pair for the
        pressure surface (x_edge, y_edge_pressure), and another for the suction
        surface of the airfoil (x_edge, y_edge_suction).

        Must run <Station>.airfoil.split_at_LE_and_TE() first.

        """
        # pressure airfoil surface --------------------------------------------
        try:
            index_right = np.nonzero(self.pressure['x']>x_edge)[0][-1]
        except AttributeError:
            raise AttributeError("Upper and pressure surface {0} coordinates\n  for station #{1} haven't been read!\n  You need to first run <Station>.airfoil.split_at_LE_and_TE().".format(self.name, self.station_num))
        index_left = index_right + 1
        f = ipl.interp1d(self.pressure[index_right:index_left+1][::-1]['x'],
                         self.pressure[index_right:index_left+1][::-1]['y'])
        y_edge_pressure = float(f(x_edge))
        # plt.plot(x_edge,y_edge_pressure,'ro')
        temp = np.append(self.pressure[:index_left],
                         np.array((x_edge,y_edge_pressure),
                                  dtype=[('x', 'f8'), ('y', 'f8')]))
        self.pressure = np.append(temp, self.pressure[index_left:])
        # suction airfoil surface ---------------------------------------------
        index_right = np.nonzero(self.suction['x']>x_edge)[0][0]
        index_left = index_right - 1
        f = ipl.interp1d(self.suction[index_left:index_right+1]['x'],
                         self.suction[index_left:index_right+1]['y'])
        y_edge_suction = float(f(x_edge))
        # plt.plot(x_edge,y_edge_suction,'gs')
        temp = np.append(self.suction[:index_right],
                         np.array((x_edge,y_edge_suction),
                                  dtype=[('x', 'f8'), ('y', 'f8')]))
        self.suction = np.append(temp, self.suction[index_right:])
        return ((x_edge,y_edge_pressure),(x_edge,y_edge_suction))


class BiplaneAirfoil(_Airfoil):
    """Define a biplane airfoil (external dimensions).

    BiplaneAirfoil
        .name
        .pitch_axis
        .twist
        .lower_name
        .lower_filename
        .lower_path
        .lower_chord
        .lower_SW_ref_pt_fraction
        .upper_name
        .upper_filename
        .upper_path
        .upper_chord
        .upper_SW_ref_pt_fraction
        .gap_to_chord_ratio
        .gap_fraction
        .stagger_to_chord_ratio
        .gap
        .stagger
        .total_chord

    """
    def __init__(self, name, name_L, filename_L, chord_L, SW_ref_pt_L, name_U,
        filename_U, chord_U, SW_ref_pt_U, pitch_axis, twist,
        gap_to_chord_ratio, gap_fraction, stagger_to_chord_ratio):
        _Airfoil.__init__(self, name, pitch_axis, twist)
        self.lower_name = name_L
        self.lower_filename = filename_L
        self.lower_path = None
        self.lower_chord = chord_L
        self.lower_SW_ref_pt_fraction = SW_ref_pt_L
        self.upper_name = name_U
        self.upper_filename = filename_U
        self.upper_path = None
        self.upper_chord = chord_U
        self.upper_SW_ref_pt_fraction = SW_ref_pt_U
        self.gap_to_chord_ratio = gap_to_chord_ratio
        self.gap_fraction = gap_fraction
        self.stagger_to_chord_ratio = stagger_to_chord_ratio
        # calculate gap, stagger, and total chord
        #   note: gap and stagger are normalized by the LOWER chord length
        self.gap = self.gap_to_chord_ratio * chord_L
        self.stagger = self.stagger_to_chord_ratio * chord_L
        self.total_chord = self.stagger + chord_L

    def __str__(self):
        return """Biplane Airfoil ---
Name:  {0}
Lower element
  Airfoil:     {1}
  Filename:    {2}
  Chord:       {3:6.4f} (meters)
Upper element
  Airfoil:     {4}
  Filename:    {5}
  Chord:       {6:6.4f} (meters)
Gap-to-chord ratio:      {7:6.4f}
Gap fraction:            {8:6.4f}
Gap:                     {9:6.4f} (normalized by lower chord)
Stagger-to-chord ratio:  {10:6.4f}
Stagger:                 {11:6.4f} (normalized by lower chord)
Total chord:             {12:6.4f}
Pitch axis:              {13:6.4f} (chord fraction)
Twist:                   {14:6.4f} (degrees)""".format(self.name, 
            self.lower_name, self.lower_filename, self.lower_chord, 
            self.upper_name, self.upper_filename, self.upper_chord, 
            self.gap_to_chord_ratio, self.gap_fraction, self.gap, 
            self.stagger_to_chord_ratio, self.stagger, self.total_chord, 
            self.pitch_axis, self.twist)

    def read_coords(self, comment_char='#'):
        """Read the airfoil coordinates into memory from a file in station_path.

        Creates a new attribute for this station: <Station>.airfoil.coords,
        which is a numpy array of airfoil coordinates.

        Note
        ----
        If the trailing edge of the airfoil is a thin feature, it must have a
        finite thickness. Check the airfoil coordinates files to make sure this
        is true.

        For "airfoils" with thick trailing edges (transition, ellipse, or
        cylinder), it is okay to have a trailing edge with zero thickness.

        """
        # lower airfoil
        try:
            self.lower_coords = np.loadtxt(self.lower_path, 
                dtype=[('x', 'f8'), ('y', 'f8')], comments=comment_char)
        except IOError:
            raise IOError("Lower airfoil file does not exist yet!\n  Run <Blade>.copy_all_airfoil_coords() first.")
        # upper airfoil
        try:
            self.upper_coords = np.loadtxt(self.upper_path, 
                dtype=[('x', 'f8'), ('y', 'f8')], comments=comment_char)
        except IOError:
            raise IOError("Upper airfoil file does not exist yet!\n  Run <Blade>.copy_all_airfoil_coords() first.")

    def scale_coords(self, upper_factor, lower_factor):
        """Scale the biplane airfoil coordinates by scale factors.

        Must run <Airfoil>.read_coords() first.

        Parameters
        ----------
        upper_factor : float, amount to scale the upper airfoil coords by
        lower_factor : float, amount to scale the lower airfoil coords by

        Usage
        -----
        # scale all airfoils in blade 'm' to the specified chord lengths
        for station in m.list_of_stations:
            station.airfoil.scale_coords(
                upper_factor=station.airfoil.upper_chord,
                lower_factor=station.airfoil.lower_chord)

        """
        # upper airfoil
        self.upper_coords['x'] = self.upper_coords['x'] * upper_factor
        self.upper_coords['y'] = self.upper_coords['y'] * upper_factor
        # lower airfoil
        self.lower_coords['x'] = self.lower_coords['x'] * lower_factor
        self.lower_coords['y'] = self.lower_coords['y'] * lower_factor

    def translate_coords_flapwise(self, upper_y, upper_direc, lower_y,
        lower_direc):
        """Translate the biplane airfoil along the flapwise direction.

        Parameters
        ----------
        upper_y : float, absolute amount to translate the upper airfoil coords
        upper_direc : str ('up' or 'down') or None, desired direction of
            translation
        lower_y : float, absolute amount to translate the lower airfoil coords
        lower_direc : str ('up' or 'down') or None, desired direction of
            translation

        Usage
        -----
        # translate all airfoils in blade 'm', so upper and lower airfoils
        # are separated by a gap, and pitch axis is at gap fraction
        for station in m.list_of_stations:
            station.airfoil.translate_coords_flapwise(
                upper_y=station.airfoil.gap_fraction*station.airfoil.gap,
                upper_direc='up',
                lower_y=(1.0-station.airfoil.gap_fraction)*station.airfoil.gap,
                lower_direc='down')

        """
        if upper_direc != 'up' and upper_direc != 'down' and upper_direc != None:
            raise ValueError("keyword 'upper_direc' must be 'up', 'down', or None.")
        if lower_direc != 'up' and lower_direc != 'down' and lower_direc != None:
            raise ValueError("keyword 'lower_direc' must be 'up', 'down', or None.")
        # upper airfoil
        if upper_direc == 'up':
            self.upper_coords['y'] = self.upper_coords['y'] + abs(upper_y)
        elif upper_direc == 'down':
            self.upper_coords['y'] = self.upper_coords['y'] - abs(upper_y)
        # lower airfoil
        if lower_direc == 'up':
            self.lower_coords['y'] = self.lower_coords['y'] + abs(lower_y)
        elif lower_direc == 'down':
            self.lower_coords['y'] = self.lower_coords['y'] - abs(lower_y)

    def translate_coords_chordwise(self, upper_x, upper_direc, lower_x,
        lower_direc):
        """Translate the biplane airfoil along the chordwise direction.

        Parameters
        ----------
        upper_x : float, absolute amount to translate the upper airfoil coords
        upper_direc : str (either 'fwd' or 'aft'), desired direction of
            translation
        lower_x : float, absolute amount to translate the lower airfoil coords
        lower_direc : str (either 'fwd' or 'aft'), desired direction of
            translation

        Usage
        -----
        # translate all airfoils in blade 'm', so pitch axis is at origin
        for station in m.list_of_stations:
            # shift lower airfoil aft by stagger
            station.airfoil.translate_coords_chordwise(
                upper_x=0.0,
                upper_direc=None,
                lower_x=station.airfoil.stagger,
                lower_direc='aft')
            # shift BOTH airfoils forward by (pitch axis frac)*(total chord)
            station.airfoil.translate_coords_chordwise(
                upper_x=station.airfoil.pitch_axis*station.airfoil.total_chord,
                upper_direc='fwd',
                lower_x=station.airfoil.pitch_axis*station.airfoil.total_chord,
                lower_direc='fwd')

        """
        if upper_direc != 'fwd' and upper_direc != 'aft' and upper_direc != None:
            raise ValueError("keyword 'upper_direc' must be 'fwd', 'aft', or None.")
        if lower_direc != 'fwd' and lower_direc != 'aft' and lower_direc != None:
            raise ValueError("keyword 'lower_direc' must be 'fwd', 'aft', or None.")
        # upper airfoil
        if upper_direc == 'fwd':
            self.upper_coords['x'] = self.upper_coords['x'] - abs(upper_x)
        elif upper_direc == 'aft':
            self.upper_coords['x'] = self.upper_coords['x'] + abs(upper_x)
        # lower airfoil
        if lower_direc == 'fwd':
            self.lower_coords['x'] = self.lower_coords['x'] - abs(lower_x)
        elif lower_direc == 'aft':
            self.lower_coords['x'] = self.lower_coords['x'] + abs(lower_x)

    def scale_and_translate_coords(self):
        """Scale and translate the biplane airfoil coords.

        First, scale up the upper and lower airfoils by their chord lengths.
        Second, translate the upper and lower airfoils apart by a gap, with the
            pitch axis at the gap fraction.
        Third, translate the lower airfoil aft by the stagger
        Fourth, translate the upper and lower airfoils forward by
            (pitch axis fraction)*(total chord)

        """
        self.scale_coords(upper_factor=self.upper_chord,
                          lower_factor=self.lower_chord)
        # shift upper and lower airfoils apart by a gap,
        # with pitch axis at gap fraction
        self.translate_coords_flapwise(
            upper_y=self.gap_fraction*self.gap,
            upper_direc='up',
            lower_y=(1.0-self.gap_fraction)*self.gap,
            lower_direc='down')
        # shift lower airfoil aft by stagger
        self.translate_coords_chordwise(
            upper_x=0.0,
            upper_direc=None,
            lower_x=self.stagger,
            lower_direc='aft')
        # shift BOTH airfoils forward by (pitch axis frac)*(total chord)
        self.translate_coords_chordwise(
            upper_x=self.pitch_axis*self.total_chord,
            upper_direc='fwd',
            lower_x=self.pitch_axis*self.total_chord,
            lower_direc='fwd')

    def rotate_coords(self):
        """Rotate the airfoil coordinates wrt the local twist angle.

        Must run <Airfoil>.read_coords() and <Airfoil>.scale_coords() first.

        """
        # rotate the lower airfoil
        for point in self.lower_coords:
            (x,y) = point
            (point['x'], point['y']) = tf.rotate_coord_pair(x, y, self.twist)
        # rotate the upper airfoil
        for point in self.upper_coords:
            (x,y) = point
            (point['x'], point['y']) = tf.rotate_coord_pair(x, y, self.twist)

    def split_at_LE_and_TE(self):
        """Split each biplane airfoil curve into suction and pressure segments.

        Must be run after <Airfoil>.scale_and_translate_coords().

        """
        # lower airfoil
        try:
            temp_list = np.nonzero(
                self.lower_coords['y']==-(1.0-self.gap_fraction)*self.gap)[0]
        except AttributeError:
            raise AttributeError("{0} lower coordinates for station #{1} haven't been read!\n  You need to first read in the lower coordinates with <Station>.airfoil.read_coords().".format(self.lower_name, self.station_num))
        else:
            # drop zeros from the list (which correspond to the TE, not the LE)
            temp_list = temp_list[np.nonzero(temp_list)[0]]
            # grab the first item in the list
            # (the last item will also correspond to the TE, not the LE)
            try:
                self.lower_LE_index = temp_list[0]
            except IndexError:
                raise IndexError("Leading edge index was not found!")
            self.lower_suction = self.lower_coords[self.lower_LE_index:]
            self.lower_pressure = self.lower_coords[:self.lower_LE_index+1]
        # upper airfoil
        try:
            temp_list = np.nonzero(
                self.upper_coords['y']==self.gap_fraction*self.gap)[0]
        except AttributeError:
            raise AttributeError("{0} upper coordinates for station #{1} haven't been read!\n  You need to first read in the upper coordinates with <Station>.airfoil.read_coords().".format(self.upper_name, self.station_num))
        else:
            # drop zeros from the list (which correspond to the TE, not the LE)
            temp_list = temp_list[np.nonzero(temp_list)[0]]
            # grab the first item in the list
            # (the last item will also correspond to the TE, not the LE)
            try:
                self.upper_LE_index = temp_list[0]
            except IndexError:
                raise IndexError("Leading edge index was not found!")
            self.upper_suction = self.upper_coords[self.upper_LE_index:]
            self.upper_pressure = self.upper_coords[:self.upper_LE_index+1]

    def plot_coords(self, axes, split_flag=False):
        """Plot the biplane airfoil coordinates of this station."""
        if split_flag:
            # lower airfoil
            try:
                axes.plot(self.lower_suction['x'], self.lower_suction['y'], 'bo-', 
                    label='lower, suction surface')
                axes.plot(self.lower_pressure['x'], self.lower_pressure['y'], 'rs-', 
                    label='lower, pressure surface')
            except AttributeError:
                raise AttributeError("Lower suction and pressure surface {0} coordinates\n  haven't been read!\n  You need to first run <Airfoil>.split_at_LE_and_TE().".format(self.lower_name))
            # upper airfoil
            try:
                axes.plot(self.upper_suction['x'], self.upper_suction['y'], 'bo-', 
                    label='upper, suction surface')
                axes.plot(self.upper_pressure['x'], self.upper_pressure['y'], 'rs-', 
                    label='upper, pressure surface')
            except AttributeError:
                raise AttributeError("Upper suction and pressure surface {0} coordinates\n  haven't been read!\n  You need to first run <Airfoil>.split_at_LE_and_TE().".format(self.upper_name))
        else:
            try:
                axes.plot(self.lower_coords['x'], self.lower_coords['y'])
            except AttributeError:
                raise AttributeError("{0} lower coordinates for station #{1} haven't been read!\n  You need to first read in the coordinates with <Station>.airfoil.read_coords().".format(self.lower_name, self.station_num))
            try:
                axes.plot(self.upper_coords['x'], self.upper_coords['y'])
            except AttributeError:
                raise AttributeError("{0} upper coordinates for station #{1} haven't been read!\n  You need to first read in the coordinates with <Station>.airfoil.read_coords().".format(self.lower_name, self.station_num))

    def find_part_edge_coords(self, x_edge, airfoil):
        """Find the airfoil coordinates at the edges of each structural part.

        Returns two coordinate pairs as tuples, one coordinate pair for the
        pressure surface (x_edge, y_edge_pressure), and another for the suction
        surface of the airfoil (x_edge, y_edge_suction).

        Must run <Station>.airfoil.split_at_LE_and_TE() first.

        Parameters
        ----------
        x_edge : float, the x-coordinate of the part edge
        airfoil : str, ('lower' or 'upper') desired airfoil to find part edge
            (x,y) coordinates on

        """
        if airfoil != 'lower' and airfoil != 'upper':
            raise ValueError("keyword 'airfoil' must be 'lower' or 'upper'.")
        # select either the lower or upper airfoil
        if airfoil == 'lower':
            pressure = self.lower_pressure
            name = self.lower_name
            suction = self.lower_suction
        elif airfoil == 'upper':
            pressure = self.upper_pressure
            name = self.upper_name
            suction = self.upper_suction
        # pressure airfoil surface ----------------------------------------
        try:
            index_right = np.nonzero(pressure['x']>x_edge)[0][-1]
        except AttributeError:
            raise AttributeError("Suction and pressure surface {0} coordinates\n  for station #{1} haven't been read!\n  You need to first run <Station>.airfoil.split_at_LE_and_TE().".format(name, self.station_num))
        index_left = index_right + 1
        f = ipl.interp1d(pressure[index_right:index_left+1][::-1]['x'],
                         pressure[index_right:index_left+1][::-1]['y'])
        y_edge_pressure = float(f(x_edge))
        # plt.plot(x_edge,y_edge_pressure,'ro')
        temp = np.append(pressure[:index_left],
                         np.array((x_edge,y_edge_pressure),
                                  dtype=[('x', 'f8'), ('y', 'f8')]))
        pressure = np.append(temp, pressure[index_left:])
        # suction airfoil surface -----------------------------------------
        index_right = np.nonzero(suction['x']>x_edge)[0][0]
        index_left = index_right - 1
        f = ipl.interp1d(suction[index_left:index_right+1]['x'],
                         suction[index_left:index_right+1]['y'])
        y_edge_suction = float(f(x_edge))
        # plt.plot(x_edge,y_edge_suction,'gs')
        temp = np.append(suction[:index_right],
                         np.array((x_edge,y_edge_suction),
                                  dtype=[('x', 'f8'), ('y', 'f8')]))
        suction = np.append(temp, suction[index_right:])
        return ((x_edge,y_edge_pressure),(x_edge,y_edge_suction))
