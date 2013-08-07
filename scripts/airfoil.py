"""A module for organizing airfoil data for a blade station.

Author: Perry Roth-Johnson
Last updated: August 7, 2013

"""


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
        self.path = ''  # initialize empty string
                        # will be assigned later by Blade.copy_airfoil_coords()
        self.chord = chord            # units [m]
    def __str__(self):
        return """Monoplane Airfoil ---
Airfoil:     {0}
Filename:    {1}
Chord:       {2:6.4f} (meters)
Pitch axis:  {3:6.4f} (chord fraction)
Twist:       {4:6.4f} (degrees)""".format(self.name, self.filename,
            self.chord, self.pitch_axis, self.twist)


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
        .upper_name
        .upper_filename
        .upper_path
        .upper_chord
        .gap_to_chord_ratio
        .gap_fraction
        .stagger_to_chord_ratio
        .gap
        .stagger
        .total_chord

    """
    def __init__(self, name, name_L, filename_L, chord_L, name_U, filename_U,
        chord_U, pitch_axis, twist, gap_to_chord_ratio, gap_fraction,
        stagger_to_chord_ratio):
        _Airfoil.__init__(self, name, pitch_axis, twist)
        self.lower_name = name_L
        self.lower_filename = filename_L
        self.lower_path = ''
        self.lower_chord = chord_L
        self.upper_name = name_U
        self.upper_filename = filename_U
        self.upper_path = ''
        self.upper_chord = chord_U
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
