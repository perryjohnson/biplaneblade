"""A module for organizing geometrical data for a blade station.

Author: Perry Roth-Johnson
Last updated: July 24, 2013

"""


import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as ipl
from math import isnan


class Coordinates:
    """Define the coordinates for a blade station."""
    def __init__(self, x1, x2, x3):
        self.x1 = x1   # units [m]
        self.x2 = x2   # units [m]
        self.x3 = x3   # units [m]
    def __str__(self):
        return """x1:  {0:7.3f} (meters)
x2:  {1:7.3f} (meters)
x3:  {2:7.3f} (meters)""".format(self.x1, self.x2, self.x3)


class Airfoil:
    """Define the airfoil and chord properties (external dimensions)."""
    def __init__(self, name, filename, pitch_axis, chord, twist):
        self.name = name
        self.filename = filename
        self.path = ''  # initialize empty string
                        # will be assigned later by Blade.copy_airfoil_coords()
        self.pitch_axis = pitch_axis  # units [-]  (chord fraction)
        self.chord = chord            # units [m]
        self.twist = twist            # units [deg]
        # self.curvatures  # implement later for biplane blade
    def __str__(self):
        return """Airfoil:     {0}
Filename:    {1}
Pitch axis:  {2:6.3f} (chord fraction)
Chord:       {3:6.3f} (meters)
Twist:       {4:6.3f} (degrees)""".format(self.name, self.filename,
            self.pitch_axis, self.chord, self.twist)


class Part:
    """Define the dimensions of a structural component."""
    def __init__(self, base, height):
        self.base = base
        self.height = height
    def __str__(self):
        return """base:    {0} (meters)
height:  {1} (meters)""".format(self.base, self.height)
    def exists(self):
        """Checks if a structural component exists at this station."""
        if isnan(self.base) and isnan(self.height):
            return False
        else:
            return True


class ShearWeb(Part):
    """Define the biax (skin) and foam (core) dimensions of a shear web."""
    def __init__(self, base_biax, base_foam, height=np.nan):
        Part.__init__(self, base=(2.0*base_biax+base_foam), height=height)
        self.base_biax = base_biax
        self.base_foam = base_foam
    def __str__(self):
        return """base:    {0:6.3f} (meters)
|-> base_biax:  {1:6.3f} (meters)
|-> base_foam:  {2:6.3f} (meters)
height:  {3} (meters)""".format(self.base, self.base_biax,
    self.base_foam, self.height)


class TE_Reinforcement(Part):
    """Define uniax and foam dimensions of a trailing edge reinforcement."""
    def __init__(self, base, height_uniax, height_foam):
        Part.__init__(self, base, height=(height_uniax+height_foam))
        self.height_uniax = height_uniax
        self.height_foam = height_foam
    def __str__(self):
        return """base:    {0:6.3f} (meters)
height:  {1:6.3f} (meters)
|-> height_uniax:  {2:6.3f} (meters)
|-> height_foam:   {3:6.3f} (meters)""".format(self.base, self.height,
    self.height_uniax, self.height_foam)


class Structure:
    """Define the laminate schedule (internal dimensions)."""
    def __init__(self, h_RB, b_SC, h_SC, b_SW1_biax, b_SW1_foam, 
                 b_SW2_biax, b_SW2_foam, b_SW3_biax, b_SW3_foam,
                 b_TE_reinf, h_TE_reinf_uniax, h_TE_reinf_foam, h_LE_panel,
                 h_aft_panel):
        self.root_buildup = Part(np.nan, h_RB)
        self.spar_cap = Part(b_SC, h_SC)
        self.shear_web_1 = ShearWeb(b_SW1_biax, b_SW1_foam)
        self.shear_web_2 = ShearWeb(b_SW2_biax, b_SW2_foam)
        self.shear_web_3 = ShearWeb(b_SW3_biax, b_SW3_foam)
        self.TE_reinforcement = TE_Reinforcement(b_TE_reinf, h_TE_reinf_uniax, 
                                                 h_TE_reinf_foam)
        self.LE_panel = Part(np.nan, h_LE_panel)
        self.aft_panel = Part(np.nan, h_aft_panel)        

    def print_summary(self):
        """Print all the internal dimensions of this structure."""
        print "--- ROOT BUILDUP ---"
        print self.root_buildup
        print "--- SPAR CAP ---"
        print self.spar_cap
        print "--- SHEAR WEB 1 ---"
        print self.shear_web
        print "--- SHEAR WEB 2 ---"
        print self.shear_web
        print "--- SHEAR WEB 3 ---"
        print self.shear_web
        print "---TE REINFORCEMENT ---"
        print self.TE_reinforcement
        print "--- LE PANEL ---"
        print self.LE_panel
        print "--- AFT PANEL ---"
        print self.aft_panel

    def get_summary(self):
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
        s += "--- AFT PANEL ---\n"
        s += str(self.aft_panel) + '\n'
        return s


class Station:
    """Define a station for a wind turbine blade.

    This module also contains methods to split the airfoil curve into separate
    segments for each structural component: spar caps, shear webs, etc.

    Usage
    -----
    import pandas as pd
    import station as stn
    df = pd.read_csv('Sandia_blade.csv', index_col=0)
    s5 = stn.Station(df.ix[5])  # import station 5

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
        .structure
            .root_buildup
                .base=np.nan
                .height : float, the root buildup height (meters)
            .spar_cap
                .base : float, the spar cap base (meters)
                .height : float, the spar cap height (meters)
            .shear_web_1
                .base : float, the shear web #1 total base (meters)
                .base_biax : float, the shear web #1 base for biax (meters)
                .base_foam : float, the shear web #1 base for foam (meters)
                .height=np.nan
            .shear_web_2
                .base : float, the shear web #2 total base (meters)
                .base_biax : float, the shear web #2 base for biax (meters)
                .base_foam : float, the shear web #2 base for foam (meters)
                .height=np.nan
            .shear_web_3
                .base : float, the shear web #3 total base (meters)
                .base_biax : float, the shear web #3 base for biax (meters)
                .base_foam : float, the shear web #3 base for foam (meters)
                .height=np.nan
            .TE_reinforcement
                .base : float, the trailing edge reinforcement base (meters)
                .height_uniax : float, the TE reinf height for uniax (meters)
                .height_foam : float, the TE reinf height for foam (meters)
                .height : float, the TE reinf total height (meters)
            .LE_panel
                .base=np.nan
                .height : the leading edge panel height (meters)
            .aft_panel
                .base=np.nan
                .height : the aft panel height (meters)

        Methods
        -------
        .structure
            .<Part>
                .exists() : bool, checks if a Part exists at this station
            .print_summary() : stdout, print all internal dimensions
            .get_summary() : str, returns all the internal dimensions

        Usage
        -----    
        Station(b._df.ix[5], 'sandia_blade')
        # this creates station #5 of the Sandia blade
        # _df is a pandas DataFrame containing properties of all blade stations
        # _df.ix[5] gets the Series object for station #5 from DataFrame _df
        # Note: Stations are usually not created directly. New Stations are
        # usually created by the Blade class.

        """
        Station.number_of_stations += 1
        self.station_num = Station.number_of_stations
        self.station_path = os.path.join(blade_path, 'stn{0:02d}'.format(self.station_num))
        try:
            os.mkdir(self.station_path)
        except WindowsError:
            print "[WindowsError] The station path '{0}' already exists!".format(os.path.split(self.station_path)[-1])
        self.logf = open(Station.logfile_name, "a")
        self.logf.write("............(Created blade station #{0})............\n".format(self.station_num))
        print " Created blade station #{0}".format(self.station_num)
        self.coords = Coordinates(stn_series['x1'], 
                                  stn_series['x2'], 
                                  stn_series['x3'])
        self.logf.write("****** COORDINATES ******\n")
        self.logf.write(str(self.coords) + '\n')
        self.airfoil = Airfoil(stn_series['airfoil'], 
                               stn_series['airfoil']+'.txt', 
                               stn_series['pitch axis'], 
                               stn_series['chord'], 
                               stn_series['twist'])
        self.logf.write("****** AIRFOIL AND CHORD PROPERTIES ******\n")
        self.logf.write(str(self.airfoil) + '\n')
        self.structure = Structure(stn_series['root buildup height'],
                                   stn_series['spar cap base'],
                                   stn_series['spar cap height'],
                                   stn_series['shear web 1 base biax'],
                                   stn_series['shear web 1 base foam'],
                                   stn_series['shear web 2 base biax'],
                                   stn_series['shear web 2 base foam'],
                                   stn_series['shear web 3 base biax'],
                                   stn_series['shear web 3 base foam'],
                                   stn_series['TE reinf base'],
                                   stn_series['TE reinf height uniax'],
                                   stn_series['TE reinf height foam'],
                                   stn_series['LE panel height'],
                                   stn_series['aft panel height'])
        self.logf.write("****** LAMINATE SCHEDULE ******\n")
        self.logf.write(self.structure.get_summary())
        self.logf.write('\n')
        self.logf.flush()
        self.logf.close()

    def __del__(self):
        Station.number_of_stations = Station.number_of_stations - 1
        print " Station deleted, and now Station.number_of_stations = {0}".format(Station.number_of_stations)

    def read_airfoil_coords(self, comment_char='#'):
        """Read the airfoil coordinates and scale wrt the airfoil dims.

        Note
        ----
        If the trailing edge of the airfoil is a thin feature, it must have a
        finite thickness. Check the airfoil coordinates files to make sure this
        is true.

        For "airfoils" with thick trailing edges (transition, ellipse, or
        cylinder), it is okay to have a trailing edge with zero thickness.

        """
        af = self.airfoil
        af.coords = np.loadtxt(af.path, dtype=[('x', 'f8'), ('y', 'f8')], 
            comments=comment_char)
        # translate the airfoil horizontally, so pitch axis is at origin
        af.coords['x'] = af.coords['x'] - af.pitch_axis
        # scale the airfoil to the specified chord length
        af.coords['x'] = af.coords['x'] * af.chord
        af.coords['y'] = af.coords['y'] * af.chord

    def plot_airfoil_coords(self, show_flag=False, savefig_flag=True):
        """Plot the airfoil coordinates of this station."""
        af = self.airfoil
        plt.figure()
        plt.title("Station #{0}, {1}".format(self.station_num,
                                             self.airfoil.name))
        plt.axes().set_aspect('equal')
        plt.plot(af.coords['x'], af.coords['y'])
        if show_flag:
            plt.show()
        if savefig_flag:
            fname = os.path.join(self.station_path, 'stn{0:02d}.png'.format(self.station_num))
            plt.savefig(fname)

    # def split_airfoil_at_LE_and_TE(self):

    # note: use <Part>.exists() method to decide whether or not to split the airfoil curve for a particular Part

    # note: clean up docs to use consistent terminology -- use "Part", not "component"

    # note: keep implementing methods from airfoil_utils.py into this Station class!!!!