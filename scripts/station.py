import os

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


class ShearWeb(Part):
    """Define the biax (skin) and foam (core) dimensions of a shear web."""
    def __init__(self, base_biax, base_foam, height=None):
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
    def __init__(self, b_RB, h_RB, b_SC, h_SC, b_SW_biax, b_SW_foam,
                 b_TE_reinf, h_TE_reinf_uniax, h_TE_reinf_foam, h_LE_panel,
                 h_aft_panel):
        self.root_buildup = Part(b_RB, h_RB)
        self.spar_cap = Part(b_SC, h_SC)
        self.shear_web = ShearWeb(b_SW_biax, b_SW_foam)
        self.TE_reinforcement = TE_Reinforcement(b_TE_reinf, h_TE_reinf_uniax, 
                                                 h_TE_reinf_foam)
        self.LE_panel = Part(None, h_LE_panel)
        self.aft_panel = Part(None, h_aft_panel)
    def print_summary(self):
        print "--- ROOT BUILDUP ---"
        print self.root_buildup
        print "--- SPAR CAP ---"
        print self.spar_cap
        print "--- SHEAR WEB ---"
        print self.shear_web
        print "---TE REINFORCEMENT ---"
        print self.TE_reinforcement
        print "--- LE PANEL ---"
        print self.LE_panel
        print "--- AFT PANEL ---"
        print self.aft_panel
    def get_summary(self):
        s = ''
        s += "--- ROOT BUILDUP ---\n"
        s += str(self.root_buildup) + '\n'
        s += "--- SPAR CAP ---\n"
        s += str(self.spar_cap) + '\n'
        s += "--- SHEAR WEB ---\n"
        s += str(self.shear_web) + '\n'
        s += "---TE REINFORCEMENT ---\n"
        s += str(self.TE_reinforcement) + '\n'
        s += "--- LE PANEL ---\n"
        s += str(self.LE_panel) + '\n'
        s += "--- AFT PANEL ---\n"
        s += str(self.aft_panel) + '\n'
        return s


class Station:
    """Define a station for a wind turbine blade.

    Usage:
    import pandas as pd
    import station as stn
    df = pd.read_csv('Sandia_blade.csv', index_col=0)
    s5 = stn.Station(df.ix[5])  # import station 5

    """
    logfile_name = 'station.log'
    number_of_stations = 0
    def __init__(self, stn_series, blade_path):
        """Create a new blade station.

        Note: Stations are usually not created directly. New Stations are
        usually created by the Blade class.

        Arguments
        ---------
        stn_series: a pandas Series object containing this station's properties
        blade_path: the local target directory for storing blade data

        Usage
        -----
        Station(b._df.ix[5], 'sandia_blade')
        # this creates station #5 of the Sandia blade
        # _df is a pandas DataFrame containing properties of all blade stations
        # _df.ix[5] gets the Series object for station #5 from DataFrame _df

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
        self.structure = Structure(stn_series['root buildup base'],
                                   stn_series['root buildup height'],
                                   stn_series['spar cap base'],
                                   stn_series['spar cap height'],
                                   stn_series['shear web base biax'],
                                   stn_series['shear web base foam'],
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
