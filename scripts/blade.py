import os
import shutil
import datetime
import pandas as pd
import station as stn
reload(stn)
import matplotlib.pyplot as plt


class Blade:
    """Define a wind turbine blade.

    Usage:
    import blade as bld
    b = bld.Blade('Sandia blade SNL100-00', 'sandia_blade')

    # access the chord length of blade station 2
    stn2 = b.list_of_stations[1]
    stn2.airfoil.chord

    """
    logfile_name = 'blade.log'
    def __init__(self, name, blade_path, defn_filename='blade_definition.csv', airfoils_path='airfoils'):
        """Create a new wind turbine blade.

        Arguments
        ---------
        name: the name of this blade
        blade_path: the local target directory for storing blade data

        Keyword Arguments
        -----------------
        defn_filename: CSV file for blade schedule of airfoils and laminates
        airfoils_path: the local directory that contains airfoil coordinates

        Usage
        -----
        Blade('Sandia blade SNL100-00', 'sandia_blade')
        Blade('Sandia blade SNL100-00', 'sandia_blade',
              defn_filename='blade_definition.csv', airfoils_path='airfoils')

        """
        # initialize the number_of_stations counter to zero
        stn.Station.number_of_stations = 0
        self.name = name
        self.logf = open(Blade.logfile_name, "a")
        self.logf.write("[{0}] Created blade: {1}\n".format(datetime.datetime.now(), self.name))
        if not os.path.exists(blade_path):
            msg = "[Error] The blade path '{0}' does not exist!\n...Exiting...".format(blade_path)
            self.logf.write("[{0}] {1}".format(datetime.datetime.now(), msg))
            print msg
        else:
            self.blade_path = os.path.join(os.getcwd(), blade_path)
            self.logf.write("[{0}] Found blade path: {1}\n".format(datetime.datetime.now(), self.blade_path))
            self.defn_filename = os.path.join(self.blade_path, defn_filename)
            self.logf.write("[{0}] Found blade definition file: {1}\n".format(datetime.datetime.now(), self.defn_filename))
            import_success = self.import_blade_definition()
            if import_success:
                self.create_all_stations()
                self.logf.write("[{0}] Created all blade stations\n".format(datetime.datetime.now()))
                self.airfoils_path = os.path.join(self.blade_path, airfoils_path)
                self.logf.write("[{0}] Found airfoils path: {1}\n".format(datetime.datetime.now(), self.airfoils_path))
        self.logf.flush()
        self.logf.close()
    def __del__(self):
        self.logf = open(Blade.logfile_name, "a")
        for station in self.list_of_stations:
            if os.path.exists(station.station_path):
                # delete the station path, even if it has contents
                shutil.rmtree(station.station_path)
        print "Deleted blade '{0}' and all its station paths.".format(self.name)
        self.logf.write("[{0}] Deleted blade '{1}' and all its station paths.\n".format(datetime.datetime.now(), self.name))
        self.logf.flush()
        self.logf.close()
    def import_blade_definition(self):
        import_result = 0
        defn_fileext = os.path.splitext(self.defn_filename)[-1]
        if defn_fileext != '.csv' and defn_fileext != '.CSV':
            # check that the blade definition file is a CSV file
            print "[Error] The blade definition file '{0}' is not a CSV file!".format(os.path.split(self.defn_filename)[-1])
            print "        The blade defintion will not be imported."
            print "        ...Exiting..."
        else:
            # import the blade definition file into a pandas DataFrame
            self._df = pd.read_csv(self.defn_filename, index_col=0)
            self.number_of_stations = len(self._df)
            import_result = 1
        return import_result
    def create_station(self, station_num):
        return stn.Station(self._df.ix[station_num], self.blade_path)
    def create_all_stations(self):
        self.list_of_stations = []
        for station in range(1, self.number_of_stations+1):
            self.list_of_stations.append(self.create_station(station))
    def copy_airfoil_coords(self, station):
        try:
            shutil.copy(os.path.join(self.airfoils_path, station.airfoil.filename), station.station_path)
            print " Copied station #{0} airfoil: {1}".format(station.station_num, station.airfoil.name)
        except IOError:
            print "[IOError] The airfoil file '{0}' for station {1} does not exist!\n          Check '{2}' for errors.".format(station.airfoil.filename, station.station_num, self.defn_filename)
    def copy_all_airfoil_coords(self):
        for station in self.list_of_stations:
            self.copy_airfoil_coords(station)
    def plot_chord_schedule(self):
        """Plot the chord vs. span."""
        plt.figure()
        plt.axes().set_aspect('equal')
        plt.plot(self._df['x1'],self._df['chord'],'bo-')
        plt.xlabel('span, x1 [m]')
        plt.ylabel('chord [m]')
        plt.show()
    def plot_twist_schedule(self):
        """Plot the twist vs. span."""
        plt.figure()
        plt.axes().set_aspect('equal')
        plt.plot(self._df['x1'],self._df['twist'],'bo-')
        plt.xlabel('span, x1 [m]')
        plt.ylabel('twist [deg]')
        plt.show()