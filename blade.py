import pandas as pd
import station as stn
import matplotlib.pyplot as plt


class Blade:
    """Define a wind turbine blade.

    Usage:
    import blade as bld
    b = bld.Blade('Sandia blade', 'Sandia_blade.csv')

    # access the chord length of blade station 2
    stn2 = b.list_of_stations[1]
    stn2.airfoil.chord

    """
    def __init__(self, name, defn_filename):
        """Initialize a new wind turbine blade."""
        self.name = name
        self.defn_filename = defn_filename
        # import the blade definition from a CSV file
        self._df = pd.read_csv(self.defn_filename, index_col=0)
        self.number_of_stations = len(self._df)
        # create all the blade stations
        self.list_of_stations = []
        for station in self._df.index:
            self.list_of_stations.append(stn.Station(self._df.ix[station]))
            print "............(Added blade station #{0})............".format(
                station)
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
    def create_station(self, station_num):
        return stn.Station(self._df.ix[station_num])