"""A module for organizing geometrical data for a blade definition.

Author: Perry Roth-Johnson
Last updated: August 6, 2013

"""


import os
import shutil
import datetime
import numpy as np
import pandas as pd
import station as stn
reload(stn)
import transformation as tf
reload(tf)
import matplotlib.pyplot as plt
from mayavi import mlab


class _Blade:
    """Define a wind turbine blade.

    The _Blade base class is not intended for use.
    Use MonoplaneBlade or BiplaneBlade instead

    Usage:
    import blade as bld
    b = bld.MonoplaneBlade('Sandia blade SNL100-00', 'sandia_blade')

    # access the chord length of blade station 2
    stn2 = b.list_of_stations[1]
    stn2.airfoil.chord

    """
    logfile_name = 'blade.log'
    def __init__(self, name, blade_path, defn_filename='blade_definition.csv',
        airfoils_path='airfoils'):
        """Create a new wind turbine blade.

        Parameters
        ----------
        name : string, the name of this blade
        blade_path : string, the local target directory for storing blade data

        defn_filename : str (for CSV file), the blade definition filename
        airfoils_path : str, local directory that contains airfoil coordinates

        Attributes
        ----------
        .airfoils_path : str, local directory that contains airfoil coords
        .blade_path : str, the local target directory for storing blade data
        .defn_filename : str, (for CSV file), the blade definition filename
        .list_of_stations : list, contains all the Stations of this blade
        .name : str, the name of this blade
        .number_of_stations : int, the total number of stations in this blade
        .logfile_name : 'blade.log', the log file for this blade

        Methods
        -------
        .copy_airfoil_coords(station) : copy airfoil coordinates from 
            airfoils_path into this station_path
        .copy_all_airfoil_coords() : copy all airfoil coordinates from
            airfoils_path into each station_path
        .create_all_stations() : create all stations for this blade
        .create_plot() : create a plot for this blade
        .create_station(station_num) : create a new station for this blade
        .get_LE_coords() : list, returns (x,y,z) coords for the blade LE
        .get_SW_cross_section_coords(sw_num) : list, returns (x,y,z) coords for
            shear web 1, 2, or 3
        .get_TE_coords() : list, returns (x,y,z) coords for the blade TE
        .import_blade_definition() : import the blade defn from a CSV file
        .plot_LE(lw) : plots the leading edge from root to tip
        .plot_TE(lw) : plots the trailing edge from root to tip
        .plot_all_SW_cross_sections(lw) : plots all shear web cross-sections
        .plot_all_SW_spanwise_lines(lw) : plots spanwise lines for all SWs
        .plot_all_airfoils(lw) : plot all the airfoils in the blade
        .plot_blade() : plots a wirerame representation of the blade
        .plot_chord_schedule() : plot the chord vs. span
        .plot_pitch_axis(lw) : plots the pitch axis from root to tip
        .plot_twist_schedule() : plot the twist vs. span
        .show_plot() : pick a nice view and show the plot

        Usage
        -----
        _Blade('Sandia blade SNL100-00', 'sandia_blade')
        _Blade('Sandia blade SNL100-00', 'sandia_blade',
              defn_filename='blade_definition.csv', airfoils_path='airfoils')

        """
        self.name = name
        self.logf = open(_Blade.logfile_name, "a")
        self.logf.write("[{0}] Created blade: {1}\n".format(datetime.datetime.now(), self.name))
        if not os.path.exists(blade_path):
            raise ValueError("The blade path '{0}' does not exist!\n  Check the 'blade_path' passed to the Blade class.".format(blade_path))
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
        """Delete a wind turbine blade.

        Also deletes all the station paths inside the blade path.

        """
        self.logf = open(_Blade.logfile_name, "a")
        for station in self.list_of_stations:
            if os.path.exists(station.station_path):
                # delete the station path, even if it has contents
                shutil.rmtree(station.station_path)
        print "Deleted blade '{0}' and all its station paths.".format(self.name)
        self.logf.write("[{0}] Deleted blade '{1}' and all its station paths.\n".format(datetime.datetime.now(), self.name))
        self.logf.flush()
        self.logf.close()

    def import_blade_definition(self):
        """Import the blade definition from a CSV file.

        Returns 1 if the import was successful, 0 otherwise.

        """
        import_result = 0
        defn_fileext = os.path.splitext(self.defn_filename)[-1]
        if defn_fileext != '.csv' and defn_fileext != '.CSV':
            # check that the blade definition file is a CSV file
            raise ValueError("blade definition file '{0}' must be of type *.csv".format(os.path.split(self.defn_filename)[-1]))
        else:
            # import the blade definition file into a pandas DataFrame
            self._df = pd.read_csv(self.defn_filename, index_col=0)
            self.number_of_stations = len(self._df)
            import_result = 1
        return import_result

    def create_all_stations(self):
        """Create all stations for this blade."""
        if stn._Station.number_of_stations != 0:
            stn._Station.number_of_stations = 0  # initialize to zero
            print " [Warning] The number of stations has been reset to zero."
            self.logf = open(_Blade.logfile_name, 'a')
            self.logf.write("[{0}] [Warning] The number of stations (_Station.number_of_stations) was reset to zero.\n".format(datetime.datetime.now()))
        self.list_of_stations = []
        for station in range(1, self.number_of_stations+1):
            self.list_of_stations.append(self.create_station(station))

    def copy_all_airfoil_coords(self):
        """Copy all airfoil coordinates from airfoils_path into each station_path."""
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

    def create_plot(self, fig_width=1000, fig_height=800):
        """Create a plot for this blade.

        Returns a handle to the figure.

        Parameters
        ----------
        fig_width : int, figure width (pixels)
        fig_height : int, figure height (pixels)

        Several settings are applied by default ---------
        Title : [blade name]
        Figure width : 1000 px
        Figure height : 800 px

        """
        mlab.figure(figure=self.name, size=(fig_width,fig_height))

    def show_plot(self, azimuth=-45.0, elevation=54.74, distance=110.0,
                  focalpoint=[60.0,0.85,0.0], print_view=False):
        """Pick a nice view and show the plot.

        Parameters
        ----------
        azimuth : float, optional. The azimuthal angle (in degrees, 0-360), 
            i.e. the angle subtended by the position vector on a sphere 
            projected on to the x-y plane with the x-axis.
        elevation : float, optional. The zenith angle (in degrees, 0-180), 
            i.e. the angle subtended by the position vector and the z-axis.
        distance : float or 'auto', optional. A positive floating point number 
            representing the distance from the focal point to place the camera.
            New in Mayavi 3.4.0: if 'auto' is passed, the distance is computed 
            to have a best fit of objects in the frame.
        focalpoint : array_like or 'auto', optional. An array of 3 floating 
            point numbers representing the focal point of the camera. New in 
            Mayavi 3.4.0: if 'auto' is passed, the focal point is positioned at
            the center of all objects in the scene.
        print_view : boolean, print/don't print the view parameters

        """
        mlab.view(azimuth, elevation, distance, focalpoint)
        if print_view:
            print "[**Mayavi mlab view parameters**]"
            print mlab.view()
        mlab.show()

    def plot_pitch_axis(self, lw):
        """Plots the pitch axis from root to tip."""
        root = self.list_of_stations[0].coords.x1
        tip = self.list_of_stations[-1].coords.x1
        mlab.plot3d([root,tip],[0,0],[0,0], color=(1,0,0), tube_radius=lw)

    def get_SW_cross_section_coords(self, sw_num, twist_flag=True):
        """Returns a list of (x,y,z) coordinates for shear web #1, #2, or #3.

        Parameters
        ----------
        sw_num = int, (1, 2, or 3) selects the desired shear web

        """
        if sw_num is not 1 and sw_num is not 2 and sw_num is not 3:
            raise ValueError("sw_num must be either 1, 2, or 3 (type: int)")
        x = []  # spanwise coordinate
        y = []  # chordwise coordinate
        z = []  # flapwise coordinate
        for station in self.list_of_stations:
            if sw_num is 1:
                if station.structure.shear_web_1.exists():
                    for point in station.structure.shear_web_1.cs_coords:
                        x.append(station.coords.x1)
                        y_temp = point[0]
                        z_temp = point[1]
                        if twist_flag:
                            # rotate the SW cross-section coords wrt twist angle
                            (y_new, z_new) = tf.rotate_coord_pair(y_temp,
                                z_temp, station.airfoil.twist)
                        else:
                            # don't rotate the SW cross-section coords
                            (y_new, z_new) = (y_temp, z_temp)
                        # append the new SW cross-section coords
                        y.append(y_new)
                        z.append(z_new)
            elif sw_num is 2:
                if station.structure.shear_web_2.exists():
                    for point in station.structure.shear_web_2.cs_coords:
                        x.append(station.coords.x1)
                        y_temp = point[0]
                        z_temp = point[1]
                        if twist_flag:
                            # rotate the SW cross-section coords wrt twist angle
                            (y_new, z_new) = tf.rotate_coord_pair(y_temp,
                                z_temp, station.airfoil.twist)
                        else:
                            # don't rotate the SW cross-section coords
                            (y_new, z_new) = (y_temp, z_temp)
                        # append the new SW cross-section coords
                        y.append(y_new)
                        z.append(z_new)
            elif sw_num is 3:
                if station.structure.shear_web_3.exists():
                    for point in station.structure.shear_web_3.cs_coords:
                        x.append(station.coords.x1)
                        y_temp = point[0]
                        z_temp = point[1]
                        if twist_flag:
                            # rotate the SW cross-section coords wrt twist angle
                            (y_new, z_new) = tf.rotate_coord_pair(y_temp,
                                z_temp, station.airfoil.twist)
                        else:
                            # don't rotate the SW cross-section coords
                            (y_new, z_new) = (y_temp, z_temp)
                        # append the new SW cross-section coords
                        y.append(y_new)
                        z.append(z_new)
        return (x,y,z)

    def plot_all_SW_cross_sections(self, lw, twist_flag=True):
        """Plots all shear web cross-sections from root to tip."""
        for sw in [1,2,3]:
            (x,y,z) = self.get_SW_cross_section_coords(sw_num=sw,
                                                       twist_flag=twist_flag)
            for i in range(0,len(x),4):
                (X,Y,Z) = (x[i:i+4], y[i:i+4], z[i:i+4])
                # add the first coord pair again, to form a closed loop
                X.append(x[i])
                Y.append(y[i])
                Z.append(z[i])
                mlab.plot3d(X,Y,Z, color=(0,0,1), tube_radius=lw)

    def plot_all_SW_spanwise_lines(self, lw, twist_flag=True):
        """Plots spanwise lines for all shear webs from root to tip."""
        for sw in [1,2,3]:
            (x,y,z) = self.get_SW_cross_section_coords(sw_num=sw, 
                                                       twist_flag=twist_flag)
            (x1, y1, z1) = ([], [], [])
            (x2, y2, z2) = ([], [], [])
            (x3, y3, z3) = ([], [], [])
            (x4, y4, z4) = ([], [], [])
            for i in range(0,len(x),4):
                x1.append(x[i])
                y1.append(y[i])
                z1.append(z[i])
                x2.append(x[i+1])
                y2.append(y[i+1])
                z2.append(z[i+1])
                x3.append(x[i+2])
                y3.append(y[i+2])
                z3.append(z[i+2])
                x4.append(x[i+3])
                y4.append(y[i+3])
                z4.append(z[i+3])
            mlab.plot3d(x1,y1,z1, color=(0,0,1), tube_radius=lw)
            mlab.plot3d(x2,y2,z2, color=(0,0,1), tube_radius=lw)
            mlab.plot3d(x3,y3,z3, color=(0,0,1), tube_radius=lw)
            mlab.plot3d(x4,y4,z4, color=(0,0,1), tube_radius=lw)

    def plot_blade(self, line_width=0.08, airfoils=True, pitch_axis=False,
        LE=True, TE=True, twist=True, SW=True):
        """Plots a wireframe representation of the blade, with Mayavi mlab.

        Parameters
        ----------
        line_width : float, line width for plotting
        airfoils : bool, plot/don't plot the airfoils at each blade station
        pitch_axis : bool, plot/don't plot the pitch axis from root to tip
        LE : bool, plot/don't plot the leading edge from root to tip
        TE : bool, plot/don't plot the trailing edge from root to tip
        twist : bool, plot/don't plot the blade with twist from root to tip
        SW : bool, plot/don't plot all shear webs along the span

        """
        self.create_plot()
        if airfoils:
            self.plot_all_airfoils(lw=line_width, twist_flag=twist)
        if pitch_axis:
            self.plot_pitch_axis(lw=line_width)
        if LE:
            self.plot_LE(lw=line_width, twist_flag=twist)
        if TE:
            self.plot_TE(lw=line_width, twist_flag=twist)
        if SW:
            self.plot_all_SW_cross_sections(lw=line_width, twist_flag=twist)
            self.plot_all_SW_spanwise_lines(lw=line_width, twist_flag=twist)
        self.show_plot()


class MonoplaneBlade(_Blade):
    """Define a monoplane (conventional) wind turbine blade."""
    def create_station(self, station_num):
        """Create a new station for this blade."""
        return stn.MonoplaneStation(self._df.ix[station_num], self.blade_path)

    def copy_airfoil_coords(self, station):
        """Copy airfoil coordinates from airfoils_path into this station_path."""
        try:
            shutil.copy(os.path.join(self.airfoils_path, station.airfoil.filename), station.station_path)
        except IOError:
            raise IOError("The airfoil file '{0}' for station {1} does not exist!\n  Check '{2}' for errors.".format(station.airfoil.filename, station.station_num, self.defn_filename))
        else:
            print " Copied station #{0} airfoil: {1}".format(station.station_num, station.airfoil.name)
            station.airfoil.path = os.path.join(station.station_path, station.airfoil.filename)
            print " ... Assigned station.airfoil.path!"
            self.logf = open(_Blade.logfile_name, 'a')
            self.logf.write("[{0}] Assigned station.airfoil.path to station #{1}: {2}\n".format(datetime.datetime.now(), station.station_num, station.airfoil.path))
            self.logf.flush()
            self.logf.close()

    def plot_all_airfoils(self, lw, twist_flag=True):
        """Plot all the airfoils in the blade with Mayavi's mlab.

        You must import the blade and its airfoil coordinates first.

        Usage
        -----
        b = bl.MonoplaneBlade('Sandia blade SNL100-00', 'sandia_blade')
        b.copy_all_airfoil_coords()
        for station in b.list_of_stations:
            station.read_airfoil_coords()
        b.plot_all_airfoils()

        """
        for station in self.list_of_stations:
            if twist_flag:
                station.rotate_airfoil_coords()  # apply twist angle to airfoil
            # assemble the airfoil coordinates for mlab
            try:
                y = station.airfoil.coords['x']  # chordwise coordinate
            except AttributeError:
                raise AttributeError("Airfoil coordinates haven't been read yet!\n Run <Station>.read_airfoil_coords() first.")
            z = station.airfoil.coords['y']  # flapwise coordinate
            l = len(y)
            x = np.ones((l,))*station.coords.x1  # spanwise coordinate
            # plot the airfoil on the screen
            mlab.plot3d(x,y,z, tube_radius=lw)

    def get_LE_coords(self, twist_flag=True):
        """Returns a list of (x,y,z) coordinates for the blade leading edge."""
        x = []  # spanwise coordinate
        y = []  # chordwise coordinate
        z = []  # flapwise coordinate
        for station in self.list_of_stations:
            x.append(station.coords.x1)
            # grab the unrotated LE coordinates
            y_temp = -(station.airfoil.chord * station.airfoil.pitch_axis)
            z_temp = 0.0
            if twist_flag:
                # rotate the LE coordinates wrt the twist angle
                (y_new, z_new) = tf.rotate_coord_pair(y_temp, z_temp,
                                                      station.airfoil.twist)
            else:
                # don't rotate the LE coordinates
                (y_new, z_new) = (y_temp, z_temp)
            # append the new LE coordinates
            y.append(y_new)
            z.append(z_new)
        return (x,y,z)

    def plot_LE(self, lw, twist_flag=True):
        """Plots the leading edge from root to tip."""
        (x,y,z) = self.get_LE_coords(twist_flag=twist_flag)
        mlab.plot3d(x,y,z, tube_radius=lw)

    def get_TE_coords(self, twist_flag=True):
        """Returns a list of (x,y,z) coordinates for the blade trailing edge."""
        x = []  # spanwise coordinate
        y = []  # chordwise coordinate
        z = []  # flapwise coordinate
        for station in self.list_of_stations:
            x.append(station.coords.x1)
            # grab the unrotated TE coordinates
            y_temp = station.airfoil.chord * (1.0-station.airfoil.pitch_axis)
            z_temp = 0.0
            if twist_flag:
                # rotate the TE coordinates wrt the twist angle
                (y_new, z_new) = tf.rotate_coord_pair(y_temp, z_temp,
                                                      station.airfoil.twist)
            else:
                # don't rotate the TE coordinates
                (y_new, z_new) = (y_temp, z_temp)
            # append the new TE coordinates
            y.append(y_new)
            z.append(z_new)
        return (x,y,z)

    def plot_TE(self, lw, twist_flag=True):
        """Plots the trailing edge from root to tip."""
        (x,y,z) = self.get_TE_coords(twist_flag=twist_flag)
        mlab.plot3d(x,y,z, tube_radius=lw)


class BiplaneBlade(_Blade):
    """Define a biplane wind turbine blade."""
    def create_station(self, station_num):
        """Create a new station for this blade."""
        if self._df.ix[station_num]['type'] == 'monoplane':
            this_stn = stn.MonoplaneStation(self._df.ix[station_num], self.blade_path)
        elif self._df.ix[station_num]['type'] == 'biplane':
            this_stn = stn.BiplaneStation(self._df.ix[station_num], self.blade_path)
        else:
            raise ValueError("Values in the 'type' column of {0} must be either 'monoplane' or 'biplane'.".format(self.defn_filename))
        return this_stn

    def copy_airfoil_coords(self, station):
        """Copy airfoil coordinates from airfoils_path into this station_path."""
        if station.type == 'monoplane':
            try:
                shutil.copy(os.path.join(self.airfoils_path, station.airfoil.filename), station.station_path)
            except IOError:
                raise IOError("The airfoil file '{0}' for station {1} does not exist!\n  Check '{2}' for errors.".format(station.airfoil.filename, station.station_num, self.defn_filename))
            else:
                print " Copied station #{0} airfoil: {1}".format(station.station_num, station.airfoil.name)
                station.airfoil.path = os.path.join(station.station_path, station.airfoil.filename)
                print " ... Assigned station.airfoil.path!"
                self.logf = open(_Blade.logfile_name, 'a')
                self.logf.write("[{0}] Assigned station.airfoil.path to station #{1}: {2}\n".format(datetime.datetime.now(), station.station_num, station.airfoil.path))
                self.logf.flush()
                self.logf.close()
        elif station.type == 'biplane':
            # lower airfoil
            try:
                shutil.copy(os.path.join(self.airfoils_path, station.airfoil.lower_filename), station.station_path)
            except IOError:
                raise IOError("The lower airfoil file '{0}' for station {1} does not exist!\n  Check '{2}' for errors.".format(station.airfoil.lower_filename, station.station_num, self.defn_filename))
            else:
                print " Copied station #{0} lower airfoil: {1}".format(station.station_num, station.airfoil.lower_name)
                station.airfoil.lower_path = os.path.join(station.station_path, station.airfoil.lower_filename)
                print " ... Assigned station.airfoil.lower_path!"
                self.logf = open(_Blade.logfile_name, 'a')
                self.logf.write("[{0}] Assigned station.airfoil.lower_path to station #{1}: {2}\n".format(datetime.datetime.now(), station.station_num, station.airfoil.lower_path))
                self.logf.flush()
                self.logf.close()
            # upper airfoil
            try:
                shutil.copy(os.path.join(self.airfoils_path, station.airfoil.upper_filename), station.station_path)
            except IOError:
                raise IOError("The upper airfoil file '{0}' for station {1} does not exist!\n  Check '{2}' for errors.".format(station.airfoil.upper_filename, station.station_num, self.defn_filename))
            else:
                print " Copied station #{0} upper airfoil: {1}".format(station.station_num, station.airfoil.upper_name)
                station.airfoil.upper_path = os.path.join(station.station_path, station.airfoil.upper_filename)
                print " ... Assigned station.airfoil.upper_path!"
                self.logf = open(_Blade.logfile_name, 'a')
                self.logf.write("[{0}] Assigned station.airfoil.upper_path to station #{1}: {2}\n".format(datetime.datetime.now(), station.station_num, station.airfoil.upper_path))
                self.logf.flush()
                self.logf.close()

    def plot_all_airfoils(self, lw, twist_flag=True):
        """Plot all the airfoils in the biplane blade with Mayavi's mlab.

        You must import the blade and its airfoil coordinates first.

        Usage
        -----
        a = bl.BiplaneBlade('biplane blade', 'biplane_blade')
        a.copy_all_airfoil_coords()
        for station in a.list_of_stations:
            station.read_airfoil_coords()
        a.plot_all_airfoils()

        """
        for station in self.list_of_stations:
            if station.type == 'monoplane':
                if twist_flag:
                    station.rotate_airfoil_coords()  # apply twist angle to airfoil
                # assemble the airfoil coordinates for mlab
                try:
                    y = station.airfoil.coords['x']  # chordwise coordinate
                except AttributeError:
                    raise AttributeError("Airfoil coordinates haven't been read yet!\n Run <Station>.read_airfoil_coords() first.")
                z = station.airfoil.coords['y']  # flapwise coordinate
                l = len(y)
                x = np.ones((l,))*station.coords.x1  # spanwise coordinate
                # plot the airfoil on the screen
                mlab.plot3d(x,y,z, tube_radius=lw)
            elif station.type == 'biplane':
                if twist_flag:
                    station.rotate_airfoil_coords()
                # assemble lower airfoil coordinates for mlab -----------------
                try:
                    y = station.airfoil.lower_coords['x']  # chordwise coordinate
                except AttributeError:
                    raise AttributeError("Lower airfoil coordinates haven't been read yet!\n  Run <Station>.read_airfoil_coords() first.")
                z = station.airfoil.lower_coords['y']  # flapwise coordinate
                l = len(y)
                x = np.ones((l,))*station.coords.x1  # spanwise coordinate
                # plot the lower airfoil on the screen
                mlab.plot3d(x,y,z, tube_radius=lw)
                # assemble upper airfoil coordinates for mlab -----------------
                try:
                    y = station.airfoil.upper_coords['x']  # chordwise coordinate
                except AttributeError:
                    raise AttributeError("Upper airfoil coordinates haven't been read yet!\n  Run <Station>.read_airfoil_coords() first.")
                z = station.airfoil.upper_coords['y']  # flapwise coordinate
                l = len(y)
                # (reuse same spanwise coordinates): x
                # plot the lower airfoil on the screen
                mlab.plot3d(x,y,z, tube_radius=lw)

    def get_LE_coords(self, twist_flag=True):
        """Returns a list of (x,y,z) coordinates for the blade leading edge."""
        xL = []  # spanwise coordinate (monoplane and lower biplane airfoils)
        yL = []  # chordwise coordinate (monoplane and lower biplane airfoils)
        zL = []  # flapwise coordinate (monoplane and lower biplane airfoils)
        xU = []  # spanwise coordinate (upper biplane airfoils)
        yU = []  # chordwise coordinate (upper biplane airfoils)
        zU = []  # flapwise coordinate (upper biplane airfoils)
        inboard_trans = False # flag to find inboard mono-biplane transition
        outboard_trans = False # flag to find outboard bi-monoplane transition
        for station in self.list_of_stations:
            if station.type == 'monoplane':
                xL.append(station.coords.x1)
                # grab the unrotated LE coordinates
                yL_temp = -(station.airfoil.chord * station.airfoil.pitch_axis)
                zL_temp = 0.0
                if twist_flag:
                    # rotate the LE coordinates wrt the twist angle
                    (yL_new, zL_new) = tf.rotate_coord_pair(
                        yL_temp, zL_temp, station.airfoil.twist)
                else:
                    # don't rotate the LE coordinates
                    (yL_new, zL_new) = (yL_temp, zL_temp)
                # append the new LE coordinates
                yL.append(yL_new)
                zL.append(zL_new)
                if inboard_trans is True and outboard_trans is False:
                    outboard_trans = True
                    # append the last LOWER coordinates that were added
                    xU.append(xL[-1])
                    yU.append(yL[-1])
                    zU.append(zL[-1])
            elif station.type == 'biplane':
                if inboard_trans is False:
                    inboard_trans = True
                    # append the last LOWER coordinates that were added
                    xU.append(xL[-1])
                    yU.append(yL[-1])
                    zU.append(zL[-1])
                xL.append(station.coords.x1)
                xU.append(station.coords.x1)
                # grab the unrotated LE coordinates
                yU_temp = -(station.airfoil.total_chord * station.airfoil.pitch_axis)
                yL_temp = yU_temp + station.airfoil.stagger
                zU_temp = station.airfoil.gap_fraction * station.airfoil.gap
                zL_temp = -(1.0-station.airfoil.gap_fraction) * station.airfoil.gap
                if twist_flag:
                    # rotate the LE coordinates wrt the twist angle
                    (yL_new, zL_new) = tf.rotate_coord_pair(
                        yL_temp, zL_temp, station.airfoil.twist)
                    (yU_new, zU_new) = tf.rotate_coord_pair(
                        yU_temp, zU_temp, station.airfoil.twist)
                else:
                    # don't rotate the LE coordinates
                    (yL_new, zL_new) = (yL_temp, zL_temp)
                    (yU_new, zU_new) = (yU_temp, zU_temp)
                # append the new LE coordinates
                yL.append(yL_new)
                zL.append(zL_new)
                yU.append(yU_new)
                zU.append(zU_new)
        return ((xL,yL,zL),(xU,yU,zU))

    def plot_LE(self, lw, twist_flag=True):
        """Plots the leading edge from root to tip."""
        ((xL,yL,zL),(xU,yU,zU)) = self.get_LE_coords(twist_flag=twist_flag)
        mlab.plot3d(xL,yL,zL, tube_radius=lw)
        mlab.plot3d(xU,yU,zU, tube_radius=lw)

    def get_TE_coords(self, twist_flag=True):
        """Returns a list of (x,y,z) coordinates for the blade trailing edge."""
        xL = []  # spanwise coordinate (monoplane and lower biplane airfoils)
        yL = []  # chordwise coordinate (monoplane and lower biplane airfoils)
        zL = []  # flapwise coordinate (monoplane and lower biplane airfoils)
        xU = []  # spanwise coordinate (upper biplane airfoils)
        yU = []  # chordwise coordinate (upper biplane airfoils)
        zU = []  # flapwise coordinate (upper biplane airfoils)
        inboard_trans = False # flag to find inboard mono-biplane transition
        outboard_trans = False # flag to find outboard bi-monoplane transition
        for station in self.list_of_stations:
            if station.type == 'monoplane':
                xL.append(station.coords.x1)
                # grab the unrotated TE coordinates
                yL_temp = station.airfoil.chord * (1.0-station.airfoil.pitch_axis)
                zL_temp = 0.0
                if twist_flag:
                    # rotate the TE coordinates wrt the twist angle
                    (yL_new, zL_new) = tf.rotate_coord_pair(
                        yL_temp, zL_temp, station.airfoil.twist)
                else:
                    # don't rotate the TE coordinates
                    (yL_new, zL_new) = (yL_temp, zL_temp)
                # append the new TE coordinates
                yL.append(yL_new)
                zL.append(zL_new)
                if inboard_trans is True and outboard_trans is False:
                    outboard_trans = True
                    # append the last LOWER coordinates that were added
                    xU.append(xL[-1])
                    yU.append(yL[-1])
                    zU.append(zL[-1])
            elif station.type == 'biplane':
                if inboard_trans is False:
                    inboard_trans = True
                    # append the last LOWER coordinates that were added
                    xU.append(xL[-1])
                    yU.append(yL[-1])
                    zU.append(zL[-1])
                xL.append(station.coords.x1)
                xU.append(station.coords.x1)
                # grab the unrotated TE coordinates
                yL_temp = station.airfoil.total_chord * (1.0-station.airfoil.pitch_axis)
                yU_temp = yL_temp - station.airfoil.stagger
                zU_temp = station.airfoil.gap_fraction * station.airfoil.gap
                zL_temp = -(1.0-station.airfoil.gap_fraction) * station.airfoil.gap
                if twist_flag:
                    # rotate the LE coordinates wrt the twist angle
                    (yL_new, zL_new) = tf.rotate_coord_pair(
                        yL_temp, zL_temp, station.airfoil.twist)
                    (yU_new, zU_new) = tf.rotate_coord_pair(
                        yU_temp, zU_temp, station.airfoil.twist)
                else:
                    # don't rotate the LE coordinates
                    (yL_new, zL_new) = (yL_temp, zL_temp)
                    (yU_new, zU_new) = (yU_temp, zU_temp)
                # append the new LE coordinates
                yL.append(yL_new)
                zL.append(zL_new)
                yU.append(yU_new)
                zU.append(zU_new)
        return ((xL,yL,zL),(xU,yU,zU))

    def plot_TE(self, lw, twist_flag=True):
        """Plots the trailing edge from root to tip."""
        ((xL,yL,zL),(xU,yU,zU)) = self.get_TE_coords(twist_flag=twist_flag)
        mlab.plot3d(xL,yL,zL, tube_radius=lw)
        mlab.plot3d(xU,yU,zU, tube_radius=lw)