"""A module for organizing geometrical data for a blade definition.

Author: Perry Roth-Johnson
Last updated: October 11, 2013

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
import material as mt
reload(mt)
import matplotlib.pyplot as plt
import matplotlib.colors as colors
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
        airfoils_path='airfoils', matl_filename='materials.csv',
        rotate_airfoil_coords=True):
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
            self.matl_filename = os.path.join(self.blade_path, matl_filename)
            self.logf.write("[{0}] Found blade definition file: {1}\n".format(datetime.datetime.now(), self.defn_filename))
            import_success = self.import_blade_definition()
            if import_success:
                self.create_all_stations()
                self.logf.write("[{0}] Created all blade stations\n".format(datetime.datetime.now()))
                self.airfoils_path = os.path.join(self.blade_path, airfoils_path)
                self.logf.write("[{0}] Found airfoils path: {1}\n".format(datetime.datetime.now(), self.airfoils_path))
                self.copy_all_airfoil_coords()
                # pre-process the airfoil coordinates and laminate schedule
                for station in self.list_of_stations:
                    station.airfoil.read_coords()
                    station.airfoil.scale_and_translate_coords()
                    station.airfoil.split_at_LE_and_TE()
                    if rotate_airfoil_coords:
                        station.airfoil.rotate_coords()
                    station.find_part_edges()
            material_import_success = self.import_material_properties()
            if material_import_success:
                self.create_all_materials()
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

        Returns True if the import was successful, False otherwise.

        """
        import_result = False
        defn_fileext = os.path.splitext(self.defn_filename)[-1]
        if defn_fileext != '.csv' and defn_fileext != '.CSV':
            # check that the blade definition file is a CSV file
            raise ValueError("blade definition file '{0}' must be of type *.csv".format(os.path.split(self.defn_filename)[-1]))
        else:
            # import the blade definition file into a pandas DataFrame
            self._df = pd.read_csv(self.defn_filename, index_col=0)
            self.number_of_stations = len(self._df)
            import_result = True
        return import_result

    def import_material_properties(self):
        """Import the material properties from a CSV file.

        Returns True if the import was successful, False otherwise.

        """
        import_result = False
        matl_fileext = os.path.splitext(self.matl_filename)[-1]
        if matl_fileext != '.csv' and matl_fileext != '.CSV':
            # check that the blade definition file is a CSV file
            raise ValueError("material properties file '{0}' must be of type *.csv".format(os.path.split(self.matl_filename)[-1]))
        else:
            # import the blade definition file into a pandas DataFrame
            self._mp = pd.read_csv(self.matl_filename, index_col=0)
            self.number_of_materials = len(self._mp)
            import_result = True
        return import_result

    def create_all_materials(self):
        """Create all materials for this blade."""
        if mt._Material.number_of_materials != 0:
            mt._Material.number_of_materials = 0  # initialize to zero
            print " [Warning] The number of materials has been reset to zero."
            self.logf = open(_Material.logfile_name, 'a')
            self.logf.write("[{0}] [Warning] The number of materials (_Material.number_of_materials) was reset to zero.\n".format(datetime.datetime.now()))
        self.dict_of_materials = {}
        for m in range(1, self.number_of_materials+1):
            if self._mp.ix[m]['type'] == 'isotropic':
                self.dict_of_materials[self._mp.ix[m]['name']] = self.create_isotropic_material(m)
            elif self._mp.ix[m]['type'] == 'orthotropic':
                self.dict_of_materials[self._mp.ix[m]['name']] = self.create_orthotropic_material(m)
            else:
                raise NotImplementedError("That material type is not implemented yet.\n  Try using either `isotropic` or `orthotropic` materials.")

    def create_isotropic_material(self, material_num):
        """Create a new isotropic material for this blade."""
        mt_series = self._mp.ix[material_num]
        return mt.IsotropicMaterial(
            name=mt_series['name'],
            E=mt_series['E1'],
            nu=mt_series['nu12'],
            rho=mt_series['rho'])

    def create_orthotropic_material(self, material_num):
        """Create a new orthotropic material for this blade."""
        mt_series = self._mp.ix[material_num]
        return mt.OrthotropicMaterial(
            name=mt_series['name'],
            E1=mt_series['E1'],
            E2=mt_series['E2'],
            E3=mt_series['E3'],
            G12=mt_series['G12'],
            G13=mt_series['G13'],
            G23=mt_series['G23'],
            nu12=mt_series['nu12'],
            nu13=mt_series['nu13'],
            nu23=mt_series['nu23'],
            rho=mt_series['rho'])

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
                  focalpoint=[60.0,0.85,0.0], print_view=False,
                  axes_triad=True):
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
        axes_triad : boolean, show/don't show a triad of orientation axes

        """
        mlab.view(azimuth, elevation, distance, focalpoint)
        if print_view:
            print "[**Mayavi mlab view parameters**]"
            print mlab.view()
        if axes_triad:
            mlab.orientation_axes(xlabel='x1', ylabel='x2', zlabel='x3')
        mlab.show()

    def plot_pitch_axis(self, lw, color='r'):
        """Plots the pitch axis from root to tip.

        Parameters
        ----------
        lw : float, line width
        color : str or RGB tuple, line color (default: 'r')

        """
        # convert the color to a RGB tuple
        if type(color) == str:
            c = colors.ColorConverter().to_rgb(color)
        else:
            c = color
        root = self.list_of_stations[0].coords.x1
        tip = self.list_of_stations[-1].coords.x1
        mlab.plot3d([root,tip],[0,0],[0,0], color=c, tube_radius=lw)

    def plot_selected_cross_sections(self, figsize=(22,12), nrows=4, ncols=3,
        selected_stations=[1,7,11,13,16,20,21,23,26,28,30,33], save_flag=True,
        plot_edges=False, plot_parts=True):
        """Plots selected cross-sections of the blade."""
        if len(selected_stations) != nrows*ncols:
            raise ValueError("The number of items in 'selected_stations' must equal nrows*ncols.")
        fig, axes = plt.subplots(figsize=figsize, nrows=nrows, ncols=ncols,
            num='Selected cross-sections for {0}'.format(self.name))
        # format each subplot
        i = 0
        for row_of_axes in axes:
            for ax in row_of_axes:
                station = self.list_of_stations[selected_stations[i]-1]
                print " Plotting station #{0}".format(station.station_num)
                st = station.structure
                ax.set_title("Station #{0}, {1}, {2}% span".format(station.station_num, station.airfoil.name, station.coords.x1))
                ax.set_aspect('equal')
                ax.grid('on')
                ax.set_xlabel('x2 [meters]')
                ax.set_ylabel('x3 [meters]')
                if plot_edges:
                    station.airfoil.plot_coords(ax, split_flag=True)
                    station.plot_part_edges(ax)
                if plot_parts:
                    station.plot_polygon(station.airfoil.polygon, ax,
                        face_color='None', edge_color='#999999', alpha=0.8)
                    (minx, miny, maxx, maxy) = station.airfoil.polygon.bounds
                    ax.set_xlim([minx*1.2,maxx*1.2])
                    ax.set_ylim([miny*1.2,maxy*1.2])
                    try:
                        if st.external_surface.exists():
                            station.plot_polygon(
                                st.external_surface.layer[0].polygon, ax,
                                face_color='#5EE54C', edge_color='#000000',
                                alpha=0.8)  # face color is purple
                            station.plot_polygon(
                                st.external_surface.layer[1].polygon, ax,
                                face_color='#5EE54C', edge_color='#000000',
                                alpha=0.8)  # face color is purple
                        if st.root_buildup.exists():
                            station.plot_polygon(
                                st.root_buildup.layer[0].polygon, ax,
                                face_color='#BE925A', edge_color='#000000',
                                alpha=0.8)  # face color is brown
                        if st.spar_cap.exists():
                            station.plot_polygon(st.spar_cap.layer[0].polygon,
                                ax, face_color='#00ACEF', edge_color='#000000',
                                alpha=0.8)  # face color is blue
                            station.plot_polygon(st.spar_cap.layer[1].polygon,
                                ax, face_color='#00ACEF', edge_color='#000000',
                                alpha=0.8)  # face color is blue
                        if st.aft_panel_1.exists():
                            station.plot_polygon(
                                st.aft_panel_1.layer[0].polygon, ax,
                                face_color='#F58612', edge_color='#000000',
                                alpha=0.8)  # face color is orange
                            station.plot_polygon(
                                st.aft_panel_1.layer[1].polygon, ax,
                                face_color='#F58612', edge_color='#000000',
                                alpha=0.8)  # face color is orange
                        if st.aft_panel_2.exists():
                            station.plot_polygon(
                                st.aft_panel_2.layer[0].polygon, ax,
                                face_color='#F58612', edge_color='#000000',
                                alpha=0.8)  # face color is orange
                            station.plot_polygon(
                                st.aft_panel_2.layer[1].polygon, ax,
                                face_color='#F58612', edge_color='#000000',
                                alpha=0.8)  # face color is orange
                        if st.LE_panel.exists():
                            station.plot_polygon(st.LE_panel.layer[0].polygon,
                                ax, face_color='#00A64F', edge_color='#000000',
                                alpha=0.8)  # face color is green
                        if st.shear_web_1.exists():
                            station.plot_polygon(
                                st.shear_web_1.layer[0].polygon, ax,
                                face_color='#FFF100', edge_color='#000000',
                                alpha=0.8)  # face color is yellow
                            station.plot_polygon(
                                st.shear_web_1.layer[1].polygon, ax,
                                face_color='#FFF100', edge_color='#000000',
                                alpha=0.8)  # face color is yellow
                            station.plot_polygon(
                                st.shear_web_1.layer[2].polygon, ax,
                                face_color='#FFF100', edge_color='#000000',
                                alpha=0.8)  # face color is yellow
                        if st.shear_web_2.exists():
                            station.plot_polygon(
                                st.shear_web_2.layer[0].polygon, ax,
                                face_color='#FFF100', edge_color='#000000',
                                alpha=0.8)  # face color is yellow
                            station.plot_polygon(
                                st.shear_web_2.layer[1].polygon, ax,
                                face_color='#FFF100', edge_color='#000000',
                                alpha=0.8)  # face color is yellow
                            station.plot_polygon(
                                st.shear_web_2.layer[2].polygon, ax,
                                face_color='#FFF100', edge_color='#000000',
                                alpha=0.8)  # face color is yellow
                        if st.shear_web_3.exists():
                            station.plot_polygon(
                                st.shear_web_3.layer[0].polygon, ax,
                                face_color='#FFF100', edge_color='#000000',
                                alpha=0.8)  # face color is yellow
                            station.plot_polygon(
                                st.shear_web_3.layer[1].polygon, ax,
                                face_color='#FFF100', edge_color='#000000',
                                alpha=0.8)  # face color is yellow
                            station.plot_polygon(
                                st.shear_web_3.layer[2].polygon, ax,
                                face_color='#FFF100', edge_color='#000000',
                                alpha=0.8)  # face color is yellow
                        if st.TE_reinforcement.exists():
                            station.plot_polygon(
                                st.TE_reinforcement.layer[0].polygon, ax,
                                face_color='#F366BA', edge_color='#000000',
                                alpha=0.8)  # face color is pink
                            try:
                                station.plot_polygon(
                                    st.TE_reinforcement.layer[1].polygon, ax,
                                    face_color='#F366BA', edge_color='#000000',
                                    alpha=0.8)  # face color is pink
                            except IndexError:  # foam region doesn't exist
                                pass
                        if st.internal_surface_1.exists():
                            station.plot_polygon(
                                st.internal_surface_1.layer[0].polygon, ax,
                                face_color='#999999', edge_color='#000000',
                                alpha=0.8)  # face color is gray
                            station.plot_polygon(
                                st.internal_surface_1.layer[1].polygon, ax,
                                face_color='#999999', edge_color='#000000',
                                alpha=0.8)  # face color is gray
                        if st.internal_surface_2.exists():
                            station.plot_polygon(
                                st.internal_surface_2.layer[0].polygon, ax,
                                face_color='#999999', edge_color='#000000',
                                alpha=0.8)  # face color is gray
                            station.plot_polygon(
                                st.internal_surface_2.layer[1].polygon, ax,
                                face_color='#999999', edge_color='#000000',
                                alpha=0.8)  # face color is gray
                        if st.internal_surface_3.exists():
                            station.plot_polygon(
                                st.internal_surface_3.layer[0].polygon, ax,
                                face_color='#999999', edge_color='#000000',
                                alpha=0.8)  # face color is gray
                            station.plot_polygon(
                                st.internal_surface_3.layer[1].polygon, ax,
                                face_color='#999999', edge_color='#000000',
                                alpha=0.8)  # face color is gray
                        if st.internal_surface_4.exists():
                            station.plot_polygon(
                                st.internal_surface_4.layer[0].polygon, ax,
                                face_color='#999999', edge_color='#000000',
                                alpha=0.8)  # face color is gray
                            station.plot_polygon(
                                st.internal_surface_4.layer[1].polygon, ax,
                                face_color='#999999', edge_color='#000000',
                                alpha=0.8)  # face color is gray
                    except AttributeError:
                        raise AttributeError("Part instance has no attribute 'polygon'.\n  Try running <station>.find_all_part_polygons() first.")
                i += 1
        fig.tight_layout()  # don't allow figure axes and labels to overlap
        plt.show()
        if save_flag:
            # Save the plot in the station path as a PNG file
            fname = os.path.join(self.blade_path, 'selected_cross-sections.png')
            fig.savefig(fname)
        return (fig, axes)

    def calculate_all_areas(self):
        """Calculate the areas of all structural parts in each station of this blade."""
        for station in self.list_of_stations:
            station.structure.calculate_area()

    def calculate_all_masses(self):
        """Calculate masses of all structural parts in each station of this blade."""
        for station in self.list_of_stations:
            station.structure.calculate_mass()

    def get_all_percent_areas(self, save_csv=True):
        """Returns a pandas.DataFrame of percent areas for each part and station.

        Parameters
        ----------
        save_csv : bool (default: True), save the DataFrame as a CSV file

        """
        self.calculate_all_areas()
        list_of_dicts = []
        for station in self.list_of_stations:
            d = station.structure.calculate_all_percent_areas()
            list_of_dicts.append(d)
        pa = pd.DataFrame(list_of_dicts, index=range(1,self.number_of_stations+1))
        pa = pa.fillna(0)  # set all NaNs to zeros
        # condense the number of columns down
        # -- spar caps --
        pa['spar caps'] = pa['spar cap (lower)'] + pa['spar cap (upper)']
        pa.pop('spar cap (lower)')
        pa.pop('spar cap (upper)')
        # -- aft panels --
        pa['aft panels'] = (pa['aft panel 1 (lower)'] + pa['aft panel 1 (upper)'] +
                            pa['aft panel 2 (lower)'] + pa['aft panel 2 (upper)'])
        pa.pop('aft panel 1 (lower)')
        pa.pop('aft panel 1 (upper)')
        pa.pop('aft panel 2 (lower)')
        pa.pop('aft panel 2 (upper)')
        # -- shear webs (biax) --
        pa['shear webs (biax)'] = (
            pa['shear web 1 (left biax)'] + pa['shear web 1 (right biax)'] + 
            pa['shear web 2 (left biax)'] + pa['shear web 2 (right biax)'] + 
            pa['shear web 3 (left biax)'] + pa['shear web 3 (right biax)'])
        pa.pop('shear web 1 (left biax)')
        pa.pop('shear web 1 (right biax)')
        pa.pop('shear web 2 (left biax)')
        pa.pop('shear web 2 (right biax)')
        pa.pop('shear web 3 (left biax)')
        pa.pop('shear web 3 (right biax)')
        # -- shear webs (foam) --
        pa['shear webs (foam)'] = (pa['shear web 1 (foam)'] + 
            pa['shear web 2 (foam)'] + pa['shear web 3 (foam)'])
        pa.pop('shear web 1 (foam)')
        pa.pop('shear web 2 (foam)')
        pa.pop('shear web 3 (foam)')
        # -- internal surface (triax) --
        pa['internal surfaces (triax)'] = (pa['internal surface 1 (triax)'] + 
            pa['internal surface 2 (triax)'] + pa['internal surface 3 (triax)'] +
            pa['internal surface 4 (triax)'])
        pa.pop('internal surface 1 (triax)')
        pa.pop('internal surface 2 (triax)')
        pa.pop('internal surface 3 (triax)')
        pa.pop('internal surface 4 (triax)')
        # -- internal surface (resin) --
        pa['internal surfaces (resin)'] = (pa['internal surface 1 (resin)'] + 
            pa['internal surface 2 (resin)'] + pa['internal surface 3 (resin)'] +
            pa['internal surface 4 (resin)'])
        pa.pop('internal surface 1 (resin)')
        pa.pop('internal surface 2 (resin)')
        pa.pop('internal surface 3 (resin)')
        pa.pop('internal surface 4 (resin)')
        if save_csv:
            # save the data to a CSV file
            pa_path = os.path.join(self.blade_path, 'percent_areas.csv')
            pa.to_csv(pa_path, index_label='blade station', cols=[
                'external surface (gelcoat)',
                'external surface (triax)',
                'root buildup',  # (triax)
                'spar caps',  # (uniax)
                'aft panels',  # (foam)
                'LE panel',  # (foam)
                'shear webs (biax)',
                'shear webs (foam)',
                'TE reinforcement (uniax)',
                'TE reinforcement (foam)',
                'internal surfaces (triax)',
                'internal surfaces (resin)'])
        return pa

    def get_all_percent_masses(self, save_csv=True):
        """Returns a pandas.DataFrame of percent masses for each part and station.

        Parameters
        ----------
        save_csv : bool (default: True), save the DataFrame as a CSV file

        """
        self.calculate_all_masses()
        list_of_dicts = []
        for station in self.list_of_stations:
            d = station.structure.calculate_all_percent_masses()
            list_of_dicts.append(d)
        pa = pd.DataFrame(list_of_dicts, index=range(1,self.number_of_stations+1))
        pa = pa.fillna(0)  # set all NaNs to zeros
        # condense the number of columns down
        # -- spar caps --
        pa['spar caps'] = pa['spar cap (lower)'] + pa['spar cap (upper)']
        pa.pop('spar cap (lower)')
        pa.pop('spar cap (upper)')
        # -- aft panels --
        pa['aft panels'] = (pa['aft panel 1 (lower)'] + pa['aft panel 1 (upper)'] +
                            pa['aft panel 2 (lower)'] + pa['aft panel 2 (upper)'])
        pa.pop('aft panel 1 (lower)')
        pa.pop('aft panel 1 (upper)')
        pa.pop('aft panel 2 (lower)')
        pa.pop('aft panel 2 (upper)')
        # -- shear webs (biax) --
        pa['shear webs (biax)'] = (
            pa['shear web 1 (left biax)'] + pa['shear web 1 (right biax)'] + 
            pa['shear web 2 (left biax)'] + pa['shear web 2 (right biax)'] + 
            pa['shear web 3 (left biax)'] + pa['shear web 3 (right biax)'])
        pa.pop('shear web 1 (left biax)')
        pa.pop('shear web 1 (right biax)')
        pa.pop('shear web 2 (left biax)')
        pa.pop('shear web 2 (right biax)')
        pa.pop('shear web 3 (left biax)')
        pa.pop('shear web 3 (right biax)')
        # -- shear webs (foam) --
        pa['shear webs (foam)'] = (pa['shear web 1 (foam)'] + 
            pa['shear web 2 (foam)'] + pa['shear web 3 (foam)'])
        pa.pop('shear web 1 (foam)')
        pa.pop('shear web 2 (foam)')
        pa.pop('shear web 3 (foam)')
        # -- internal surface (triax) --
        pa['internal surfaces (triax)'] = (pa['internal surface 1 (triax)'] + 
            pa['internal surface 2 (triax)'] + pa['internal surface 3 (triax)'] +
            pa['internal surface 4 (triax)'])
        pa.pop('internal surface 1 (triax)')
        pa.pop('internal surface 2 (triax)')
        pa.pop('internal surface 3 (triax)')
        pa.pop('internal surface 4 (triax)')
        # -- internal surface (resin) --
        pa['internal surfaces (resin)'] = (pa['internal surface 1 (resin)'] + 
            pa['internal surface 2 (resin)'] + pa['internal surface 3 (resin)'] +
            pa['internal surface 4 (resin)'])
        pa.pop('internal surface 1 (resin)')
        pa.pop('internal surface 2 (resin)')
        pa.pop('internal surface 3 (resin)')
        pa.pop('internal surface 4 (resin)')
        if save_csv:
            # save the data to a CSV file
            pa_path = os.path.join(self.blade_path, 'percent_masses.csv')
            pa.to_csv(pa_path, index_label='blade station', cols=[
                'external surface (gelcoat)',
                'external surface (triax)',
                'root buildup',  # (triax)
                'spar caps',  # (uniax)
                'aft panels',  # (foam)
                'LE panel',  # (foam)
                'shear webs (biax)',
                'shear webs (foam)',
                'TE reinforcement (uniax)',
                'TE reinforcement (foam)',
                'internal surfaces (triax)',
                'internal surfaces (resin)'])
        return pa

    def plot_percent_areas(self):
        """Plot the percent areas as a 100% stacked bar plot.

        ref: http://matplotlib.org/examples/pylab_examples/bar_stacked.html

        """
        pa = self.get_all_percent_areas()
        plt.figure(figsize=(22,12))
        ind = np.arange(self.number_of_stations)  # the x locations for each blade station
        width = 0.45                           # the width of the bars
        p1 = plt.bar(ind, pa['external surface (gelcoat)'], width, color='#5EE54C')
        p2 = plt.bar(ind, pa['external surface (triax)'], width, color='#5EE54C', hatch='.',
            bottom=pa['external surface (gelcoat)'])
        p3 = plt.bar(ind, pa['root buildup'], width, color='#BE925A', 
            bottom=(pa['external surface (gelcoat)'] + pa['external surface (triax)']))
        p4 = plt.bar(ind, pa['spar caps'], width, color='#00ACEF', 
            bottom=(pa['external surface (gelcoat)'] + pa['external surface (triax)'] +
                pa['root buildup']))
        p5 = plt.bar(ind, pa['aft panels'], width, color='#F58612', 
            bottom=(pa['external surface (gelcoat)'] + pa['external surface (triax)'] +
                pa['root buildup'] + pa['spar caps']))
        p6 = plt.bar(ind, pa['LE panel'], width, color='#00A64F', 
            bottom=(pa['external surface (gelcoat)'] + pa['external surface (triax)'] +
                pa['root buildup'] + pa['spar caps'] + pa['aft panels']))
        p7 = plt.bar(ind, pa['shear webs (biax)'], width, color='#FFF100', hatch='x',
            bottom=(pa['external surface (gelcoat)'] + pa['external surface (triax)'] +
                pa['root buildup'] + pa['spar caps'] + pa['aft panels'] +
                pa['LE panel']))
        p8 = plt.bar(ind, pa['shear webs (foam)'], width, color='#FFF100', 
            bottom=(pa['external surface (gelcoat)'] + pa['external surface (triax)'] +
                pa['root buildup'] + pa['spar caps'] + pa['aft panels'] +
                pa['LE panel'] + pa['shear webs (biax)']))
        p9 = plt.bar(ind, pa['TE reinforcement (uniax)'], width, color='#F366BA', hatch='/',
            bottom=(pa['external surface (gelcoat)'] + pa['external surface (triax)'] +
                pa['root buildup'] + pa['spar caps'] + pa['aft panels'] +
                pa['LE panel'] + pa['shear webs (biax)'] + pa['shear webs (foam)']))
        p10 = plt.bar(ind, pa['TE reinforcement (foam)'], width, color='#F366BA', 
            bottom=(pa['external surface (gelcoat)'] + pa['external surface (triax)'] +
                pa['root buildup'] + pa['spar caps'] + pa['aft panels'] +
                pa['LE panel'] + pa['shear webs (biax)'] + pa['shear webs (foam)'] +
                pa['TE reinforcement (uniax)']))
        p11 = plt.bar(ind, pa['internal surfaces (triax)'], width, color='#999999', hatch='.',
            bottom=(pa['external surface (gelcoat)'] + pa['external surface (triax)'] +
                pa['root buildup'] + pa['spar caps'] + pa['aft panels'] +
                pa['LE panel'] + pa['shear webs (biax)'] + pa['shear webs (foam)'] +
                pa['TE reinforcement (uniax)'] + pa['TE reinforcement (foam)']))
        p12 = plt.bar(ind, pa['internal surfaces (resin)'], width, color='#999999', 
            bottom=(pa['external surface (gelcoat)'] + pa['external surface (triax)'] +
                pa['root buildup'] + pa['spar caps'] + pa['aft panels'] +
                pa['LE panel'] + pa['shear webs (biax)'] + pa['shear webs (foam)'] +
                pa['TE reinforcement (uniax)'] + pa['TE reinforcement (foam)'] +
                pa['internal surfaces (triax)']))
        plt.ylabel('percent area')
        plt.xlabel('blade station')
        plt.xticks(ind+width/2.0, ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34'))
        plt.yticks( np.arange(0.0,1.1,0.1), ('0%', '10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', '100%'))
        plt.ylim([0.0,1.0])
        plt.xlim([-width/2.0, self.number_of_stations-0.5+width/2.0])
        plt.legend( 
            (
            p12[0],
            p11[0], 
            p10[0], 
            p9[0], 
            p8[0], 
            p7[0], 
            p6[0], 
            p5[0], 
            p4[0], 
            p3[0], 
            p2[0], 
            p1[0]
            ), 
            (
            'internal surfaces (resin)',
            'internal surfaces (triax)',
            'TE reinforcement (foam)',
            'TE reinforcement (uniax)',
            'shear webs (foam)',
            'shear webs (biax)',
            'LE panel (foam)',
            'aft panels (foam)',
            'spar caps (uniax)',
            'root buildup (triax)',
            'external surface (triax)',
            'external surface (gelcoat)'
            ), 
            bbox_to_anchor=(1.02, 0.5),
            loc='center left',
            borderaxespad=0.0)
        plt.subplots_adjust(left=0.05, bottom=0.05, right=0.82, top=0.95)
        plt.grid(axis='y')
        plt.show()

    def plot_percent_masses(self):
        """Plot the percent masses as a 100% stacked bar plot.

        ref: http://matplotlib.org/examples/pylab_examples/bar_stacked.html

        """
        pa = self.get_all_percent_masses()
        plt.figure(figsize=(22,12))
        ind = np.arange(self.number_of_stations)  # the x locations for each blade station
        width = 0.45                           # the width of the bars
        p1 = plt.bar(ind, pa['external surface (gelcoat)'], width, color='#5EE54C')
        p2 = plt.bar(ind, pa['external surface (triax)'], width, color='#5EE54C', hatch='.',
            bottom=pa['external surface (gelcoat)'])
        p3 = plt.bar(ind, pa['root buildup'], width, color='#BE925A', 
            bottom=(pa['external surface (gelcoat)'] + pa['external surface (triax)']))
        p4 = plt.bar(ind, pa['spar caps'], width, color='#00ACEF', 
            bottom=(pa['external surface (gelcoat)'] + pa['external surface (triax)'] +
                pa['root buildup']))
        p5 = plt.bar(ind, pa['aft panels'], width, color='#F58612', 
            bottom=(pa['external surface (gelcoat)'] + pa['external surface (triax)'] +
                pa['root buildup'] + pa['spar caps']))
        p6 = plt.bar(ind, pa['LE panel'], width, color='#00A64F', 
            bottom=(pa['external surface (gelcoat)'] + pa['external surface (triax)'] +
                pa['root buildup'] + pa['spar caps'] + pa['aft panels']))
        p7 = plt.bar(ind, pa['shear webs (biax)'], width, color='#FFF100', hatch='x',
            bottom=(pa['external surface (gelcoat)'] + pa['external surface (triax)'] +
                pa['root buildup'] + pa['spar caps'] + pa['aft panels'] +
                pa['LE panel']))
        p8 = plt.bar(ind, pa['shear webs (foam)'], width, color='#FFF100', 
            bottom=(pa['external surface (gelcoat)'] + pa['external surface (triax)'] +
                pa['root buildup'] + pa['spar caps'] + pa['aft panels'] +
                pa['LE panel'] + pa['shear webs (biax)']))
        p9 = plt.bar(ind, pa['TE reinforcement (uniax)'], width, color='#F366BA', hatch='/',
            bottom=(pa['external surface (gelcoat)'] + pa['external surface (triax)'] +
                pa['root buildup'] + pa['spar caps'] + pa['aft panels'] +
                pa['LE panel'] + pa['shear webs (biax)'] + pa['shear webs (foam)']))
        p10 = plt.bar(ind, pa['TE reinforcement (foam)'], width, color='#F366BA', 
            bottom=(pa['external surface (gelcoat)'] + pa['external surface (triax)'] +
                pa['root buildup'] + pa['spar caps'] + pa['aft panels'] +
                pa['LE panel'] + pa['shear webs (biax)'] + pa['shear webs (foam)'] +
                pa['TE reinforcement (uniax)']))
        p11 = plt.bar(ind, pa['internal surfaces (triax)'], width, color='#999999', hatch='.',
            bottom=(pa['external surface (gelcoat)'] + pa['external surface (triax)'] +
                pa['root buildup'] + pa['spar caps'] + pa['aft panels'] +
                pa['LE panel'] + pa['shear webs (biax)'] + pa['shear webs (foam)'] +
                pa['TE reinforcement (uniax)'] + pa['TE reinforcement (foam)']))
        p12 = plt.bar(ind, pa['internal surfaces (resin)'], width, color='#999999', 
            bottom=(pa['external surface (gelcoat)'] + pa['external surface (triax)'] +
                pa['root buildup'] + pa['spar caps'] + pa['aft panels'] +
                pa['LE panel'] + pa['shear webs (biax)'] + pa['shear webs (foam)'] +
                pa['TE reinforcement (uniax)'] + pa['TE reinforcement (foam)'] +
                pa['internal surfaces (triax)']))
        plt.ylabel('percent mass (per unit length)')
        plt.xlabel('blade station')
        plt.xticks(ind+width/2.0, ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34'))
        plt.yticks( np.arange(0.0,1.1,0.1), ('0%', '10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', '100%'))
        plt.ylim([0.0,1.0])
        plt.xlim([-width/2.0, self.number_of_stations-0.5+width/2.0])
        plt.legend( 
            (
            p12[0],
            p11[0], 
            p10[0], 
            p9[0], 
            p8[0], 
            p7[0], 
            p6[0], 
            p5[0], 
            p4[0], 
            p3[0], 
            p2[0], 
            p1[0]
            ), 
            (
            'internal surfaces (resin)',
            'internal surfaces (triax)',
            'TE reinforcement (foam)',
            'TE reinforcement (uniax)',
            'shear webs (foam)',
            'shear webs (biax)',
            'LE panel (foam)',
            'aft panels (foam)',
            'spar caps (uniax)',
            'root buildup (triax)',
            'external surface (triax)',
            'external surface (gelcoat)'
            ), 
            bbox_to_anchor=(1.02, 0.5),
            loc='center left',
            borderaxespad=0.0)
        plt.subplots_adjust(left=0.05, bottom=0.05, right=0.82, top=0.95)
        plt.grid(axis='y')
        plt.show()


class MonoplaneBlade(_Blade):
    """Define a monoplane (conventional) wind turbine blade."""
    def create_station(self, station_num):
        """Create a new station for this blade."""
        return stn.MonoplaneStation(self._df.ix[station_num], self.blade_path,
            parent_blade=self)

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
            # self.logf.flush()
            # self.logf.close()

    def plot_all_airfoils(self, lw, color='k', twist_flag=True,
        export_flag=True):
        """Plot all the airfoils in the blade with Mayavi's mlab.

        You must import the blade and its airfoil coordinates first.

        Parameters
        ----------
        lw : float, line width
        color : str or RGB tuple, line color (default: 'k')
        twist_flag : bool, do/don't twist the airfoils about the pitch axis
        export_flag : bool, do/don't write airfoil coords to a text file in the
            station path, which can be imported into SolidWorks as an XYZ curve

        Usage
        -----
        b = bl.MonoplaneBlade('Sandia blade SNL100-00', 'sandia_blade')
        b.copy_all_airfoil_coords()
        for station in b.list_of_stations:
            station.airfoil.read_coords()
        b.plot_all_airfoils()

        """
        # convert the color to a RGB tuple
        if type(color) == str:
            c = colors.ColorConverter().to_rgb(color)
        else:
            c = color
        for station in self.list_of_stations:
            if twist_flag:
                station.airfoil.rotate_coords()  # apply twist angle to airfoil
            # assemble the airfoil coordinates for mlab
            try:
                y = station.airfoil.coords['x']  # chordwise coordinate
            except AttributeError:
                raise AttributeError("Airfoil coordinates haven't been read yet!\n Run <Station>.airfoil.read_coords() first.")
            z = station.airfoil.coords['y']  # flapwise coordinate
            l = len(y)
            x = np.ones((l,))*station.coords.x1  # spanwise coordinate
            # plot the airfoil on the screen
            mlab.plot3d(x,y,z, color=c, tube_radius=lw)
            if export_flag:
                filename = os.path.join(station.station_path,
                    'stn{0:02d}_coords.txt'.format(station.station_num))
                f = open(filename, 'w')
                for i in range(l):
                    f.write('{0:12.9f}\t{1:12.9f}\t{2:12.9f}\n'.format(x[i],y[i],z[i]))
                f.close()
                print ' Wrote airfoil coordinates to {0}'.format(filename)
                self.logf = open(_Blade.logfile_name, "a")
                self.logf.write("[{0}] Wrote airfoil coordinates to: {1}\n".format(datetime.datetime.now(), filename))
                self.logf.flush()
                self.logf.close()

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

    def plot_LE(self, lw, color='k', twist_flag=True):
        """Plots the leading edge from root to tip.

        Parameters
        ----------
        lw : float, line width
        color : str or RGB tuple, line color (default: 'k')
        twist_flag : bool, do/don't twist the airfoils about the pitch axis

        """
        # convert the color to a RGB tuple
        if type(color) == str:
            c = colors.ColorConverter().to_rgb(color)
        else:
            c = color
        (x,y,z) = self.get_LE_coords(twist_flag=twist_flag)
        mlab.plot3d(x,y,z, color=c, tube_radius=lw)

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

    def plot_TE(self, lw, color='k', twist_flag=True):
        """Plots the trailing edge from root to tip.

        Parameters
        ----------
        lw : float, line width
        color : str or RGB tuple, line color (default: 'k')
        twist_flag : bool, do/don't twist the airfoils about the pitch axis

        """
        # convert the color to a RGB tuple
        if type(color) == str:
            c = colors.ColorConverter().to_rgb(color)
        else:
            c = color
        (x,y,z) = self.get_TE_coords(twist_flag=twist_flag)
        mlab.plot3d(x,y,z, color=c, tube_radius=lw)

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
            SW_cs_coords = np.array([])
            if sw_num is 1 and station.structure.shear_web_1.exists():
                try:
                    SW_cs_coords = station.structure.shear_web_1.cs_coords
                except AttributeError:
                    raise AttributeError("Shear web #1's cross-section coordinates have not been read yet!\n  Try running <Station>.find_all_part_cs_coords() first.")
            elif sw_num is 2 and station.structure.shear_web_2.exists():
                try:
                    SW_cs_coords = station.structure.shear_web_2.cs_coords
                except AttributeError:
                    raise AttributeError("Shear web #2's cross-section coordinates have not been read yet!\n  Try running <Station>.find_all_part_cs_coords() first.")
            elif sw_num is 3 and station.structure.shear_web_3.exists():
                try:
                    SW_cs_coords = station.structure.shear_web_3.cs_coords
                except AttributeError:
                    raise AttributeError("Shear web #3's cross-section coordinates have not been read yet!\n  Try running <Station>.find_all_part_cs_coords() first.")
            for point in SW_cs_coords:
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

    def plot_all_SW_cross_sections(self, lw, color='g', twist_flag=True):
        """Plots all shear web cross-sections from root to tip.

        Parameters
        ----------
        lw : float, line width
        color : str or RGB tuple, line color (default: 'g')
        twist_flag : bool, do/don't twist the airfoils about the pitch axis
        
        """
        # convert the color to a RGB tuple
        if type(color) == str:
            c = colors.ColorConverter().to_rgb(color)
        else:
            c = color
        for sw in [1,2,3]:
            (x,y,z) = self.get_SW_cross_section_coords(sw_num=sw,
                                                       twist_flag=twist_flag)
            for i in range(0,len(x),4):
                (X,Y,Z) = (x[i:i+4], y[i:i+4], z[i:i+4])
                # add the first coord pair again, to form a closed loop
                X.append(x[i])
                Y.append(y[i])
                Z.append(z[i])
                mlab.plot3d(X,Y,Z, color=c, tube_radius=lw)

    def plot_all_SW_spanwise_lines(self, lw, color='g', twist_flag=True):
        """Plots spanwise lines for all shear webs from root to tip.

        Parameters
        ----------
        lw : float, line width
        color : str or RGB tuple, line color (default: 'g')
        twist_flag : bool, do/don't twist the airfoils about the pitch axis

        """
        # convert the color to a RGB tuple
        if type(color) == str:
            c = colors.ColorConverter().to_rgb(color)
        else:
            c = color
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
            mlab.plot3d(x1,y1,z1, color=c, tube_radius=lw)
            mlab.plot3d(x2,y2,z2, color=c, tube_radius=lw)
            mlab.plot3d(x3,y3,z3, color=c, tube_radius=lw)
            mlab.plot3d(x4,y4,z4, color=c, tube_radius=lw)

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
            SW_cs_coords = np.array([])
            if sw_num is 1 and station.structure.shear_web_1.exists():
                try:
                    SW_cs_coords = station.structure.shear_web_1.cs_coords
                except AttributeError:
                    raise AttributeError("Shear web #1's cross-section coordinates have not been read yet!\n  Try running <Station>.find_all_part_cs_coords() first.")
            elif sw_num is 2 and station.structure.shear_web_2.exists():
                try:
                    SW_cs_coords = station.structure.shear_web_2.cs_coords
                except AttributeError:
                    raise AttributeError("Shear web #2's cross-section coordinates have not been read yet!\n  Try running <Station>.find_all_part_cs_coords() first.")
            elif sw_num is 3 and station.structure.shear_web_3.exists():
                try:
                    SW_cs_coords = station.structure.shear_web_3.cs_coords
                except AttributeError:
                    raise AttributeError("Shear web #3's cross-section coordinates have not been read yet!\n  Try running <Station>.find_all_part_cs_coords() first.")
            for point in SW_cs_coords:
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

    def plot_station_nums(self):
        """Plots station numbers next to each cross-section."""
        for station in self.list_of_stations:
            mlab.text3d(
                x=station.coords.x1,
                y=(station.coords.x2 - 
                   station.airfoil.chord*station.airfoil.pitch_axis - 4.0),
                z=station.coords.x3,
                text='{0}'.format(station.station_num),
                scale=0.5,
                color=(120.0/256.0,0.0,146.0/256.0))

    def plot_blade(self, line_width=0.08, airfoils=True, pitch_axis=False,
        LE=True, TE=True, twist=True, SW=True, color_airfoils='0.1',
        color_pitch_axis='r', color_LE='0.1', color_TE='0.1',
        color_SW=(37.0/256.0,197.0/256.0,85.0/256.0), stn_nums=False,
        export=True):
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
        color_airfoils : str or RGB tuple, line color of airfoils
            (default: dark grey)
        color_pitch_axis : str or RGB tuple, line color of pitch axis
            (default: red)
        color_LE : str or RGB tuple, line color of leading edge
            (default: dark grey)
        color_TE : str or RGB tuple, line color of trailing edge
            (default: dark grey)
        color_SW : str or RGB tuple, line color of shear webs
            (default: light green)
        stn_nums : boolean, plot/don't plot station numbers
        export : bool, do/don't write airfoil coords to a text file in the
            station path, which can be imported into SolidWorks as an XYZ curve

        """
        self.create_plot()
        if airfoils:
            self.plot_all_airfoils(lw=line_width, color=color_airfoils,
                twist_flag=twist, export_flag=export)
        if pitch_axis:
            self.plot_pitch_axis(lw=line_width, color=color_pitch_axis)
        if LE:
            self.plot_LE(lw=line_width, color=color_LE, twist_flag=twist)
        if TE:
            self.plot_TE(lw=line_width, color=color_TE, twist_flag=twist)
        if SW:
            self.plot_all_SW_cross_sections(lw=line_width, color=color_SW,
                twist_flag=twist)
            self.plot_all_SW_spanwise_lines(lw=line_width, color=color_SW,
                twist_flag=twist)
        if stn_nums:
            self.plot_station_nums()
        self.show_plot()


class BiplaneBlade(_Blade):
    """Define a biplane wind turbine blade."""
    def __init__(self, name, blade_path, defn_filename='blade_definition.csv',
        airfoils_path='airfoils'):
        """Create a new biplane wind turbine blade."""
        _Blade.__init__(self, name, blade_path, defn_filename, airfoils_path)
        (self.root_joint_station,
            self.midblade_joint_station) = self.assign_joint_stations()
        print " Root joint found at station #{0}".format(
            self.root_joint_station)
        print " Mid-blade joint found at station #{0}".format(
            self.midblade_joint_station)
        self.logf = open(_Blade.logfile_name, "a")
        self.logf.write("[{0}] Root joint found at station #{1}\n".format(
            datetime.datetime.now(), self.root_joint_station))
        self.logf.write("[{0}] Mid-blade joint found at station #{1}\n".format(
            datetime.datetime.now(), self.midblade_joint_station))
        self.logf.flush()
        self.logf.close()

    def create_station(self, station_num):
        """Create a new station for this blade."""
        if self._df.ix[station_num]['type'] == 'monoplane':
            this_stn = stn.MonoplaneStation(self._df.ix[station_num], self.blade_path)
        elif self._df.ix[station_num]['type'] == 'biplane':
            this_stn = stn.BiplaneStation(self._df.ix[station_num], self.blade_path)
        else:
            raise ValueError("Values in the 'type' column of {0} must be either 'monoplane' or 'biplane'.".format(self.defn_filename))
        return this_stn

    def assign_joint_stations(self):
        """Identify and mark the root and mid-blade joint stations.

        Creates a '<_Station>.joint' attribute, which is 'root',
        'midblade', or None.

        """
        root_joint = 0
        midblade_joint = 0
        # identify joint stations
        for station in self.list_of_stations:
            if station.type == 'biplane' and root_joint == 0:
                root_joint = station.station_num - 1
            if (station.type == 'monoplane' and root_joint != 0 and 
                midblade_joint == 0):
                midblade_joint = station.station_num
        # mark joint stations, by creating the '.joint' attribute
        for station in self.list_of_stations:
            if station.station_num == root_joint:
                station.joint = 'root'
            elif station.station_num == midblade_joint:
                station.joint = 'midblade'
            else:
                station.joint = None
        return (root_joint, midblade_joint)

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
                # self.logf.flush()
                # self.logf.close()
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
                # self.logf.flush()
                # self.logf.close()
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
                # self.logf.flush()
                # self.logf.close()

    def plot_all_airfoils(self, lw, color='k', twist_flag=True,
        export_flag=True):
        """Plot all the airfoils in the biplane blade with Mayavi's mlab.

        You must import the blade and its airfoil coordinates first.

        Parameters
        ----------
        lw : float, line width
        color : str or RGB tuple, line color (default: 'k')
        twist_flag : bool, do/don't twist the airfoils about the pitch axis
        export_flag : bool, do/don't write airfoil coords to a text file in the
            station path, which can be imported into SolidWorks as an XYZ curve

        Usage
        -----
        a = bl.BiplaneBlade('biplane blade', 'biplane_blade')
        a.copy_all_airfoil_coords()
        for station in a.list_of_stations:
            station.airfoil.read_coords()
        a.plot_all_airfoils()

        """
        # convert the color to a RGB tuple
        if type(color) == str:
            c = colors.ColorConverter().to_rgb(color)
        else:
            c = color
        for station in self.list_of_stations:
            if station.type == 'monoplane':
                if twist_flag:
                    station.airfoil.rotate_coords()  # apply twist angle to airfoil
                # assemble the airfoil coordinates for mlab
                try:
                    y = station.airfoil.coords['x']  # chordwise coordinate
                except AttributeError:
                    raise AttributeError("Airfoil coordinates haven't been read yet!\n Run <Station>.airfoil.read_coords() first.")
                z = station.airfoil.coords['y']  # flapwise coordinate
                l = len(y)
                x = np.ones((l,))*station.coords.x1  # spanwise coordinate
                # plot the airfoil on the screen
                mlab.plot3d(x,y,z, color=c, tube_radius=lw)
                if export_flag:
                    filename = os.path.join(station.station_path,
                        'stn{0:02d}_coords.txt'.format(station.station_num))
                    f = open(filename, 'w')
                    for i in range(l):
                        f.write('{0:12.9f}\t{1:12.9f}\t{2:12.9f}\n'.format(x[i],y[i],z[i]))
                    f.close()
                    print ' Wrote airfoil coordinates to {0}'.format(filename)
                    self.logf = open(_Blade.logfile_name, "a")
                    self.logf.write("[{0}] Wrote airfoil coordinates to: {1}\n".format(datetime.datetime.now(), filename))
                    self.logf.flush()
                    self.logf.close()
            elif station.type == 'biplane':
                if twist_flag:
                    station.airfoil.rotate_coords()
                # assemble lower airfoil coordinates for mlab -----------------
                try:
                    y = station.airfoil.lower_coords['x']  # chordwise coordinate
                except AttributeError:
                    raise AttributeError("Lower airfoil coordinates haven't been read yet!\n  Run <Station>.airfoil.read_coords() first.")
                z = station.airfoil.lower_coords['y']  # flapwise coordinate
                l = len(y)
                x = np.ones((l,))*station.coords.x1  # spanwise coordinate
                # plot the lower airfoil on the screen
                mlab.plot3d(x,y,z, color=c, tube_radius=lw)
                if export_flag:
                    filename = os.path.join(station.station_path,
                        'stn{0:02d}_lower_coords.txt'.format(station.station_num))
                    f = open(filename, 'w')
                    for i in range(l):
                        f.write('{0:12.9f}\t{1:12.9f}\t{2:12.9f}\n'.format(x[i],y[i],z[i]))
                    f.close()
                    print ' Wrote lower airfoil coordinates to {0}'.format(filename)
                    self.logf = open(_Blade.logfile_name, "a")
                    self.logf.write("[{0}] Wrote lower airfoil coordinates to: {1}\n".format(datetime.datetime.now(), filename))
                    self.logf.flush()
                    self.logf.close()
                # assemble upper airfoil coordinates for mlab -----------------
                try:
                    y = station.airfoil.upper_coords['x']  # chordwise coordinate
                except AttributeError:
                    raise AttributeError("Upper airfoil coordinates haven't been read yet!\n  Run <Station>.airfoil.read_coords() first.")
                z = station.airfoil.upper_coords['y']  # flapwise coordinate
                l = len(y)
                # (reuse same spanwise coordinates): x
                # plot the lower airfoil on the screen
                mlab.plot3d(x,y,z, color=c, tube_radius=lw)
                if export_flag:
                    filename = os.path.join(station.station_path,
                        'stn{0:02d}_upper_coords.txt'.format(station.station_num))
                    f = open(filename, 'w')
                    for i in range(l):
                        f.write('{0:12.9f}\t{1:12.9f}\t{2:12.9f}\n'.format(x[i],y[i],z[i]))
                    f.close()
                    print ' Wrote upper airfoil coordinates to {0}'.format(filename)
                    self.logf = open(_Blade.logfile_name, "a")
                    self.logf.write("[{0}] Wrote upper airfoil coordinates to: {1}\n".format(datetime.datetime.now(), filename))
                    self.logf.flush()
                    self.logf.close()

    def get_LE_coords(self, twist_flag=True):
        """Returns a list of (x,y,z) coordinates for the blade leading edge."""
        xL = []  # spanwise coordinate (monoplane and lower biplane airfoils)
        yL = []  # chordwise coordinate (monoplane and lower biplane airfoils)
        zL = []  # flapwise coordinate (monoplane and lower biplane airfoils)
        xU = []  # spanwise coordinate (upper biplane airfoils)
        yU = []  # chordwise coordinate (upper biplane airfoils)
        zU = []  # flapwise coordinate (upper biplane airfoils)
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
                if station.joint is not None:
                    # append the last LOWER coordinates to the UPPER coords
                    xU.append(xL[-1])
                    yU.append(yL[-1])
                    zU.append(zL[-1])
            elif station.type == 'biplane':
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

    def plot_LE(self, lw, color='k', twist_flag=True):
        """Plots the leading edge from root to tip.

        Parameters
        ----------
        lw : float, line width
        color : str or RGB tuple, line color (default: 'k')
        twist_flag : bool, do/don't twist the airfoils about the pitch axis

        """
        # convert the color to a RGB tuple
        if type(color) == str:
            c = colors.ColorConverter().to_rgb(color)
        else:
            c = color
        ((xL,yL,zL),(xU,yU,zU)) = self.get_LE_coords(twist_flag=twist_flag)
        mlab.plot3d(xL,yL,zL, color=c, tube_radius=lw)
        mlab.plot3d(xU,yU,zU, color=c, tube_radius=lw)

    def get_TE_coords(self, twist_flag=True):
        """Returns a list of (x,y,z) coordinates for the blade trailing edge."""
        xL = []  # spanwise coordinate (monoplane and lower biplane airfoils)
        yL = []  # chordwise coordinate (monoplane and lower biplane airfoils)
        zL = []  # flapwise coordinate (monoplane and lower biplane airfoils)
        xU = []  # spanwise coordinate (upper biplane airfoils)
        yU = []  # chordwise coordinate (upper biplane airfoils)
        zU = []  # flapwise coordinate (upper biplane airfoils)
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
                if station.joint is not None:
                    # append the last LOWER coordinates to the UPPER coords
                    xU.append(xL[-1])
                    yU.append(yL[-1])
                    zU.append(zL[-1])
            elif station.type == 'biplane':
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

    def plot_TE(self, lw, color='k', twist_flag=True):
        """Plots the trailing edge from root to tip.

        Parameters
        ----------
        lw : float, line width
        color : str or RGB tuple, line color (default: 'k')
        twist_flag : bool, do/don't twist the airfoils about the pitch axis

        """
        # convert the color to a RGB tuple
        if type(color) == str:
            c = colors.ColorConverter().to_rgb(color)
        else:
            c = color
        ((xL,yL,zL),(xU,yU,zU)) = self.get_TE_coords(twist_flag=twist_flag)
        mlab.plot3d(xL,yL,zL, color=c, tube_radius=lw)
        mlab.plot3d(xU,yU,zU, color=c, tube_radius=lw)

    def get_SW_cross_section_coords(self, sw_num, twist_flag=True):
        """Returns a list of (x,y,z) coordinates for shear web #1, #2, or #3.

        Parameters
        ----------
        sw_num = int, (1, 2, or 3) selects the desired shear web

        """
        if sw_num is not 1 and sw_num is not 2 and sw_num is not 3:
            raise ValueError("sw_num must be either 1, 2, or 3 (type: int)")
        xL = []  # spanwise coordinate (monoplane and lower biplane airfoils)
        yL = []  # chordwise coordinate (monoplane and lower biplane airfoils)
        zL = []  # flapwise coordinate (monoplane and lower biplane airfoils)
        xU = []  # spanwise coordinate (upper biplane airfoils)
        yU = []  # chordwise coordinate (upper biplane airfoils)
        zU = []  # flapwise coordinate (upper biplane airfoils)
        for station in self.list_of_stations:
            # shear web CS coords for monoplane and lower biplane airfoils
            SW_cs_coords_L = np.array([])
            # shear web CS coords for upper biplane airfoils
            SW_cs_coords_U = np.array([])
            if station.type == 'monoplane':
                if sw_num is 1 and station.structure.shear_web_1.exists():
                    try:
                        SW_cs_coords_L = station.structure.shear_web_1.cs_coords
                    except AttributeError:
                        raise AttributeError("Shear web #1's cross-section coordinates have not been read yet!\n  Try running <Station>.find_all_part_cs_coords() first.")
                elif sw_num is 2 and station.structure.shear_web_2.exists():
                    try:
                        SW_cs_coords_L = station.structure.shear_web_2.cs_coords
                    except AttributeError:
                        raise AttributeError("Shear web #2's cross-section coordinates have not been read yet!\n  Try running <Station>.find_all_part_cs_coords() first.")
                elif sw_num is 3 and station.structure.shear_web_3.exists():
                    try:
                        SW_cs_coords_L = station.structure.shear_web_3.cs_coords
                    except AttributeError:
                        raise AttributeError("Shear web #3's cross-section coordinates have not been read yet!\n  Try running <Station>.find_all_part_cs_coords() first.")
            elif station.type == 'biplane':
                if sw_num is 1:
                    if station.structure.lower_shear_web_1.exists():
                        try:
                            SW_cs_coords_L = station.structure.lower_shear_web_1.cs_coords
                        except AttributeError:
                            raise AttributeError("Lower shear web #1's cross-section coordinates have not been read yet!\n  Try running <Station>.find_all_part_cs_coords() first.")
                    if station.structure.upper_shear_web_1.exists():
                        try:
                            SW_cs_coords_U = station.structure.upper_shear_web_1.cs_coords
                        except AttributeError:
                            raise AttributeError("Upper shear web #1's cross-section coordinates have not been read yet!\n  Try running <Station>.find_all_part_cs_coords() first.")
                elif sw_num is 2:
                    if station.structure.lower_shear_web_2.exists():
                        try:
                            SW_cs_coords_L = station.structure.lower_shear_web_2.cs_coords
                        except AttributeError:
                            raise AttributeError("Lower shear web #2's cross-section coordinates have not been read yet!\n  Try running <Station>.find_all_part_cs_coords() first.")
                    if station.structure.upper_shear_web_2.exists():
                        try:
                            SW_cs_coords_U = station.structure.upper_shear_web_2.cs_coords
                        except AttributeError:
                            raise AttributeError("Upper shear web #2's cross-section coordinates have not been read yet!\n  Try running <Station>.find_all_part_cs_coords() first.")
                elif sw_num is 3:
                    if station.structure.lower_shear_web_3.exists():
                        try:
                            SW_cs_coords_L = station.structure.lower_shear_web_3.cs_coords
                        except AttributeError:
                            raise AttributeError("Lower shear web #3's cross-section coordinates have not been read yet!\n  Try running <Station>.find_all_part_cs_coords() first.")
                    if station.structure.upper_shear_web_3.exists():
                        try:
                            SW_cs_coords_U = station.structure.upper_shear_web_3.cs_coords
                        except AttributeError:
                            raise AttributeError("Upper shear web #3's cross-section coordinates have not been read yet!\n  Try running <Station>.find_all_part_cs_coords() first.")
            # assemble the lower coordinates
            for point in SW_cs_coords_L:
                xL.append(station.coords.x1)
                yL_temp = point[0]
                zL_temp = point[1]
                if twist_flag:
                    # rotate the lower SW cross-section coords wrt twist angle
                    (yL_new, zL_new) = tf.rotate_coord_pair(yL_temp,
                        zL_temp, station.airfoil.twist)
                else:
                    # don't rotate the SW cross-section coords
                    (yL_new, zL_new) = (yL_temp, zL_temp)
                # append the new SW cross-section coords
                yL.append(yL_new)
                zL.append(zL_new)
            # assemble the upper coordinates
            for point in SW_cs_coords_U:
                xU.append(station.coords.x1)
                yU_temp = point[0]
                zU_temp = point[1]
                if twist_flag:
                    # rotate the lower SW cross-section coords wrt twist angle
                    (yU_new, zU_new) = tf.rotate_coord_pair(yU_temp,
                        zU_temp, station.airfoil.twist)
                else:
                    # don't rotate the SW cross-section coords
                    (yU_new, zU_new) = (yU_temp, zU_temp)
                # append the new SW cross-section coords
                yU.append(yU_new)
                zU.append(zU_new)
        return ((xL,yL,zL),(xU,yU,zU))

    def get_SW_cross_section_coords_at_joint(self, joint, sw_num,
        twist_flag=True):
        """Returns (x,y,z) coordinates for shear web #1, #2, or #3 at joints.

        Parameters
        ----------
        joint : str, ('root' or 'midblade') selects the desired joint
        sw_num : int, (1, 2, or 3) selects the desired shear web

        """
        if sw_num is not 1 and sw_num is not 2 and sw_num is not 3:
            raise ValueError("sw_num must be either 1, 2, or 3 (type: int)")
        x = []  # spanwise coordinate
        y = []  # chordwise coordinate
        z = []  # flapwise coordinate
        if joint != 'root' and joint != 'midblade':
            raise ValueError("joint must be 'root' or 'midblade' (type: str)")
        if joint == 'root':
            stn_num = self.root_joint_station
        elif joint == 'midblade':
            stn_num = self.midblade_joint_station
        station = self.list_of_stations[stn_num-1]
        SW_cs_coords = np.array([])
        if sw_num is 1 and station.structure.shear_web_1.exists():
            try:
                SW_cs_coords = station.structure.shear_web_1.cs_coords
            except AttributeError:
                raise AttributeError("Shear web #1's cross-section coordinates have not been read yet!\n  Try running <Station>.find_all_part_cs_coords() first.")
        elif sw_num is 2 and station.structure.shear_web_2.exists():
            try:
                SW_cs_coords = station.structure.shear_web_2.cs_coords
            except AttributeError:
                raise AttributeError("Shear web #2's cross-section coordinates have not been read yet!\n  Try running <Station>.find_all_part_cs_coords() first.")
        elif sw_num is 3 and station.structure.shear_web_3.exists():
            try:
                SW_cs_coords = station.structure.shear_web_3.cs_coords
            except AttributeError:
                raise AttributeError("Shear web #3's cross-section coordinates have not been read yet!\n  Try running <Station>.find_all_part_cs_coords() first.")
        for point in SW_cs_coords:
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

    def plot_all_SW_cross_sections(self, lw, color='g', twist_flag=True,
        lower=True, upper=True):
        """Plots all shear web cross-sections from root to tip.

        Parameters
        ----------
        lw : float, line width
        color : str or RGB tuple, line color (default: 'g')
        twist_flag : bool, do/don't twist the airfoils about the pitch axis

        """
        # convert the color to a RGB tuple
        if type(color) == str:
            c = colors.ColorConverter().to_rgb(color)
        else:
            c = color
        for sw in [1,2,3]:
            ((xL,yL,zL),(xU,yU,zU)) = self.get_SW_cross_section_coords(
                sw_num=sw, twist_flag=twist_flag)
            if lower:
                # plot the lower shear web
                for i in range(0,len(xL),4):
                    (X,Y,Z) = (xL[i:i+4], yL[i:i+4], zL[i:i+4])
                    # add the first coord pair again, to form a closed loop
                    X.append(xL[i])
                    Y.append(yL[i])
                    Z.append(zL[i])
                    mlab.plot3d(X,Y,Z, color=c, tube_radius=lw)
            if upper:
                # plot the upper shear web
                for i in range(0,len(xU),4):
                    (X,Y,Z) = (xU[i:i+4], yU[i:i+4], zU[i:i+4])
                    # add the first coord pair again, to form a closed loop
                    X.append(xU[i])
                    Y.append(yU[i])
                    Z.append(zU[i])
                    mlab.plot3d(X,Y,Z, color=c, tube_radius=lw)

    def plot_all_SW_spanwise_lines(self, lw, color='g', twist_flag=True,
        lower=True, upper=True, sw_list=[1,2,3]):
        """Plots spanwise lines for all shear webs from root to tip.

        Parameters
        ----------
        lw : float, line width
        color : str or RGB tuple, line color (default: 'g')
        twist_flag : bool, do/don't twist the airfoils about the pitch axis

        """
        # convert the color to a RGB tuple
        if type(color) == str:
            c = colors.ColorConverter().to_rgb(color)
        else:
            c = color
        for sw in sw_list:
            ((xL,yL,zL),(xU,yU,zU)) = self.get_SW_cross_section_coords(
                sw_num=sw, twist_flag=twist_flag)
            if lower:
                # plot the lower shear web
                (x1, y1, z1) = ([], [], [])
                (x2, y2, z2) = ([], [], [])
                (x3, y3, z3) = ([], [], [])
                (x4, y4, z4) = ([], [], [])
                for i in range(0,len(xL),4):
                    x1.append(xL[i])
                    y1.append(yL[i])
                    z1.append(zL[i])
                    x2.append(xL[i+1])
                    y2.append(yL[i+1])
                    z2.append(zL[i+1])
                    x3.append(xL[i+2])
                    y3.append(yL[i+2])
                    z3.append(zL[i+2])
                    x4.append(xL[i+3])
                    y4.append(yL[i+3])
                    z4.append(zL[i+3])
                mlab.plot3d(x1,y1,z1, color=c, tube_radius=lw)
                mlab.plot3d(x2,y2,z2, color=c, tube_radius=lw)
                mlab.plot3d(x3,y3,z3, color=c, tube_radius=lw)
                mlab.plot3d(x4,y4,z4, color=c, tube_radius=lw)
            if upper:
                # plot the upper shear web
                (x1, y1, z1) = ([], [], [])
                (x2, y2, z2) = ([], [], [])
                (x3, y3, z3) = ([], [], [])
                (x4, y4, z4) = ([], [], [])
                # append coords for shear web (if it exists) at the root joint
                (xr,yr,zr) = self.get_SW_cross_section_coords_at_joint(
                    joint='root', sw_num=sw, twist_flag=twist_flag)
                if len(xr) > 0:
                    x1.append(xr[0])
                    y1.append(yr[0])
                    z1.append(zr[0])
                    x2.append(xr[1])
                    y2.append(yr[1])
                    z2.append(zr[1])
                    x3.append(xr[2])
                    y3.append(yr[2])
                    z3.append(zr[2])
                    x4.append(xr[3])
                    y4.append(yr[3])
                    z4.append(zr[3])
                for i in range(0,len(xU),4):
                    x1.append(xU[i])
                    y1.append(yU[i])
                    z1.append(zU[i])
                    x2.append(xU[i+1])
                    y2.append(yU[i+1])
                    z2.append(zU[i+1])
                    x3.append(xU[i+2])
                    y3.append(yU[i+2])
                    z3.append(zU[i+2])
                    x4.append(xU[i+3])
                    y4.append(yU[i+3])
                    z4.append(zU[i+3])
                # append coords for shear web (if it exists) at the midblade joint
                (xm,ym,zm) = self.get_SW_cross_section_coords_at_joint(
                    joint='midblade', sw_num=sw, twist_flag=twist_flag)
                if len(xm) > 0:
                    x1.append(xm[0])
                    y1.append(ym[0])
                    z1.append(zm[0])
                    x2.append(xm[1])
                    y2.append(ym[1])
                    z2.append(zm[1])
                    x3.append(xm[2])
                    y3.append(ym[2])
                    z3.append(zm[2])
                    x4.append(xm[3])
                    y4.append(ym[3])
                    z4.append(zm[3])
                mlab.plot3d(x1,y1,z1, color=c, tube_radius=lw)
                mlab.plot3d(x2,y2,z2, color=c, tube_radius=lw)
                mlab.plot3d(x3,y3,z3, color=c, tube_radius=lw)
                mlab.plot3d(x4,y4,z4, color=c, tube_radius=lw)

    def plot_station_nums(self):
        """Plots station numbers next to each cross-section."""
        for station in self.list_of_stations:
            if station.type == 'monoplane':
                y_coord = (station.coords.x2 - 
                       station.airfoil.chord*station.airfoil.pitch_axis - 4.0)
            elif station.type == 'biplane':
                y_coord = (station.coords.x2 - 4.0 -
                   station.airfoil.total_chord*station.airfoil.pitch_axis)
            mlab.text3d(
                x=station.coords.x1,
                y=y_coord,
                z=station.coords.x3,
                text='{0}'.format(station.station_num),
                scale=0.5,
                color=(120.0/256.0,0.0,146.0/256.0))

    def plot_blade(self, line_width=0.08, airfoils=True, pitch_axis=False,
        LE=True, TE=True, twist=True, SW=True, color_airfoils='0.1',
        color_pitch_axis='r', color_LE='0.1', color_TE='0.1',
        color_SW=(37.0/256.0,197.0/256.0,85.0/256.0), stn_nums=False,
        export=True):
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
        color_airfoils : str or RGB tuple, line color of airfoils
            (default: dark grey)
        color_pitch_axis : str or RGB tuple, line color of pitch axis
            (default: red)
        color_LE : str or RGB tuple, line color of leading edge
            (default: dark grey)
        color_TE : str or RGB tuple, line color of trailing edge
            (default: dark grey)
        color_SW : str or RGB tuple, line color of shear webs
            (default: light green)
        stn_nums : boolean, plot/don't plot station numbers
        export : bool, do/don't write airfoil coords to a text file in the
            station path, which can be imported into SolidWorks as an XYZ curve

        """
        self.create_plot()
        if airfoils:
            self.plot_all_airfoils(lw=line_width, color=color_airfoils,
                twist_flag=twist, export_flag=export)
        if pitch_axis:
            self.plot_pitch_axis(lw=line_width, color=color_pitch_axis)
        if LE:
            self.plot_LE(lw=line_width, color=color_LE, twist_flag=twist)
        if TE:
            self.plot_TE(lw=line_width, color=color_TE, twist_flag=twist)
        if SW:
            self.plot_all_SW_cross_sections(lw=line_width, color=color_SW, 
                twist_flag=twist)
            self.plot_all_SW_spanwise_lines(lw=line_width, color=color_SW,
                twist_flag=twist, lower=True, upper=True, sw_list=[1,2,3])
        if stn_nums:
            self.plot_station_nums()
        self.show_plot()
            