"""Plot properties of different blades on the same plot.

Author: Perry Roth-Johnson
Last updated: April 17, 2014

"""


import matplotlib.pyplot as plt


def plot_chord_schedule(list_of_blades):
    """Plot the chord vs. span."""
    plt.figure(figsize=(7,4))
    # plt.axes().set_aspect('equal')
    for blade in list_of_blades:
        plt.plot(blade._df['x1'], blade._df['chord'], blade.linestyle,
            label=blade.name, lw=blade.lw, mfc=blade.mfc, mec=blade.mec,
            mew=blade.mew, ms=blade.ms)
    plt.xlabel('span, x1 [m]')
    plt.ylabel('chord [m]')
    plt.legend(loc='lower left')
    plt.grid('on')
    plt.show()

def plot_twist_schedule(list_of_blades):
    """Plot the twist vs. span."""
    plt.figure(figsize=(7,4))
    # plt.axes().set_aspect('equal')
    for blade in list_of_blades:
        plt.plot(blade._df['x1'], blade._df['twist'], blade.linestyle,
            label=blade.name, lw=blade.lw, mfc=blade.mfc, mec=blade.mec,
            mew=blade.mew, ms=blade.ms)
    plt.xlabel('span, x1 [m]')
    plt.ylabel('twist [deg]')
    plt.legend(loc='lower left')
    plt.grid('on')
    plt.show()

def plot_thickness_schedule(list_of_blades):
    """Plot the thickness vs. span."""
    plt.figure(figsize=(7,4))
    # plt.axes().set_aspect('equal')
    for blade in list_of_blades:
        plt.plot(blade._df['x1'], blade._df['thickness-to-chord ratio']*blade._df['chord'], blade.linestyle,
            label=blade.name, lw=blade.lw, mfc=blade.mfc, mec=blade.mec,
            mew=blade.mew, ms=blade.ms)
    plt.xlabel('span, x1 [m]')
    plt.ylabel('thickness [m]')
    plt.legend(loc='upper right')
    plt.grid('on')
    plt.show()

def plot_gap_schedule(list_of_blades):
    """Plot the gap vs. span."""
    plt.figure(figsize=(7,4))
    # plt.axes().set_aspect('equal')
    for blade in list_of_blades:
        plt.plot(blade._df['x1'], blade._df['gap-to-chord ratio']*blade._df['chord'], blade.linestyle,
            label=blade.name, lw=blade.lw, mfc=blade.mfc, mec=blade.mec,
            mew=blade.mew, ms=blade.ms)
    plt.xlabel('span, x1 [m]')
    plt.ylabel('gap [m]')
    plt.xlim((0,100))
    plt.legend()
    plt.grid('on')
    plt.show()

def plot_thickness_to_chord_schedule(list_of_blades):
    """Plot the thickness-to-chord ratio vs. span."""
    plt.figure(figsize=(7,4))
    # plt.axes().set_aspect('equal')
    for blade in list_of_blades:
        plt.plot(blade._df['x1'], blade._df['thickness-to-chord ratio'], blade.linestyle,
            label=blade.name, lw=blade.lw, mfc=blade.mfc, mec=blade.mec,
            mew=blade.mew, ms=blade.ms)
    plt.xlabel('span, x1 [m]')
    plt.ylabel('thickness-to-chord ratio [-]')
    plt.legend(loc='upper right')
    plt.grid('on')
    plt.show()

def plot_gap_to_chord_schedule(list_of_blades):
    """Plot the gap-to-chord ratio vs. span."""
    plt.figure(figsize=(7,4))
    # plt.axes().set_aspect('equal')
    for blade in list_of_blades:
        plt.plot(blade._df['x1'], blade._df['gap-to-chord ratio'], blade.linestyle,
            label=blade.name, lw=blade.lw, mfc=blade.mfc, mec=blade.mec,
            mew=blade.mew, ms=blade.ms)
    plt.xlabel('span, x1 [m]')
    plt.ylabel('gap-to-chord ratio [-]')
    plt.xlim((0,100))
    plt.legend()
    plt.grid('on')
    plt.show()

def plot_mass_schedule(blade1,blade2):
    """Plot the mass vs. span."""
    blade1.calculate_all_masses()
    masses1 = []
    for station in blade1.list_of_stations:
        masses1.append(station.structure.mass)
    blade2.calculate_all_masses()
    masses2 = []
    for station in blade2.list_of_stations:
        masses2.append(station.structure.mass)
    plt.figure(figsize=(12,8))
    plt.plot(blade1._df['x1'],masses1,'bo-',label=blade1.name)
    plt.plot(blade2._df['x1'],masses2,'rx-',label=blade2.name)
    plt.xlabel('span, x1 [m]')
    plt.ylabel('mass [kg/m]')
    plt.legend()
    plt.show()
