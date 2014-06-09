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

def plot_mass_schedule(blade1, blade2, show_stn_nums=False, 
    blade1_stn_nums=[], blade2_stn_nums=[], blade1_label='', blade2_label='',
    print_flag=False):
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
    if show_stn_nums:
        plt.text(blade1._df['x1'].ix[blade1_stn_nums[0]-1], 
            (masses1[blade1_stn_nums[0]-1]+500), (blade1_label+' station #:'),
            va='bottom')
        for i in range(blade1_stn_nums[0]-1, blade1_stn_nums[1]):
            plt.text(blade1._df['x1'].ix[i+1], masses1[i]+100, str(i+1),
                va='bottom', ha='center')
        plt.text(blade2._df['x1'].ix[blade2_stn_nums[0]-1], 
            (masses2[blade2_stn_nums[0]-1]-1000), (blade2_label+' station #:'),
            va='top')
        for i in range(blade2_stn_nums[0]-1, blade2_stn_nums[1]):
            plt.text(blade2._df['x1'].ix[i+1], masses2[i]-100, str(i+1),
                va='top', ha='center')
    if print_flag:
        print "stn#", "x1", "mass1"
        for i in range(1,len(masses1)):
            print str(i+1), blade1._df['x1'].ix[i+1], masses1[i]
        print ""
        print "stn#", "x1", "mass2"
        for i in range(1,len(masses2)):
            print str(i+1), blade2._df['x1'].ix[i+1], masses2[i]
    plt.xlabel('span, x1 [m]')
    plt.ylabel('mass [kg/m]')
    plt.legend()
    plt.show()
