"""Copy the monoplane stations shared w/ the Sandia blade to the biplane blade.

This script must be run from the project root directory (spardesign2), not this
file's directory (biplane_blade_lib).

Author: Perry Roth-Johnson
Last modified: April 23, 2014

"""

import os
from shutil import copy

def copy_vabs_input_file(Sandia_blade_station, biplane_blade_station):
    fmt_s = 'sandia_blade/stn{0:02d}/mesh_stn{0:02d}.vabs'.format(
        Sandia_blade_station)
    fmt_b = 'biplane_blade/stn{0:02d}/mesh_stn{0:02d}.vabs'.format(
        biplane_blade_station)
    copy(fmt_s, fmt_b)
    p = " Copied VABS input file: Sandia stn #{0:02d} --> biplane stn #{1:02d}"
    pstr = p.format(Sandia_blade_station, biplane_blade_station)
    print pstr

# make sure we're in the project root directory
f = os.getcwd().split(os.sep)[-1]
if f != 'spardesign2':
    raise Warning("Current working directory is not `spardesign2`!\n  Change directories and try running this script again.")

# if we're in the correct directory, proceed.
copy_vabs_input_file(Sandia_blade_station=1, biplane_blade_station=1)
copy_vabs_input_file(Sandia_blade_station=2, biplane_blade_station=2)
copy_vabs_input_file(Sandia_blade_station=3, biplane_blade_station=3)
copy_vabs_input_file(Sandia_blade_station=4, biplane_blade_station=4)
copy_vabs_input_file(Sandia_blade_station=5, biplane_blade_station=5)
copy_vabs_input_file(Sandia_blade_station=6, biplane_blade_station=6)
copy_vabs_input_file(Sandia_blade_station=7, biplane_blade_station=7)
copy_vabs_input_file(Sandia_blade_station=8, biplane_blade_station=8)
copy_vabs_input_file(Sandia_blade_station=9, biplane_blade_station=9)
copy_vabs_input_file(Sandia_blade_station=21, biplane_blade_station=25)
copy_vabs_input_file(Sandia_blade_station=22, biplane_blade_station=28)
copy_vabs_input_file(Sandia_blade_station=23, biplane_blade_station=29)
copy_vabs_input_file(Sandia_blade_station=24, biplane_blade_station=30)
copy_vabs_input_file(Sandia_blade_station=25, biplane_blade_station=31)
copy_vabs_input_file(Sandia_blade_station=26, biplane_blade_station=32)
copy_vabs_input_file(Sandia_blade_station=27, biplane_blade_station=33)
copy_vabs_input_file(Sandia_blade_station=28, biplane_blade_station=34)
copy_vabs_input_file(Sandia_blade_station=29, biplane_blade_station=35)
copy_vabs_input_file(Sandia_blade_station=30, biplane_blade_station=36)
copy_vabs_input_file(Sandia_blade_station=31, biplane_blade_station=37)
copy_vabs_input_file(Sandia_blade_station=32, biplane_blade_station=38)
copy_vabs_input_file(Sandia_blade_station=33, biplane_blade_station=39)
copy_vabs_input_file(Sandia_blade_station=34, biplane_blade_station=40)
