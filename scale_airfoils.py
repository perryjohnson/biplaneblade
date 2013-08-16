"""Create new airfoil profiles by scaling existing profiles to new thicknesses.

Usage
-----
use insert_station.py to find the thickness-to-chord ratio at a new station
update the parameters in the 'SET THESE PARAMETERS' section below
from an IPython terminal, run this script:
    |> %run scale_airfoils

Author: Perry Roth-Johnson
Last updated: August 16, 2013

"""


import scripts.scale_airfoil_thickness as sat


airfoil_list = []

### SET THESE PARAMETERS ------------------------------------------------------
# read and write airfoil coordinates to this directory
airfoil_dir = 'biplane_flap-sym_no-stagger/airfoils'
# create entries for each new airfoil that we need to make
airfoil_list.append({
    'original_path' : 'DU91-W2-250_26.txt',
    'original_tc_ratio' : 0.26,
    'new_path' : 'DU91-W2-250_24pt3.txt',
    'new_tc_ratio' : 0.243
    })
### ---------------------------------------------------------------------------

for airfoil in airfoil_list:
    sat.write_and_plot(
        airfoil_dir,
        airfoil['original_path'],
        airfoil['original_tc_ratio'],
        airfoil['new_path'],
        airfoil['new_tc_ratio']
        )
