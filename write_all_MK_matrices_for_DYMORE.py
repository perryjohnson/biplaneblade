"""Write all mass and stiffness matrices for a DYMORE input file.

Usage
-----
start an IPython (qt)console with the pylab flag:
$ ipython qtconsole --pylab
or
$ ipython --pylab
Then, from the prompt, run this script:
|> %run write_all_MK_matrices_for_DYMORE
Finally, open the 'MK_matrices.txt' file, and manually copy the contents into
your DYMORE input file.

Author: Perry Roth-Johnson
Last updated: April 9, 2014

"""


import lib.blade as bl
import lib.dymore_utils as du

# load the Sandia blade
m = bl.MonoplaneBlade('Sandia blade SNL100-00', 'sandia_blade',
    rotate_airfoil_coords=False)

# open a new text file
MKfile = 'MK_matrices.txt'
f = open(MKfile, 'w')

print " Writing [M] and [K] matrices for DYMORE for:"
for station in m.list_of_stations:
    print "   Station #{0}...".format(station.station_num)
    # get the VABS output filename
    vabsMK = 'sandia_blade/stn{0:02d}/mesh_stn{0:02d}.vabs.K'.format(
        station.station_num)
    # write the mass and stiffness matrices in the DYMORE input file format
    du.writeMKmatrices(f, vabsMK,
        {'eta': station.coords.x1/m.list_of_stations[-1].coords.x1},
        CoordType='ETA_COORDINATE', debug_flag=False)

# close the text file
f.close()
print " See '{0}' for the DYMORE-formatted [M] and [K] matrices.".format(MKfile)
