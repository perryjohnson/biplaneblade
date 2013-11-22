from scipy import interpolate
import numpy as np

def write_curve_with_new_points(old_curve_file, new_curve_file, xnew):
    """Writes the coordinates of a new curve with interpolated points.

    Parameters:
    -----------
    old_curve_file : str, text file with original curve coordinates (to read)
    new_curve_file : str, text file with new curve coordinates (to write)
    xnew : list of floats, the x-coords of new interpolated points to be added

    """
    nc = add_points(old_curve_file, xnew)
    write_curve_file(new_curve_file, nc)
    print "Read original curve coordinates from:", old_curve_file
    print "Added coordinates at new x-coordinates..."
    for point in xnew:
        print " ", point
    print "Wrote new curve coordinates to:", new_curve_file

def add_points(curve_file, xnew):
    """Returns a new curve with interpolated points added at new x-coordinates.

    Parameters:
    -----------
    curve_file : str, the name of the text file with the original curve points
    xnew : list of floats, the x-coords of new interpolated points to be added

    Usage:
    ------
    add_points('curve110.txt', [-0.750, 0.750, 0.753, 0.833, 0.836])

    """
    # Load a curve text file into an array; only keep the xy-coordinates
    curve = np.loadtxt(curve_file)[:,:2]
    # check if the curve is ordered from xmax-->xmin
    flip_flag = False
    if curve[0,0] > curve[-1,0]:
        # curve is ordered from xmax-->xmin, so flip it to be xmin-->xmax
        curve = curve[::-1,:]
        flip_flag = True
    x = curve[:,0]
    y = curve[:,1]
    f = interpolate.interp1d(x,y)
    xnew = np.array(xnew)
    # get y-coordinates for new interpolated points (at each element of xnew)
    ynew = f(xnew)
    temp = np.array([xnew, ynew])
    temp = temp.T
    new_curve = np.vstack((curve,temp))
    # sort the coordinates from xmin-->xmax
    #   note: new_curve.sort(axis=0) doesn't work!
    new_curve = new_curve[new_curve[:,0].argsort()]
    if flip_flag:
        # if the original curve was ordered from xmax-->xmin,
        # flip the new curve, so it will also be ordered from xmax-->xmin
        new_curve = new_curve[::-1,:]
    return new_curve

def write_curve_file(new_curve_file, new_curve):
    """Write the coordinates of new_curve to a text file: new_curve_file.

    Parameters:
    -----------
    new_curve_file : str, the name of the text file that will be written
    new_curve : np.array, coordinates that will be written to new_curve_file

    Usage:
    ------
    nc = add_points('curve110.txt', [-0.750, 0.750, 0.753, 0.833, 0.836])
    write_curve_file('curve110_new.txt', nc)

    """
    f = open(new_curve_file, 'w')
    # f.write('curd 110 lp3\n')
    for cd_pair in new_curve:
        f.write('{0: .8f}  {1: .8f}  0.0\n'.format(cd_pair[0], cd_pair[1]))
    f.write(';;\n\n')
    f.close()
