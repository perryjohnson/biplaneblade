import numpy as np
import matplotlib.pyplot as plt

class Layup:
    """
    The Layup class defines a monoplane or biplane blade layup file.

    Initialization
    --------------
    Layup(filename, biplane_flag)
        filename : <string>
            The path and filename of the layup file.

        biplane_flag : <logical>
            True for biplane, False for monoplane.

    
    Public Variables
    ----------------
    filename : <string>
        The path and filename of the layup file.

    biplane_flag : <logical>
        True for biplane, False for monoplane.

    dict : <dict>
        A dictionary linking the data names and column numbers.

    data : <np array>
        Raw data loaded from the layup file.

    num_stations : <int>
        Number of stations in layup file.

    Public Methods
    --------------
    extractStationData(stn, print_flag=False)
        Parameters
            stn : <int>
                Station number.

            print_flag : <logical>
                True to print output.
        
        Returns
            stationData : <dict>
                Dictionary containing data and headings for specified spar station.
    """
    
    def __init__(self, filename, biplane_flag):
        """ Initializes the Layup class and loads the data into array """
        self.filename = filename
        self.biplane_flag = biplane_flag

        # load dictionary
        self.dict = self.loadDict()

        # read layup file into array
        self.data = self.readLayupFile(filename)

        # calculate number of stations in layup file
        self.num_stations = self.data.shape[0]


    def readLayupFile(self, fname):
        data = np.loadtxt(fname)  # read layup file and store all data in one array
        return data

    def loadDict(self):
        """
        Load a dictionary for columns in the 2D array of data from the layup file

        Parameters
        ----------
        None

        Returns
        -------
        dataDict : <dictionary>
            A dictionary of keywords that correspond to different columns of the 
            data in the layup file.
        """

        dataDict = {'spar station'                :  0,  # station numbers for the spar (begins at 1)
                    'x1'                          :  1,  # distances from root (in meters)
                    'eta'                         :  2,  # distances from root (dimensionless, normalized by total spar length)
                    'spar fraction'               :  3,  # percentage distances from spar root (dimensionless, normalized by total spar length and multiplied by 100)
                    'spar cap base'               :  4,  # spar cap base lengths (in meters)
                    'spar cap height'             :  5,  # spar cap height lengths (in meters)
                    'root buildup base'           :  6,  # root buildup base lengths (in meters)
                    'root buildup height'         :  7,  # root buildup height lengths (in meters)
                    'shear web base'              :  8,  # shear web base lengths (in meters)
                    'shear web height'            :  9,  # shear web height lengths (in meters)
                    'shear web foam base'         : 10,  # shear web foam base lengths (in meters)
                    'shear web biaxial GFRP base' : 11,  # shear web biaxial GFRP base lengths (in meters)
                    'twist degrees'               : 12,  # twist angles (in degrees)
                    'twist radians'               : 13,  # twist angles (in radians)
                    'k1'                          : 14,  # twist rates (in radians/meter)
                    'blade station'               : 15,  # station numbers for the blade (begins at 7)
                    'blade fraction'              : 16 } # percentage distances from blade root (dimensionless, normalized by total blade length and multiplied by 100)

        if self.biplane_flag:
            dataDict['chord']                     = 17   # chord length of airfoil that goes over this spar cross-section (in meters)
            dataDict['gap-to-chord ratio']        = 18   # gap-to-chord ratio of this cross-section for this layup (dimensionless, normalized by chord length)
            dataDict['x3']                        = 19   # vertical distance from pitch axis (in meters)
            dataDict['k2']                        = 20   # curvature about the x2 (edgewise) axis (in radians/meter)

        return dataDict 


    def extractStationData(self, stn, print_flag=False):
        """
        Extract a specific row from the layup data

        Parameters
        -----------
        stn : <int>
            spanwise station number
        
        biplane_flag : <logical>
        
        print_flag : <logical}

        Returns
        -------                
        stationData : <np.array>
            array of the data for the desired station
        """

        row_data = self.data[stn-1,:]
        if int(row_data[0]) != stn:
            print "ERROR: wrong station data pulled!"
        else:
            if print_flag:
                print "correct station data pulled!  :)"

        stationData = { 'spar station'                : int(row_data[self.dict['spar station']]),
                        'x1'                          : row_data[self.dict['x1']],
                        'eta'                         : row_data[self.dict['eta']],
                        'spar fraction'               : row_data[self.dict['spar fraction']],
                        'spar cap base'               : row_data[self.dict['spar cap base']],
                        'spar cap height'             : row_data[self.dict['spar cap height']],
                        'root buildup base'           : row_data[self.dict['root buildup base']],
                        'root buildup height'         : row_data[self.dict['root buildup height']],
                        'shear web base'              : row_data[self.dict['shear web base']],
                        'shear web height'            : row_data[self.dict['shear web height']],
                        'shear web foam base'         : row_data[self.dict['shear web foam base']],
                        'shear web biaxial GFRP base' : row_data[self.dict['shear web biaxial GFRP base']],
                        'twist degrees'               : row_data[self.dict['twist degrees']],
                        'twist radians'               : row_data[self.dict['twist radians']],
                        'k1'                          : row_data[self.dict['k1']],
                        'blade station'               : int(row_data[self.dict['blade station']]),
                        'blade fraction'              : row_data[self.dict['blade fraction']] }
        if self.biplane_flag:
            stationData['chord']                      = row_data[self.dict['chord']]
            stationData['gap-to-chord ratio']         = row_data[self.dict['gap-to-chord ratio']]
            stationData['x3']                         = row_data[self.dict['x3']]

        return stationData