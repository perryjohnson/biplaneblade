"""Delete leftover station paths in ./sandia_blade and in ./biplane_blade"""


import os
import shutil


# --- sandia blade ------------------------------------------------------------
# generate list of station paths
list_of_station_paths = []
root_path = r".\sandia_blade"
for station in range(1,34+1):
    d = 'stn{0:02d}'.format(station)
    p = os.path.join(root_path, d)
    list_of_station_paths.append(p)

# delete each station path
for station_path in list_of_station_paths:
    if os.path.exists(station_path):
        # delete the station path, even if it has contents
        shutil.rmtree(station_path)
        print " [Deleted station path] {0}".format(station_path)

# --- biplane blade -----------------------------------------------------------
# generate list of station paths
list_of_station_paths = []
root_path = r".\biplane_blade"
for station in range(1,34+1):
    d = 'stn{0:02d}'.format(station)
    p = os.path.join(root_path, d)
    list_of_station_paths.append(p)

# delete each station path
for station_path in list_of_station_paths:
    if os.path.exists(station_path):
        # delete the station path, even if it has contents
        shutil.rmtree(station_path)
        print " [Deleted station path] {0}".format(station_path)
