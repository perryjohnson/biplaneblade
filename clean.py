"""Delete leftover files and paths in each blade_path"""


import os
import shutil


def clean_blade_path(blade_path,
    except_list=['blade_definition.csv',
                 'layers.csv',
                 'materials.csv',
                 'airfoils',
                 'percent_areas.csv',
                 'percent_masses.csv',
                 'config_name.txt',
                 'blade_props_from_Sandia.csv',
                 'blade_props_from_VABS.csv']):
    """Removes all files except blade defintion and airfoil path."""
    clean_list = os.listdir(blade_path)
    for item in clean_list:
        if item not in except_list:
            item = os.path.join(blade_path,item)
            if os.path.exists(item):
                if os.path.isfile(item):
                    os.remove(item)
                    print " [Deleted file] {0}".format(item)
                else:
                    shutil.rmtree(item)
                    print " [Deleted path] {0}".format(item)

clean_blade_path('sandia_blade')
clean_blade_path('biplane_blade')
