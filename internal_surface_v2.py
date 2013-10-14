"""A script to merge all the structural parts into one polygon at each station.

Author: Perry Roth-Johnson
Last updated: October 11, 2013

"""


import lib.blade as bl
reload(bl)
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

# --- sandia blade ------------------------------------------------------------
m = bl.MonoplaneBlade('Sandia blade SNL100-00', 'sandia_blade',
    rotate_airfoil_coords=False)

# pre-process the airfoil coordinates
for station in m.list_of_stations:
    # station.find_SW_cs_coords()
    station.airfoil.create_polygon()
    try:
        station.structure.create_all_layers()
    except Warning:
        print "  ...skipping Station #{0}************".format(station.station_num)
        print "length of list_of_polygons =", len(station.structure.list_of_polygons)
        pass
    # station.write_all_part_polygons()
    station.structure.calculate_area()
    # print "Station #{0}, area = {1:6.4f} m^2".format(station.station_num, station.structure.area)

# m.plot_selected_cross_sections(plot_edges=False, plot_parts=True,
    # selected_stations=[7,11,14,16,18,19,23,26,30,31,32,33])

# assemble a pandas DataFrame of percent areas for each structural part
list_of_dicts = []
for station in m.list_of_stations:
    d = station.structure.calculate_all_percent_areas()
    list_of_dicts.append(d)

pa = pd.DataFrame(list_of_dicts, index=range(1,m.number_of_stations+1))
pa = pa.fillna(0)
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
# save the data to a CSV file
pa_path = os.path.join(m.blade_path, 'percent_areas.csv')
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

# make a stacked bar plot of all the percent areas for each structural part
# ref: http://matplotlib.org/examples/pylab_examples/bar_stacked.html
plt.figure(figsize=(22,12))
ind = np.arange(m.number_of_stations)  # the x locations for each blade station
width = 0.35                           # the width of the bars
p1 = plt.bar(ind, pa['external surface (gelcoat)'], width, color='#4000FF')
p2 = plt.bar(ind, pa['external surface (triax)'], width, color='#4000AA', 
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
p7 = plt.bar(ind, pa['shear webs (biax)'], width, color='#AAF100', 
    bottom=(pa['external surface (gelcoat)'] + pa['external surface (triax)'] +
        pa['root buildup'] + pa['spar caps'] + pa['aft panels'] +
        pa['LE panel']))
p8 = plt.bar(ind, pa['shear webs (foam)'], width, color='#FFF100', 
    bottom=(pa['external surface (gelcoat)'] + pa['external surface (triax)'] +
        pa['root buildup'] + pa['spar caps'] + pa['aft panels'] +
        pa['LE panel'] + pa['shear webs (biax)']))
p9 = plt.bar(ind, pa['TE reinforcement (uniax)'], width, color='#F366BA', 
    bottom=(pa['external surface (gelcoat)'] + pa['external surface (triax)'] +
        pa['root buildup'] + pa['spar caps'] + pa['aft panels'] +
        pa['LE panel'] + pa['shear webs (biax)'] + pa['shear webs (foam)']))
p10 = plt.bar(ind, pa['TE reinforcement (foam)'], width, color='#F300BA', 
    bottom=(pa['external surface (gelcoat)'] + pa['external surface (triax)'] +
        pa['root buildup'] + pa['spar caps'] + pa['aft panels'] +
        pa['LE panel'] + pa['shear webs (biax)'] + pa['shear webs (foam)'] +
        pa['TE reinforcement (uniax)']))
p11 = plt.bar(ind, pa['internal surfaces (triax)'], width, color='#999999', 
    bottom=(pa['external surface (gelcoat)'] + pa['external surface (triax)'] +
        pa['root buildup'] + pa['spar caps'] + pa['aft panels'] +
        pa['LE panel'] + pa['shear webs (biax)'] + pa['shear webs (foam)'] +
        pa['TE reinforcement (uniax)'] + pa['TE reinforcement (foam)']))
p12 = plt.bar(ind, pa['internal surfaces (resin)'], width, color='#999944', 
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
plt.xlim([-width/2.0, m.number_of_stations-0.5+width/2.0])
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
    'LE panel',
    'aft panels',
    'spar caps',
    'root buildup',
    'external surface (triax)',
    'external surface (gelcoat)'
    ), 
    bbox_to_anchor=(1.02, 0.5),
    loc='center left',
    borderaxespad=0.0)

plt.subplots_adjust(left=0.05, bottom=0.05, right=0.82, top=0.95)
plt.grid(axis='y')
plt.show()
