# -*- coding: latin-1 -*-

"""
Construct a model equivalent to Distribution_CN_postflo_V2.mbk and run it.
The outputs are saved in 'cnwheat_output.csv'

The package cnwheat must be correctly installed before running this example.

CSV files must contain only ASCII characters and ',' as separator.
"""

import os
import time
import pandas as pd

t0 = time.time()

from cnwheat import cnwheat
from cnwheat import organ

DATA_DIRPATH = 'data'

CNWHEAT_OUTPUT_FILENAME = 'cnwheat_output.csv'


def read_t_data(curr_data_dirpath, data_filename):
    data_filepath = os.path.join(curr_data_dirpath, data_filename)
    return pd.read_csv(data_filepath, index_col='t')


# create the chaff
name='chaff'
PAR_df = read_t_data(DATA_DIRPATH, 'PAR_%s.csv' % name)
chaff = organ.Chaff(Area=0.00075, Mstruct=0.21, PAR=PAR_df.PAR,
                    STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0, name=name)

# create the internodes
number_of_internodes = 3
Areas = [(0.0012, 0.0003), 0.0004, 0.00025]
Mstructs = [(0.148, 0.04), 0.18, 0.154]
internodes = []
for i in xrange(number_of_internodes):
    internode_index = i + 1
    name = 'internode%d' % internode_index
    if internode_index == 1:
        current_names = (name + '_enclosed', name + '_exposed')
        current_Areas = Areas[i]
        current_Mstructs = Mstructs[i]
    else:
        current_names = (name,)
        current_Areas = (Areas[i],)
        current_Mstructs = (Mstructs[i],)
    for j in xrange(len(current_names)):
        PAR_df = read_t_data(DATA_DIRPATH, 'PAR_%s.csv' % current_names[j])
        internodes.append(organ.Internode(Area=current_Areas[j], Mstruct=current_Mstructs[j],
                                          PAR=PAR_df.PAR, FRUCTAN_0=0, STORAGE_0=0, SUCROSE_0=0,
                                          TRIOSESP_0=0, name=current_names[j]))

# create the laminae
number_of_laminae = 3
Areas = [0.00346, 0.0034, 0.00228]
Mstructs = [0.14, 0.09, 0.05]
laminae = []
for i in xrange(number_of_laminae):
    name = 'lamina%d' % (i + 1)
    PAR_df = read_t_data(DATA_DIRPATH, 'PAR_%s.csv' % name)
    laminae.append(organ.Lamina(Area=Areas[i], Mstruct=Mstructs[i],
                                PAR=PAR_df.PAR,
                                STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0,
                                name=name))

# create the peduncles
peduncles = []
names = ['peduncle_enclosed', 'peduncle_exposed']
Areas = [0.00155, 0.00085]
Mstructs = [0.168, 0.089]
for i in xrange(len(names)):
    PAR_df = read_t_data(DATA_DIRPATH, 'PAR_%s.csv' % names[i])
    peduncles.append(organ.Peduncle(Area=Areas[i], Mstruct=Mstructs[i], PAR=PAR_df.PAR,
                                    FRUCTAN_0=0, STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0,
                                    name=names[i]))

# create the sheaths
number_of_sheaths = 3
sheaths = []
Areas = [0.0006, 0.0005, 0.0004]
Mstructs = [0.103, 0.069, 0.043]
for i in xrange(number_of_sheaths):
    name = 'sheath%d' % (i + 1)
    PAR_df = read_t_data(DATA_DIRPATH, 'PAR_%s.csv' % name)
    sheaths.append(organ.Sheath(Area=Areas[i], Mstruct=Mstructs[i], PAR=PAR_df.PAR,
                                FRUCTAN_0=0, STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0, name=name))

# create the grains
grains = organ.Grains(STORAGE_0=0, STRUCTURE_0=10850, name='grains')

# create the roots
roots = organ.Roots(Mstruct=0.504, Sucrose_0=0, name='roots')

# create the phloem
phloem = organ.Phloem(SUCROSE_0=0, name='phloem')

# get meteo data
meteo_df = read_t_data(DATA_DIRPATH, 'meteo.csv')

# run the model
cnwheat_output_df = cnwheat.run(start_time=0, stop_time=960, number_of_output_steps=241,
                                phloem=phloem,
                                organs=[chaff]+internodes+laminae+peduncles+sheaths+[grains]+[roots],
                                meteo=meteo_df,
                                photosynthesis_computation_interval=2)

cnwheat_output_df.to_csv(CNWHEAT_OUTPUT_FILENAME, na_rep='NA', index=False)


print 'Model executed in ', time.time()-t0

# Graphs output
cnwheat_output_df = pd.read_csv(r'C:\Users\rbarillot\Documents\PostDoc_Grignon\Modeles\Distribution_CN\CN-Wheat_Python\trunk\example\complex_system\cnwheat_output.csv')

Areas={'lamina1':0.00346, 'lamina2':0.0034, 'lamina3':0.00228, 'internode1_enclosed':0.0012, 'internode1_exposed':0.0003, 'internode2': 0.0004, 'internode3': 0.00025,
        'sheath1':0.0006, 'sheath2':0.0005, 'sheath3':0.0004, 'peduncle_enclosed':0.00155, 'peduncle_exposed': 0.00085, 'chaff':0.00075}

kwargs_lamina, kwargs_sheath, kwargs_internode, kwargs_peduncle_chaff ={},{},{},{}

for col in cnwheat_output_df.columns:
    data = col.split('_')
    if data[0]=='Photosynthesis':
        org = data[1]
        if len(data)>2:
            loc = data[2]
            org_pos = (org + '_' + loc)
            Area = Areas[org_pos]
        else:
            Area = Areas[(org)]
            org_pos = org
        new_header = (col + '_Surfacic_rate' )
        cnwheat_output_df[new_header] = cnwheat_output_df[col] / (3600 * Area)
        if org[:-1] =='lamina':
            kwargs_lamina[new_header]={'label': org_pos}
        elif org[:-1] =='sheath':
            kwargs_sheath[new_header]={'label': org_pos}
        elif org[:-1] =='internode':
            kwargs_internode[new_header]={'label': org_pos}
        elif org =='peduncle' or org=='chaff':
            kwargs_peduncle_chaff[new_header]={'label': org_pos}

import imp
y_label_graphs = u'Net photosynthesis (µmol m$^{-2}$ s$^{-1}$)'
plot_columns = imp.load_source('plot_columns', r'C:\Users\rbarillot\Documents\PostDoc_Grignon\Modeles\Distribution_CN\CN-Wheat_Python\trunk\tools\plot_columns.py')
plot_columns.plot_dataframe(cnwheat_output_df, y_label= y_label_graphs, column_to_matplotlib_kwargs = kwargs_lamina, plot_filepath = r'C:\Users\rbarillot\Documents\PostDoc_Grignon\Modeles\Distribution_CN\CN-Wheat_Python\trunk\example\complex_system\Graphs\Photosynthesis_lamina.PNG')
plot_columns.plot_dataframe(cnwheat_output_df, y_label= y_label_graphs, column_to_matplotlib_kwargs = kwargs_sheath, plot_filepath = r'C:\Users\rbarillot\Documents\PostDoc_Grignon\Modeles\Distribution_CN\CN-Wheat_Python\trunk\example\complex_system\Graphs\Photosynthesis_sheath.PNG')
plot_columns.plot_dataframe(cnwheat_output_df, y_label= y_label_graphs, column_to_matplotlib_kwargs = kwargs_internode, plot_filepath = r'C:\Users\rbarillot\Documents\PostDoc_Grignon\Modeles\Distribution_CN\CN-Wheat_Python\trunk\example\complex_system\Graphs\Photosynthesis_internode.PNG')
plot_columns.plot_dataframe(cnwheat_output_df, y_label= y_label_graphs, column_to_matplotlib_kwargs = kwargs_peduncle_chaff, plot_filepath = r'C:\Users\rbarillot\Documents\PostDoc_Grignon\Modeles\Distribution_CN\CN-Wheat_Python\trunk\example\complex_system\Graphs\Photosynthesis_peduncle_chaff.PNG')
