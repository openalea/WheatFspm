# -*- coding: latin-1 -*-

"""
Construct a model equivalent to Distribution_CN_postflo_V2.mbk and run it.
The outputs are saved in 'cnwheat_output.csv'

The package cnwheat must be correctly installed before running this example.

CSV files must contain only ASCII characters and ',' as separator.
"""

"""
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
"""

import os
import time, datetime
import pandas as pd

import logging
import logging.config
import json


t0 = time.time()

from cnwheat import simulation
from cnwheat import model

DATA_DIRPATH = 'data'

CNWHEAT_OUTPUT_FILENAME = 'cnwheat_output.csv'

LOGGING_CONFIG_FILEPATH = os.path.join('..', 'logging.json')

LOGGING_LEVEL = logging.INFO # change to logging.DEBUG for debugging

def setup_logging(config_filepath='logging.json', level=logging.INFO):
    """Setup logging configuration
    """
    if os.path.exists(config_filepath):
        with open(config_filepath, 'r') as f:
            config = json.load(f)
        logging.config.dictConfig(config)
    else:
        logging.basicConfig()
    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    logging.getLogger('cnwheat.compartments').disabled = True
    logging.getLogger('cnwheat.derivatives').disabled = True

setup_logging(LOGGING_CONFIG_FILEPATH, LOGGING_LEVEL)

def read_t_data(curr_data_dirpath, data_filename):
    data_filepath = os.path.join(curr_data_dirpath, data_filename)
    return pd.read_csv(data_filepath, sep=None, index_col='t', engine = 'python')

organs = []
# create the chaff
name='chaff'
diameter = 0.02 # Approximated diamter for transpiration computations
PAR_df = read_t_data(DATA_DIRPATH, 'PAR_%s.csv' % name)
chaff = model.Chaff(area=0.00075, mstruct=0.21, width=diameter, height= 0.7, PAR=PAR_df.PAR,
                    starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0, amino_acids=0, proteins=1800, name=name)
organs.append(chaff)

# create the internodes
number_of_internodes = 3
areas = [(0.0012, 0.0003), 0.0004, 0.00025]
mstructs = [(0.148, 0.04), 0.18, 0.154]
diameters = [(0.042, 0.042), 0.043, 0.04] # Diameters of internodes, approximated from Ljutovac thesis. Diameters used as widths for transpiration computations
heights = [(0.4, 0.3), 0.18, 0.8]
prot = [(605, 500), 140, 0]
internodes = []
for i in xrange(number_of_internodes):
    internode_index = i + 1
    name = 'internode%d' % internode_index
    if internode_index == 1:
        current_names = (name + '_enclosed', name + '_exposed')
        current_areas = areas[i]
        current_mstructs = mstructs[i]
        current_diameters = diameters[i]
        current_heights = heights[i]
        current_proteins = prot[i]
    else:
        current_names = (name,)
        current_areas = (areas[i],)
        current_mstructs = (mstructs[i],)
        current_diameters = (diameters[i],)
        current_heights = (heights[i],)
        current_proteins = (prot[i],)

    for j in xrange(len(current_names)):
        PAR_df = read_t_data(DATA_DIRPATH, 'PAR_%s.csv' % current_names[j])
        internode = model.Internode(area=current_areas[j], mstruct=current_mstructs[j], width=current_diameters[j], height=current_heights[j], PAR=PAR_df.PAR,
                                          starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0, amino_acids=0, proteins=current_proteins[j], name=current_names[j])
        internodes.append(internode)
        organs.append(internode)

# create the laminae
number_of_laminae = 3
areas = [0.00346, 0.0034, 0.00228]
mstructs = [0.14, 0.09, 0.05]
widths = [0.018, 0.014, 0.0125]  # Widths of laminae, approximated from Ljutovac thesis. Widths used for transpiration computations
heights = [0.6, 0.38, 0.24]
prot = [2600, 1420, 590]
laminae = []
for i in xrange(number_of_laminae):
    name = 'lamina%d' % (i + 1)
    PAR_df = read_t_data(DATA_DIRPATH, 'PAR_%s.csv' % name)
    lamina = model.Lamina(area=areas[i], mstruct=mstructs[i], width= widths[i], height=heights[i], PAR=PAR_df.PAR,
                                starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0, amino_acids=0, proteins=prot[i], name=name)
    laminae.append(lamina)
    organs.append(lamina)

# create the peduncles
peduncles = []
names = ['peduncle_enclosed', 'peduncle_exposed']
areas = [0.00155, 0.00085]
mstructs = [0.168, 0.089]
diameters = [0.031, 0.031] # Diameters of peduncles, approximated from Ljutovac thesis. Diameters used as widths for transpiration computations
heights = [0.65, 0.5]
prot = [1240, 1000]
for i in xrange(len(names)):
    PAR_df = read_t_data(DATA_DIRPATH, 'PAR_%s.csv' % names[i])
    peduncle = model.Peduncle(area=areas[i], mstruct=mstructs[i], width= diameters[i], height=heights[i], PAR=PAR_df.PAR,
                                    starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0, amino_acids=0, proteins=prot[i], name=names[i])
    peduncles.append(peduncle)
    organs.append(peduncle)

# create the sheaths
number_of_sheaths = 3
sheaths = []
areas = [0.0006, 0.0005, 0.0004]
mstructs = [0.103, 0.069, 0.043]
diameters = [0.042, 0.043, 0.04] # Same as diameters of internodes, approximated from Ljutovac thesis. Diameters used as widths for transpiration computations
heights =[0.5, 0.3, 0.18]
prot = [325, 90, 10]
for i in xrange(number_of_sheaths):
    name = 'sheath%d' % (i + 1)
    PAR_df = read_t_data(DATA_DIRPATH, 'PAR_%s.csv' % name)
    sheath = model.Sheath(area=areas[i], mstruct=mstructs[i], width=diameters[i], height=heights[i], PAR=PAR_df.PAR,
                                starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0 , amino_acids=0, proteins=prot[i], name=name)
    sheaths.append(sheath)
    organs.append(sheath)

# create the grains
grains = model.Grains(starch=0, structure=10850, proteins=1180, name='grains')
organs.append(grains)

# create the roots
roots = model.Roots(mstruct=0.504, sucrose=0, nitrates=0, amino_acids=0, name='roots')
organs.append(roots)

# create the phloem
phloem = model.Phloem(sucrose=0, amino_acids=0, name='phloem')
organs.append(phloem)

# get meteo data
meteo_df = read_t_data(DATA_DIRPATH, 'meteo.csv')

# initialize the simulator
cnwheat_ = simulation.CNWheat(organs=organs, meteo=meteo_df)

# run the model
cnwheat_output_df = cnwheat_.run(start_time=0, stop_time=500, number_of_output_steps=125,
                                photosynthesis_computation_interval=2, show_progressbar=True)
##cnwheat_output_df = cnwheat_.run(start_time=0, stop_time=960, number_of_output_steps=241,
##                                photosynthesis_computation_interval=2, show_progressbar=True)

cnwheat_output_df.to_csv(CNWHEAT_OUTPUT_FILENAME, na_rep='NA', index=False)
execution_time = int(time.time()-t0)
print '\n', 'Model executed in ', str(datetime.timedelta(seconds=execution_time))


##POST-PROCESSING##
import imp
plot_columns = imp.load_source('plot_columns', r'C:\Users\rbarillot\Documents\PostDoc_Grignon\Modeles\Distribution_CN\CN-Wheat_Python\trunk\tools\plot_columns.py')

cnwheat_output_df = pd.read_csv(r'C:\Users\rbarillot\Documents\PostDoc_Grignon\Modeles\Distribution_CN\CN-Wheat_Python\trunk\test\system\cnwheat_output.csv')

# Computes Surfacic rate of photosynthesis
path_graphs = r'C:\Users\rbarillot\Documents\PostDoc_Grignon\Modeles\Distribution_CN\CN-Wheat_Python\trunk\test\system\Graphs'

graph_variables = {'An_': u'Net photosynthesis (µmol m$^{-2}$ s$^{-1}$)', 'Conc_TriosesP_': u'[TriosesP] (µmol g$^{-1}$ mstruct)', 'Conc_Starch_':u'[Starch] (µmol g$^{-1}$ mstruct)',
                   'Conc_Sucrose_':u'[Sucrose] (µmol g$^{-1}$ mstruct)', 'Conc_Fructan_':u'[Fructan] (µmol g$^{-1}$ mstruct)', 'Loading_Sucrose_': u'Loading Sucrose (µmol C sucrose g$^{-1}$ mstruct over delta_t)', 'Loading_Amino_Acids_': u'Loading Amino acids (µmol N amino acids g$^{-1}$ mstruct over delta_t)',
                   'Conc_Nitrates_': u'[Nitrates] (µmol g$^{-1}$ mstruct)', 'Conc_Amino_Acids_': u'[Amino_Acids] (µmol g$^{-1}$ mstruct)', 'Conc_Proteins_': u'[Proteins] (g g$^{-1}$ mstruct)', 'S_Proteins_': u'[Rate of protein synthesis] (µmol N g$^{-1}$ mstruct h$^{-1}$)',
                   'Nitrates_import_': u'Total nitrates imported (µmol h$^{-1}$)', 'Transpiration_':u'Organ transpiration (mm H$_{2}$0 h$^{-1}$)', 'Amino_Acids_import_': u'Total amino acids imported (µmol N over delta_t)',
                   'Conc_Amino_Acids_phloem':u'[Amino Acids phloem] (µmol g$^{-1}$ mstruct)', 'Conc_Sucrose_phloem':u'[Sucrose phloem] (µmol g$^{-1}$ mstruct)',
                   'Conc_Nitrates_roots': u'[Nitrates] (µmol g$^{-1}$ mstruct)', 'Conc_Amino_Acids_roots': u'[Amino_Acids] (µmol g$^{-1}$ mstruct)',
                   'S_Amino_Acids_': u'[Rate of amino acids synthesis] (µmol N g$^{-1}$ mstruct h$^{-1}$)',
                   'Uptake_Nitrates_roots':u'Nitrates uptake (µmol over delta_t)','Conc_Sucrose_roots': u'[Sucrose] (µmol g$^{-1}$ mstruct)','Conc_Nitrates_Soil_roots':u'[Nitrates] (µmol m$^{-3}$)',
                   'Unloading_Sucrose_roots': u'Unloaded_Sucrose (µmol C sucrose g$^{-1}$ mstruct h$^{-1}$)', 'Unloading_Amino_Acids_roots':u'Unloaded Amino Acids (µmol N AA g$^{-1}$ mstruct h$^{-1}$)',
                   'S_Amino_Acids_roots': u'[Rate of amino acids synthesis] (µmol N g$^{-1}$ mstruct h$^{-1}$)',
                   'Dry_Mass_grains':'Dry mass grains (g)', 'Proteins_N_Mass_grains': u'[N Proteins] (g)', 'Unloading_Sucrose_grains': u'Unloading_Sucrose (µmol sucrose over delta_t)'}

for var in graph_variables.keys():
    kwargs = {}
    if var in ('Dry_Mass_grains', 'Unloading_Amino_Acids_roots', 'Conc_Sucrose_phloem', 'Conc_Amino_Acids_phloem', 'Conc_Nitrates_roots', 'Conc_Amino_Acids_roots', 'Proteins_N_Mass_grains', 'Uptake_Nitrates_roots',
    'Conc_Sucrose_roots','Conc_Nitrates_Soil_roots', 'Unloading_Sucrose_roots', 'S_Amino_Acids_roots', 'Unloading_Sucrose_grains'):
##        if var=='Uptake_Nitrates_roots':
##            kwargs[var]={'label':'Actual NO3- uptake'}
##            plot_columns.plot_dataframe(cnwheat_output_df, y_label= graph_variables[var], column_to_matplotlib_kwargs = kwargs, plot_filepath = path_graphs + '\\' + var + '.PNG')
##            kwargs['Potential_Uptake_Nitrates_roots']={'label':'Potential NO3- uptake'}
##            plot_columns.plot_dataframe(cnwheat_output_df, y_label= graph_variables[var], column_to_matplotlib_kwargs = kwargs, plot_filepath = path_graphs + '\\' + var + '.PNG')
        kwargs={var:{}}
        plot_columns.plot_dataframe(cnwheat_output_df, y_label= graph_variables[var], column_to_matplotlib_kwargs = kwargs, plot_filepath = path_graphs + '\\' + var + '.PNG')

    else:
        kwargs_lamina, kwargs_sheath, kwargs_internode, kwargs_peduncle_chaff ={},{},{},{}
        for lam in laminae:
            header = var + lam.name
            kwargs_lamina[header]={'label':lam.name}
            plot_columns.plot_dataframe(cnwheat_output_df, y_label= graph_variables[var], column_to_matplotlib_kwargs = kwargs_lamina, plot_filepath = path_graphs + '\\' + var + 'lamina.PNG')
        for sh in sheaths:
            header = var + sh.name
            kwargs_sheath[header]={'label':sh.name}
            plot_columns.plot_dataframe(cnwheat_output_df, y_label= graph_variables[var], column_to_matplotlib_kwargs = kwargs_sheath, plot_filepath = path_graphs + '\\' + var + 'sheath.PNG')
        for inte in internodes:
            header = var + inte.name
            kwargs_internode[header]={'label':inte.name}
            plot_columns.plot_dataframe(cnwheat_output_df, y_label= graph_variables[var], column_to_matplotlib_kwargs = kwargs_internode, plot_filepath = path_graphs + '\\' + var + 'internode.PNG')
        for ped in peduncles:
            header = var + ped.name
            kwargs_peduncle_chaff[header]={'label':ped.name}
        header = var + chaff.name
        kwargs_peduncle_chaff[header]={'label':chaff.name}
        plot_columns.plot_dataframe(cnwheat_output_df, y_label= graph_variables[var], column_to_matplotlib_kwargs = kwargs_peduncle_chaff, plot_filepath = path_graphs + '\\' + var + 'peduncle_chaff.PNG')