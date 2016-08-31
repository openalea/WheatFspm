# -*- coding: latin-1 -*-

'''
    main
    ~~~~

    An example to show how to couple models CN-Wheat, Farquhar-Wheat and Senesc-Wheat using a static topology from Adel-Wheat.
    This example uses the format MTG to exchange data between the models.

    You must first install :mod:`alinea.adel`, :mod:`cnwheat`, :mod:`farquharwheat` and :mod:`senescwheat` (and add them to your PYTHONPATH)
    before running this script with the command `python`.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2015.
'''

'''
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
'''

import os
import time, datetime
import profile, pstats

import logging

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from alinea.adel.astk_interface import AdelWheat

from cnwheat import simulation as cnwheat_simulation, model as cnwheat_model, parameters as cnwheat_parameters, converter as cnwheat_converter, tools as cnwheat_tools, run_caribu
from farquharwheat import simulation as farquharwheat_simulation, model as farquharwheat_model, converter as farquharwheat_converter
from senescwheat import simulation as senescwheat_simulation, model as senescwheat_model, converter as senescwheat_converter
from growthwheat import simulation as growthwheat_simulation, model as growthwheat_model, converter as growthwheat_converter, interface as  growthwheat_interface

INPUTS_DIRPATH = 'inputs'
GRAPHS_DIRPATH = 'graphs'

# adelwheat inputs at t0
ADELWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'adelwheat') # the directory adelwheat must contain files 'adel0000.pckl' and 'scene0000.bgeom'

# cnwheat inputs at t0
CNWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'cnwheat')
CNWHEAT_PLANTS_INPUTS_FILEPATH = os.path.join(CNWHEAT_INPUTS_DIRPATH, 'plants_inputs.csv')
CNWHEAT_AXES_INPUTS_FILEPATH = os.path.join(CNWHEAT_INPUTS_DIRPATH, 'axes_inputs.csv')
CNWHEAT_METAMERS_INPUTS_FILEPATH = os.path.join(CNWHEAT_INPUTS_DIRPATH, 'metamers_inputs.csv')
CNWHEAT_ORGANS_INPUTS_FILEPATH = os.path.join(CNWHEAT_INPUTS_DIRPATH, 'organs_inputs.csv')
CNWHEAT_HGZ_INPUTS_FILEPATH = os.path.join(CNWHEAT_INPUTS_DIRPATH, 'hgzs_inputs.csv')
CNWHEAT_ELEMENTS_INPUTS_FILEPATH = os.path.join(CNWHEAT_INPUTS_DIRPATH, 'elements_inputs.csv')
CNWHEAT_SOILS_INPUTS_FILEPATH = os.path.join(CNWHEAT_INPUTS_DIRPATH, 'soils_inputs.csv')

# farquharwheat inputs at t0
FARQUHARWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'farquharwheat')
FARQUHARWHEAT_INPUTS_FILEPATH = os.path.join(FARQUHARWHEAT_INPUTS_DIRPATH, 'inputs.csv')
METEO_FILEPATH = os.path.join(FARQUHARWHEAT_INPUTS_DIRPATH, 'meteo_Clermont_rebuild.csv')

# senescwheat inputs at t0
SENESCWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'senescwheat')
SENESCWHEAT_ROOTS_INPUTS_FILEPATH = os.path.join(SENESCWHEAT_INPUTS_DIRPATH, 'roots_inputs.csv')
SENESCWHEAT_ELEMENTS_INPUTS_FILEPATH = os.path.join(SENESCWHEAT_INPUTS_DIRPATH, 'elements_inputs.csv')

# growthwheat inputs at t0
GROWTHWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'growthwheat')
GROWTHWHEAT_HGZ_INPUTS_FILEPATH = os.path.join(GROWTHWHEAT_INPUTS_DIRPATH, 'hgzs_inputs.csv')
GROWTHWHEAT_EXPOSED_ELEMENT_INPUTS_FILEPATH = os.path.join(GROWTHWHEAT_INPUTS_DIRPATH, 'organs_inputs.csv')

# the path of the CSV files where to save the states of the modeled system at each step
OUTPUTS_DIRPATH = 'outputs'
AXES_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'axes_states.csv')
ORGANS_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'organs_states.csv')
HGZS_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'hgzs_states.csv')
ELEMENTS_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'elements_states.csv')
SOILS_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'soils_states.csv')

INPUTS_OUTPUTS_PRECISION = 10

LOGGING_CONFIG_FILEPATH = os.path.join('..', '..', 'logging.json')

LOGGING_LEVEL = logging.INFO # can be one of: DEBUG, INFO, WARNING, ERROR, CRITICAL

cnwheat_tools.setup_logging(LOGGING_CONFIG_FILEPATH, LOGGING_LEVEL, log_model=False, log_compartments=False, log_derivatives=False)


meteo = pd.read_csv(METEO_FILEPATH, index_col='t')

current_time_of_the_system = time.time()

# define the time step in hours for each simulator
caribu_ts = 2
senescwheat_ts = 2
farquharwheat_ts = 2
growthwheat_ts = 1
cnwheat_ts = 1

hour_to_second_conversion_factor = 3600
# create the simulators
senescwheat_simulation_ = senescwheat_simulation.Simulation(senescwheat_ts * hour_to_second_conversion_factor)
farquharwheat_simulation_ = farquharwheat_simulation.Simulation()
cnwheat_simulation_ = cnwheat_simulation.Simulation(cnwheat_ts * hour_to_second_conversion_factor)

# read adelwheat inputs at t0
adel_wheat = AdelWheat(seed=1234)
g = adel_wheat.load(dir=ADELWHEAT_INPUTS_DIRPATH)[0]

# read cnwheat inputs at t0
cnwheat_plants_inputs_t0 = pd.read_csv(CNWHEAT_PLANTS_INPUTS_FILEPATH)
cnwheat_axes_inputs_t0 = pd.read_csv(CNWHEAT_AXES_INPUTS_FILEPATH)
cnwheat_metamers_inputs_t0 = pd.read_csv(CNWHEAT_METAMERS_INPUTS_FILEPATH)
cnwheat_organs_inputs_t0 = pd.read_csv(CNWHEAT_ORGANS_INPUTS_FILEPATH)
cnwheat_hgzs_inputs_t0 = pd.read_csv(CNWHEAT_HGZ_INPUTS_FILEPATH)
cnwheat_elements_inputs_t0 = pd.read_csv(CNWHEAT_ELEMENTS_INPUTS_FILEPATH)
cnwheat_soils_inputs_t0 = pd.read_csv(CNWHEAT_SOILS_INPUTS_FILEPATH)

# read farquharwheat inputs at t0
farquharwheat_elements_inputs_t0 = pd.read_csv(FARQUHARWHEAT_INPUTS_FILEPATH)

# read senescwheat inputs at t0
senescwheat_roots_inputs_t0 = pd.read_csv(SENESCWHEAT_ROOTS_INPUTS_FILEPATH)
senescwheat_elements_inputs_t0 = pd.read_csv(SENESCWHEAT_ELEMENTS_INPUTS_FILEPATH)

# read growthwheat inputs at t0
growthwheat_hgzs_inputs_t0 = pd.read_csv(GROWTHWHEAT_HGZ_INPUTS_FILEPATH)
growthwheat_exposed_element_inputs_t0 = pd.read_csv(GROWTHWHEAT_EXPOSED_ELEMENT_INPUTS_FILEPATH)

# define the start and the end of the whole simulation (in hours)
start_time = 0
stop_time = 4

# define lists of dataframes to store the state of the system at each step.
axes_all_data_list = []
organs_all_data_list = [] # organs which belong to axes: roots, phloem, grains
hgzs_all_data_list = []
elements_all_data_list = []
soils_all_data_list = []

# initialize dataframes to share data between the models
# organs
cnwheat_organs_inputs_t0_reindexed = pd.DataFrame(cnwheat_organs_inputs_t0.values,
                                                  index=sorted(cnwheat_organs_inputs_t0.groupby(cnwheat_simulation.Simulation.ORGANS_INPUTS_INDEXES).groups.keys()),
                                                  columns=cnwheat_organs_inputs_t0.columns)
senescwheat_roots_inputs_t0_with_organ_column = senescwheat_roots_inputs_t0.copy()
senescwheat_roots_inputs_t0_with_organ_column.loc[:, 'organ'] = 'roots'
senescwheat_roots_inputs_t0_with_organ_column_reindexed = pd.DataFrame(senescwheat_roots_inputs_t0_with_organ_column.values,
                                                                       index=sorted(senescwheat_roots_inputs_t0_with_organ_column.groupby(senescwheat_converter.ROOTS_TOPOLOGY_COLUMNS + ['organ']).groups.keys()),
                                                                       columns=senescwheat_roots_inputs_t0_with_organ_column.columns)
organs_inputs_t0 = cnwheat_organs_inputs_t0_reindexed.combine_first(senescwheat_roots_inputs_t0_with_organ_column_reindexed)
organs_states = organs_inputs_t0.reindex_axis(cnwheat_simulation.Simulation.ORGANS_INPUTS_INDEXES + sorted(set(cnwheat_simulation.Simulation.ORGANS_STATE  + senescwheat_converter.SENESCWHEAT_ROOTS_INPUTS_OUTPUTS)), axis=1)
# hidden growing zones
cnwheat_hgzs_inputs_t0_reindexed = pd.DataFrame(cnwheat_hgzs_inputs_t0.values,
                                                index=sorted(cnwheat_hgzs_inputs_t0.groupby(cnwheat_simulation.Simulation.HIDDENGROWINGZONE_INPUTS_INDEXES).groups.keys()),
                                                columns=cnwheat_hgzs_inputs_t0.columns)
growthwheat_hgzs_inputs_t0_reindexed = pd.DataFrame(growthwheat_hgzs_inputs_t0.values.tolist(),
                                                    index=sorted(growthwheat_hgzs_inputs_t0.groupby(growthwheat_converter.HGZ_TOPOLOGY_COLUMNS).groups.keys()),
                                                    columns=growthwheat_hgzs_inputs_t0.columns)
hgzs_inputs_t0 = cnwheat_hgzs_inputs_t0_reindexed.combine_first(growthwheat_hgzs_inputs_t0_reindexed)
dtypes = growthwheat_hgzs_inputs_t0_reindexed.dtypes.combine_first(cnwheat_hgzs_inputs_t0_reindexed.dtypes)
for k, v in dtypes.iteritems():
    hgzs_inputs_t0[k] = hgzs_inputs_t0[k].astype(v)

hgzs_states = hgzs_inputs_t0.reindex_axis(cnwheat_simulation.Simulation.HIDDENGROWINGZONE_INPUTS_INDEXES + sorted(set(cnwheat_simulation.Simulation.HIDDENGROWINGZONE_STATE + growthwheat_simulation.HGZ_INPUTS_OUTPUTS)), axis=1)
# elements
cnwheat_elements_inputs_t0_reindexed = pd.DataFrame(cnwheat_elements_inputs_t0.values,
                                                    index=sorted(cnwheat_elements_inputs_t0.groupby(cnwheat_simulation.Simulation.ELEMENTS_INPUTS_INDEXES).groups.keys()),
                                                    columns=cnwheat_elements_inputs_t0.columns)
farquharwheat_elements_inputs_t0_reindexed = pd.DataFrame(farquharwheat_elements_inputs_t0.values,
                                                          index=sorted(farquharwheat_elements_inputs_t0.groupby(farquharwheat_converter.DATAFRAME_TOPOLOGY_COLUMNS).groups.keys()),
                                                          columns=farquharwheat_elements_inputs_t0.columns)
growthwheat_elements_inputs_t0_with_element_column = growthwheat_exposed_element_inputs_t0.copy()
growthwheat_elements_inputs_t0_with_element_column.loc[growthwheat_elements_inputs_t0_with_element_column.organ == 'blade', 'element'] = 'LeafElement1'
growthwheat_elements_inputs_t0_with_element_column.loc[growthwheat_elements_inputs_t0_with_element_column.organ != 'blade', 'element'] = 'StemElement'
growthwheat_exposed_element_inputs_t0_with_element_column_reindexed = pd.DataFrame(growthwheat_elements_inputs_t0_with_element_column.values,
                                                                                   index=sorted(growthwheat_elements_inputs_t0_with_element_column.groupby(growthwheat_converter.ORGAN_TOPOLOGY_COLUMNS + ['element']).groups.keys()),
                                                                                   columns=growthwheat_elements_inputs_t0_with_element_column.columns)
senescwheat_elements_inputs_t0_reindexed = pd.DataFrame(senescwheat_elements_inputs_t0.values,
                                                        index=sorted(senescwheat_elements_inputs_t0.groupby(senescwheat_converter.ELEMENTS_TOPOLOGY_COLUMNS).groups.keys()),
                                                        columns=senescwheat_elements_inputs_t0.columns)
elements_inputs_t0 = cnwheat_elements_inputs_t0_reindexed.combine_first(farquharwheat_elements_inputs_t0_reindexed).combine_first(growthwheat_exposed_element_inputs_t0_with_element_column_reindexed).combine_first(senescwheat_elements_inputs_t0_reindexed)
elements_states = elements_inputs_t0.reindex_axis(cnwheat_simulation.Simulation.ELEMENTS_INPUTS_INDEXES + sorted(set(cnwheat_simulation.Simulation.ELEMENTS_STATE + farquharwheat_converter.FARQUHARWHEAT_INPUTS_OUTPUTS + growthwheat_simulation.ORGAN_INPUTS_OUTPUTS + senescwheat_converter.SENESCWHEAT_ELEMENTS_INPUTS_OUTPUTS)), axis=1)
# soils
cnwheat_soils_inputs_t0_reindexed = pd.DataFrame(cnwheat_soils_inputs_t0.values,
                                                 index=sorted(cnwheat_soils_inputs_t0.groupby(cnwheat_simulation.Simulation.SOILS_INPUTS_INDEXES).groups.keys()),
                                                 columns=cnwheat_soils_inputs_t0.columns)
soils_inputs_t0 = cnwheat_soils_inputs_t0_reindexed
soils_states = soils_inputs_t0.reindex_axis(cnwheat_converter.SOILS_STATE_VARIABLES, axis=1)

# Initialise simulations
growthwheat_hgzs_inputs = hgzs_states.loc[:, growthwheat_converter.HGZ_TOPOLOGY_COLUMNS + growthwheat_simulation.HGZ_INPUTS].reset_index(drop=True)
growthwheat_elements_inputs = elements_states.loc[(elements_states.element == 'LeafElement1') | (elements_states.element == 'StemElement'),
                                                 growthwheat_converter.ORGAN_TOPOLOGY_COLUMNS + growthwheat_simulation.ORGAN_INPUTS].reset_index(drop=True)
growthwheat_interface.initialize(g, {'hgz_inputs':growthwheat_hgzs_inputs, 'organ_inputs':growthwheat_elements_inputs}, adel_wheat)

all_simulation_steps = [] # to store the steps of the simulation

# run the simulators
current_time_of_the_system = time.time()

for t_caribu in xrange(start_time, stop_time, caribu_ts):
    #run_caribu.run_caribu(g, adel_wheat)
    for t_senescwheat in xrange(t_caribu, t_caribu + caribu_ts, senescwheat_ts):
        # initialize and run senescwheat
        senescwheat_roots_inputs = organs_states.loc[organs_states.organ == 'roots', senescwheat_converter.ROOTS_TOPOLOGY_COLUMNS + senescwheat_converter.SENESCWHEAT_ROOTS_INPUTS].reset_index(drop=True)
        senescwheat_elements_inputs = elements_states.loc[:, senescwheat_converter.ELEMENTS_TOPOLOGY_COLUMNS + senescwheat_converter.SENESCWHEAT_ELEMENTS_INPUTS].reset_index(drop=True)
        senescwheat_simulation_.initialize(senescwheat_converter.from_MTG(g, senescwheat_roots_inputs, senescwheat_elements_inputs))
        senescwheat_simulation_.run()
        senescwheat_roots_outputs, senescwheat_elements_outputs = senescwheat_converter.to_dataframes(senescwheat_simulation_.outputs)
        senescwheat_converter.update_MTG(senescwheat_simulation_.inputs, senescwheat_simulation_.outputs, g)
        # update the shared data
        senescwheat_roots_outputs_with_organ_column = senescwheat_roots_outputs.copy()
        senescwheat_roots_outputs_with_organ_column.loc[:, 'organ'] = 'roots'
        senescwheat_roots_outputs_with_organ_column_reindexed = pd.DataFrame(senescwheat_roots_outputs_with_organ_column.values,
                                                                             index=sorted(senescwheat_roots_outputs_with_organ_column.groupby(senescwheat_converter.ROOTS_TOPOLOGY_COLUMNS + ['organ']).groups.keys()),
                                                                             columns=senescwheat_roots_outputs_with_organ_column.columns)
        organs_states.update(senescwheat_roots_outputs_with_organ_column_reindexed)
        senescwheat_elements_outputs_reindexed = pd.DataFrame(senescwheat_elements_outputs.values,
                                                              index=sorted(senescwheat_elements_outputs.groupby(senescwheat_converter.ELEMENTS_TOPOLOGY_COLUMNS).groups.keys()),
                                                              columns=senescwheat_elements_outputs.columns)
        elements_states.update(senescwheat_elements_outputs_reindexed)

        for t_farquharwheat in xrange(t_senescwheat, t_senescwheat + senescwheat_ts, farquharwheat_ts):
            # get the meteo of the current step
            Ta, ambient_CO2, RH, Ur, PARi = meteo.loc[t_farquharwheat, ['air_temperature', 'ambient_CO2', 'humidity', 'Wind', 'PARi']]
            # initialize and run farquharwheat
            farquharwheat_elements_inputs = elements_states.loc[:, farquharwheat_converter.DATAFRAME_TOPOLOGY_COLUMNS + farquharwheat_converter.FARQUHARWHEAT_INPUTS].reset_index(drop=True)
            farquharwheat_simulation_.initialize(farquharwheat_converter.from_MTG(g, farquharwheat_elements_inputs))
            farquharwheat_simulation_.run(Ta, ambient_CO2, RH, Ur, PARi)
            farquharwheat_outputs = farquharwheat_converter.to_dataframe(farquharwheat_simulation_.outputs)
            farquharwheat_converter.update_MTG(farquharwheat_simulation_.inputs, farquharwheat_simulation_.outputs, g)
            # update the shared data
            farquharwheat_outputs_reindexed = pd.DataFrame(farquharwheat_outputs.values,
                                                           index=sorted(farquharwheat_outputs.groupby(farquharwheat_converter.DATAFRAME_TOPOLOGY_COLUMNS).groups.keys()),
                                                           columns=farquharwheat_outputs.columns)
            elements_states.update(farquharwheat_outputs_reindexed)

            for t_growthwheat in xrange(t_farquharwheat, t_farquharwheat + farquharwheat_ts, growthwheat_ts):
                # initialize and run growthwheat
                _, growthwheat_hgzs_outputs, growthwheat_elements_outputs = growthwheat_interface.run(g, growthwheat_ts * hour_to_second_conversion_factor, adel_wheat)
                # update the shared data
                growthwheat_hgzs_outputs_reindexed = pd.DataFrame(growthwheat_hgzs_outputs.values,
                                                                  index=sorted(growthwheat_hgzs_outputs.groupby(growthwheat_converter.HGZ_TOPOLOGY_COLUMNS).groups.keys()),
                                                                  columns=growthwheat_hgzs_outputs.columns)
                hgzs_states.update(growthwheat_hgzs_outputs_reindexed)
                growthwheat_elements_outputs_with_element_column = growthwheat_elements_outputs.copy()
                growthwheat_elements_outputs_with_element_column.loc[growthwheat_elements_outputs_with_element_column.organ == 'blade', 'element'] = 'LeafElement1'
                growthwheat_elements_outputs_with_element_column.loc[growthwheat_elements_outputs_with_element_column.organ != 'blade', 'element'] = 'StemElement'
                growthwheat_elements_outputs_with_element_column_reindexed = pd.DataFrame(growthwheat_elements_outputs_with_element_column.values,
                                                                                          index=sorted(growthwheat_elements_outputs_with_element_column.groupby(growthwheat_converter.ORGAN_TOPOLOGY_COLUMNS + ['element']).groups.keys()),
                                                                                          columns=growthwheat_elements_outputs_with_element_column.columns)
                elements_states.update(growthwheat_elements_outputs_with_element_column_reindexed)

                for t_cnwheat in xrange(t_growthwheat, t_growthwheat + growthwheat_ts, cnwheat_ts):
                    # initialize and run cnwheat
                    cnwheat_organs_inputs = organs_states.loc[:, cnwheat_converter.ORGANS_STATE_VARIABLES].reset_index(drop=True)
                    cnwheat_hgzs_inputs = hgzs_states.loc[:, cnwheat_converter.HGZS_STATE_VARIABLES].reset_index(drop=True)
                    cnwheat_elements_inputs = elements_states.loc[:, cnwheat_converter.ELEMENTS_STATE_VARIABLES].reset_index(drop=True)
                    cnwheat_soils_inputs = soils_states.reset_index(drop=True)
                    population = cnwheat_converter.from_MTG(g, organs_inputs=cnwheat_organs_inputs, hgzs_inputs=cnwheat_hgzs_inputs,
                                                                          elements_inputs=cnwheat_elements_inputs)
                    cnwheat_simulation_.initialize(population, cnwheat_converter.from_dataframes(soils_inputs=cnwheat_soils_inputs))
                    cnwheat_simulation_.run(start_time=t_cnwheat, stop_time=t_cnwheat+cnwheat_ts, number_of_output_steps=cnwheat_ts+1)
                    cnwheat_converter.update_MTG(population, g)

                    (cnwheat_plants_all_data,
                     cnwheat_axes_all_data,
                     cnwheat_metamers_all_data,
                     cnwheat_organs_all_data,
                     cnwheat_hgzs_all_data,
                     cnwheat_elements_all_data,
                     cnwheat_soils_all_data) = cnwheat_simulation_.postprocessings()

                    # update the shared data
                    cnwheat_organs_outputs_reindexed = pd.DataFrame(cnwheat_organs_all_data.values,
                                                                    index=sorted(cnwheat_organs_all_data.groupby(cnwheat_simulation.Simulation.ORGANS_OUTPUTS_INDEXES).groups.keys()),
                                                                    columns=cnwheat_organs_all_data.columns)
                    organs_states.update(cnwheat_organs_outputs_reindexed)
                    cnwheat_hgzs_outputs_reindexed = pd.DataFrame(cnwheat_hgzs_all_data.values,
                                                                  index=sorted(cnwheat_hgzs_all_data.groupby(cnwheat_simulation.Simulation.HIDDENGROWINGZONE_OUTPUTS_INDEXES).groups.keys()),
                                                                  columns=cnwheat_hgzs_all_data.columns)
                    hgzs_states.update(cnwheat_hgzs_outputs_reindexed)
                    cnwheat_elements_outputs_reindexed = pd.DataFrame(cnwheat_elements_all_data.values,
                                                                      index=sorted(cnwheat_elements_all_data.groupby(cnwheat_simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES).groups.keys()),
                                                                      columns=cnwheat_elements_all_data.columns)
                    elements_states.update(cnwheat_elements_outputs_reindexed)
                    cnwheat_soils_outputs_reindexed = pd.DataFrame(cnwheat_soils_all_data.values,
                                                                   index=sorted(cnwheat_soils_all_data.groupby(cnwheat_simulation.Simulation.SOILS_OUTPUTS_INDEXES).groups.keys()),
                                                                   columns=cnwheat_soils_all_data.columns)
                    soils_states.update(cnwheat_soils_outputs_reindexed)

                    # append the computed states to global list of states
                    all_simulation_steps.append(t_cnwheat)
                    axes_all_data_list.append(cnwheat_axes_all_data.loc[cnwheat_axes_all_data.t == t_cnwheat])
                    organs_all_data_list.append(organs_states.copy())
                    hgzs_all_data_list.append(hgzs_states.copy())
                    elements_all_data_list.append(elements_states.copy())
                    soils_all_data_list.append(soils_states.copy())

execution_time = int(time.time() - current_time_of_the_system)
print '\n', 'Simulation run in ', str(datetime.timedelta(seconds=execution_time))

# write all the computed states to CSV files
all_axes_states = pd.concat(axes_all_data_list)
all_axes_states.reset_index(inplace=True, drop=True)
all_axes_states.to_csv(AXES_STATES_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

all_organs_states = pd.concat(organs_all_data_list, keys=all_simulation_steps)
all_organs_states.reset_index(0, inplace=True)
all_organs_states.rename_axis({'level_0': 't'}, axis=1, inplace=True)
all_organs_states.to_csv(ORGANS_STATES_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

all_hgzs_states = pd.concat(hgzs_all_data_list, keys=all_simulation_steps)
all_hgzs_states.reset_index(0, inplace=True)
all_hgzs_states.rename_axis({'level_0': 't'}, axis=1, inplace=True)
all_hgzs_states.to_csv(HGZS_STATES_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

all_elements_states = pd.concat(elements_all_data_list, keys=all_simulation_steps)
all_elements_states.reset_index(0, inplace=True)
all_elements_states.rename_axis({'level_0': 't'}, axis=1, inplace=True)
all_elements_states.to_csv(ELEMENTS_STATES_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

all_soils_states = pd.concat(soils_all_data_list, keys=all_simulation_steps)
all_soils_states.reset_index(0, inplace=True)
all_soils_states.rename_axis({'level_0': 't'}, axis=1, inplace=True)
all_soils_states.to_csv(SOILS_STATES_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))


##from cnwheat import parameters
##from cnwheat import tools
##x_name = 't'
##x_label='Time (Hour)'
##
### 1) Photosynthetic organs
##ph_elements_output_df = pd.read_csv(ELEMENTS_INPUTS_OUTPUTS_FILEPATH)
##
##graph_variables_ph_elements = {'Ag': u'Gross photosynthesis (µmol m$^{-2}$ s$^{-1}$)','An': u'Net photosynthesis (µmol m$^{-2}$ s$^{-1}$)', 'Tr':u'Organ surfacic transpiration rate (mmol H$_{2}$0 m$^{-2}$ s$^{-1}$)', 'Transpiration':u'Organ transpiration rate (mmol H$_{2}$0 s$^{-1}$)', 'Rd': u'Mitochondrial respiration rate of organ in light (µmol C h$^{-1}$)', 'Ts': u'Temperature surface (°C)', 'gs': u'Conductance stomatique (mol m$^{-2}$ s$^{-1}$)',
##                   'Conc_TriosesP': u'[TriosesP] (µmol g$^{-1}$ mstruct)', 'Conc_Starch':u'[Starch] (µmol g$^{-1}$ mstruct)', 'Conc_Sucrose':u'[Sucrose] (µmol g$^{-1}$ mstruct)', 'Conc_Fructan':u'[Fructan] (µmol g$^{-1}$ mstruct)',
##                   'Conc_Nitrates': u'[Nitrates] (µmol g$^{-1}$ mstruct)', 'Conc_Amino_Acids': u'[Amino_Acids] (µmol g$^{-1}$ mstruct)', 'Conc_Proteins': u'[Proteins] (g g$^{-1}$ mstruct)',
##                   'Nitrates_import': u'Total nitrates imported (µmol h$^{-1}$)', 'Amino_Acids_import': u'Total amino acids imported (µmol N h$^{-1}$)',
##                   'S_Amino_Acids': u'[Rate of amino acids synthesis] (µmol N g$^{-1}$ mstruct h$^{-1}$)', 'S_Proteins': u'Rate of protein synthesis (µmol N g$^{-1}$ mstruct h$^{-1}$)', 'D_Proteins': u'Rate of protein degradation (µmol N g$^{-1}$ mstruct h$^{-1}$)', 'k_proteins': u'Relative rate of protein degradation (s$^{-1}$)',
##                   'Loading_Sucrose': u'Loading Sucrose (µmol C sucrose h$^{-1}$)', 'Loading_Amino_Acids': u'Loading Amino acids (µmol N amino acids h$^{-1}$)',
##                   'green_area': u'Green area (m$^{2}$)', 'R_phloem_loading': u'Respiration phloem loading (µmol C h$^{-1}$)', 'R_Nnit_red': u'Respiration nitrate reduction (µmol C h$^{-1}$)', 'R_residual': u'Respiration residual (µmol C h$^{-1}$)', 'R_maintenance': u'Respiration residual (µmol C h$^{-1}$)',
##                   'mstruct': u'Structural mass (g)', 'Nstruct': u'Structural N mass (g)',
##                   'Conc_cytokinins':u'[cytokinins] (UA g$^{-1}$ mstruct)', 'D_cytokinins':u'Cytokinin degradation (UA g$^{-1}$ mstruct)', 'cytokinins_import':u'Cytokinin import (UA)'}
##
##
##for org_ph in (['blade'], ['sheath'], ['internode'], ['peduncle', 'ear']):
##    for variable_name, variable_label in graph_variables_ph_elements.iteritems():
##        graph_name = variable_name + '_' + '_'.join(org_ph) + '.PNG'
##        tools.plot_cnwheat_ouputs(ph_elements_output_df,
##                      x_name = x_name,
##                      y_name = variable_name,
##                      x_label=x_label,
##                      y_label=variable_label,
##                      filters={'organ': org_ph},
##                      plot_filepath=os.path.join(GRAPHS_DIRPATH, graph_name),
##                      explicit_label=False)