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
from elongwheat import simulation as elongwheat_simulation, model as elongwheat_model, converter as elongwheat_converter, interface as elongwheat_interface
from growthwheat import simulation as growthwheat_simulation, model as growthwheat_model, converter as growthwheat_converter, interface as growthwheat_interface

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
CNWHEAT_HIDDENZONE_INPUTS_FILEPATH = os.path.join(CNWHEAT_INPUTS_DIRPATH, 'hiddenzones_inputs.csv')
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

# elongwheat inputs at t0
ELONGWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'elongwheat')
ELONGWHEAT_HIDDENZONE_INPUTS_FILEPATH = os.path.join(ELONGWHEAT_INPUTS_DIRPATH, 'hiddenzones_inputs.csv')
ELONGWHEAT_ORGANS_INPUTS_FILEPATH = os.path.join(ELONGWHEAT_INPUTS_DIRPATH, 'organs_inputs.csv')

# growthwheat inputs at t0
GROWTHWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'growthwheat')
GROWTHWHEAT_HIDDENZONE_INPUTS_FILEPATH = os.path.join(GROWTHWHEAT_INPUTS_DIRPATH, 'hiddenzones_inputs.csv')
GROWTHWHEAT_ORGANS_INPUTS_FILEPATH = os.path.join(GROWTHWHEAT_INPUTS_DIRPATH, 'organs_inputs.csv')
GROWTHWHEAT_ROOTS_INPUTS_FILEPATH = os.path.join(GROWTHWHEAT_INPUTS_DIRPATH, 'roots_inputs.csv')

# the path of the CSV files where to save the states of the modeled system at each step
OUTPUTS_DIRPATH = 'outputs'
AXES_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'axes_states.csv')
ORGANS_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'organs_states.csv')
HIDDENZONES_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'hiddenzones_states.csv')
ELEMENTS_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'elements_states.csv')
SOILS_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'soils_states.csv')

INPUTS_OUTPUTS_PRECISION = 10

LOGGING_CONFIG_FILEPATH = os.path.join('..', '..', 'logging.json')

LOGGING_LEVEL = logging.INFO # can be one of: DEBUG, INFO, WARNING, ERROR, CRITICAL

cnwheat_tools.setup_logging(LOGGING_CONFIG_FILEPATH, LOGGING_LEVEL, log_model=False, log_compartments=False, log_derivatives=False)

def main(stop_time, run_simu=True, make_graphs=True):
    if run_simu:
        meteo = pd.read_csv(METEO_FILEPATH, index_col='t')

        current_time_of_the_system = time.time()

        # define the time step in hours for each simulator
        caribu_ts = 4
        senescwheat_ts = 2
        farquharwheat_ts = 2
        elongwheat_ts = 1
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
        properties_to_convert = {'lengths': ['shape_mature_length', 'shape_max_width', 'length', 'visible_length', 'width'], 'areas': ['green_area']}
        adel_wheat.convert_to_SI_units(g, properties_to_convert)

        # read cnwheat inputs at t0
        cnwheat_plants_inputs_t0 = pd.read_csv(CNWHEAT_PLANTS_INPUTS_FILEPATH)
        cnwheat_axes_inputs_t0 = pd.read_csv(CNWHEAT_AXES_INPUTS_FILEPATH)
        cnwheat_metamers_inputs_t0 = pd.read_csv(CNWHEAT_METAMERS_INPUTS_FILEPATH)
        cnwheat_organs_inputs_t0 = pd.read_csv(CNWHEAT_ORGANS_INPUTS_FILEPATH)
        cnwheat_hiddenzones_inputs_t0 = pd.read_csv(CNWHEAT_HIDDENZONE_INPUTS_FILEPATH)
        cnwheat_elements_inputs_t0 = pd.read_csv(CNWHEAT_ELEMENTS_INPUTS_FILEPATH)
        cnwheat_soils_inputs_t0 = pd.read_csv(CNWHEAT_SOILS_INPUTS_FILEPATH)

        # read farquharwheat inputs at t0
        farquharwheat_elements_inputs_t0 = pd.read_csv(FARQUHARWHEAT_INPUTS_FILEPATH)

        # read senescwheat inputs at t0
        senescwheat_roots_inputs_t0 = pd.read_csv(SENESCWHEAT_ROOTS_INPUTS_FILEPATH)
        senescwheat_elements_inputs_t0 = pd.read_csv(SENESCWHEAT_ELEMENTS_INPUTS_FILEPATH)

        # read elongwheat inputs at t0
        elongwheat_hiddenzones_inputs_t0 = pd.read_csv(ELONGWHEAT_HIDDENZONE_INPUTS_FILEPATH)
        elongwheat_organ_inputs_t0 = pd.read_csv(ELONGWHEAT_ORGANS_INPUTS_FILEPATH)

        # read growthwheat inputs at t0
        growthwheat_hiddenzones_inputs_t0 = pd.read_csv(GROWTHWHEAT_HIDDENZONE_INPUTS_FILEPATH)
        growthwheat_organ_inputs_t0 = pd.read_csv(GROWTHWHEAT_ORGANS_INPUTS_FILEPATH)
        growthwheat_root_inputs_t0 = pd.read_csv(GROWTHWHEAT_ROOTS_INPUTS_FILEPATH)

        # Initialise simulations
        elongwheat_interface.initialize(g, {'hiddenzone_inputs':elongwheat_hiddenzones_inputs_t0, 'organ_inputs':elongwheat_organ_inputs_t0}, adel_wheat)
        growthwheat_interface.initialize(g, {'hiddenzone_inputs':growthwheat_hiddenzones_inputs_t0, 'organ_inputs':growthwheat_organ_inputs_t0, 'root_inputs':growthwheat_root_inputs_t0})

        # Update geometry
        adel_wheat.update_geometry(g, SI_units=True, properties_to_convert=properties_to_convert) # Return mtg with non-SI units
        #adel_wheat.plot(g)
        adel_wheat.convert_to_SI_units(g, properties_to_convert)

        # define the start and the end of the whole simulation (in hours)
        start_time = 0
        stop_time = stop_time

        # define lists of dataframes to store the state of the system at each step.
        axes_all_data_list = []
        organs_all_data_list = [] # organs which belong to axes: roots, phloem, grains
        hiddenzones_all_data_list = []
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
        organs_inputs_outputs = organs_inputs_t0.reindex_axis(cnwheat_simulation.Simulation.ORGANS_INPUTS_INDEXES + sorted(set(cnwheat_simulation.Simulation.ORGANS_INPUTS_OUTPUTS + senescwheat_converter.SENESCWHEAT_ROOTS_INPUTS_OUTPUTS)), axis=1)

        # hidden zones
        cnwheat_hiddenzones_inputs_t0_reindexed = pd.DataFrame(cnwheat_hiddenzones_inputs_t0.values,
                                                        index=sorted(cnwheat_hiddenzones_inputs_t0.groupby(cnwheat_simulation.Simulation.HIDDENZONE_INPUTS_INDEXES).groups.keys()),
                                                        columns=cnwheat_hiddenzones_inputs_t0.columns)
        elongwheat_hiddenzones_inputs_t0_reindexed = pd.DataFrame(elongwheat_hiddenzones_inputs_t0.values.tolist(),
                                                            index=sorted(elongwheat_hiddenzones_inputs_t0.groupby(elongwheat_converter.HIDDENZONE_TOPOLOGY_COLUMNS).groups.keys()),
                                                            columns=elongwheat_hiddenzones_inputs_t0.columns)
        hiddenzones_inputs_t0 = cnwheat_hiddenzones_inputs_t0_reindexed.combine_first(elongwheat_hiddenzones_inputs_t0_reindexed)
        dtypes = elongwheat_hiddenzones_inputs_t0_reindexed.dtypes.combine_first(cnwheat_hiddenzones_inputs_t0_reindexed.dtypes)
        for k, v in dtypes.iteritems():
            hiddenzones_inputs_t0[k] = hiddenzones_inputs_t0[k].astype(v)

        hiddenzones_inputs_outputs = hiddenzones_inputs_t0.reindex_axis(cnwheat_simulation.Simulation.HIDDENZONE_INPUTS_INDEXES + sorted(set(cnwheat_simulation.Simulation.HIDDENZONE_INPUTS_OUTPUTS + elongwheat_simulation.HIDDENZONE_INPUTS_OUTPUTS + growthwheat_simulation.HIDDENZONE_INPUTS_OUTPUTS)), axis=1)
        # elements
        cnwheat_elements_inputs_t0_reindexed = pd.DataFrame(cnwheat_elements_inputs_t0.values,
                                                            index=sorted(cnwheat_elements_inputs_t0.groupby(cnwheat_simulation.Simulation.ELEMENTS_INPUTS_INDEXES).groups.keys()),
                                                            columns=cnwheat_elements_inputs_t0.columns)
        farquharwheat_elements_inputs_t0_reindexed = pd.DataFrame(farquharwheat_elements_inputs_t0.values,
                                                                  index=sorted(farquharwheat_elements_inputs_t0.groupby(farquharwheat_converter.DATAFRAME_TOPOLOGY_COLUMNS).groups.keys()),
                                                                  columns=farquharwheat_elements_inputs_t0.columns)
        elongwheat_elements_inputs_t0_with_element_column = elongwheat_organ_inputs_t0.copy()
        elongwheat_elements_inputs_t0_with_element_column.loc[elongwheat_elements_inputs_t0_with_element_column.organ == 'blade', 'element'] = 'LeafElement1'
        elongwheat_elements_inputs_t0_with_element_column.loc[elongwheat_elements_inputs_t0_with_element_column.organ != 'blade', 'element'] = 'StemElement'
        elongwheat_organ_inputs_t0_with_element_column_reindexed = pd.DataFrame(elongwheat_elements_inputs_t0_with_element_column.values,
                                                                                           index=sorted(elongwheat_elements_inputs_t0_with_element_column.groupby(elongwheat_converter.ORGAN_TOPOLOGY_COLUMNS + ['element']).groups.keys()),
                                                                                           columns=elongwheat_elements_inputs_t0_with_element_column.columns)
        senescwheat_elements_inputs_t0_reindexed = pd.DataFrame(senescwheat_elements_inputs_t0.values,
                                                                index=sorted(senescwheat_elements_inputs_t0.groupby(senescwheat_converter.ELEMENTS_TOPOLOGY_COLUMNS).groups.keys()),
                                                                columns=senescwheat_elements_inputs_t0.columns)
        elements_inputs_t0 = cnwheat_elements_inputs_t0_reindexed.combine_first(farquharwheat_elements_inputs_t0_reindexed).combine_first(elongwheat_organ_inputs_t0_with_element_column_reindexed).combine_first(senescwheat_elements_inputs_t0_reindexed)
        elements_inputs_outputs = elements_inputs_t0.reindex_axis(cnwheat_simulation.Simulation.ELEMENTS_INPUTS_INDEXES + sorted(set(cnwheat_simulation.Simulation.ELEMENTS_INPUTS_OUTPUTS + farquharwheat_converter.FARQUHARWHEAT_INPUTS_OUTPUTS + elongwheat_simulation.ORGAN_INPUTS_OUTPUTS + senescwheat_converter.SENESCWHEAT_ELEMENTS_INPUTS_OUTPUTS)), axis=1)
        # soils
        cnwheat_soils_inputs_t0_reindexed = pd.DataFrame(cnwheat_soils_inputs_t0.values,
                                                         index=sorted(cnwheat_soils_inputs_t0.groupby(cnwheat_simulation.Simulation.SOILS_INPUTS_INDEXES).groups.keys()),
                                                         columns=cnwheat_soils_inputs_t0.columns)
        soils_inputs_t0 = cnwheat_soils_inputs_t0_reindexed
        soils_inputs_outputs = soils_inputs_t0.reindex_axis(cnwheat_simulation.Simulation.SOILS_INPUTS_INDEXES + sorted(set(cnwheat_simulation.Simulation.SOILS_INPUTS_OUTPUTS)), axis=1)

        all_simulation_steps = [] # to store the steps of the simulation

        # run the simulators
        current_time_of_the_system = time.time()

        for t_caribu in xrange(start_time, stop_time, caribu_ts):
            caribu_outputs = run_caribu.run_caribu(g, adel_wheat)
            # update the shared data
            caribu_outputs_reindexed = pd.DataFrame(caribu_outputs.values,
                                                           index=sorted(caribu_outputs.groupby(run_caribu.DATAFRAME_TOPOLOGY_COLUMNS).groups.keys()),
                                                           columns=caribu_outputs.columns)
            elements_inputs_outputs.update(caribu_outputs_reindexed)

            for t_senescwheat in xrange(t_caribu, t_caribu + caribu_ts, senescwheat_ts):
                # initialize and run senescwheat
                senescwheat_roots_inputs = organs_inputs_outputs.loc[organs_inputs_outputs.organ == 'roots', senescwheat_converter.ROOTS_TOPOLOGY_COLUMNS + senescwheat_converter.SENESCWHEAT_ROOTS_INPUTS].reset_index(drop=True)
                senescwheat_elements_inputs = elements_inputs_outputs.loc[:, senescwheat_converter.ELEMENTS_TOPOLOGY_COLUMNS + senescwheat_converter.SENESCWHEAT_ELEMENTS_INPUTS].reset_index(drop=True)
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
                organs_inputs_outputs.update(senescwheat_roots_outputs_with_organ_column_reindexed)
                senescwheat_elements_outputs_reindexed = pd.DataFrame(senescwheat_elements_outputs.values,
                                                                      index=sorted(senescwheat_elements_outputs.groupby(senescwheat_converter.ELEMENTS_TOPOLOGY_COLUMNS).groups.keys()),
                                                                      columns=senescwheat_elements_outputs.columns)
                elements_inputs_outputs.update(senescwheat_elements_outputs_reindexed)

                for t_farquharwheat in xrange(t_senescwheat, t_senescwheat + senescwheat_ts, farquharwheat_ts):
                    # get the meteo of the current step
                    Ta, ambient_CO2, RH, Ur, PARi = meteo.loc[t_farquharwheat, ['air_temperature', 'ambient_CO2', 'humidity', 'Wind', 'PARi']]
                    # initialize and run farquharwheat
                    farquharwheat_elements_inputs = elements_inputs_outputs.loc[:, farquharwheat_converter.DATAFRAME_TOPOLOGY_COLUMNS + farquharwheat_converter.FARQUHARWHEAT_INPUTS].reset_index(drop=True)
                    farquharwheat_simulation_.initialize(farquharwheat_converter.from_MTG(g, farquharwheat_elements_inputs))
                    farquharwheat_simulation_.run(Ta, ambient_CO2, RH, Ur, PARi)
                    farquharwheat_outputs = farquharwheat_converter.to_dataframe(farquharwheat_simulation_.outputs)
                    farquharwheat_converter.update_MTG(farquharwheat_simulation_.inputs, farquharwheat_simulation_.outputs, g)
                    # update the shared data
                    farquharwheat_outputs_reindexed = pd.DataFrame(farquharwheat_outputs.values,
                                                                   index=sorted(farquharwheat_outputs.groupby(farquharwheat_converter.DATAFRAME_TOPOLOGY_COLUMNS).groups.keys()),
                                                                   columns=farquharwheat_outputs.columns)
                    elements_inputs_outputs.update(farquharwheat_outputs_reindexed)

                    for t_elongwheat in xrange(t_farquharwheat, t_farquharwheat + farquharwheat_ts, elongwheat_ts):
                        # Run elongwheat
                        _, elongwheat_hiddenzones_outputs, elongwheat_organs_outputs = elongwheat_interface.run(g, elongwheat_ts * hour_to_second_conversion_factor, adel_wheat)
                        # update the shared data
                        elongwheat_hiddenzones_outputs_reindexed = pd.DataFrame(elongwheat_hiddenzones_outputs.values,
                                                                          index=sorted(elongwheat_hiddenzones_outputs.groupby(elongwheat_converter.HIDDENZONE_TOPOLOGY_COLUMNS).groups.keys()),
                                                                          columns=elongwheat_hiddenzones_outputs.columns)
                        hiddenzones_inputs_outputs.update(elongwheat_hiddenzones_outputs_reindexed)
                        elongwheat_organs_outputs_with_organ_column = elongwheat_organs_outputs.copy()
                        elongwheat_organs_outputs_with_organ_column.loc[elongwheat_organs_outputs_with_organ_column.organ == 'blade', 'element'] = 'LeafElement1'
                        elongwheat_organs_outputs_with_organ_column.loc[elongwheat_organs_outputs_with_organ_column.organ != 'blade', 'element'] = 'StemElement'
                        elongwheat_organs_outputs_with_organ_column_reindexed = pd.DataFrame(elongwheat_organs_outputs_with_organ_column.values,
                                                                                                  index=sorted(elongwheat_organs_outputs_with_organ_column.groupby(elongwheat_converter.ORGAN_TOPOLOGY_COLUMNS + ['element']).groups.keys()),
                                                                                                  columns=elongwheat_organs_outputs_with_organ_column.columns)
                        elements_inputs_outputs.update(elongwheat_organs_outputs_with_organ_column_reindexed)
                        # Update geometry
                        adel_wheat.update_geometry(g, SI_units=True, properties_to_convert=properties_to_convert) # Return mtg with non-SI units
                        #adel_wheat.plot(g)
                        adel_wheat.convert_to_SI_units(g, properties_to_convert)

                        for t_growthwheat in xrange(t_elongwheat, t_elongwheat + elongwheat_ts, growthwheat_ts):
                            # Run growthwheat
                            _, growthwheat_hiddenzones_outputs, growthwheat_organs_outputs, growthwheat_roots_outputs = growthwheat_interface.run(g, growthwheat_ts * hour_to_second_conversion_factor)
                            # update the shared data
                            growthwheat_hiddenzones_outputs_reindexed = pd.DataFrame(growthwheat_hiddenzones_outputs.values,
                                                                              index=sorted(growthwheat_hiddenzones_outputs.groupby(growthwheat_converter.HIDDENZONE_TOPOLOGY_COLUMNS).groups.keys()),
                                                                              columns=growthwheat_hiddenzones_outputs.columns)
                            hiddenzones_inputs_outputs.update(growthwheat_hiddenzones_outputs_reindexed)
                            growthwheat_organs_outputs_with_organ_column = growthwheat_organs_outputs.copy()
                            growthwheat_organs_outputs_with_organ_column.loc[growthwheat_organs_outputs_with_organ_column.organ == 'blade', 'element'] = 'LeafElement1'
                            growthwheat_organs_outputs_with_organ_column.loc[growthwheat_organs_outputs_with_organ_column.organ != 'blade', 'element'] = 'StemElement'
                            growthwheat_organs_outputs_with_organ_column_reindexed = pd.DataFrame(growthwheat_organs_outputs_with_organ_column.values,
                                                                                                      index=sorted(growthwheat_organs_outputs_with_organ_column.groupby(growthwheat_converter.ORGAN_TOPOLOGY_COLUMNS + ['element']).groups.keys()),
                                                                                                      columns=growthwheat_organs_outputs_with_organ_column.columns)
                            elements_inputs_outputs.update(growthwheat_organs_outputs_with_organ_column_reindexed)
                            growthwheat_roots_outputs_reindexed = pd.DataFrame(growthwheat_roots_outputs.values,
                                                                              index=sorted(growthwheat_roots_outputs.groupby(growthwheat_converter.ROOT_TOPOLOGY_COLUMNS).groups.keys()),
                                                                              columns=growthwheat_roots_outputs.columns)
                            organs_inputs_outputs.update(growthwheat_roots_outputs_reindexed)

                            for t_cnwheat in xrange(t_growthwheat, t_growthwheat + growthwheat_ts, cnwheat_ts):
                                # initialize and run cnwheat
                                cnwheat_organs_inputs = organs_inputs_outputs.loc[:, cnwheat_converter.ORGANS_STATE_VARIABLES].reset_index(drop=True)
                                cnwheat_hiddenzones_inputs = hiddenzones_inputs_outputs.loc[:, cnwheat_converter.HIDDENZONES_STATE_VARIABLES].reset_index(drop=True)
                                cnwheat_elements_inputs = elements_inputs_outputs.loc[:, cnwheat_converter.ELEMENTS_STATE_VARIABLES].reset_index(drop=True)
                                cnwheat_soils_inputs = soils_inputs_outputs.reset_index(drop=True)
                                population = cnwheat_converter.from_MTG(g, organs_inputs=cnwheat_organs_inputs, hiddenzones_inputs=cnwheat_hiddenzones_inputs,
                                                                                      elements_inputs=cnwheat_elements_inputs)
                                cnwheat_simulation_.initialize(population, cnwheat_converter.from_dataframes(soils_inputs=cnwheat_soils_inputs))
                                cnwheat_simulation_.run(start_time=t_cnwheat, stop_time=t_cnwheat+cnwheat_ts, number_of_output_steps=cnwheat_ts+1)
                                cnwheat_converter.update_MTG(population, g)

                                (cnwheat_plants_all_data,
                                 cnwheat_axes_all_data,
                                 cnwheat_metamers_all_data,
                                 cnwheat_organs_all_data,
                                 cnwheat_hiddenzones_all_data,
                                 cnwheat_elements_all_data,
                                 cnwheat_soils_all_data) = cnwheat_simulation_.postprocessings()

                                # update the shared data
                                ## Organs
                                cnwheat_organs_all_data = cnwheat_organs_all_data.loc[cnwheat_organs_all_data.t == t_cnwheat+1, :].reset_index(drop=True)
                                cnwheat_organs_outputs_reindexed = pd.DataFrame(cnwheat_organs_all_data.values,
                                                                                index=sorted(cnwheat_organs_all_data.groupby(cnwheat_simulation.Simulation.ORGANS_INPUTS_INDEXES).groups.keys()),
                                                                                columns=cnwheat_organs_all_data.columns)
                                organs_inputs_outputs.update(cnwheat_organs_outputs_reindexed)
                                #organs_inputs_outputs = organs_inputs_outputs.combine_first(cnwheat_organs_outputs_reindexed)
                                ## Hidden zones
                                cnwheat_hiddenzones_all_data = cnwheat_hiddenzones_all_data.loc[cnwheat_hiddenzones_all_data.t == t_cnwheat+1, :].reset_index(drop=True)
                                cnwheat_hiddenzones_outputs_reindexed = pd.DataFrame(cnwheat_hiddenzones_all_data.values,
                                                                              index=sorted(cnwheat_hiddenzones_all_data.groupby(cnwheat_simulation.Simulation.HIDDENZONE_INPUTS_INDEXES).groups.keys()),
                                                                              columns=cnwheat_hiddenzones_all_data.columns)
                                hiddenzones_inputs_outputs.update(cnwheat_hiddenzones_outputs_reindexed)
                                #hiddenzones_inputs_outputs = hiddenzones_inputs_outputs.combine_first(cnwheat_hiddenzones_outputs_reindexed)
                                ## Elements
                                cnwheat_elements_all_data = cnwheat_elements_all_data.loc[cnwheat_elements_all_data.t == t_cnwheat+1, :].reset_index(drop=True)
                                cnwheat_elements_outputs_reindexed = pd.DataFrame(cnwheat_elements_all_data.values,
                                                                                  index=sorted(cnwheat_elements_all_data.groupby(cnwheat_simulation.Simulation.ELEMENTS_INPUTS_INDEXES).groups.keys()),
                                                                                  columns=cnwheat_elements_all_data.columns)
                                elements_inputs_outputs.update(cnwheat_elements_outputs_reindexed)
                                #elements_inputs_outputs = elements_inputs_outputs.combine_first(cnwheat_elements_outputs_reindexed)
                                ## Soil
                                cnwheat_soils_all_data = cnwheat_soils_all_data.loc[cnwheat_soils_all_data.t == t_cnwheat+1, :].reset_index(drop=True)
                                cnwheat_soils_outputs_reindexed = pd.DataFrame(cnwheat_soils_all_data.values,
                                                                               index=sorted(cnwheat_soils_all_data.groupby(cnwheat_simulation.Simulation.SOILS_INPUTS_INDEXES).groups.keys()),
                                                                               columns=cnwheat_soils_all_data.columns)
                                soils_inputs_outputs.update(cnwheat_soils_outputs_reindexed)
                                #soils_inputs_outputs = soils_inputs_outputs.combine_first(cnwheat_soils_outputs_reindexed)

                                # append the computed states to global list of states
                                all_simulation_steps.append(t_cnwheat)
                                axes_all_data_list.append(cnwheat_axes_all_data.loc[cnwheat_axes_all_data.t == t_cnwheat+1])
                                organs_all_data_list.append(organs_inputs_outputs.copy())
                                hiddenzones_all_data_list.append(hiddenzones_inputs_outputs.copy())
                                elements_all_data_list.append(elements_inputs_outputs.copy())
                                soils_all_data_list.append(soils_inputs_outputs.copy())

        execution_time = int(time.time() - current_time_of_the_system)
        print '\n', 'Simulation run in ', str(datetime.timedelta(seconds=execution_time))

        # write all the computed states to CSV files
        all_axes_inputs_outputs = pd.concat(axes_all_data_list)
        all_axes_inputs_outputs.reset_index(inplace=True, drop=True)
        all_axes_inputs_outputs.to_csv(AXES_STATES_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

        all_organs_inputs_outputs = pd.concat(organs_all_data_list, keys=all_simulation_steps)
        all_organs_inputs_outputs.reset_index(0, inplace=True)
        all_organs_inputs_outputs.rename_axis({'level_0': 't'}, axis=1, inplace=True)
        all_organs_inputs_outputs.to_csv(ORGANS_STATES_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

        all_hiddenzones_inputs_outputs = pd.concat(hiddenzones_all_data_list, keys=all_simulation_steps)
        all_hiddenzones_inputs_outputs.reset_index(0, inplace=True)
        all_hiddenzones_inputs_outputs.rename_axis({'level_0': 't'}, axis=1, inplace=True)
        all_hiddenzones_inputs_outputs.to_csv(HIDDENZONES_STATES_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

        all_elements_inputs_outputs = pd.concat(elements_all_data_list, keys=all_simulation_steps)
        all_elements_inputs_outputs.reset_index(0, inplace=True)
        all_elements_inputs_outputs.rename_axis({'level_0': 't'}, axis=1, inplace=True)
        all_elements_inputs_outputs.to_csv(ELEMENTS_STATES_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

        all_soils_inputs_outputs = pd.concat(soils_all_data_list, keys=all_simulation_steps)
        all_soils_inputs_outputs.reset_index(0, inplace=True)
        all_soils_inputs_outputs.rename_axis({'level_0': 't'}, axis=1, inplace=True)
        all_soils_inputs_outputs.to_csv(SOILS_STATES_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))


    ########POST-PROCESSING##
    if make_graphs:
        from cnwheat import parameters
        from cnwheat import tools
        x_name = 't'
        x_label='Time (Hour)'

        # 1) Photosynthetic organs
        ph_elements_output_df = pd.read_csv(ELEMENTS_STATES_FILEPATH)

        graph_variables_ph_elements = {'Eabsm2': u'Absorbed PAR (µmol m$^{-2}$ s$^{-1}$)', 'Ag': u'Gross photosynthesis (µmol m$^{-2}$ s$^{-1}$)','An': u'Net photosynthesis (µmol m$^{-2}$ s$^{-1}$)', 'Tr':u'Organ surfacic transpiration rate (mmol H$_{2}$0 m$^{-2}$ s$^{-1}$)', 'Transpiration':u'Organ transpiration rate (mmol H$_{2}$0 s$^{-1}$)', 'Rd': u'Mitochondrial respiration rate of organ in light (µmol C h$^{-1}$)', 'Ts': u'Temperature surface (°C)', 'gs': u'Conductance stomatique (mol m$^{-2}$ s$^{-1}$)',
                           'Conc_TriosesP': u'[TriosesP] (µmol g$^{-1}$ mstruct)', 'Conc_Starch':u'[Starch] (µmol g$^{-1}$ mstruct)', 'Conc_Sucrose':u'[Sucrose] (µmol g$^{-1}$ mstruct)', 'Conc_Fructan':u'[Fructan] (µmol g$^{-1}$ mstruct)',
                           'Conc_Nitrates': u'[Nitrates] (µmol g$^{-1}$ mstruct)', 'Conc_Amino_Acids': u'[Amino_Acids] (µmol g$^{-1}$ mstruct)', 'Conc_Proteins': u'[Proteins] (g g$^{-1}$ mstruct)',
                           'Nitrates_import': u'Total nitrates imported (µmol h$^{-1}$)', 'Amino_Acids_import': u'Total amino acids imported (µmol N h$^{-1}$)',
                           'S_Amino_Acids': u'[Rate of amino acids synthesis] (µmol N g$^{-1}$ mstruct h$^{-1}$)', 'S_Proteins': u'Rate of protein synthesis (µmol N g$^{-1}$ mstruct h$^{-1}$)', 'D_Proteins': u'Rate of protein degradation (µmol N g$^{-1}$ mstruct h$^{-1}$)', 'k_proteins': u'Relative rate of protein degradation (s$^{-1}$)',
                           'Loading_Sucrose': u'Loading Sucrose (µmol C sucrose h$^{-1}$)', 'Loading_Amino_Acids': u'Loading Amino acids (µmol N amino acids h$^{-1}$)',
                           'green_area': u'Green area (m$^{2}$)', 'R_phloem_loading': u'Respiration phloem loading (µmol C h$^{-1}$)', 'R_Nnit_red': u'Respiration nitrate reduction (µmol C h$^{-1}$)', 'R_residual': u'Respiration residual (µmol C h$^{-1}$)', 'R_maintenance': u'Respiration residual (µmol C h$^{-1}$)',
                           'mstruct': u'Structural mass (g)', 'Nstruct': u'Structural N mass (g)',
                           'Conc_cytokinins':u'[cytokinins] (UA g$^{-1}$ mstruct)', 'D_cytokinins':u'Cytokinin degradation (UA g$^{-1}$ mstruct)', 'cytokinins_import':u'Cytokinin import (UA)'}


        for org_ph in (['blade'], ['sheath'], ['internode'], ['peduncle', 'ear']):
            for variable_name, variable_label in graph_variables_ph_elements.iteritems():
                graph_name = variable_name + '_' + '_'.join(org_ph) + '.PNG'
                tools.plot_cnwheat_ouputs(ph_elements_output_df,
                              x_name = x_name,
                              y_name = variable_name,
                              x_label=x_label,
                              y_label=variable_label,
                              filters={'organ': org_ph},
                              plot_filepath=os.path.join(GRAPHS_DIRPATH, graph_name),
                              explicit_label=False)

        # 2) Roots, grains and phloem
        organs_output_df = pd.read_csv(ORGANS_STATES_FILEPATH)

        graph_variables_organs = {'Conc_Sucrose':u'[Sucrose] (µmol g$^{-1}$ mstruct)', 'Dry_Mass':'Dry mass (g)',
                            'Conc_Nitrates': u'[Nitrates] (µmol g$^{-1}$ mstruct)', 'Conc_Amino_Acids':u'[Amino Acids] (µmol g$^{-1}$ mstruct)', 'Proteins_N_Mass': u'[N Proteins] (g)',
                            'Uptake_Nitrates':u'Nitrates uptake (µmol h$^{-1}$)', 'Unloading_Sucrose':u'Unloaded sucrose (µmol C g$^{-1}$ mstruct h$^{-1}$)', 'Unloading_Amino_Acids':u'Unloaded Amino Acids (µmol N AA g$^{-1}$ mstruct h$^{-1}$)',
                            'S_Amino_Acids': u'Rate of amino acids synthesis (µmol N g$^{-1}$ mstruct h$^{-1}$)', 'S_Proteins': u'Rate of protein synthesis (µmol N h$^{-1}$)', 'Export_Nitrates': u'Total export of nitrates (µmol N h$^{-1}$)', 'Export_Amino_Acids': u'Total export of Amino acids (µmol N h$^{-1}$)',
                            'R_Nnit_upt': u'Respiration nitrates uptake (µmol C h$^{-1}$)', 'R_Nnit_red': u'Respiration nitrate reduction (µmol C h$^{-1}$)', 'R_residual': u'Respiration residual (µmol C h$^{-1}$)', 'R_maintenance': u'Respiration residual (µmol C h$^{-1}$)',
                            'R_grain_growth_struct': u'Respiration grain structural growth (µmol C h$^{-1}$)', 'R_grain_growth_starch': u'Respiration grain starch growth (µmol C h$^{-1}$)',
                            'R_growth': u'Growth respiration of roots (µmol C h$^{-1}$)', 'mstruct': u'Structural mass (g)',
                            'C_exudation': u'Carbon lost by root exudation (µmol C g$^{-1}$ mstruct h$^{-1}$', 'N_exudation': u'Nitrogen lost by root exudation (µmol N g$^{-1}$ mstruct h$^{-1}$',
                            'Conc_cytokinins':u'[cytokinins] (UA g$^{-1}$ mstruct)', 'S_cytokinins':u'Rate of cytokinins synthesis (UA g$^{-1}$ mstruct)', 'Export_cytokinins': 'Export of cytokinins from roots (UA h$^{-1}$)',
                            'HATS_LATS': u'Potential uptake (µmol h$^{-1}$)' , 'regul_transpiration':'Regulating transpiration function'}

        for org in (['roots'], ['grains'], ['phloem']):
            for variable_name, variable_label in graph_variables_organs.iteritems():
                graph_name = variable_name + '_' + '_'.join(org) + '.PNG'
                tools.plot_cnwheat_ouputs(organs_output_df,
                              x_name = x_name,
                              y_name = variable_name,
                              x_label=x_label,
                              y_label=variable_label,
                              filters={'organ': org},
                              plot_filepath=os.path.join(GRAPHS_DIRPATH, graph_name),
                              explicit_label=False)

        # 3) Soil
        soil_output_df = pd.read_csv(SOILS_STATES_FILEPATH)

        fig, (ax1) = plt.subplots(1)
        conc_nitrates_soil = soil_output_df['Conc_Nitrates_Soil']*14E-6
        ax1.plot(soil_output_df['t'], conc_nitrates_soil)
        ax1.set_ylabel(u'[Nitrates] (g m$^{-3}$)')
        ax1.set_xlabel('Time from flowering (hour)')
        ax1.set_title = 'Conc Nitrates Soil'
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Conc_Nitrates_Soil.PNG'), format='PNG', bbox_inches='tight')
        plt.close()

        # 4) Hidden zones
        all_hiddenzones_inputs_outputs_df = pd.read_csv(HIDDENZONES_STATES_FILEPATH)
        graph_variables_hiddenzones = {'hiddenzone_L': u'Hidden zone length (m)','leaf_L': u'Leaf length (m)', 'delta_leaf_L':u'Delta leaf length (m)',
                                       'Conc_Sucrose':u'[Sucrose] (µmol g$^{-1}$ mstruct)', 'Conc_Amino_Acids':u'[Amino Acids] (µmol g$^{-1}$ mstruct)', 'Conc_Proteins': u'[Proteins] (g g$^{-1}$ mstruct)', 'Conc_Fructan':u'[Fructan] (µmol g$^{-1}$ mstruct)',
                                        'Unloading_Sucrose':u'Sucrose unloading (µmol C)', 'Unloading_Amino_Acids':u'Amino_acids unloading (µmol N)', 'mstruct': u'Structural mass (g)', 'Respi_growth': u'Growth respiration (µmol C)', 'sucrose_consumption_mstruct': u'Consumption of sucrose for growth (µmol C)'}

        for variable_name, variable_label in graph_variables_hiddenzones.iteritems():
            graph_name = variable_name + '_hz' + '.PNG'
            cnwheat_tools.plot_cnwheat_ouputs(all_hiddenzones_inputs_outputs_df,
                          x_name = x_name,
                          y_name = variable_name,
                          x_label = x_label,
                          y_label = variable_label,
                          filters={'plant': 1, 'axis': 'MS'},
                          plot_filepath=os.path.join(GRAPHS_DIRPATH, graph_name),
                          explicit_label=False)

        # 5) Organs
        all_organs_inputs_outputs_df = pd.read_csv(ELEMENTS_STATES_FILEPATH)
        graph_variables_organs = {'visible_length': u'Length (m)'}

        for variable_name, variable_label in graph_variables_organs.iteritems():
            graph_name = variable_name + '.PNG'
            cnwheat_tools.plot_cnwheat_ouputs(all_organs_inputs_outputs_df,
                          x_name = x_name,
                          y_name = variable_name,
                          x_label = x_label,
                          y_label = variable_label,
                          filters={'plant': 1, 'axis': 'MS'},
                          plot_filepath=os.path.join(GRAPHS_DIRPATH, graph_name),
                          explicit_label=False)

if __name__ == '__main__':
    main(200, run_simu=True, make_graphs=False)