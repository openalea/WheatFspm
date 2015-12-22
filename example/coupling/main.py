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

from alinea.adel import astk_interface

from cnwheat import simulation as cnwheat_simulation, model as cnwheat_model, parameters as cnwheat_parameters, converter as cnwheat_converter, tools as cnwheat_tools
from farquharwheat import simulation as farquharwheat_simulation, model as farquharwheat_model, converter as farquharwheat_converter
from senescwheat import simulation as senescwheat_simulation, model as senescwheat_model, converter as senescwheat_converter

INPUTS_DIRPATH = 'inputs'

# adelwheat inputs at t0
ADELWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'adelwheat') # the directory adelwheat must contain files 'adel0000.pckl' and 'scene0000.bgeom'

# cnwheat inputs at t0
CNWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'cnwheat')
CNWHEAT_ORGANS_INPUTS_FILEPATH = os.path.join(CNWHEAT_INPUTS_DIRPATH, 'organs_inputs.csv')
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

# inputs and outputs of all models for each step
OUTPUTS_DIRPATH = 'outputs'
AXES_INPUTS_OUTPUTS_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'axes_inputs_outputs.csv')
ORGANS_INPUTS_OUTPUTS_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'organs_inputs_outputs.csv')
ELEMENTS_INPUTS_OUTPUTS_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'elements_inputs_outputs.csv')
SOILS_INPUTS_OUTPUTS_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'soils_inputs_outputs.csv')

INPUTS_OUTPUTS_PRECISION = 6

LOGGING_CONFIG_FILEPATH = os.path.join('..', 'logging.json')

LOGGING_LEVEL = logging.INFO # can be one of: DEBUG, INFO, WARNING, ERROR, CRITICAL

cnwheat_tools.setup_logging(LOGGING_CONFIG_FILEPATH, LOGGING_LEVEL, log_model=False, log_compartments=False, log_derivatives=False)


meteo_df = pd.read_csv(METEO_FILEPATH, index_col='t')
    
current_time_of_the_system = time.time()

# define the time step in hours for each simulator 
senescwheat_ts = 2
farquharwheat_ts = 2
cnwheat_ts = 1

hour_to_second_conversion_factor = 3600
# create the simulators
senescwheat_simulation_ = senescwheat_simulation.Simulation(senescwheat_ts * hour_to_second_conversion_factor)
farquharwheat_simulation_ = farquharwheat_simulation.Simulation()
cnwheat_simulation_ = cnwheat_simulation.Simulation(cnwheat_ts * hour_to_second_conversion_factor)

# read adelwheat inputs at t0
adel_wheat = astk_interface.AdelWheat(seed=1234)
g = adel_wheat.load(dir=ADELWHEAT_INPUTS_DIRPATH)[0]

# read cnwheat inputs at t0
cnwheat_organs_inputs_df = pd.read_csv(CNWHEAT_ORGANS_INPUTS_FILEPATH)
cnwheat_elements_inputs_df = pd.read_csv(CNWHEAT_ELEMENTS_INPUTS_FILEPATH)
cnwheat_soils_inputs_df = pd.read_csv(CNWHEAT_SOILS_INPUTS_FILEPATH)

# read farquharwheat inputs at t0
farquharwheat_inputs_df = pd.read_csv(FARQUHARWHEAT_INPUTS_FILEPATH)

# read senescwheat inputs at t0
senescwheat_roots_inputs_df = pd.read_csv(SENESCWHEAT_ROOTS_INPUTS_FILEPATH)
senescwheat_elements_inputs_df = pd.read_csv(SENESCWHEAT_ELEMENTS_INPUTS_FILEPATH)

# define the start and the end of the whole simulation (in hours)
start_time = 0
stop_time = 8

# define lists of dataframes to store the inputs and outputs of cnwheat, farquharwheat and senescwheat at each step, and use them for postprocessings
axes_inputs_ouputs_df_list = []
organs_inputs_ouputs_df_list = []
elements_inputs_ouputs_df_list = []
soils_inputs_ouputs_df_list = []

# Initialize dataframes to store the inputs and outputs of cnwheat, farquharwheat and senescwheat at a given step.
# Since the topology is static, these dataframes keep the same shape during all the simulation.
organs_inputs_ouputs_df = pd.DataFrame(index=range(len(cnwheat_organs_inputs_df.index)), 
                                       columns=cnwheat_simulation.Simulation.ORGANS_OUTPUTS_INDEXES + sorted(set(cnwheat_simulation.Simulation.ORGANS_STATE + cnwheat_simulation.Simulation.ORGANS_INTERMEDIATE_VARIABLES + cnwheat_simulation.Simulation.ORGANS_FLUXES + cnwheat_simulation.Simulation.ORGANS_INTEGRATIVE_VARIABLES + cnwheat_simulation.Simulation.ORGANS_POSTPROCESSING_VARIABLES + senescwheat_converter.SENESCWHEAT_ROOTS_INPUTS_OUTPUTS)))
organs_inputs_indexes = cnwheat_organs_inputs_df.loc[:, cnwheat_simulation.Simulation.ORGANS_INPUTS_INDEXES]
organs_inputs_indexes.sort_index(by=cnwheat_simulation.Simulation.ORGANS_INPUTS_INDEXES, inplace=True)
organs_inputs_indexes.reset_index(drop=True, inplace=True)
organs_inputs_ouputs_df.loc[:, cnwheat_simulation.Simulation.ORGANS_INPUTS_INDEXES] = organs_inputs_indexes.loc[:, cnwheat_simulation.Simulation.ORGANS_INPUTS_INDEXES]

elements_inputs_ouputs_df = pd.DataFrame(index=cnwheat_elements_inputs_df.index, 
                                         columns=cnwheat_simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES + sorted(set(cnwheat_simulation.Simulation.ELEMENTS_STATE + cnwheat_simulation.Simulation.ELEMENTS_INTERMEDIATE_VARIABLES + cnwheat_simulation.Simulation.ELEMENTS_FLUXES + cnwheat_simulation.Simulation.ELEMENTS_INTEGRATIVE_VARIABLES + cnwheat_simulation.Simulation.ELEMENTS_POSTPROCESSING_VARIABLES + farquharwheat_converter.FARQUHARWHEAT_INPUTS_OUTPUTS + senescwheat_converter.SENESCWHEAT_ELEMENTS_INPUTS_OUTPUTS)))
elements_inputs_ouputs_df.loc[:, cnwheat_simulation.Simulation.ELEMENTS_INPUTS_INDEXES] = cnwheat_elements_inputs_df.loc[:, cnwheat_simulation.Simulation.ELEMENTS_INPUTS_INDEXES]

# run the simulators
for t_senescwheat in xrange(start_time, stop_time, senescwheat_ts):
    print 'run senescwheat at', t_senescwheat
    # initialize and run senescwheat, and update the global mtg
    senescwheat_simulation_.initialize(senescwheat_converter.from_MTG(g, senescwheat_roots_inputs_df, senescwheat_elements_inputs_df))
    senescwheat_simulation_.run()
    senescwheat_converter.update_MTG(senescwheat_simulation_.inputs, senescwheat_simulation_.outputs, g)
    # fill the global dataframes for post-processings
    senescwheat_roots_inputs_df, senescwheat_elements_inputs_df = senescwheat_converter.to_dataframes(senescwheat_simulation_.inputs)
    senescwheat_roots_outputs_df, senescwheat_elements_outputs_df = senescwheat_converter.to_dataframes(senescwheat_simulation_.outputs)
    senescwheat_roots_inputs_outputs_df = senescwheat_roots_outputs_df.combine_first(senescwheat_roots_inputs_df)
    senescwheat_elements_inputs_outputs_df = senescwheat_elements_outputs_df.combine_first(senescwheat_elements_inputs_df)
    organs_inputs_ouputs_df.loc[organs_inputs_ouputs_df.organ == 'roots', senescwheat_converter.SENESCWHEAT_ROOTS_INPUTS_OUTPUTS] = senescwheat_roots_inputs_outputs_df.loc[:, senescwheat_converter.SENESCWHEAT_ROOTS_INPUTS_OUTPUTS].values
    elements_inputs_ouputs_df.loc[:, senescwheat_converter.SENESCWHEAT_ELEMENTS_INPUTS_OUTPUTS] = senescwheat_elements_inputs_outputs_df.loc[:, senescwheat_converter.SENESCWHEAT_ELEMENTS_INPUTS_OUTPUTS].values
    
    for t_farquharwheat in xrange(t_senescwheat, t_senescwheat + senescwheat_ts, farquharwheat_ts):
        
        print 'run farquharwheat at', t_farquharwheat
        # get the meteo of the current step
        Ta, ambient_CO2, RH, Ur, PARi = meteo_df.loc[t_farquharwheat, ['air_temperature', 'ambient_CO2', 'humidity', 'Wind', 'PARi']]
        # initialize and run farquharwheat, and update the global mtg
        farquharwheat_simulation_.initialize(farquharwheat_converter.from_MTG(g, farquharwheat_inputs_df))
        farquharwheat_simulation_.run(Ta, ambient_CO2, RH, Ur, PARi)
        farquharwheat_converter.update_MTG(farquharwheat_simulation_.inputs, farquharwheat_simulation_.outputs, g)
        # fill the global dataframes for post-processings
        farquharwheat_inputs_df = farquharwheat_converter.to_dataframe(farquharwheat_simulation_.inputs)
        farquharwheat_outputs_df = farquharwheat_converter.to_dataframe(farquharwheat_simulation_.outputs)
        farquharwheat_inputs_outputs_df = farquharwheat_outputs_df.combine_first(farquharwheat_inputs_df)
        elements_inputs_ouputs_df.loc[:, farquharwheat_converter.FARQUHARWHEAT_INPUTS_OUTPUTS] = farquharwheat_inputs_outputs_df.loc[:, farquharwheat_converter.FARQUHARWHEAT_INPUTS_OUTPUTS].values
        
        for t_cnwheat in xrange(t_farquharwheat, t_farquharwheat + farquharwheat_ts, cnwheat_ts):
            
            print 'run cnwheat at', t_cnwheat
            # initialize and run cnwheat, and update the global mtg
            cnwheat_simulation_.initialize(cnwheat_converter.from_MTG(g, cnwheat_organs_inputs_df, cnwheat_elements_inputs_df),
                                           cnwheat_converter.from_dataframes(soils_inputs=cnwheat_soils_inputs_df))
            cnwheat_simulation_.run(start_time=t_cnwheat, stop_time=t_cnwheat+cnwheat_ts, number_of_output_steps=cnwheat_ts+1)
            cnwheat_converter.update_MTG(cnwheat_simulation_.population, g)
            cnwheat_soils_inputs_df = cnwheat_converter.to_dataframes(soils=cnwheat_simulation_.soils)
            # run cnwheat post-processings
            (_, 
             cnwheat_axes_postprocessing_df,
             _,
             cnwheat_organs_postprocessing_df,
             cnwheat_elements_postprocessing_df,
             cnwheat_soils_postprocessing_df) = cnwheat_simulation_.postprocessings()
            # fill the global dataframes for post-processings
            axes_inputs_ouputs_df = cnwheat_axes_postprocessing_df.loc[cnwheat_axes_postprocessing_df.t == t_cnwheat, :].reset_index(drop=True)
            cnwheat_organs_postprocessing_df = cnwheat_organs_postprocessing_df.loc[cnwheat_organs_postprocessing_df.t == t_cnwheat, :].reset_index(drop=True)
            organs_inputs_ouputs_df.loc[organs_inputs_ouputs_df.organ.isin(cnwheat_converter.MTG_TO_CNWHEAT_AXES_ORGANS_MAPPING), cnwheat_simulation.Simulation.ORGANS_ALL_VARIABLES] = cnwheat_organs_postprocessing_df.loc[:, cnwheat_simulation.Simulation.ORGANS_ALL_VARIABLES].values
            organs_inputs_ouputs_df.loc[:, cnwheat_simulation.Simulation.T_INDEX] = t_cnwheat
            cnwheat_elements_postprocessing_df = cnwheat_elements_postprocessing_df.loc[cnwheat_elements_postprocessing_df.t == t_cnwheat, :].reset_index(drop=True)
            elements_inputs_ouputs_df = cnwheat_elements_postprocessing_df.combine_first(elements_inputs_ouputs_df)
            soils_inputs_ouputs_df = cnwheat_soils_postprocessing_df.loc[cnwheat_soils_postprocessing_df.t == t_cnwheat, :].reset_index(drop=True)
            
            axes_inputs_ouputs_df_list.append(axes_inputs_ouputs_df)
            organs_inputs_ouputs_df_list.append(organs_inputs_ouputs_df)
            elements_inputs_ouputs_df_list.append(elements_inputs_ouputs_df)
            soils_inputs_ouputs_df_list.append(soils_inputs_ouputs_df)
            
            organs_inputs_ouputs_df = organs_inputs_ouputs_df.copy()
            elements_inputs_ouputs_df = elements_inputs_ouputs_df.copy()

global_axes_df = pd.concat(axes_inputs_ouputs_df_list, ignore_index=True)
axes_inputs_outputs_variable_names = sorted([column for column in global_axes_df.columns if column not in cnwheat_simulation.Simulation.AXES_OUTPUTS_INDEXES])
global_axes_df = global_axes_df.reindex_axis(cnwheat_simulation.Simulation.AXES_OUTPUTS_INDEXES + axes_inputs_outputs_variable_names, axis=1, copy=False)
global_axes_df.to_csv(AXES_INPUTS_OUTPUTS_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

global_organs_df = pd.concat(organs_inputs_ouputs_df_list, ignore_index=True)
organs_inputs_outputs_variable_names = sorted([column for column in global_organs_df.columns if column not in cnwheat_simulation.Simulation.ORGANS_OUTPUTS_INDEXES])
global_organs_df = global_organs_df.reindex_axis(cnwheat_simulation.Simulation.ORGANS_OUTPUTS_INDEXES + organs_inputs_outputs_variable_names, axis=1, copy=False)
global_organs_df.to_csv(ORGANS_INPUTS_OUTPUTS_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

global_elements_df = pd.concat(elements_inputs_ouputs_df_list, ignore_index=True)
elements_inputs_outputs_variable_names = [column for column in global_elements_df.columns if column not in cnwheat_simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES]
global_elements_df = global_elements_df.reindex_axis(cnwheat_simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES + elements_inputs_outputs_variable_names, axis=1, copy=False)
global_elements_df.to_csv(ELEMENTS_INPUTS_OUTPUTS_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

global_soils_df = pd.concat(soils_inputs_ouputs_df_list, ignore_index=True)
soils_inputs_outputs_variable_names = [column for column in global_soils_df.columns if column not in cnwheat_simulation.Simulation.SOILS_OUTPUTS_INDEXES]
global_soils_df = global_soils_df.reindex_axis(cnwheat_simulation.Simulation.SOILS_OUTPUTS_INDEXES + soils_inputs_outputs_variable_names, axis=1, copy=False)
global_soils_df.to_csv(SOILS_INPUTS_OUTPUTS_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

execution_time = int(time.time() - current_time_of_the_system)
print '\n', 'Model executed in ', str(datetime.timedelta(seconds=execution_time))

