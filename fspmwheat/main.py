# -*- coding: latin-1 -*-

'''
    main
    ~~~~

    An example to show how to couple models CN-Wheat, Farquhar-Wheat and Senesc-Wheat using a static topology from Adel-Wheat.
    This example uses the format MTG to exchange data between the models.

    You must first install :mod:`alinea.adel`, :mod:`cnwheat`, :mod:`farquharwheat` and :mod:`senescwheat` (and add them to your PYTHONPATH)
    before running this script with the command `python`.

    :copyright: Copyright 2014-2016 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2016.
'''

'''
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
'''

import datetime
import logging
import os
import random
import time

import numpy as np
import pandas as pd

import caribu_facade
import cnwheat_facade
import elongwheat_facade
import farquharwheat_facade
import growthwheat_facade
import senescwheat_facade

from alinea.adel.adel_dynamic import AdelWheatDyn

random.seed(1234)
np.random.seed(1234)

# number of seconds in 1 hour
HOUR_TO_SECOND_CONVERSION_FACTOR = 3600

INPUTS_DIRPATH = 'inputs'

# adelwheat inputs at t0
ADELWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'adelwheat') # the directory adelwheat must contain files 'adel0000.pckl' and 'scene0000.bgeom'

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
FARQUHARWHEAT_ELEMENTS_INPUTS_FILEPATH = os.path.join(FARQUHARWHEAT_INPUTS_DIRPATH, 'elements_inputs.csv')
FARQUHARWHEAT_SAMS_INPUTS_FILEPATH = os.path.join(FARQUHARWHEAT_INPUTS_DIRPATH, 'SAMs_inputs.csv')
METEO_FILEPATH = os.path.join(FARQUHARWHEAT_INPUTS_DIRPATH, 'meteo_Clermont_rebuild.csv')

# senescwheat inputs at t0
SENESCWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'senescwheat')
SENESCWHEAT_ROOTS_INPUTS_FILEPATH = os.path.join(SENESCWHEAT_INPUTS_DIRPATH, 'roots_inputs.csv')
SENESCWHEAT_ELEMENTS_INPUTS_FILEPATH = os.path.join(SENESCWHEAT_INPUTS_DIRPATH, 'elements_inputs.csv')

# elongwheat inputs at t0
ELONGWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'elongwheat')
ELONGWHEAT_HIDDENZONE_INPUTS_FILEPATH = os.path.join(ELONGWHEAT_INPUTS_DIRPATH, 'hiddenzones_inputs.csv')
ELONGWHEAT_ELEMENTS_INPUTS_FILEPATH = os.path.join(ELONGWHEAT_INPUTS_DIRPATH, 'elements_inputs.csv')
ELONGWHEAT_SAM_INPUTS_FILEPATH = os.path.join(ELONGWHEAT_INPUTS_DIRPATH, 'SAM_inputs.csv')

# growthwheat inputs at t0
GROWTHWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'growthwheat')
GROWTHWHEAT_HIDDENZONE_INPUTS_FILEPATH = os.path.join(GROWTHWHEAT_INPUTS_DIRPATH, 'hiddenzones_inputs.csv')
GROWTHWHEAT_ELEMENTS_INPUTS_FILEPATH = os.path.join(GROWTHWHEAT_INPUTS_DIRPATH, 'elements_inputs.csv')
GROWTHWHEAT_ROOTS_INPUTS_FILEPATH = os.path.join(GROWTHWHEAT_INPUTS_DIRPATH, 'roots_inputs.csv')

# the path of the CSV files where to save the states of the modeled system at each step
OUTPUTS_DIRPATH = 'outputs'
AXES_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'axes_states.csv')
SAM_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'SAM_states.csv')
ORGANS_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'organs_states.csv')
HIDDENZONES_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'hiddenzones_states.csv')
ELEMENTS_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'elements_states.csv')
SOILS_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'soils_states.csv')

# post-processing directory path
POSTPROCESSING_DIRPATH = 'postprocessing'
AXES_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'axes_postprocessing.csv')
ORGANS_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'organs_postprocessing.csv')
HIDDENZONES_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'hiddenzones_postprocessing.csv')
ELEMENTS_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'elements_postprocessing.csv')
SOILS_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'soils_postprocessing.csv')

GRAPHS_DIRPATH = 'graphs'

AXES_INDEX_COLUMNS = ['t','plant','axis']
ELEMENTS_INDEX_COLUMNS = ['t','plant','axis', 'metamer', 'organ', 'element']
HIDDENZONES_INDEX_COLUMNS = ['t','plant','axis', 'metamer']
ORGANS_INDEX_COLUMNS = ['t','plant','axis', 'organ']
SAM_INDEX_COLUMNS = ['t','plant','axis']
SOILS_INDEX_COLUMNS = ['t','plant','axis']

# Define culm density (culm m-2)
CULM_DENSITY = {1:410}

INPUTS_OUTPUTS_PRECISION = 5 # 10

LOGGING_CONFIG_FILEPATH = os.path.join('..', '..', 'logging.json')

LOGGING_LEVEL = logging.INFO # can be one of: DEBUG, INFO, WARNING, ERROR, CRITICAL

def main(stop_time, run_simu=True, run_postprocessing=True, generate_graphs=True):
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

        # read adelwheat inputs at t0
        adel_wheat = AdelWheatDyn(seed=1234, convUnit=1)
        g = adel_wheat.load(dir=ADELWHEAT_INPUTS_DIRPATH)
        adel_wheat.domain = g.get_vertex_property(0)['meta']['domain'] # temp (until Christian's commit)
        adel_wheat.nplants = g.get_vertex_property(0)['meta']['nplants'] # temp (until Christian's commit)


        # create empty dataframes to shared data between the models
        shared_axes_inputs_outputs_df = pd.DataFrame()
        shared_SAM_inputs_outputs_df = pd.DataFrame()
        shared_organs_inputs_outputs_df = pd.DataFrame()
        shared_hiddenzones_inputs_outputs_df = pd.DataFrame()
        shared_elements_inputs_outputs_df = pd.DataFrame()
        shared_soils_inputs_outputs_df = pd.DataFrame()

        # read the inputs at t0 and create the facades
        # caribu
        caribu_facade_ = caribu_facade.CaribuFacade(g,
                                                    shared_elements_inputs_outputs_df,
                                                    adel_wheat)

        # senescwheat
        senescwheat_roots_inputs_t0 = pd.read_csv(SENESCWHEAT_ROOTS_INPUTS_FILEPATH)
        senescwheat_elements_inputs_t0 = pd.read_csv(SENESCWHEAT_ELEMENTS_INPUTS_FILEPATH)
        senescwheat_facade_ = senescwheat_facade.SenescWheatFacade(g,
                                                                   senescwheat_ts * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                                   senescwheat_roots_inputs_t0,
                                                                   senescwheat_elements_inputs_t0,
                                                                   shared_organs_inputs_outputs_df,
                                                                   shared_elements_inputs_outputs_df)

        # farquharwheat
        farquharwheat_elements_inputs_t0 = pd.read_csv(FARQUHARWHEAT_ELEMENTS_INPUTS_FILEPATH)
        farquharwheat_SAMs_inputs_t0 = pd.read_csv(FARQUHARWHEAT_SAMS_INPUTS_FILEPATH)
        farquharwheat_facade_ = farquharwheat_facade.FarquharWheatFacade(g,
                                                                         farquharwheat_elements_inputs_t0,
                                                                         farquharwheat_SAMs_inputs_t0,
                                                                         shared_elements_inputs_outputs_df)

        # elongwheat
        elongwheat_hiddenzones_inputs_t0 = pd.read_csv(ELONGWHEAT_HIDDENZONE_INPUTS_FILEPATH)
        elongwheat_element_inputs_t0 = pd.read_csv(ELONGWHEAT_ELEMENTS_INPUTS_FILEPATH)
        elongwheat_SAM_inputs_t0 = pd.read_csv(ELONGWHEAT_SAM_INPUTS_FILEPATH)
        elongwheat_facade_ = elongwheat_facade.ElongWheatFacade(g,
                                                                elongwheat_ts * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                                elongwheat_SAM_inputs_t0,
                                                                elongwheat_hiddenzones_inputs_t0,
                                                                elongwheat_element_inputs_t0,
                                                                shared_SAM_inputs_outputs_df,
                                                                shared_hiddenzones_inputs_outputs_df,
                                                                shared_elements_inputs_outputs_df,
                                                                adel_wheat)

        # growthwheat
        growthwheat_hiddenzones_inputs_t0 = pd.read_csv(GROWTHWHEAT_HIDDENZONE_INPUTS_FILEPATH)
        growthwheat_element_inputs_t0 = pd.read_csv(GROWTHWHEAT_ELEMENTS_INPUTS_FILEPATH)
        growthwheat_root_inputs_t0 = pd.read_csv(GROWTHWHEAT_ROOTS_INPUTS_FILEPATH)
        growthwheat_facade_ = growthwheat_facade.GrowthWheatFacade(g,
                                                                   growthwheat_ts * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                                   growthwheat_hiddenzones_inputs_t0,
                                                                   growthwheat_element_inputs_t0,
                                                                   growthwheat_root_inputs_t0,
                                                                   shared_organs_inputs_outputs_df,
                                                                   shared_hiddenzones_inputs_outputs_df,
                                                                   shared_elements_inputs_outputs_df)

        # cnwheat
        cnwheat_organs_inputs_t0 = pd.read_csv(CNWHEAT_ORGANS_INPUTS_FILEPATH)
        cnwheat_hiddenzones_inputs_t0 = pd.read_csv(CNWHEAT_HIDDENZONE_INPUTS_FILEPATH)
        cnwheat_elements_inputs_t0 = pd.read_csv(CNWHEAT_ELEMENTS_INPUTS_FILEPATH)
        cnwheat_soils_inputs_t0 = pd.read_csv(CNWHEAT_SOILS_INPUTS_FILEPATH)
        cnwheat_facade_ = cnwheat_facade.CNWheatFacade(g,
                                                       cnwheat_ts * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                       CULM_DENSITY,
                                                       cnwheat_organs_inputs_t0,
                                                       cnwheat_hiddenzones_inputs_t0,
                                                       cnwheat_elements_inputs_t0,
                                                       cnwheat_soils_inputs_t0,
                                                       shared_axes_inputs_outputs_df,
                                                       shared_organs_inputs_outputs_df,
                                                       shared_hiddenzones_inputs_outputs_df,
                                                       shared_elements_inputs_outputs_df,
                                                       shared_soils_inputs_outputs_df)

        ##        # Update geometry
        adel_wheat.update_geometry(g)#, SI_units=True, properties_to_convert=properties_to_convert) # Returns mtg with non-SI units
        ####        adel_wheat.convert_to_SI_units(g, properties_to_convert)
        adel_wheat.plot(g)

        # define the start and the end of the whole simulation (in hours)
        start_time = 0
        stop_time = stop_time

        # define lists of dataframes to store the inputs and the outputs of the models at each step.
        axes_all_data_list = []
        organs_all_data_list = [] # organs which belong to axes: roots, phloem, grains
        hiddenzones_all_data_list = []
        elements_all_data_list = []
        SAM_all_data_list = []
        soils_all_data_list = []

        all_simulation_steps = [] # to store the steps of the simulation

        print('Simulation starts')

        # run the simulators
        current_time_of_the_system = time.time()

        for t_caribu in xrange(start_time, stop_time, caribu_ts):

            # run Caribu
            PARi = meteo.loc[t_caribu, ['PARi']].iloc[0]
            caribu_facade_.run(energy=PARi)
            print('t caribu is {}'.format(t_caribu))

            for t_senescwheat in xrange(t_caribu, t_caribu + caribu_ts, senescwheat_ts):

                # run SenescWheat
                print('t senescwheat is {}'.format(t_senescwheat))
                senescwheat_facade_.run()

                for t_farquharwheat in xrange(t_senescwheat, t_senescwheat + senescwheat_ts, farquharwheat_ts):
                    # get the meteo of the current step
                    Ta, ambient_CO2, RH, Ur = meteo.loc[t_farquharwheat, ['air_temperature', 'ambient_CO2', 'humidity', 'Wind']]

                    # run FarquharWheat
                    farquharwheat_facade_.run(Ta, ambient_CO2, RH, Ur)

                    for t_elongwheat in xrange(t_farquharwheat, t_farquharwheat + farquharwheat_ts, elongwheat_ts):
                        print('t elongwheat is {}'.format(t_elongwheat))

                        # run ElongWheat
                        Ta, Tsol = meteo.loc[t_elongwheat, ['air_temperature', 'air_temperature']] # TODO: Add soil temperature in the weather input file
                        elongwheat_facade_.run(Ta, Tsol)

                        # Update geometry
                        adel_wheat.update_geometry(g)#, SI_units=True, properties_to_convert=properties_to_convert) # Return mtg with non-SI units
                        adel_wheat.plot(g)
                        ##                        adel_wheat.convert_to_SI_units(g, properties_to_convert)

                        for t_growthwheat in xrange(t_elongwheat, t_elongwheat + elongwheat_ts, growthwheat_ts):

                            # run GrowthWheat
                            print('t growthwheat is {}'.format(t_growthwheat))
                            growthwheat_facade_.run()

                            for t_cnwheat in xrange(t_growthwheat, t_growthwheat + growthwheat_ts, cnwheat_ts):
                                if t_cnwheat > 0:

                                    # run CNWheat
                                    print('t cnwheat is {}'.format(t_cnwheat))
                                    cnwheat_facade_.run()

                                # append the inputs and outputs at current step to global lists
                                all_simulation_steps.append(t_cnwheat)
                                axes_all_data_list.append(shared_axes_inputs_outputs_df.copy())
                                organs_all_data_list.append(shared_organs_inputs_outputs_df.copy())
                                hiddenzones_all_data_list.append(shared_hiddenzones_inputs_outputs_df.copy())
                                elements_all_data_list.append(shared_elements_inputs_outputs_df.copy())
                                SAM_all_data_list.append(shared_SAM_inputs_outputs_df.copy())
                                soils_all_data_list.append(shared_soils_inputs_outputs_df.copy())


        execution_time = int(time.time() - current_time_of_the_system)
        print '\n', 'Simulation run in ', str(datetime.timedelta(seconds=execution_time))

        # write all inputs and outputs to CSV files
        all_axes_inputs_outputs = pd.concat(axes_all_data_list, keys=all_simulation_steps)
        all_axes_inputs_outputs.reset_index(0, inplace=True)
        all_axes_inputs_outputs.rename_axis({'level_0': 't'}, axis=1, inplace=True)
        all_axes_inputs_outputs = all_axes_inputs_outputs.reindex_axis(AXES_INDEX_COLUMNS+all_axes_inputs_outputs.columns.difference(AXES_INDEX_COLUMNS).tolist(), axis=1, copy=False)
        all_axes_inputs_outputs.to_csv(AXES_STATES_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

        all_organs_inputs_outputs = pd.concat(organs_all_data_list, keys=all_simulation_steps)
        all_organs_inputs_outputs.reset_index(0, inplace=True)
        all_organs_inputs_outputs.rename_axis({'level_0': 't'}, axis=1, inplace=True)
        all_organs_inputs_outputs = all_organs_inputs_outputs.reindex_axis(ORGANS_INDEX_COLUMNS+all_organs_inputs_outputs.columns.difference(ORGANS_INDEX_COLUMNS).tolist(), axis=1, copy=False)
        all_organs_inputs_outputs.to_csv(ORGANS_STATES_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

        all_hiddenzones_inputs_outputs = pd.concat(hiddenzones_all_data_list, keys=all_simulation_steps)
        all_hiddenzones_inputs_outputs.reset_index(0, inplace=True)
        all_hiddenzones_inputs_outputs.rename_axis({'level_0': 't'}, axis=1, inplace=True)
        all_hiddenzones_inputs_outputs = all_hiddenzones_inputs_outputs.reindex_axis(HIDDENZONES_INDEX_COLUMNS+all_hiddenzones_inputs_outputs.columns.difference(HIDDENZONES_INDEX_COLUMNS).tolist(), axis=1, copy=False)
        all_hiddenzones_inputs_outputs.to_csv(HIDDENZONES_STATES_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

        all_SAM_inputs_outputs = pd.concat(SAM_all_data_list, keys=all_simulation_steps)
        all_SAM_inputs_outputs.reset_index(0, inplace=True)
        all_SAM_inputs_outputs.rename_axis({'level_0': 't'}, axis=1, inplace=True)
        all_SAM_inputs_outputs = all_SAM_inputs_outputs.reindex_axis(SAM_INDEX_COLUMNS+all_SAM_inputs_outputs.columns.difference(SAM_INDEX_COLUMNS).tolist(), axis=1, copy=False)
        all_SAM_inputs_outputs.to_csv(SAM_STATES_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

        all_elements_inputs_outputs = pd.concat(elements_all_data_list, keys=all_simulation_steps)
        all_elements_inputs_outputs = all_elements_inputs_outputs.loc[(all_elements_inputs_outputs.plant == 1) & # TODO: temporary ; to remove when there will be default input values for each element in the mtg
                                                                      (all_elements_inputs_outputs.axis == 'MS')]
        all_elements_inputs_outputs.reset_index(0, inplace=True)
        all_elements_inputs_outputs.rename_axis({'level_0': 't'}, axis=1, inplace=True)
        all_elements_inputs_outputs = all_elements_inputs_outputs.reindex_axis(ELEMENTS_INDEX_COLUMNS+all_elements_inputs_outputs.columns.difference(ELEMENTS_INDEX_COLUMNS).tolist(), axis=1, copy=False)
        all_elements_inputs_outputs.to_csv(ELEMENTS_STATES_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

        all_soils_inputs_outputs = pd.concat(soils_all_data_list, keys=all_simulation_steps)
        all_soils_inputs_outputs.reset_index(0, inplace=True)
        all_soils_inputs_outputs.rename_axis({'level_0': 't'}, axis=1, inplace=True)
        all_soils_inputs_outputs = all_soils_inputs_outputs.reindex_axis(SOILS_INDEX_COLUMNS+all_soils_inputs_outputs.columns.difference(SOILS_INDEX_COLUMNS).tolist(), axis=1, copy=False)
        all_soils_inputs_outputs.to_csv(SOILS_STATES_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))


    if run_postprocessing:

        # cnwheat postprocessing only

        # Retrieve outputs dataframes from precedent simulation run
        states_df_dict = {}
        for states_filepath in (AXES_STATES_FILEPATH,
                                ORGANS_STATES_FILEPATH,
                                HIDDENZONES_STATES_FILEPATH,
                                ELEMENTS_STATES_FILEPATH,
                                SOILS_STATES_FILEPATH):
            states_df = pd.read_csv(states_filepath)
            states_file_basename = os.path.basename(states_filepath).split('.')[0]
            states_df_dict[states_file_basename] = states_df
        time_grid = states_df_dict.values()[0].t
        delta_t = (time_grid.loc[1] - time_grid.loc[0]) * HOUR_TO_SECOND_CONVERSION_FACTOR

        # run the postprocessing
        axes_postprocessing_file_basename = os.path.basename(AXES_POSTPROCESSING_FILEPATH).split('.')[0]
        hiddenzones_postprocessing_file_basename = os.path.basename(HIDDENZONES_POSTPROCESSING_FILEPATH).split('.')[0]
        organs_postprocessing_file_basename = os.path.basename(ORGANS_POSTPROCESSING_FILEPATH).split('.')[0]
        elements_postprocessing_file_basename = os.path.basename(ELEMENTS_POSTPROCESSING_FILEPATH).split('.')[0]
        soils_postprocessing_file_basename = os.path.basename(SOILS_POSTPROCESSING_FILEPATH).split('.')[0]
        postprocessing_df_dict = {}
        (postprocessing_df_dict[axes_postprocessing_file_basename],
         postprocessing_df_dict[hiddenzones_postprocessing_file_basename],
         postprocessing_df_dict[organs_postprocessing_file_basename],
         postprocessing_df_dict[elements_postprocessing_file_basename],
         postprocessing_df_dict[soils_postprocessing_file_basename]) \
            = cnwheat_facade.CNWheatFacade.postprocessing(axes_outputs_df=states_df_dict[os.path.basename(AXES_STATES_FILEPATH).split('.')[0]],
                                                          hiddenzone_outputs_df=states_df_dict[os.path.basename(HIDDENZONES_STATES_FILEPATH).split('.')[0]],
                                                          organs_outputs_df=states_df_dict[os.path.basename(ORGANS_STATES_FILEPATH).split('.')[0]],
                                                          elements_outputs_df=states_df_dict[os.path.basename(ELEMENTS_STATES_FILEPATH).split('.')[0]],
                                                          soils_outputs_df=states_df_dict[os.path.basename(SOILS_STATES_FILEPATH).split('.')[0]],
                                                          delta_t=delta_t)

        # save the postprocessing to disk
        for postprocessing_file_basename, postprocessing_filepath in ((axes_postprocessing_file_basename, AXES_POSTPROCESSING_FILEPATH),
                                                                      (hiddenzones_postprocessing_file_basename, HIDDENZONES_POSTPROCESSING_FILEPATH),
                                                                      (organs_postprocessing_file_basename, ORGANS_POSTPROCESSING_FILEPATH),
                                                                      (elements_postprocessing_file_basename, ELEMENTS_POSTPROCESSING_FILEPATH),
                                                                      (soils_postprocessing_file_basename, SOILS_POSTPROCESSING_FILEPATH)):
            postprocessing_df_dict[postprocessing_file_basename].to_csv(postprocessing_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))


    if generate_graphs:

        # Retrieve last computed post-processing dataframes
        organs_postprocessing_file_basename = os.path.basename(ORGANS_POSTPROCESSING_FILEPATH).split('.')[0]
        hiddenzones_postprocessing_file_basename = os.path.basename(HIDDENZONES_POSTPROCESSING_FILEPATH).split('.')[0]
        elements_postprocessing_file_basename = os.path.basename(ELEMENTS_POSTPROCESSING_FILEPATH).split('.')[0]
        soils_postprocessing_file_basename = os.path.basename(SOILS_POSTPROCESSING_FILEPATH).split('.')[0]
        postprocessing_df_dict = {}
        for (postprocessing_filepath, postprocessing_file_basename) in ((ORGANS_POSTPROCESSING_FILEPATH, organs_postprocessing_file_basename),
                                                                        (HIDDENZONES_POSTPROCESSING_FILEPATH, hiddenzones_postprocessing_file_basename),
                                                                        (ELEMENTS_POSTPROCESSING_FILEPATH, elements_postprocessing_file_basename),
                                                                        (SOILS_POSTPROCESSING_FILEPATH, soils_postprocessing_file_basename)):
            postprocessing_df = pd.read_csv(postprocessing_filepath)
            postprocessing_df_dict[postprocessing_file_basename] = postprocessing_df

        # Generate graphs
        cnwheat_facade.CNWheatFacade.graphs(hiddenzones_postprocessing_df=postprocessing_df_dict[hiddenzones_postprocessing_file_basename],
                                            organs_postprocessing_df=postprocessing_df_dict[organs_postprocessing_file_basename],
                                            elements_postprocessing_df=postprocessing_df_dict[elements_postprocessing_file_basename],
                                            soils_postprocessing_df=postprocessing_df_dict[soils_postprocessing_file_basename],
                                            graphs_dirpath=GRAPHS_DIRPATH)


if __name__ == '__main__':
    main(500, run_simu=True, run_postprocessing=True, generate_graphs=True)
