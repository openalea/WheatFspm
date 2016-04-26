# -*- coding: latin-1 -*-

'''
    main
    ~~~~

    An example to show how to run the model CN-Wheat.

    You must first install :mod:`cnwheat` (and add it to your PYTHONPATH)
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

import pandas as pd
import numpy as np

from cnwheat import simulation, parameters, converter, model

# inputs paths
INPUTS_DIRPATH = 'inputs'

PLANTS_INPUTS_FILENAME = 'plants_inputs.csv'
AXES_INPUTS_FILENAME = 'axes_inputs.csv'
METAMERS_INPUTS_FILENAME = 'metamers_inputs.csv'
ORGANS_INPUTS_FILENAME = 'organs_inputs.csv'
HGZS_INPUTS_FILENAME = 'hgzs_inputs.csv'
ELEMENTS_INPUTS_FILENAME = 'elements_inputs.csv'
SOILS_INPUTS_FILENAME = 'soils_inputs.csv'

PHOTOSYNTHESIS_ELEMENTS_DATA_FILENAME = 'photosynthesis_elements_data.csv'
SENESCENCE_ROOTS_DATA_FILENAME = 'senescence_roots_data.csv'
SENESCENCE_ELEMENTS_DATA_FILENAME = 'senescence_elements_data.csv'

# outputs paths
OUTPUTS_DIRPATH = 'outputs'

PLANTS_OUTPUTS_FILENAME = 'plants_outputs.csv'
AXES_OUTPUTS_FILENAME = 'axes_outputs.csv'
METAMERS_OUTPUTS_FILENAME = 'metamers_outputs.csv'
ORGANS_OUTPUTS_FILENAME = 'organs_outputs.csv'
HGZS_OUTPUTS_FILENAME = 'hgzs_outputs.csv'
ELEMENTS_OUTPUTS_FILENAME = 'elements_outputs.csv'
SOILS_OUTPUTS_FILENAME = 'soils_outputs.csv'

PLANTS_OUTPUTS_FILEPATH = os.path.join(OUTPUTS_DIRPATH, PLANTS_OUTPUTS_FILENAME)
AXES_OUTPUTS_FILEPATH = os.path.join(OUTPUTS_DIRPATH, AXES_OUTPUTS_FILENAME)
METAMERS_OUTPUTS_FILEPATH = os.path.join(OUTPUTS_DIRPATH, METAMERS_OUTPUTS_FILENAME)
ORGANS_OUTPUTS_FILEPATH = os.path.join(OUTPUTS_DIRPATH, ORGANS_OUTPUTS_FILENAME)
HGZS_OUTPUTS_FILEPATH = os.path.join(OUTPUTS_DIRPATH, HGZS_OUTPUTS_FILENAME)
ELEMENTS_OUTPUTS_FILEPATH = os.path.join(OUTPUTS_DIRPATH, ELEMENTS_OUTPUTS_FILENAME)
SOILS_OUTPUTS_FILEPATH = os.path.join(OUTPUTS_DIRPATH, SOILS_OUTPUTS_FILENAME)

# precision of floats in the output CSV files
OUTPUTS_PRECISION = 6


### 1. Initialize from CSV files

# Read CN exchange inputs
inputs_dataframes = {}
for inputs_filename in (PLANTS_INPUTS_FILENAME, AXES_INPUTS_FILENAME, METAMERS_INPUTS_FILENAME, ORGANS_INPUTS_FILENAME, HGZS_INPUTS_FILENAME, ELEMENTS_INPUTS_FILENAME, SOILS_INPUTS_FILENAME):
    inputs_dataframes[inputs_filename] = pd.read_csv(os.path.join(INPUTS_DIRPATH, inputs_filename))

# Initialize a simulation from CN exchange inputs
simulation_ = simulation.Simulation(delta_t=3600)
population, soils = converter.from_dataframes(inputs_dataframes[PLANTS_INPUTS_FILENAME],
                                              inputs_dataframes[AXES_INPUTS_FILENAME],
                                              inputs_dataframes[METAMERS_INPUTS_FILENAME],
                                              inputs_dataframes[ORGANS_INPUTS_FILENAME],
                                              inputs_dataframes[HGZS_INPUTS_FILENAME],
                                              inputs_dataframes[ELEMENTS_INPUTS_FILENAME],
                                              inputs_dataframes[SOILS_INPUTS_FILENAME])
simulation_.initialize(population, soils)

# Read photosynthesis inputs
photosynthesis_elements_data_filepath = os.path.join(INPUTS_DIRPATH, PHOTOSYNTHESIS_ELEMENTS_DATA_FILENAME)
photosynthesis_elements_data_df = pd.read_csv(photosynthesis_elements_data_filepath)
photosynthesis_elements_data_grouped = photosynthesis_elements_data_df.groupby(simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES)

# Read senescence and growth inputs
senescence_roots_data_filepath = os.path.join(INPUTS_DIRPATH, SENESCENCE_ROOTS_DATA_FILENAME)
senescence_roots_data_df = pd.read_csv(senescence_roots_data_filepath)
senescence_roots_data_grouped = senescence_roots_data_df.groupby(simulation.Simulation.AXES_OUTPUTS_INDEXES)
senescence_elements_data_filepath = os.path.join(INPUTS_DIRPATH, SENESCENCE_ELEMENTS_DATA_FILENAME)
senescence_elements_data_df = pd.read_csv(senescence_elements_data_filepath)
senescence_elements_data_grouped = senescence_elements_data_df.groupby(simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES)


### 2. Run the simulation

# Define the time grid
start_time = 0
stop_time = 2
time_step = 1

# Initialize the lists of outputs. There is one list per type of outputs.
all_plants_df_list = []
all_axes_df_list = []
all_metamers_df_list = []
all_organs_df_list = []
all_hgzs_df_list = []
all_elements_df_list = []
all_soils_df_list = []

# Start the execution loop
for t in xrange(start_time, stop_time, time_step):
    # force the senescence and photosynthesis of the population
    for plant in simulation_.population.plants:
        for axis in plant.axes:
            group = senescence_roots_data_grouped.get_group((t, plant.index, axis.label))
            senescence_data_to_use = group.loc[group.first_valid_index(), simulation.Simulation.ORGANS_STATE].dropna().to_dict()
            axis.roots.__dict__.update(senescence_data_to_use)
            for phytomer in axis.phytomers:
                for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                    if organ is None:
                        continue
                    for element in (organ.exposed_element, organ.enclosed_element):
                        if element is None:
                            continue
                        group_senesc = senescence_elements_data_grouped.get_group((t, plant.index, axis.label, phytomer.index, organ.label, element.label))
                        senescence_data_to_use = group_senesc.loc[group_senesc.first_valid_index(), simulation.Simulation.ELEMENTS_STATE].dropna().to_dict()
                        element.__dict__.update(senescence_data_to_use)
                        group_photo = photosynthesis_elements_data_grouped.get_group((t, plant.index, axis.label, phytomer.index, organ.label, element.label))
                        photosynthesis_elements_data_to_use = group_photo.loc[group_photo.first_valid_index(), simulation.Simulation.ELEMENTS_STATE].dropna().to_dict()
                        element.__dict__.update(photosynthesis_elements_data_to_use)

    # run the model of CN exchanges ; the population is internally updated by the model
    simulation_.run(start_time=t, stop_time=t+time_step, number_of_output_steps=time_step+1)

    # run post-processings
    all_plants_df, all_axes_df, all_metamers_df, all_organs_df, all_hgzs_df, all_elements_df, all_soils_df = simulation_.postprocessings()

    all_plants_df_list.append(all_plants_df)
    all_axes_df_list.append(all_axes_df)
    all_metamers_df_list.append(all_metamers_df)
    all_organs_df_list.append(all_organs_df)
    all_hgzs_df_list.append(all_hgzs_df)
    all_elements_df_list.append(all_elements_df)
    all_soils_df_list.append(all_soils_df)


### 3. Write the outputs to CSV files

global_plants_df = pd.concat(all_plants_df_list, ignore_index=True)
global_plants_df.drop_duplicates(subset=simulation.Simulation.PLANTS_OUTPUTS_INDEXES, inplace=True)
global_plants_df.to_csv(PLANTS_OUTPUTS_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))

global_axes_df = pd.concat(all_axes_df_list, ignore_index=True)
global_axes_df.drop_duplicates(subset=simulation.Simulation.AXES_OUTPUTS_INDEXES, inplace=True)
global_axes_df.to_csv(AXES_OUTPUTS_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))

global_metamers_df = pd.concat(all_metamers_df_list, ignore_index=True)
global_metamers_df.drop_duplicates(subset=simulation.Simulation.PHYTOMERS_OUTPUTS_INDEXES, inplace=True)
global_metamers_df.to_csv(METAMERS_OUTPUTS_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))

global_organs_df = pd.concat(all_organs_df_list, ignore_index=True)
global_organs_df.drop_duplicates(subset=simulation.Simulation.ORGANS_OUTPUTS_INDEXES, inplace=True)
global_organs_df.to_csv(ORGANS_OUTPUTS_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))

global_hgzs_df = pd.concat(all_hgzs_df_list, ignore_index=True)
global_hgzs_df.drop_duplicates(subset=simulation.Simulation.HIDDENGROWINGZONE_OUTPUTS_INDEXES, inplace=True)
global_hgzs_df.to_csv(HGZS_OUTPUTS_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))

global_elements_df = pd.concat(all_elements_df_list, ignore_index=True)
global_elements_df.drop_duplicates(subset=simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES, inplace=True)
global_elements_df.to_csv(ELEMENTS_OUTPUTS_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))

global_soils_df = pd.concat(all_soils_df_list, ignore_index=True)
global_soils_df.drop_duplicates(subset=simulation.Simulation.SOILS_OUTPUTS_INDEXES, inplace=True)
global_soils_df.to_csv(SOILS_OUTPUTS_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))