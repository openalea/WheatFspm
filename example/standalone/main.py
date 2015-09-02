# -*- coding: latin-1 -*-

'''
    main
    ~~~~

    An example to show how to run the model CN-Wheat.

    You must first install :mod:`cnwheat` (and add it to your PYTHONPATH)
    before running this script with the command `python`.

    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2014.
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

from cnwheat import simulation, parameters, converter

# inputs paths
INPUTS_DIRPATH = 'inputs'

PLANTS_INPUTS_FILENAME = 'plants_inputs.csv'
AXES_INPUTS_FILENAME = 'axes_inputs.csv'
PHYTOMERS_INPUTS_FILENAME = 'phytomers_inputs.csv'
ORGANS_INPUTS_FILENAME = 'organs_inputs.csv'
ELEMENTS_INPUTS_FILENAME = 'elements_inputs.csv'

PHOTOSYNTHESIS_DATA_FILENAME = 'photosynthesis_data.csv'
SENESCENCE_DATA_FILENAME = 'senescence_data.csv'

# outputs paths
OUTPUTS_DIRPATH = 'outputs'

PLANTS_OUTPUTS_FILENAME = 'plants_outputs.csv'
AXES_OUTPUTS_FILENAME = 'axes_outputs.csv'
PHYTOMERS_OUTPUTS_FILENAME = 'phytomers_outputs.csv'
ORGANS_OUTPUTS_FILENAME = 'organs_outputs.csv'
ELEMENTS_OUTPUTS_FILENAME = 'elements_outputs.csv'

plants_outputs_filepath = os.path.join(OUTPUTS_DIRPATH, PLANTS_OUTPUTS_FILENAME)
axes_outputs_filepath = os.path.join(OUTPUTS_DIRPATH, AXES_OUTPUTS_FILENAME)
phytomers_outputs_filepath = os.path.join(OUTPUTS_DIRPATH, PHYTOMERS_OUTPUTS_FILENAME)
organs_outputs_filepath = os.path.join(OUTPUTS_DIRPATH, ORGANS_OUTPUTS_FILENAME)
elements_outputs_filepath = os.path.join(OUTPUTS_DIRPATH, ELEMENTS_OUTPUTS_FILENAME)

# precision of floats in the output CSV files
OUTPUTS_PRECISION = 6


### 1. Initialize from CSV files

# Read CN exchange inputs
inputs_dataframes = {}
for inputs_filename in (PLANTS_INPUTS_FILENAME, AXES_INPUTS_FILENAME, PHYTOMERS_INPUTS_FILENAME, ORGANS_INPUTS_FILENAME, ELEMENTS_INPUTS_FILENAME):
    inputs_dataframes[inputs_filename] = pd.read_csv(os.path.join(INPUTS_DIRPATH, inputs_filename))
    
# Initialize a simulation from CN exchange inputs
simulation_ = simulation.Simulation()
population = converter.from_dataframes(inputs_dataframes[PLANTS_INPUTS_FILENAME], 
                                       inputs_dataframes[AXES_INPUTS_FILENAME],
                                       inputs_dataframes[PHYTOMERS_INPUTS_FILENAME],
                                       inputs_dataframes[ORGANS_INPUTS_FILENAME],
                                       inputs_dataframes[ELEMENTS_INPUTS_FILENAME])
simulation_.initialize(population)

# Read photosynthesis inputs
photosynthesis_data_filepath = os.path.join(INPUTS_DIRPATH, PHOTOSYNTHESIS_DATA_FILENAME)
photosynthesis_data_df = pd.read_csv(photosynthesis_data_filepath)
photosynthesis_data_grouped = photosynthesis_data_df.groupby(simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES)

# Read senescence and growth inputs
senescence_data_filepath = os.path.join(INPUTS_DIRPATH, SENESCENCE_DATA_FILENAME)
senescence_data_df = pd.read_csv(senescence_data_filepath)
senescence_data_grouped = senescence_data_df.groupby(simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES)


### 2. Run the simulation

# Define the time grid
start_time = 0
stop_time = 48
timestep = 1

# Initialize the lists of outputs. There is one list per type of outputs.
all_plants_df_list = []
all_axes_df_list = []
all_phytomers_df_list = []
all_organs_df_list = []
all_elements_df_list = []

# Start the execution loop
for t in xrange(start_time, stop_time, timestep):
    # Update the population
    for plant in simulation_.population.plants:
        plant_index = plant.index
        for axis in plant.axes:
            axis_id = axis.id

            # Update the state of roots from senescence and growth data
            group = senescence_data_grouped.get_group((t, plant_index, axis_id, 0, 'Roots', 'enclosed'))
            senescence_data_to_use = group.loc[group.first_valid_index(), simulation.Simulation.ORGANS_STATE].dropna().to_dict()
            axis.roots.__dict__.update(senescence_data_to_use)

            # Update the state of each photosynthetic element from senescence and photosynthesis data
            for phytomer in axis.phytomers:
                phytomer_index = phytomer.index
                for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                    if organ is None:
                        continue
                    organ_type = organ.__class__.__name__
                    for element, element_type in ((organ.exposed_element, 'exposed'), (organ.enclosed_element, 'enclosed')):
                        if element is None:
                            continue
                        # Senescence
                        group_senesc = senescence_data_grouped.get_group((t, plant_index, axis_id, phytomer_index, organ_type, element_type))
                        senescence_data_to_use = group_senesc.loc[group_senesc.first_valid_index(), simulation.Simulation.ELEMENTS_STATE].dropna().to_dict()
                        element.__dict__.update(senescence_data_to_use)
                        # Photosynthesis
                        group_photo = photosynthesis_data_grouped.get_group((t, plant_index, axis_id, phytomer_index, organ_type, element_type))
                        photosynthesis_data_to_use = group_photo.loc[group_photo.first_valid_index(), simulation.Simulation.ELEMENTS_STATE].dropna().to_dict()
                        element.__dict__.update(photosynthesis_data_to_use)

    # Run the model ; the population is internally updated by the model
    simulation_.run(start_time=t, stop_time=t+timestep, number_of_output_steps=timestep+1)
    
    # Run post-processings
    all_plants_df, all_axes_df, all_phytomers_df, all_organs_df, all_elements_df = simulation_.postprocessings()
    
    # Update the lists of outputs 
    all_plants_df_list.append(all_plants_df)
    all_axes_df_list.append(all_axes_df)
    all_phytomers_df_list.append(all_phytomers_df)
    all_organs_df_list.append(all_organs_df)
    all_elements_df_list.append(all_elements_df)


### 3. Write the outputs to CSV files

global_plants_df = pd.concat(all_plants_df_list, ignore_index=True)
global_plants_df.drop_duplicates(subset=simulation.Simulation.PLANTS_OUTPUTS_INDEXES, inplace=True)
global_plants_df.to_csv(plants_outputs_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))

global_axes_df = pd.concat(all_axes_df_list, ignore_index=True)
global_axes_df.drop_duplicates(subset=simulation.Simulation.AXES_OUTPUTS_INDEXES, inplace=True)
global_axes_df.to_csv(axes_outputs_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))

global_phytomers_df = pd.concat(all_phytomers_df_list, ignore_index=True)
global_phytomers_df.drop_duplicates(subset=simulation.Simulation.PHYTOMERS_OUTPUTS_INDEXES, inplace=True)
global_phytomers_df.to_csv(phytomers_outputs_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))

global_organs_df = pd.concat(all_organs_df_list, ignore_index=True)
global_organs_df.drop_duplicates(subset=simulation.Simulation.ORGANS_OUTPUTS_INDEXES, inplace=True)
global_organs_df.to_csv(organs_outputs_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))

global_elements_df = pd.concat(all_elements_df_list, ignore_index=True)
global_elements_df.drop_duplicates(subset=simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES, inplace=True)
global_elements_df.to_csv(elements_outputs_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))

