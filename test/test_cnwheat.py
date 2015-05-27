# -*- coding: latin-1 -*-
"""
    test_cnwheat
    ~~~~~~~~~~~~

    Test the CN-Wheat model.

    You must first install :mod:`cnwheat` (and add it to your PYTHONPATH)
    before running this script with the command `python`.

    CSV files must contain only ASCII characters and ',' as separator.

    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.
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

import numpy as np
import pandas as pd

from cnwheat import simulation
from cnwheat import model as cnwheat_model
from cnwheat import tools

import time, datetime

t0 = time.time()

INPUTS_DIRPATH = 'inputs'

PLANTS_INPUTS_FILENAME = 'plants_inputs.csv'
AXES_INPUTS_FILENAME = 'axes_inputs.csv'
PHYTOMERS_INPUTS_FILENAME = 'phytomers_inputs.csv'
ORGANS_INPUTS_FILENAME = 'organs_inputs.csv'
ELEMENTS_INPUTS_FILENAME = 'elements_inputs.csv'

PHOTOSYNTHESIS_DATA_FILENAME = 'photosynthesis_data.csv'
SENESCENCE_DATA_FILENAME = 'senescence_data.csv'

OUTPUTS_DIRPATH = 'outputs'

DESIRED_PLANTS_OUTPUTS_FILENAME = 'desired_plants_outputs.csv'
DESIRED_AXES_OUTPUTS_FILENAME = 'desired_axes_outputs.csv'
DESIRED_PHYTOMERS_OUTPUTS_FILENAME = 'desired_phytomers_outputs.csv'
DESIRED_ORGANS_OUTPUTS_FILENAME = 'desired_organs_outputs.csv'
DESIRED_ELEMENTS_OUTPUTS_FILENAME = 'desired_elements_outputs.csv'

ACTUAL_PLANTS_OUTPUTS_FILENAME = 'actual_plants_outputs.csv'
ACTUAL_AXES_OUTPUTS_FILENAME = 'actual_axes_outputs.csv'
ACTUAL_PHYTOMERS_OUTPUTS_FILENAME = 'actual_phytomers_outputs.csv'
ACTUAL_ORGANS_OUTPUTS_FILENAME = 'actual_organs_outputs.csv'
ACTUAL_ELEMENTS_OUTPUTS_FILENAME = 'actual_elements_outputs.csv'


def test_run():
    
    # create the simulation
    simulation_ = simulation.Simulation()
    # read inputs from Pandas dataframes
    inputs_dataframes = {}
    for inputs_filename in (PLANTS_INPUTS_FILENAME, AXES_INPUTS_FILENAME, PHYTOMERS_INPUTS_FILENAME, ORGANS_INPUTS_FILENAME, ELEMENTS_INPUTS_FILENAME):
        inputs_dataframes[inputs_filename] = pd.read_csv(os.path.join(INPUTS_DIRPATH, inputs_filename))
    # initialize the simulation with the inputs
    simulation_.initialize(plants_inputs=inputs_dataframes[PLANTS_INPUTS_FILENAME], 
                           axes_inputs=inputs_dataframes[AXES_INPUTS_FILENAME],
                           phytomers_inputs=inputs_dataframes[PHYTOMERS_INPUTS_FILENAME],
                           organs_inputs=inputs_dataframes[ORGANS_INPUTS_FILENAME],
                           elements_inputs=inputs_dataframes[ELEMENTS_INPUTS_FILENAME])
    # format the inputs to Pandas dataframes
    formatted_inputs_dataframes = {}
    formatted_inputs_dataframes[PLANTS_INPUTS_FILENAME], \
    formatted_inputs_dataframes[AXES_INPUTS_FILENAME], \
    formatted_inputs_dataframes[PHYTOMERS_INPUTS_FILENAME], \
    formatted_inputs_dataframes[ORGANS_INPUTS_FILENAME], \
    formatted_inputs_dataframes[ELEMENTS_INPUTS_FILENAME] = simulation_.format_inputs()
    # compare inputs
    for inputs_filename, inputs_df in formatted_inputs_dataframes.iteritems():
        tools.compare_actual_to_desired(INPUTS_DIRPATH, inputs_df, inputs_filename)
        
    # Get photosynthesis data
    photosynthesis_data_filepath = os.path.join(INPUTS_DIRPATH, PHOTOSYNTHESIS_DATA_FILENAME)
    photosynthesis_data_df = pd.read_csv(photosynthesis_data_filepath)
    photosynthesis_data_grouped = photosynthesis_data_df.groupby(simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES)
 
    # Get senescence and growth data
    senescence_data_filepath = os.path.join(INPUTS_DIRPATH, SENESCENCE_DATA_FILENAME)
    senescence_data_df = pd.read_csv(senescence_data_filepath)
    senescence_data_grouped = senescence_data_df.groupby(simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES)

    start_time = 0
    stop_time = 48
    time_step = 4

    all_plants_df_list = []
    all_axes_df_list = []
    all_phytomers_df_list = []
    all_organs_df_list = []
    all_elements_df_list = []

    for t in xrange(start_time, stop_time, time_step):
        # update the population
        for plant in simulation_.population.plants:
            plant_index = plant.index
            for axis in plant.axes:
                axis_id = axis.id

                # Root growth and senescence
                group = senescence_data_grouped.get_group((t, plant_index, axis_id, 0, 'Roots', 'enclosed'))
                row_index = group.first_valid_index()
                axis.roots.mstruct_C_growth = group.mstruct_growth[row_index]
                axis.roots.Nstruct_N_growth = group.Nstruct_N_growth[row_index]
                axis.roots.mstruct_senescence = group.mstruct_senescence[row_index]
                axis.roots.mstruct = group.mstruct[row_index]
                axis.roots.Nstruct = group.Nstruct[row_index]

                for phytomer in axis.phytomers:
                    phytomer_index = phytomer.index
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ is None:
                            continue
                        organ_type = organ.__class__.__name__
                        for element, element_type in ((organ.exposed_element, 'exposed'), (organ.enclosed_element, 'enclosed')):
                            if element is None:
                                continue
                            group_photo = photosynthesis_data_grouped.get_group((t, plant_index, axis_id, phytomer_index, organ_type, element_type))
                            group_senesc = senescence_data_grouped.get_group((t, plant_index, axis_id, phytomer_index, organ_type, element_type))
                            row_index_photo = group_photo.first_valid_index()
                            row_index_sensc = group_senesc.first_valid_index()

                            # Senescence
                            element.green_area = group_senesc.green_area[row_index_sensc]
                            element.relative_delta_green_area = group_senesc.relative_delta_green_area[row_index_sensc]
                            element.mstruct = group_senesc.mstruct[row_index_sensc]
                            element.Nstruct = group_senesc.Nstruct[row_index_sensc]
                            element.surfacic_nitrogen = group_senesc.SLN[row_index_sensc]

                            element.Ag = group_photo.Ag[row_index_photo]
                            element.An = group_photo.An[row_index_photo]
                            element.Rd = group_photo.Rd[row_index_photo]
                            element.Tr = group_photo.Tr[row_index_photo]
                            element.Ts = group_photo.Ts[row_index_photo]
                            element.gs = group_photo.gs[row_index_photo]

        # run the model of CN exchanges ; the population is internally updated by the model of CN exchanges
        simulation_.run(start_time=t, stop_time=t+time_step, number_of_output_steps=time_step+1)
        
        # format the outputs to Pandas dataframes
        all_plants_df, all_axes_df, all_phytomers_df, all_organs_df, all_elements_df = simulation_.format_outputs()
        
        all_plants_df_list.append(all_plants_df)
        all_axes_df_list.append(all_axes_df)
        all_phytomers_df_list.append(all_phytomers_df)
        all_organs_df_list.append(all_organs_df)
        all_elements_df_list.append(all_elements_df)

    execution_time = int(time.time()-t0)
    print '\n', 'Model executed in ', str(datetime.timedelta(seconds=execution_time))

    global_plants_df = pd.concat(all_plants_df_list, ignore_index=True)
    global_plants_df.drop_duplicates(subset=simulation.Simulation.PLANTS_OUTPUTS_INDEXES, inplace=True)
    tools.compare_actual_to_desired(OUTPUTS_DIRPATH, global_plants_df, DESIRED_PLANTS_OUTPUTS_FILENAME, ACTUAL_PLANTS_OUTPUTS_FILENAME)

    global_axes_df = pd.concat(all_axes_df_list, ignore_index=True)
    global_axes_df.drop_duplicates(subset=simulation.Simulation.AXES_OUTPUTS_INDEXES, inplace=True)
    tools.compare_actual_to_desired(OUTPUTS_DIRPATH, global_axes_df, DESIRED_AXES_OUTPUTS_FILENAME, ACTUAL_AXES_OUTPUTS_FILENAME)

    global_phytomers_df = pd.concat(all_phytomers_df_list, ignore_index=True)
    global_phytomers_df.drop_duplicates(subset=simulation.Simulation.PHYTOMERS_OUTPUTS_INDEXES, inplace=True)
    tools.compare_actual_to_desired(OUTPUTS_DIRPATH, global_phytomers_df, DESIRED_PHYTOMERS_OUTPUTS_FILENAME, ACTUAL_PHYTOMERS_OUTPUTS_FILENAME)

    global_organs_df = pd.concat(all_organs_df_list, ignore_index=True)
    global_organs_df.drop_duplicates(subset=simulation.Simulation.ORGANS_OUTPUTS_INDEXES, inplace=True)
    tools.compare_actual_to_desired(OUTPUTS_DIRPATH, global_organs_df, DESIRED_ORGANS_OUTPUTS_FILENAME, ACTUAL_ORGANS_OUTPUTS_FILENAME)

    global_elements_df = pd.concat(all_elements_df_list, ignore_index=True)
    global_elements_df.drop_duplicates(subset=simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES, inplace=True)
    tools.compare_actual_to_desired(OUTPUTS_DIRPATH, global_elements_df, DESIRED_ELEMENTS_OUTPUTS_FILENAME, ACTUAL_ELEMENTS_OUTPUTS_FILENAME)

if __name__ == '__main__':
    test_run()
