# -*- coding: latin-1 -*-
"""
    test_cnwheat
    ~~~~~~~~~~~~

    Test the CN-Wheat model.

    You must first install :mod:`cnwheat` (and add it to your PYTHONPATH)
    before running this script with the command `python`.

    To get a coverage report of :mod:`cnwheat`, use the command:
    `nosetests --with-coverage --cover-package=cnwheat test_cnwheat.py`.

    CSV files must contain only ASCII characters and ',' as separator.

    :copyright: Copyright 2014-2016 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2016.
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

import pandas as pd

from cnwheat import simulation, tools, converter, model

INPUTS_DIRPATH = 'inputs'

PLANTS_INPUTS_FILENAME = 'plants_inputs.csv'
AXES_INPUTS_FILENAME = 'axes_inputs.csv'
METAMERS_INPUTS_FILENAME = 'metamers_inputs.csv'
ORGANS_INPUTS_FILENAME = 'organs_inputs.csv'
HIDDENZONES_INPUTS_FILENAME = 'hiddenzones_inputs.csv'
ELEMENTS_INPUTS_FILENAME = 'elements_inputs.csv'
SOILS_INPUTS_FILENAME = 'soils_inputs.csv'

PHOTOSYNTHESIS_ELEMENTS_DATA_FILENAME = 'photosynthesis_elements_data.csv'
SENESCENCE_ROOTS_DATA_FILENAME = 'senescence_roots_data.csv'
SENESCENCE_ELEMENTS_DATA_FILENAME = 'senescence_elements_data.csv'

OUTPUTS_DIRPATH = 'outputs'

DESIRED_PLANTS_OUTPUTS_FILENAME = 'desired_plants_outputs.csv'
DESIRED_AXES_OUTPUTS_FILENAME = 'desired_axes_outputs.csv'
DESIRED_METAMERS_OUTPUTS_FILENAME = 'desired_metamers_outputs.csv'
DESIRED_ORGANS_OUTPUTS_FILENAME = 'desired_organs_outputs.csv'
DESIRED_HIDDENZONES_OUTPUTS_FILENAME = 'desired_hiddenzones_outputs.csv'
DESIRED_ELEMENTS_OUTPUTS_FILENAME = 'desired_elements_outputs.csv'
DESIRED_SOILS_OUTPUTS_FILENAME = 'desired_soils_outputs.csv'

ACTUAL_PLANTS_OUTPUTS_FILENAME = 'actual_plants_outputs.csv'
ACTUAL_AXES_OUTPUTS_FILENAME = 'actual_axes_outputs.csv'
ACTUAL_METAMERS_OUTPUTS_FILENAME = 'actual_metamers_outputs.csv'
ACTUAL_ORGANS_OUTPUTS_FILENAME = 'actual_organs_outputs.csv'
ACTUAL_HIDDENZONES_OUTPUTS_FILENAME = 'actual_hiddenzones_outputs.csv'
ACTUAL_ELEMENTS_OUTPUTS_FILENAME = 'actual_elements_outputs.csv'
ACTUAL_SOILS_OUTPUTS_FILENAME = 'actual_soils_outputs.csv'


def test_run():

    # create the simulation
    simulation_ = simulation.Simulation(delta_t=3600)
    # read inputs from Pandas dataframes
    inputs_dataframes = {}
    for inputs_filename in (PLANTS_INPUTS_FILENAME, AXES_INPUTS_FILENAME, METAMERS_INPUTS_FILENAME, ORGANS_INPUTS_FILENAME, HIDDENZONES_INPUTS_FILENAME, ELEMENTS_INPUTS_FILENAME, SOILS_INPUTS_FILENAME):
        inputs_dataframes[inputs_filename] = pd.read_csv(os.path.join(INPUTS_DIRPATH, inputs_filename))

    # convert inputs to a population of plants and a dictionary of soils
    population, soils = converter.from_dataframes(inputs_dataframes[ORGANS_INPUTS_FILENAME],
                                                  inputs_dataframes[HIDDENZONES_INPUTS_FILENAME],
                                                  inputs_dataframes[ELEMENTS_INPUTS_FILENAME],
                                                  inputs_dataframes[SOILS_INPUTS_FILENAME])

    # initialize the simulation from the population and the soils
    simulation_.initialize(population, soils)
    # convert the population and the soils to Pandas dataframes
    formatted_inputs_dataframes = {}
    formatted_inputs_dataframes[AXES_INPUTS_FILENAME], \
    formatted_inputs_dataframes[ORGANS_INPUTS_FILENAME], \
    formatted_inputs_dataframes[HIDDENZONES_INPUTS_FILENAME], \
    formatted_inputs_dataframes[ELEMENTS_INPUTS_FILENAME], \
    formatted_inputs_dataframes[SOILS_INPUTS_FILENAME] = converter.to_dataframes(population, soils)
##    # compare inputs
##    for inputs_filename, inputs_df in formatted_inputs_dataframes.iteritems():
##        tools.compare_actual_to_desired(INPUTS_DIRPATH, inputs_df, inputs_filename)

    # Get photosynthesis data
    photosynthesis_elements_data_filepath = os.path.join(INPUTS_DIRPATH, PHOTOSYNTHESIS_ELEMENTS_DATA_FILENAME)
    photosynthesis_elements_data_df = pd.read_csv(photosynthesis_elements_data_filepath)
    photosynthesis_elements_data_grouped = photosynthesis_elements_data_df.groupby(simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES)

    # Get senescence and growth data
    senescence_roots_data_filepath = os.path.join(INPUTS_DIRPATH, SENESCENCE_ROOTS_DATA_FILENAME)
    senescence_roots_data_df = pd.read_csv(senescence_roots_data_filepath)
    senescence_roots_data_grouped = senescence_roots_data_df.groupby(simulation.Simulation.AXES_OUTPUTS_INDEXES)
    senescence_elements_data_filepath = os.path.join(INPUTS_DIRPATH, SENESCENCE_ELEMENTS_DATA_FILENAME)
    senescence_elements_data_df = pd.read_csv(senescence_elements_data_filepath)
    senescence_elements_data_grouped = senescence_elements_data_df.groupby(simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES)

    start_time = 0
    stop_time = 48
    time_step = 4

    all_plants_df_list = []
    all_axes_df_list = []
    all_metamers_df_list = []
    all_organs_df_list = []
    all_hiddenzones_df_list = []
    all_elements_df_list = []
    all_soils_df_list = []

    for t in xrange(start_time, stop_time, time_step):
        # force the senescence and photosynthesis of the population
        for plant in simulation_.population.plants:
            for axis in plant.axes:
                # Root growth and senescence
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
                            # Element senescence
                            group_senesc = senescence_elements_data_grouped.get_group((t, plant.index, axis.label, phytomer.index, organ.label, element.label))
                            senescence_data_to_use = group_senesc.loc[group_senesc.first_valid_index(), simulation.Simulation.ELEMENTS_STATE].dropna().to_dict()
                            element.__dict__.update(senescence_data_to_use)
                            # Element photosynthesis
                            group_photo = photosynthesis_elements_data_grouped.get_group((t, plant.index, axis.label, phytomer.index, organ.label, element.label))
                            photosynthesis_elements_data_to_use = group_photo.loc[group_photo.first_valid_index(), simulation.Simulation.ELEMENTS_STATE].dropna().to_dict()
                            element.__dict__.update(photosynthesis_elements_data_to_use)

        # run the model of CN exchanges ; the population is internally updated by the model
        simulation_.run(start_time=t, stop_time=t+time_step, number_of_output_steps=time_step+1)

        # run post-processings
        all_plants_df, all_axes_df, all_metamers_df, all_organs_df, all_hiddenzones_df, all_elements_df, all_soils_df = simulation_.postprocessings()

        all_plants_df_list.append(all_plants_df)
        all_axes_df_list.append(all_axes_df)
        all_metamers_df_list.append(all_metamers_df)
        all_organs_df_list.append(all_organs_df)
        all_hiddenzones_df_list.append(all_hiddenzones_df)
        all_elements_df_list.append(all_elements_df)
        all_soils_df_list.append(all_soils_df)

    global_plants_df = pd.concat(all_plants_df_list, ignore_index=True)
    global_plants_df.drop_duplicates(subset=simulation.Simulation.PLANTS_OUTPUTS_INDEXES, inplace=True)
    tools.compare_actual_to_desired(OUTPUTS_DIRPATH, global_plants_df, DESIRED_PLANTS_OUTPUTS_FILENAME, ACTUAL_PLANTS_OUTPUTS_FILENAME)

    global_axes_df = pd.concat(all_axes_df_list, ignore_index=True)
    global_axes_df.drop_duplicates(subset=simulation.Simulation.AXES_OUTPUTS_INDEXES, inplace=True)
    tools.compare_actual_to_desired(OUTPUTS_DIRPATH, global_axes_df, DESIRED_AXES_OUTPUTS_FILENAME, ACTUAL_AXES_OUTPUTS_FILENAME)

    global_metamers_df = pd.concat(all_metamers_df_list, ignore_index=True)
    global_metamers_df.drop_duplicates(subset=simulation.Simulation.PHYTOMERS_OUTPUTS_INDEXES, inplace=True)
    tools.compare_actual_to_desired(OUTPUTS_DIRPATH, global_metamers_df, DESIRED_METAMERS_OUTPUTS_FILENAME, ACTUAL_METAMERS_OUTPUTS_FILENAME)

    global_organs_df = pd.concat(all_organs_df_list, ignore_index=True)
    global_organs_df.drop_duplicates(subset=simulation.Simulation.ORGANS_OUTPUTS_INDEXES, inplace=True)
    tools.compare_actual_to_desired(OUTPUTS_DIRPATH, global_organs_df, DESIRED_ORGANS_OUTPUTS_FILENAME, ACTUAL_ORGANS_OUTPUTS_FILENAME)

    global_hiddenzones_df = pd.concat(all_hiddenzones_df_list, ignore_index=True)
    global_hiddenzones_df.drop_duplicates(subset=simulation.Simulation.HIDDENZONE_OUTPUTS_INDEXES, inplace=True)
    tools.compare_actual_to_desired(OUTPUTS_DIRPATH, global_hiddenzones_df, DESIRED_HIDDENZONES_OUTPUTS_FILENAME, ACTUAL_HIDDENZONES_OUTPUTS_FILENAME)

    global_elements_df = pd.concat(all_elements_df_list, ignore_index=True)
    global_elements_df.drop_duplicates(subset=simulation.Simulation.ELEMENTS_OUTPUTS_INDEXES, inplace=True)
    tools.compare_actual_to_desired(OUTPUTS_DIRPATH, global_elements_df, DESIRED_ELEMENTS_OUTPUTS_FILENAME, ACTUAL_ELEMENTS_OUTPUTS_FILENAME)

    global_soils_df = pd.concat(all_soils_df_list, ignore_index=True)
    global_soils_df.drop_duplicates(subset=simulation.Simulation.SOILS_OUTPUTS_INDEXES, inplace=True)
    tools.compare_actual_to_desired(OUTPUTS_DIRPATH, global_soils_df, DESIRED_SOILS_OUTPUTS_FILENAME, ACTUAL_SOILS_OUTPUTS_FILENAME)

if __name__ == '__main__':
    test_run()