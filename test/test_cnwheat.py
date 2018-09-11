# -*- coding: latin-1 -*-
"""
    test_cnwheat
    ~~~~~~~~~~~~

    Test the model CN-Wheat.

    You must first install :mod:`cnwheat` (and add it to your PYTHONPATH)
    before running this script with the command `python`.

    To get a coverage report of :mod:`cnwheat`, use the command:
    `nosetests --with-coverage --cover-package=cnwheat test_cnwheat.py`.

    CSV files must contain only ASCII characters and ',' as separator.

    :copyright: Copyright 2014-2017 INRA-ECOSYS, see AUTHORS.
    :license: CeCILL-C, see LICENSE for details.
    
    **Acknowledgments**: The research leading these results has received funding through the 
    Investment for the Future programme managed by the Research National Agency 
    (BreedWheat project ANR-10-BTBR-03).
    
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

from respiwheat import model as respiwheat_model

from cnwheat import simulation as cnwheat_simulation, tools as cnwheat_tools, converter as cnwheat_converter

# inputs directory path
INPUTS_DIRPATH = 'inputs'

# the file names of the inputs
ORGANS_INPUTS_FILENAME = 'organs_inputs.csv'
HIDDENZONES_INPUTS_FILENAME = 'hiddenzones_inputs.csv'
ELEMENTS_INPUTS_FILENAME = 'elements_inputs.csv'
SOILS_INPUTS_FILENAME = 'soils_inputs.csv'

# the file names of the data used to force photosynthesis and senescence parameters
PHOTOSYNTHESIS_ELEMENTS_DATA_FILENAME = 'photosynthesis_elements_data.csv'
SENESCENCE_ROOTS_DATA_FILENAME = 'senescence_roots_data.csv'
SENESCENCE_ELEMENTS_DATA_FILENAME = 'senescence_elements_data.csv'

# outputs directory path
OUTPUTS_DIRPATH = 'outputs'

# desired outputs filenames
DESIRED_AXES_OUTPUTS_FILENAME = 'desired_axes_outputs.csv'
DESIRED_ORGANS_OUTPUTS_FILENAME = 'desired_organs_outputs.csv'
DESIRED_HIDDENZONES_OUTPUTS_FILENAME = 'desired_hiddenzones_outputs.csv'
DESIRED_ELEMENTS_OUTPUTS_FILENAME = 'desired_elements_outputs.csv'
DESIRED_SOILS_OUTPUTS_FILENAME = 'desired_soils_outputs.csv'

# actual outputs filenames
ACTUAL_AXES_OUTPUTS_FILENAME = 'actual_axes_outputs.csv'
ACTUAL_ORGANS_OUTPUTS_FILENAME = 'actual_organs_outputs.csv'
ACTUAL_HIDDENZONES_OUTPUTS_FILENAME = 'actual_hiddenzones_outputs.csv'
ACTUAL_ELEMENTS_OUTPUTS_FILENAME = 'actual_elements_outputs.csv'
ACTUAL_SOILS_OUTPUTS_FILENAME = 'actual_soils_outputs.csv'

# culm density (culm m-2)
CULM_DENSITY = {1:410}

# number of seconds in 1 hour  
HOUR_TO_SECOND_CONVERSION_FACTOR = 3600


def force_senescence_and_photosynthesis(t, population, senescence_roots_data_grouped, senescence_elements_data_grouped, photosynthesis_elements_data_grouped):
        '''Force the senescence and photosynthesis data of the population at `t` from input grouped dataframes'''
        for plant in population.plants:
            for axis in plant.axes:
                # Root growth and senescence
                group = senescence_roots_data_grouped.get_group((t, plant.index, axis.label))
                senescence_data_to_use = group.loc[group.first_valid_index(), cnwheat_simulation.Simulation.ORGANS_STATE].dropna().to_dict()
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
                            senescence_data_to_use = group_senesc.loc[group_senesc.first_valid_index(), cnwheat_simulation.Simulation.ELEMENTS_STATE].dropna().to_dict()
                            element.__dict__.update(senescence_data_to_use)
                            # Element photosynthesis
                            group_photo = photosynthesis_elements_data_grouped.get_group((t, plant.index, axis.label, phytomer.index, organ.label, element.label))
                            photosynthesis_elements_data_to_use = group_photo.loc[group_photo.first_valid_index(), cnwheat_simulation.Simulation.ELEMENTS_STATE].dropna().to_dict()
                            element.__dict__.update(photosynthesis_elements_data_to_use)


def test_run():

    time_step_hours = 1
    time_step_seconds = time_step_hours * HOUR_TO_SECOND_CONVERSION_FACTOR

    # create the simulation
    simulation_ = cnwheat_simulation.Simulation(respiration_model=respiwheat_model, delta_t=time_step_seconds, culm_density=CULM_DENSITY)
    
    # read inputs from Pandas dataframes
    inputs_dataframes = {}
    for inputs_filename in (ORGANS_INPUTS_FILENAME, HIDDENZONES_INPUTS_FILENAME, ELEMENTS_INPUTS_FILENAME, SOILS_INPUTS_FILENAME):
        inputs_dataframes[inputs_filename] = pd.read_csv(os.path.join(INPUTS_DIRPATH, inputs_filename))

    # convert inputs to a population of plants and a dictionary of soils
    population, soils = cnwheat_converter.from_dataframes(inputs_dataframes[ORGANS_INPUTS_FILENAME],
                                                          inputs_dataframes[HIDDENZONES_INPUTS_FILENAME],
                                                          inputs_dataframes[ELEMENTS_INPUTS_FILENAME],
                                                          inputs_dataframes[SOILS_INPUTS_FILENAME])

    # initialize the simulation from the population and the soils
    simulation_.initialize(population, soils)

    # get photosynthesis data
    photosynthesis_elements_data_filepath = os.path.join(INPUTS_DIRPATH, PHOTOSYNTHESIS_ELEMENTS_DATA_FILENAME)
    photosynthesis_elements_data_df = pd.read_csv(photosynthesis_elements_data_filepath)
    photosynthesis_elements_data_grouped = photosynthesis_elements_data_df.groupby(cnwheat_simulation.Simulation.ELEMENTS_T_INDEXES)

    # get senescence and growth data
    senescence_roots_data_filepath = os.path.join(INPUTS_DIRPATH, SENESCENCE_ROOTS_DATA_FILENAME)
    senescence_roots_data_df = pd.read_csv(senescence_roots_data_filepath)
    senescence_roots_data_grouped = senescence_roots_data_df.groupby(cnwheat_simulation.Simulation.AXES_T_INDEXES)
    senescence_elements_data_filepath = os.path.join(INPUTS_DIRPATH, SENESCENCE_ELEMENTS_DATA_FILENAME)
    senescence_elements_data_df = pd.read_csv(senescence_elements_data_filepath)
    senescence_elements_data_grouped = senescence_elements_data_df.groupby(cnwheat_simulation.Simulation.ELEMENTS_T_INDEXES)
    
    # create empty lists of dataframes to store the outputs at each step
    axes_outputs_df_list = []
    organs_outputs_df_list = []
    hiddenzones_outputs_df_list = []
    elements_outputs_df_list = []
    soils_outputs_df_list = []
    
    # define the time grid to run the model on
    start_time = 0
    stop_time = 2
    time_grid = xrange(start_time, stop_time+time_step_hours, time_step_hours)
    
    # force the senescence and photosynthesis of the population
    force_senescence_and_photosynthesis(0, population, senescence_roots_data_grouped, senescence_elements_data_grouped, photosynthesis_elements_data_grouped)
    
    # reinitialize the simulation from forced population and soils
    simulation_.initialize(population, soils)

    # run the model on the time grid
    for t in time_grid:

        if t > 0: # do not run the model at t = 0
            # run the model of CN exchanges ; the population is internally updated by the model
            simulation_.run()
        
        # convert outputs to dataframes
        _, axes_outputs_df, _, organs_outputs_df, hiddenzones_outputs_df, elements_outputs_df, soils_outputs_df = cnwheat_converter.to_dataframes(simulation_.population, simulation_.soils)
        
        # append the outputs at current t to the lists of dataframes
        for df, list_ in ((axes_outputs_df, axes_outputs_df_list), (organs_outputs_df, organs_outputs_df_list), 
                          (hiddenzones_outputs_df, hiddenzones_outputs_df_list), (elements_outputs_df, elements_outputs_df_list), 
                          (soils_outputs_df, soils_outputs_df_list)):
            df.insert(0, 't', t)
            list_.append(df)
        
        if t > 0 and t < stop_time:
            
            # force the senescence and photosynthesis of the population
            force_senescence_and_photosynthesis(t, population, senescence_roots_data_grouped, senescence_elements_data_grouped, photosynthesis_elements_data_grouped)
            # reinitialize the simulation from forced population and soils
            simulation_.initialize(population, soils)
    
    # compare actual to desired outputs at each scale level (an exception is raised if the test failed)  
    for (outputs_df_list, 
         desired_outputs_filename, 
         actual_outputs_filename,
         state_variables_names) \
         in ((axes_outputs_df_list, DESIRED_AXES_OUTPUTS_FILENAME, ACTUAL_AXES_OUTPUTS_FILENAME, cnwheat_simulation.Simulation.AXES_T_INDEXES + cnwheat_simulation.Simulation.AXES_STATE),
             (organs_outputs_df_list, DESIRED_ORGANS_OUTPUTS_FILENAME, ACTUAL_ORGANS_OUTPUTS_FILENAME, cnwheat_simulation.Simulation.ORGANS_T_INDEXES + cnwheat_simulation.Simulation.ORGANS_STATE),
             (hiddenzones_outputs_df_list, DESIRED_HIDDENZONES_OUTPUTS_FILENAME, ACTUAL_HIDDENZONES_OUTPUTS_FILENAME, cnwheat_simulation.Simulation.HIDDENZONE_T_INDEXES + cnwheat_simulation.Simulation.HIDDENZONE_STATE),
             (elements_outputs_df_list, DESIRED_ELEMENTS_OUTPUTS_FILENAME, ACTUAL_ELEMENTS_OUTPUTS_FILENAME, cnwheat_simulation.Simulation.ELEMENTS_T_INDEXES + cnwheat_simulation.Simulation.ELEMENTS_STATE),
             (soils_outputs_df_list, DESIRED_SOILS_OUTPUTS_FILENAME, ACTUAL_SOILS_OUTPUTS_FILENAME, cnwheat_simulation.Simulation.SOILS_T_INDEXES + cnwheat_simulation.Simulation.SOILS_STATE)):
        outputs_df = pd.concat(outputs_df_list, ignore_index=True)
        outputs_df = outputs_df.loc[:, state_variables_names] # compare only the values of the compartments
        print 'Compare', actual_outputs_filename, 'to', desired_outputs_filename 
        cnwheat_tools.compare_actual_to_desired(OUTPUTS_DIRPATH, outputs_df, desired_outputs_filename, actual_outputs_filename)
        print actual_outputs_filename, 'OK!' 
        

if __name__ == '__main__':
    test_run()
