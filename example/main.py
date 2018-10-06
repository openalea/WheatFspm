# -*- coding: latin-1 -*-

"""
    main
    ~~~~

    An example to show how to run the model CN-Wheat, compute the post-processing, and generate the plots for validation.

    You must first install :mod:`cnwheat` (and add it to your PYTHONPATH)
    before running this script with the command `python`.

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


###############################################
####### CONFIGURATION OF THE SIMULATION #######
###############################################


### INPUTS CONFIGURATION ###

# Path of the directory which contains the inputs of the model
INPUTS_DIRPATH = 'inputs'

# Name of the CSV files which describe the initial state of the system
ORGANS_INITIAL_STATE_FILENAME = 'organs_initial_state.csv'
HIDDENZONES_INITIAL_STATE_FILENAME = 'hiddenzones_initial_state.csv'
ELEMENTS_INITIAL_STATE_FILENAME = 'elements_initial_state.csv'
SOILS_INITIAL_STATE_FILENAME = 'soils_initial_state.csv'

# Name of the CSV files which contain the photosynthesis and senescence forcings
ELEMENTS_PHOTOSYNTHESIS_FORCINGS_FILENAME = 'elements_photosynthesis_forcings.csv'
ROOTS_SENESCENCE_FORCINGS_FILENAME = 'roots_senescence_forcings.csv'
ELEMENTS_SENESCENCE_FORCINGS_FILENAME = 'elements_senescence_forcings.csv'


### OUTPUTS CONFIGURATION ###

# Path of the directory where to write the outputs of the model
OUTPUTS_DIRPATH = 'outputs'

# Name of the CSV files which will contain the outputs of the model
AXES_OUTPUTS_FILENAME = 'axes_outputs.csv'
ORGANS_OUTPUTS_FILENAME = 'organs_outputs.csv'
HIDDENZONES_OUTPUTS_FILENAME = 'hiddenzones_outputs.csv'
ELEMENTS_OUTPUTS_FILENAME = 'elements_outputs.csv'
SOILS_OUTPUTS_FILENAME = 'soils_outputs.csv'


### POSTPROCESSING CONFIGURATION ###

# Path of the directory where to write the postprocessing of the model
POSTPROCESSING_DIRPATH = 'postprocessing'

# Name of the CSV files which will contain the postprocessing of the model
AXES_POSTPROCESSING_FILENAME = 'axes_postprocessing.csv'
ORGANS_POSTPROCESSING_FILENAME = 'organs_postprocessing.csv'
HIDDENZONES_POSTPROCESSING_FILENAME = 'hiddenzones_postprocessing.csv'
ELEMENTS_POSTPROCESSING_FILENAME = 'elements_postprocessing.csv'
SOILS_POSTPROCESSING_FILENAME = 'soils_postprocessing.csv'


### GRAPHS CONFIGURATION ###

# Path of the directory where to save the generated graphs
GRAPHS_DIRPATH = 'graphs'


### SIMULATION PARAMETERS ###

# Start time of the simulation
START_TIME = 0

# Length of the simulation (in hours)
SIMULATION_LENGTH = 48

# Time step of the simulation (in hours)
TIME_STEP = 4

# Do run the simulation?
RUN_SIMU = True

# Do run the postprocessing?
RUN_POSTPROCESSING = True

# Do generate the graphs?
GENERATE_GRAPHS = True

# Do log the execution?
LOG_EXECUTION = False

# Config file path for logging
LOGGING_CONFIG_FILEPATH = 'logging.json'

# Do interpolate the forcings?
INTERPOLATE_FORCINGS = False

# Culm density (culm m-2)
CULM_DENSITY = {1:410}



###############################################
#######      RUN OF THE SIMULATION      #######
###############################################

# Warning: you should not modify the following code excepting you really know what you are doing

import os
import logging
import datetime

import pandas as pd

from respiwheat import model as respiwheat_model
from cnwheat import simulation as cnwheat_simulation, converter as cnwheat_converter, \
    tools as cnwheat_tools, postprocessing as cnwheat_postprocessing

# Number of seconds in 1 hour
HOUR_TO_SECOND_CONVERSION_FACTOR = 3600

# Precision of floats used to write and format the output CSV files
OUTPUTS_PRECISION = 6


def force_senescence_and_photosynthesis(t, population, senescence_roots_data_grouped, senescence_elements_data_grouped, photosynthesis_elements_data_grouped):
        """Force the senescence and photosynthesis data of the population at `t` from input grouped dataframes"""
        for plant in population.plants:
            for axis in plant.axes:
                # Root growth and senescence
                group = senescence_roots_data_grouped.get_group((t, plant.index, axis.label))
                senescence_data_to_use = group.loc[group.first_valid_index(), group.columns.intersection(cnwheat_simulation.Simulation.ORGANS_STATE)].dropna().to_dict()
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
                            senescence_data_to_use = group_senesc.loc[group_senesc.first_valid_index(), group_senesc.columns.intersection(cnwheat_simulation.Simulation.ELEMENTS_STATE)].dropna().to_dict()
                            element.__dict__.update(senescence_data_to_use)
                            # Element photosynthesis
                            group_photo = photosynthesis_elements_data_grouped.get_group((t, plant.index, axis.label, phytomer.index, organ.label, element.label))
                            photosynthesis_elements_data_to_use = group_photo.loc[group_photo.first_valid_index(), group_photo.columns.intersection(cnwheat_simulation.Simulation.ELEMENTS_STATE)].dropna().to_dict()
                            element.__dict__.update(photosynthesis_elements_data_to_use)


if RUN_SIMU:

    print 'Prepare the simulation...'

    time_step_seconds = TIME_STEP * HOUR_TO_SECOND_CONVERSION_FACTOR

    if LOG_EXECUTION:
        # Setup the logging
        cnwheat_tools.setup_logging(config_filepath=LOGGING_CONFIG_FILEPATH, level=logging.DEBUG,
                  log_model=True, log_compartments=True, log_derivatives=True, remove_old_logs=True)

    # Read the inputs from CSV files and create inputs dataframes
    inputs_dataframes = {}
    for inputs_filename in (
    ORGANS_INITIAL_STATE_FILENAME, HIDDENZONES_INITIAL_STATE_FILENAME, ELEMENTS_INITIAL_STATE_FILENAME,
    SOILS_INITIAL_STATE_FILENAME):
        inputs_dataframes[inputs_filename] = pd.read_csv(os.path.join(INPUTS_DIRPATH, inputs_filename))

    # Convert the inputs dataframes to a population of plants and a dictionary of soils
    population, soils = cnwheat_converter.from_dataframes(inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME],
                                                          inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME],
                                                          inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME],
                                                          inputs_dataframes[SOILS_INITIAL_STATE_FILENAME])

    # Set simulation parameters for interpolation
    if INTERPOLATE_FORCINGS:
        senescence_forcings_delta_t = time_step_seconds
        photosynthesis_forcings_delta_t = time_step_seconds
    else:
        senescence_forcings_delta_t = None
        photosynthesis_forcings_delta_t = None

    # Create the simulation
    simulation_ = cnwheat_simulation.Simulation(respiration_model=respiwheat_model, delta_t=time_step_seconds, culm_density=CULM_DENSITY,
                                                interpolate_forcings=INTERPOLATE_FORCINGS, senescence_forcings_delta_t=senescence_forcings_delta_t,
                                                photosynthesis_forcings_delta_t=photosynthesis_forcings_delta_t)

    # Initialize the simulation from the population of plants and the dictionary of soils created previously
    simulation_.initialize(population, soils)

    # Read photosynthesis and senescence forcings from CSV files, create dataframes, and group the dataframes by object index
    photosynthesis_elements_data_filepath = os.path.join(INPUTS_DIRPATH, ELEMENTS_PHOTOSYNTHESIS_FORCINGS_FILENAME)
    photosynthesis_elements_data_df = pd.read_csv(photosynthesis_elements_data_filepath)
    photosynthesis_elements_data_grouped = photosynthesis_elements_data_df.groupby(cnwheat_simulation.Simulation.ELEMENTS_T_INDEXES)
    senescence_roots_data_filepath = os.path.join(INPUTS_DIRPATH, ROOTS_SENESCENCE_FORCINGS_FILENAME)
    senescence_roots_data_df = pd.read_csv(senescence_roots_data_filepath)
    senescence_roots_data_grouped = senescence_roots_data_df.groupby(cnwheat_simulation.Simulation.AXES_T_INDEXES)
    senescence_elements_data_filepath = os.path.join(INPUTS_DIRPATH, ELEMENTS_SENESCENCE_FORCINGS_FILENAME)
    senescence_elements_data_df = pd.read_csv(senescence_elements_data_filepath)
    senescence_elements_data_grouped = senescence_elements_data_df.groupby(cnwheat_simulation.Simulation.ELEMENTS_T_INDEXES)

    # Force the senescence and photosynthesis of the population
    force_senescence_and_photosynthesis(0, population, senescence_roots_data_grouped, senescence_elements_data_grouped, photosynthesis_elements_data_grouped)

    # Reinitialize the simulation from forced population and soils
    simulation_.initialize(population, soils)

    # Define the time grid of the simulation
    time_grid = xrange(START_TIME, SIMULATION_LENGTH + TIME_STEP, TIME_STEP)

    # Create empty lists of dataframes to store the outputs at each step of the simulation
    axes_outputs_df_list = []
    organs_outputs_df_list = []
    hiddenzones_outputs_df_list = []
    elements_outputs_df_list = []
    soils_outputs_df_list = []

    print 'Prepare the simulation... DONE!'

    print 'Run the simulation...'
    current_time_of_the_system = datetime.datetime.now()

    for t in time_grid:

        if t > 0:
            # Run the model of CN exchanges ; the population is internally updated by the model
            print '\tt =', t
            simulation_.run()

        # Convert the model outputs to dataframes
        _, axes_outputs_df, _, organs_outputs_df, hiddenzones_outputs_df, elements_outputs_df, soils_outputs_df = cnwheat_converter.to_dataframes(simulation_.population, simulation_.soils)

        # Append the outputs dataframes at current t to the global lists of dataframes
        for df, list_ in ((axes_outputs_df, axes_outputs_df_list), (organs_outputs_df, organs_outputs_df_list),
                          (hiddenzones_outputs_df, hiddenzones_outputs_df_list), (elements_outputs_df, elements_outputs_df_list),
                          (soils_outputs_df, soils_outputs_df_list)):
            df.insert(0, 't', t)
            list_.append(df)

        if t > 0 and t < SIMULATION_LENGTH:

            # Force the senescence and photosynthesis of the population
            force_senescence_and_photosynthesis(t, population, senescence_roots_data_grouped, senescence_elements_data_grouped, photosynthesis_elements_data_grouped)
            # Reinitialize the simulation from forced population and soils
            simulation_.initialize(population, soils)

    print 'Run the simulation... DONE!'

    execution_time = datetime.datetime.now() - current_time_of_the_system
    print  'Simulation run in ', execution_time

    print 'Total RHS evaluations: ', simulation_.nfev_total

    print 'Write the outputs to CSV files...'

    outputs_df_dict = {}
    for outputs_df_list, outputs_filename in ((axes_outputs_df_list, AXES_OUTPUTS_FILENAME),
                                              (organs_outputs_df_list, ORGANS_OUTPUTS_FILENAME),
                                              (hiddenzones_outputs_df_list, HIDDENZONES_OUTPUTS_FILENAME),
                                              (elements_outputs_df_list, ELEMENTS_OUTPUTS_FILENAME),
                                              (soils_outputs_df_list, SOILS_OUTPUTS_FILENAME)):
        outputs_filepath = os.path.join(OUTPUTS_DIRPATH, outputs_filename)
        outputs_df = pd.concat(outputs_df_list, ignore_index=True)
        outputs_df.to_csv(outputs_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))
        outputs_file_basename = outputs_filename.split('.')[0]
        outputs_df_dict[outputs_file_basename] = outputs_df

    print 'Write the outputs to CSV files... DONE!'


if RUN_POSTPROCESSING:

    if not RUN_SIMU:

        print 'Retrieve outputs dataframes from precedent simulation run...'

        outputs_df_dict = {}

        for outputs_filename in (AXES_OUTPUTS_FILENAME,
                                 ORGANS_OUTPUTS_FILENAME,
                                 HIDDENZONES_OUTPUTS_FILENAME,
                                 ELEMENTS_OUTPUTS_FILENAME,
                                 SOILS_OUTPUTS_FILENAME):
            outputs_filepath = os.path.join(OUTPUTS_DIRPATH, outputs_filename)
            outputs_df = pd.read_csv(outputs_filepath)
            outputs_file_basename = outputs_filename.split('.')[0]
            outputs_df_dict[outputs_file_basename] = outputs_df

        time_grid = outputs_df_dict.values()[0].t
        delta_t = (time_grid.loc[1] - time_grid.loc[0]) * HOUR_TO_SECOND_CONVERSION_FACTOR

        print 'Retrieve outputs dataframes from precedent simulation run... DONE!'
    else:
        delta_t = simulation_.delta_t

    print 'Compute the post-processing...'

    axes_postprocessing_file_basename = AXES_POSTPROCESSING_FILENAME.split('.')[0]
    hiddenzones_postprocessing_file_basename = HIDDENZONES_POSTPROCESSING_FILENAME.split('.')[0]
    organs_postprocessing_file_basename = ORGANS_POSTPROCESSING_FILENAME.split('.')[0]
    elements_postprocessing_file_basename = ELEMENTS_POSTPROCESSING_FILENAME.split('.')[0]
    soils_postprocessing_file_basename = SOILS_POSTPROCESSING_FILENAME.split('.')[0]

    postprocessing_df_dict = {}

    (postprocessing_df_dict[axes_postprocessing_file_basename],
     postprocessing_df_dict[hiddenzones_postprocessing_file_basename],
     postprocessing_df_dict[organs_postprocessing_file_basename],
     postprocessing_df_dict[elements_postprocessing_file_basename],
     postprocessing_df_dict[soils_postprocessing_file_basename]) \
        = cnwheat_postprocessing.postprocessing(axes_df=outputs_df_dict[AXES_OUTPUTS_FILENAME.split('.')[0]],
                                                  hiddenzones_df=outputs_df_dict[HIDDENZONES_OUTPUTS_FILENAME.split('.')[0]],
                                                  organs_df=outputs_df_dict[ORGANS_OUTPUTS_FILENAME.split('.')[0]],
                                                  elements_df=outputs_df_dict[ELEMENTS_OUTPUTS_FILENAME.split('.')[0]],
                                                  soils_df=outputs_df_dict[SOILS_OUTPUTS_FILENAME.split('.')[0]],
                                                  delta_t=delta_t)

    print 'Compute the post-processing... DONE!'

    print 'Write the postprocessing to CSV files...'

    for postprocessing_file_basename, postprocessing_filename in ((axes_postprocessing_file_basename, AXES_POSTPROCESSING_FILENAME),
                                                                    (hiddenzones_postprocessing_file_basename, HIDDENZONES_POSTPROCESSING_FILENAME),
                                                                    (organs_postprocessing_file_basename, ORGANS_POSTPROCESSING_FILENAME),
                                                                    (elements_postprocessing_file_basename, ELEMENTS_POSTPROCESSING_FILENAME),
                                                                    (soils_postprocessing_file_basename, SOILS_POSTPROCESSING_FILENAME)):

        postprocessing_filepath = os.path.join(POSTPROCESSING_DIRPATH, postprocessing_filename)
        postprocessing_df_dict[postprocessing_file_basename].to_csv(postprocessing_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))

    print 'Write the postprocessing to CSV files... DONE!'

if GENERATE_GRAPHS:

    if not RUN_POSTPROCESSING:

        print 'Retrieve last computed post-processing dataframes...'

        postprocessing_df_dict = {}

        for postprocessing_filename in (ORGANS_POSTPROCESSING_FILENAME,
                                         HIDDENZONES_POSTPROCESSING_FILENAME,
                                         ELEMENTS_POSTPROCESSING_FILENAME,
                                         SOILS_POSTPROCESSING_FILENAME):

            postprocessing_filepath = os.path.join(POSTPROCESSING_DIRPATH, postprocessing_filename)
            postprocessing_df = pd.read_csv(postprocessing_filepath)
            postprocessing_file_basename = postprocessing_filename.split('.')[0]
            postprocessing_df_dict[postprocessing_file_basename] = postprocessing_df

        print 'Retrieve last computed post-processing dataframes... DONE!'


    print 'Generate graphs for validation...'

    cnwheat_postprocessing.generate_graphs(hiddenzones_df=postprocessing_df_dict[HIDDENZONES_POSTPROCESSING_FILENAME.split('.')[0]],
                                            organs_df=postprocessing_df_dict[ORGANS_POSTPROCESSING_FILENAME.split('.')[0]],
                                            elements_df=postprocessing_df_dict[ELEMENTS_POSTPROCESSING_FILENAME.split('.')[0]],
                                            soils_df=postprocessing_df_dict[SOILS_POSTPROCESSING_FILENAME.split('.')[0]],
                                            graphs_dirpath=GRAPHS_DIRPATH)

    print 'Generate graphs for validation... DONE!'