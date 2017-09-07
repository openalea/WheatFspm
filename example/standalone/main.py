# -*- coding: latin-1 -*-

'''
    main
    ~~~~

    An example to show how to run the model CN-Wheat, compute the post-processings, and generate the plots for validation.

    You must first install :mod:`cnwheat` (and add it to your PYTHONPATH)
    before running this script with the command `python`.

    :copyright: Copyright 2014-2017 INRA-ECOSYS, see AUTHORS.
    :license: CeCILL-C, see LICENSE for details.
    
    **Acknowledgments**: The research leading these results has received funding through the 
    Investment for the Future programme managed by the Research National Agency 
    (BreedWheat project ANR-10-BTBR-03).
    
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

import os

import pandas as pd

from respiwheat import model as respiwheat_model

from cnwheat import simulation as cnwheat_simulation, converter as cnwheat_converter, tools as cnwheat_tools, postprocessings as cnwheat_postprocessings


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

# CSV file paths to save the outputs of the model in
AXES_OUTPUTS_FILENAME = 'axes_outputs.csv'
ORGANS_OUTPUTS_FILENAME = 'organs_outputs.csv'
HIDDENZONES_OUTPUTS_FILENAME = 'hiddenzones_outputs.csv'
ELEMENTS_OUTPUTS_FILENAME = 'elements_outputs.csv'
SOILS_OUTPUTS_FILENAME = 'soils_outputs.csv'

# CSV file paths to save the post-processings of the model in
AXES_POSTPROCESSINGS_FILENAME = 'axes_postprocessings.csv'
ORGANS_POSTPROCESSINGS_FILENAME = 'organs_postprocessings.csv'
HIDDENZONES_POSTPROCESSINGS_FILENAME = 'hiddenzones_postprocessings.csv'
ELEMENTS_POSTPROCESSINGS_FILENAME = 'elements_postprocessings.csv'
SOILS_POSTPROCESSINGS_FILENAME = 'soils_postprocessings.csv'

# the path of the directory to save the generated graphs in
GRAPHS_DIRPATH = 'graphs'

# culm density (culm m-2)
CULM_DENSITY = {1:410}

# precision of floats used to write and format the output CSV files
OUTPUTS_PRECISION = 6

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
                            


def main(stop_time, run_simu=True, run_postprocessings=True, generate_graphs=True):
    
    if run_simu:


        print 'Prepare the simulation...'
        
        # create the simulation
        simulation_ = cnwheat_simulation.Simulation(respiration_model=respiwheat_model, delta_t=3600, culm_density=CULM_DENSITY)
        
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
        time_step = 4
        time_grid = xrange(start_time, stop_time+time_step, time_step)
        
        # force the senescence and photosynthesis of the population
        force_senescence_and_photosynthesis(0, population, senescence_roots_data_grouped, senescence_elements_data_grouped, photosynthesis_elements_data_grouped)
    
        #TODO: shouldn't we call "simulation_.initialize(population, soils)" here, to reinitialize the simulation from forced population and soils?
        
        print 'Prepare the simulation... DONE!'
        
        
        print 'Run the simulation...'
        
        for t in time_grid:
            
            if t > 0:
                # run the model of CN exchanges ; the population is internally updated by the model
                print '\tt =', t
                simulation_.run()
        
            # convert model outputs to dataframes
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
                #TODO: shouldn't we call "simulation_.initialize(population, soils)" here, to reinitialize the simulation from forced population and soils?
            
        print 'Run the simulation... DONE!'
        
        
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



    if run_postprocessings:
        
        if not run_simu:
            
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
                
        print 'Compute the post-processings...'

        axes_postprocessings_file_basename = AXES_POSTPROCESSINGS_FILENAME.split('.')[0]
        hiddenzones_postprocessings_file_basename = HIDDENZONES_POSTPROCESSINGS_FILENAME.split('.')[0]
        organs_postprocessings_file_basename = ORGANS_POSTPROCESSINGS_FILENAME.split('.')[0]
        elements_postprocessings_file_basename = ELEMENTS_POSTPROCESSINGS_FILENAME.split('.')[0]
        soils_postprocessings_file_basename = SOILS_POSTPROCESSINGS_FILENAME.split('.')[0]

        postprocessings_df_dict = {}
        
        (postprocessings_df_dict[axes_postprocessings_file_basename], 
         postprocessings_df_dict[hiddenzones_postprocessings_file_basename], 
         postprocessings_df_dict[organs_postprocessings_file_basename], 
         postprocessings_df_dict[elements_postprocessings_file_basename], 
         postprocessings_df_dict[soils_postprocessings_file_basename]) \
            = cnwheat_postprocessings.postprocessings(axes_df=outputs_df_dict[AXES_OUTPUTS_FILENAME.split('.')[0]], 
                                                      hiddenzones_df=outputs_df_dict[HIDDENZONES_OUTPUTS_FILENAME.split('.')[0]], 
                                                      organs_df=outputs_df_dict[ORGANS_OUTPUTS_FILENAME.split('.')[0]], 
                                                      elements_df=outputs_df_dict[ELEMENTS_OUTPUTS_FILENAME.split('.')[0]], 
                                                      soils_df=outputs_df_dict[SOILS_OUTPUTS_FILENAME.split('.')[0]], 
                                                      delta_t=delta_t)
            
        print 'Compute the post-processings... DONE!'
        
        
        print 'Write the postprocessings to CSV files...'
        
        for postprocessings_file_basename, postprocessings_filename in ((axes_postprocessings_file_basename, AXES_POSTPROCESSINGS_FILENAME),
                                                                        (hiddenzones_postprocessings_file_basename, HIDDENZONES_POSTPROCESSINGS_FILENAME),
                                                                        (organs_postprocessings_file_basename, ORGANS_POSTPROCESSINGS_FILENAME),
                                                                        (elements_postprocessings_file_basename, ELEMENTS_POSTPROCESSINGS_FILENAME),
                                                                        (soils_postprocessings_file_basename, SOILS_POSTPROCESSINGS_FILENAME)):
            
            postprocessings_filepath = os.path.join(OUTPUTS_DIRPATH, postprocessings_filename)
            postprocessings_df_dict[postprocessings_file_basename].to_csv(postprocessings_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))
            
        print 'Write the postprocessings to CSV files... DONE!'



    if generate_graphs:
        
        if not run_postprocessings:
            
            print 'Retrieve last computed post-processings dataframes...'
            
            postprocessings_df_dict = {}
            
            for postprocessings_filename in (ORGANS_POSTPROCESSINGS_FILENAME,
                                             HIDDENZONES_POSTPROCESSINGS_FILENAME,
                                             ELEMENTS_POSTPROCESSINGS_FILENAME,
                                             SOILS_POSTPROCESSINGS_FILENAME):
                
                postprocessings_filepath = os.path.join(OUTPUTS_DIRPATH, postprocessings_filename)
                postprocessings_df = pd.read_csv(postprocessings_filepath)
                postprocessings_file_basename = postprocessings_filename.split('.')[0]
                postprocessings_df_dict[postprocessings_file_basename] = postprocessings_df
                
            print 'Retrieve last computed post-processings dataframes... DONE!'    
        
        
        print 'Generate graphs for validation...'
        
        cnwheat_postprocessings.generate_graphs(hiddenzones_df=postprocessings_df_dict[HIDDENZONES_POSTPROCESSINGS_FILENAME.split('.')[0]],
                                                organs_df=postprocessings_df_dict[ORGANS_POSTPROCESSINGS_FILENAME.split('.')[0]],
                                                elements_df=postprocessings_df_dict[ELEMENTS_POSTPROCESSINGS_FILENAME.split('.')[0]],
                                                soils_df=postprocessings_df_dict[SOILS_POSTPROCESSINGS_FILENAME.split('.')[0]],
                                                graphs_dirpath=GRAPHS_DIRPATH)
        
        print 'Generate graphs for validation... DONE!'
        

if __name__ == '__main__':
    main(48, run_simu=True, run_postprocessings=True, generate_graphs=True)
    
