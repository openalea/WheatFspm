# -*- coding: latin-1 -*-

'''
    main
    ~~~~

    An example to show how to initialize and run the model CN-Wheat.
    
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

import logging
import logging.config
import json

import pandas as pd

from cnwheat import simulation
from cnwheat import model

DATA_DIRPATH = 'data'

PLANTS_OUTPUTS_FILENAME = 'plants_outputs.csv'
AXES_OUTPUTS_FILENAME = 'axes_outputs.csv'
PHYTOMERS_OUTPUTS_FILENAME = 'phytomers_outputs.csv'
ORGANS_OUTPUTS_FILENAME = 'organs_outputs.csv'

OUTPUT_PRECISION = 2

LOGGING_CONFIG_FILEPATH = os.path.join('..', 'logging.json')

LOGGING_LEVEL = logging.INFO # set to logging.DEBUG for debugging

def setup_logging(config_filepath='logging.json', level=logging.INFO, 
                  log_compartments=False, log_derivatives=False, log_model=False):
    """Setup logging configuration.
    
    :Parameters:

        - `config_filepath` (:class:`str`) - the file path of the logging 
          configuration.
          
        - `level` (:class:`int`) - the global level of the logging. Use either 
          `logging.DEBUG`, `logging.INFO`, `logging.WARNING`, `logging.ERROR` or 
          `logging.CRITICAL`.
        
        - `log_compartments` (:class:`bool`) - if `True`, log the values of the compartments.
          `False` otherwise.
        
        - `log_derivatives` (:class:`bool`) - if `True`, log the values of the derivatives.
          `False` otherwise.
          
        - `log_model` (:class:`bool`) - if `True`, log the messages from :mod:`cnwheat.model`.
          `False` otherwise.
         
    """
    if os.path.exists(config_filepath):
        with open(config_filepath, 'r') as f:
            config = json.load(f)
        logging.config.dictConfig(config)
    else:
        logging.basicConfig()
    root_logger = logging.getLogger()
    root_logger.setLevel(level)
    
    logging.getLogger('cnwheat.compartments').disabled = not log_compartments # set to False to log the compartments
    logging.getLogger('cnwheat.derivatives').disabled = not log_derivatives # set to False to log the derivatives
    logging.getLogger('cnwheat.model').disabled = not log_model # set to False to log the derivatives
        
setup_logging(LOGGING_CONFIG_FILEPATH, LOGGING_LEVEL)


def read_t_data(curr_data_dirpath, data_filename):
    data_filepath = os.path.join(curr_data_dirpath, data_filename)
    return pd.read_csv(data_filepath, sep=None, index_col='t', engine = 'python')


if __name__ == '__main__':
    t = 0
    
    population = model.Population(t)
    
    plant = model.Plant()
    population.plants.append(plant)
    
    axis = model.Axis()
    plant.axes.append(axis)
    
    axis.grains = model.Grains(starch=0, structure=10850, proteins=0)

    axis.roots = model.Roots(mstruct=0.504, sucrose=0, nitrates=250, amino_acids=0)

    axis.phloem = model.Phloem(sucrose=0, amino_acids=0)
    
    phytomer1 = model.Phytomer(index=1)
    
    phytomer1.lamina = model.Lamina(area=0.00346, mstruct=0.14, width= 0.018, height=0.6, 
                                    starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0, 
                                    amino_acids=0, proteins=0)
    
    phytomer1.internode = model.Internode(area=0.0012, mstruct=0.148, width=0.042, height=0.4, 
                                          starch=0, sucrose=0, triosesP=0, fructan=0, 
                                          nitrates=0, amino_acids=0, proteins=0)
    
    phytomer1.sheath = model.Sheath(area=0.0006, mstruct=0.103, width=0.042, height=0.5, 
                                    starch=0, sucrose=0, triosesP=0, fructan=0, 
                                    nitrates=0 , amino_acids=0, proteins=0)
    axis.phytomers.append(phytomer1)
    
    phytomer2 = model.Phytomer(index=2)
    
    phytomer2.peduncle = model.Peduncle(area=0.00155, mstruct=0.168, width= 0.031, height=0.65, 
                                        starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0, 
                                        amino_acids=0, proteins=0)
    axis.phytomers.append(phytomer2)
    
    phytomer3 = model.Phytomer(index=3)
    
    phytomer3.chaff = model.Chaff(area=0.00075, mstruct=0.21, width=0.02, height= 0.7, starch=0, 
                                  sucrose=0, triosesP=0, fructan=0, nitrates=0, amino_acids=0, 
                                  proteins=0)
    axis.phytomers.append(phytomer3)
    
    An_Tr_dict = {}
    for organ_name in ('Chaff', 'Internode', 'Lamina', 'Peduncle', 'Sheath'):
        An_Tr_dict[organ_name] = read_t_data(DATA_DIRPATH, 'An_Tr_%s.csv' % organ_name)

    # initialize the model
    cnwheat_ = simulation.CNWheat(population=population)

    start_time = 0
    stop_time = 48
    time_step = 4
     
    all_plants_df_list = []
    all_axes_df_list = []
    all_phytomers_df_list = []
    all_organs_df_list = []
     
    for t in xrange(start_time, stop_time, time_step):
        # update the population
        population.t = t
        for plant in population.plants:
            for axis in plant.axes:
                for phytomer in axis.phytomers:
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ is not None: 
                            organ.An = An_Tr_dict[organ.__class__.__name__]['An'][t]
                            organ.Tr = An_Tr_dict[organ.__class__.__name__]['Tr'][t]
        # run the model
        all_plants_df, all_axes_df, all_phytomers_df, all_organs_df = cnwheat_.run(start_time=t, stop_time=t+time_step, number_of_output_steps=time_step+1)
        all_plants_df_list.append(all_plants_df)
        all_axes_df_list.append(all_axes_df)
        all_phytomers_df_list.append(all_phytomers_df)
        all_organs_df_list.append(all_organs_df)
    
    global_plants_df = pd.concat(all_plants_df_list, ignore_index=True)
    global_plants_df.drop_duplicates(subset=simulation.CNWheat.PLANTS_INDEXES, inplace=True)
    global_plants_df.to_csv(PLANTS_OUTPUTS_FILENAME, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUT_PRECISION))
    
    global_axes_df = pd.concat(all_axes_df_list, ignore_index=True)
    global_axes_df.drop_duplicates(subset=simulation.CNWheat.AXES_INDEXES, inplace=True)
    global_axes_df.to_csv(AXES_OUTPUTS_FILENAME, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUT_PRECISION))
    
    global_phytomers_df = pd.concat(all_phytomers_df_list, ignore_index=True)
    global_phytomers_df.drop_duplicates(subset=simulation.CNWheat.PHYTOMERS_INDEXES, inplace=True)
    global_phytomers_df.to_csv(PHYTOMERS_OUTPUTS_FILENAME, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUT_PRECISION))
    
    global_organs_df = pd.concat(all_organs_df_list, ignore_index=True)
    global_organs_df.drop_duplicates(subset=simulation.CNWheat.ORGANS_INDEXES, inplace=True)
    global_organs_df.to_csv(ORGANS_OUTPUTS_FILENAME, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUT_PRECISION))
 

