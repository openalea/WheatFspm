# -*- coding: latin-1 -*-

'''
    main
    ~~~~

    An example to show how to couple CN-Wheat and Farquhar-Wheat.
    
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

import numpy as np
import pandas as pd

from cnwheat import simulation
from cnwheat import model as cnwheat_model
from cnwheat.farquharwheat import model as photosynthesis_model

DATA_DIRPATH = 'data'

OUTPUT_FILEPATH = 'cnwheat_output.csv'

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
    organs = []
    t = 0
    PAR_dict = {}
    # create the chaff
    name='Chaff'
    PAR_dict[name] = read_t_data(DATA_DIRPATH, 'PAR_%s.csv' % name).PAR
    chaff = cnwheat_model.Chaff(t=t, area=0.00075, mstruct=0.21, width=0.02, height= 0.7, 
                                starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0, 
                                amino_acids=0, proteins=0, name=name)
    organs.append(chaff)

    # create the internodes
    name='Internode'
    PAR_dict[name] = read_t_data(DATA_DIRPATH, 'PAR_%s.csv' % name).PAR
    internode = cnwheat_model.Internode(t=t, area=0.0012, mstruct=0.148, width=0.042, height=0.4, 
                                        starch=0, sucrose=0, triosesP=0, fructan=0, 
                                        nitrates=0, amino_acids=0, proteins=0, name=name)
    organs.append(internode)

    # create the laminae
    name='Lamina'
    PAR_dict[name] = read_t_data(DATA_DIRPATH, 'PAR_%s.csv' % name).PAR
    lamina = cnwheat_model.Lamina(t=t, area=0.00346, mstruct=0.14, width= 0.018, height=0.6, 
                                  starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0, 
                                  amino_acids=0, proteins=0, name=name)
    organs.append(lamina)

    # create the peduncles
    name = 'Peduncle'
    PAR_dict[name] = read_t_data(DATA_DIRPATH, 'PAR_%s.csv' % name).PAR
    peduncle = cnwheat_model.Peduncle(t=t, area=0.00155, mstruct=0.168, width= 0.031, height=0.65, 
                                      starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0, 
                                      amino_acids=0, proteins=0, name=name)
    organs.append(peduncle)

    # create the sheaths
    name = 'Sheath'
    PAR_dict[name] = read_t_data(DATA_DIRPATH, 'PAR_%s.csv' % name).PAR
    sheath = cnwheat_model.Sheath(t=t, area=0.0006, mstruct=0.103, width=0.042, height=0.5, 
                                  starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0, 
                                  amino_acids=0, proteins=0, name=name)
    organs.append(sheath)

    # create the grains
    grains = cnwheat_model.Grains(t=t, starch=0, structure=10850, proteins=0, name='Grains')
    organs.append(grains)

    # create the roots
    roots = cnwheat_model.Roots(t=t, mstruct=0.504, sucrose=0, nitrates=250, amino_acids=0, name='Roots')
    organs.append(roots)

    # create the phloem
    phloem = cnwheat_model.Phloem(t=t, sucrose=0, amino_acids=0, name='Phloem')
    organs.append(phloem)

    # get meteo data
    meteo_df = read_t_data(DATA_DIRPATH, 'meteo.csv')
    
    # initialize the model of CN exchanges
    cnwheat_ = simulation.CNWheat(organs=organs)

    # run the models
    start_time = 0
    stop_time = 48
    time_step = 4
    
    output_df_list = []
    
    for t in xrange(start_time, stop_time, time_step):
        # run the model of photosynthesis
        for organ in organs:
            organ.t = t
            if isinstance(organ, cnwheat_model.PhotosyntheticOrgan):
                PAR = PAR_dict[organ.name][t]
                An, Tr = photosynthesis_model.PhotosynthesisModel.calculate_An(t, organ.width, 
                    organ.height, PAR, meteo_df['air_temperature'][t],
                    meteo_df['ambient_CO2'][t], meteo_df['humidity'][t],
                    meteo_df['Wind'][t])
                organ.An = An
                organ.Tr = Tr
        # run the model of CN exchanges
        output_df = cnwheat_.run(start_time=t, stop_time=t+time_step, number_of_output_steps=time_step+1)
        output_df_list.append(output_df)
            
    global_output_df = pd.concat(output_df_list, ignore_index=True)
    global_output_df.drop_duplicates(subset='t', inplace=True)

    global_output_df.to_csv(OUTPUT_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUT_PRECISION))
    

