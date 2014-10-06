# -*- coding: latin-1 -*-
"""
    test_cn_wheat
    ~~~~~~~~~~~~~
    
    Test the cnwheat.cnwheat module.  

    CSV files must contain only ASCII characters and ',' as separator.

    Be sure to add the 'cnwheat' directory to your PYTHONPATH before running this script.
    
    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.
"""

import os

import numpy as np
import pandas as pd

from cnwheat import cnwheat
from cnwheat import organ

DATA_DIRPATH = 'data'
DESIRED_OUTPUT_FILENAME = 'desired_output.csv'
ACTUAL_OUTPUT_FILENAME = 'actual_output.csv'


RELATIVE_TOLERANCE = 10e-3
ABSOLUTE_TOLERANCE = 10e-3


def read_t_data(DATA_DIRPATH, data_filename):
    data_filepath = os.path.join(DATA_DIRPATH, data_filename)
    return pd.read_csv(data_filepath, index_col='t')


def compare_actual_to_desired(DATA_DIRPATH, actual_output_df, save_actual_output=False):
    # read desired output
    desired_output_filepath = os.path.join(DATA_DIRPATH, DESIRED_OUTPUT_FILENAME)
    desired_output_df = pd.read_csv(desired_output_filepath)
      
    # keep only the columns to test
    actual_output_df = actual_output_df[desired_output_df.columns]
      
    if save_actual_output:
        actual_output_filepath = os.path.join(DATA_DIRPATH, ACTUAL_OUTPUT_FILENAME)
        actual_output_df.to_csv(actual_output_filepath, na_rep='NA', index=False)
      
    # compare to the desired output
    np.testing.assert_allclose(actual_output_df.values, desired_output_df.values, RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE)
    
    
def test_run():
    
    organs = []
    # create the chaff
    name='chaff'
    PAR_df = read_t_data(DATA_DIRPATH, 'PAR_%s.csv' % name)
    chaff = organ.Chaff(Area=0.00075, Mstruct=0.21, PAR=PAR_df.PAR, STORAGE_0=0, 
                        SUCROSE_0=0, TRIOSESP_0=0, FRUCTAN_0=0, name=name)
    organs.append(chaff)
     
    # create the internode
    name='internode'
    PAR_df = read_t_data(DATA_DIRPATH, 'PAR_%s.csv' % name)
    internode = organ.Internode(Area=0.0004, Mstruct=0.18, PAR=PAR_df.PAR, 
                                STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0, 
                                FRUCTAN_0=0, name=name)
    organs.append(internode)
    
    # create the lamina
    name = 'lamina'
    PAR_df = read_t_data(DATA_DIRPATH, 'PAR_%s.csv' % name)
    lamina = organ.Lamina(Area=0.0034, Mstruct=0.09, PAR=PAR_df.PAR, 
                          STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0, FRUCTAN_0=0,
                          name=name)
    organs.append(lamina)
         
    # create the peduncle
    name = 'peduncle'
    PAR_df = read_t_data(DATA_DIRPATH, 'PAR_%s.csv' % name)
    peduncle = organ.Peduncle(Area=0.00155, Mstruct=0.168, PAR=PAR_df.PAR, 
                              STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0, 
                              FRUCTAN_0=0, name=name)
    organs.append(peduncle)
    
    # create the sheath
    name = 'sheath'
    PAR_df = read_t_data(DATA_DIRPATH, 'PAR_%s.csv' % name)
    sheath = organ.Sheath(Area=0.0005, Mstruct=0.069, PAR=PAR_df.PAR, 
                          STORAGE_0=0, SUCROSE_0=0, TRIOSESP_0=0, 
                          FRUCTAN_0=0, name=name)
    organs.append(sheath)
         
    # create the grains
    grains = organ.Grains(STORAGE_0=0, STRUCTURE_0=10850, name='grains')
    organs.append(grains)
     
    # create the roots
    roots = organ.Roots(Mstruct=0.504, Sucrose_0=0, name='roots')
    organs.append(roots)
     
    # create the phloem
    phloem = organ.Phloem(SUCROSE_0=0, name='phloem')
    
    # get meteo data
    meteo_df = read_t_data(DATA_DIRPATH, 'meteo.csv')
     
    # run the model
    actual_output_df = cnwheat.run(start_time=0, stop_time=48, number_of_output_steps=7,
                                    phloem=phloem,
                                    organs=organs, 
                                    meteo=meteo_df,
                                    photosynthesis_computation_interval=4)
    
    compare_actual_to_desired(DATA_DIRPATH, actual_output_df)
    

if __name__ == '__main__':
    test_run()


