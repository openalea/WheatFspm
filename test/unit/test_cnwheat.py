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

PRECISION = 2
RELATIVE_TOLERANCE = 10**-PRECISION
ABSOLUTE_TOLERANCE = RELATIVE_TOLERANCE


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
        actual_output_df.to_csv(actual_output_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(PRECISION))
      
    # compare to the desired output
    np.testing.assert_allclose(actual_output_df.values, desired_output_df.values, RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE)
    
    
def test_run():
    
    organs = []
    # create the chaff
    name='Chaff'
    PAR_df = read_t_data(DATA_DIRPATH, 'PAR_%s.csv' % name)
    chaff = organ.Chaff(area=0.00075, mstruct=0.21, PAR=PAR_df.PAR, storage_0=0, 
                        sucrose_0=0, triosesP_0=0, fructan_0=0, name=name)
    organs.append(chaff)
     
    # create the internode
    name='Internode'
    PAR_df = read_t_data(DATA_DIRPATH, 'PAR_%s.csv' % name)
    internode = organ.Internode(area=0.0004, mstruct=0.18, PAR=PAR_df.PAR, 
                                storage_0=0, sucrose_0=0, triosesP_0=0, 
                                fructan_0=0, name=name)
    organs.append(internode)
    
    # create the lamina
    name = 'Lamina'
    PAR_df = read_t_data(DATA_DIRPATH, 'PAR_%s.csv' % name)
    lamina = organ.Lamina(area=0.0034, mstruct=0.09, PAR=PAR_df.PAR, 
                          storage_0=0, sucrose_0=0, triosesP_0=0, fructan_0=0,
                          name=name)
    organs.append(lamina)
         
    # create the peduncle
    name = 'Peduncle'
    PAR_df = read_t_data(DATA_DIRPATH, 'PAR_%s.csv' % name)
    peduncle = organ.Peduncle(area=0.00155, mstruct=0.168, PAR=PAR_df.PAR, 
                              storage_0=0, sucrose_0=0, triosesP_0=0, 
                              fructan_0=0, name=name)
    organs.append(peduncle)
    
    # create the sheath
    name = 'Sheath'
    PAR_df = read_t_data(DATA_DIRPATH, 'PAR_%s.csv' % name)
    sheath = organ.Sheath(area=0.0005, mstruct=0.069, PAR=PAR_df.PAR, 
                          storage_0=0, sucrose_0=0, triosesP_0=0, 
                          fructan_0=0, name=name)
    organs.append(sheath)
         
    # create the grains
    grains = organ.Grains(storage_0=0, structure_0=10850, name='Grains')
    organs.append(grains)
     
    # create the roots
    roots = organ.Roots(mstruct=0.504, sucrose_0=0, name='Roots')
    organs.append(roots)
     
    # create the phloem
    phloem = organ.Phloem(sucrose_0=0, name='Phloem')
    organs.append(phloem)
    
    # get meteo data
    meteo_df = read_t_data(DATA_DIRPATH, 'meteo.csv')
     
    # run the model
    actual_output_df = cnwheat.run(start_time=0, stop_time=48, number_of_output_steps=7,
                                    organs=organs, 
                                    meteo=meteo_df,
                                    photosynthesis_computation_interval=4)
    
    compare_actual_to_desired(DATA_DIRPATH, actual_output_df)
    

if __name__ == '__main__':
    test_run()


