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
from cnwheat import model

DATA_DIRPATH = 'data'
DESIRED_OUTPUT_FILENAME = 'desired_output.csv'
ACTUAL_OUTPUT_FILENAME = 'actual_output.csv'

PRECISION = 2
RELATIVE_TOLERANCE = 10**-PRECISION
ABSOLUTE_TOLERANCE = RELATIVE_TOLERANCE


def read_t_data(curr_data_dirpath, data_filename):
    data_filepath = os.path.join(curr_data_dirpath, data_filename)
    return pd.read_csv(data_filepath, sep=None, index_col='t', engine = 'python')


def compare_actual_to_desired(DATA_DIRPATH, actual_output_df, save_actual_output=False):
    # read desired output
    desired_output_filepath = os.path.join(DATA_DIRPATH, DESIRED_OUTPUT_FILENAME)
    desired_output_df = pd.read_csv(desired_output_filepath)

    # keep only the rows to test
    actual_output_df = actual_output_df[actual_output_df['t'].isin(desired_output_df['t'])]

    # keep only the columns to test
    actual_output_df = actual_output_df[desired_output_df.columns]
    
    if save_actual_output:
        actual_output_filepath = os.path.join(DATA_DIRPATH, ACTUAL_OUTPUT_FILENAME)
        actual_output_df.to_csv(actual_output_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(PRECISION))

    # compare to the desired output
    np.testing.assert_allclose(actual_output_df.values, desired_output_df.values, RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE)


def test_run():

    organs = []
    t = 0
    An_Tr_dict = {}
    # create the chaff
    name='Chaff'
    An_Tr_dict[name] = read_t_data(DATA_DIRPATH, 'An_Tr_%s.csv' % name)
    chaff = model.Chaff(t=t, area=0.00075, mstruct=0.21, width=0.02, height= 0.7, starch=0, 
                        sucrose=0, triosesP=0, fructan=0, nitrates=0, amino_acids=0, 
                        proteins=0, name=name)
    organs.append(chaff)

    # create the internode
    name='Internode'
    An_Tr_dict[name] = read_t_data(DATA_DIRPATH, 'An_Tr_%s.csv' % name)
    internode = model.Internode(t=t, area=0.0012, mstruct=0.148, width=0.042, height=0.4, 
                                starch=0, sucrose=0, triosesP=0, fructan=0, 
                                nitrates=0, amino_acids=0, proteins=0, name=name)
    organs.append(internode)

    # create the lamina
    name='Lamina'
    An_Tr_dict[name] = read_t_data(DATA_DIRPATH, 'An_Tr_%s.csv' % name)
    lamina = model.Lamina(t=t, area=0.00346, mstruct=0.14, width= 0.018, height=0.6, 
                          starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0, 
                          amino_acids=0, proteins=0, name=name)
    organs.append(lamina)

    # create the peduncle
    name = 'Peduncle'
    An_Tr_dict[name] = read_t_data(DATA_DIRPATH, 'An_Tr_%s.csv' % name)
    peduncle = model.Peduncle(t=t, area=0.00155, mstruct=0.168, width= 0.031, height=0.65, 
                              starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0, 
                              amino_acids=0, proteins=0, name=name)
    organs.append(peduncle)

    # create the sheath
    name = 'Sheath'
    An_Tr_dict[name] = read_t_data(DATA_DIRPATH, 'An_Tr_%s.csv' % name)
    sheath = model.Sheath(t=t, area=0.0006, mstruct=0.103, width=0.042, height=0.5, 
                          starch=0, sucrose=0, triosesP=0, fructan=0, 
                          nitrates=0 , amino_acids=0, proteins=0, name=name)
    organs.append(sheath)

    # create the grains
    grains = model.Grains(t=t, starch=0, structure=10850, proteins=0, name='Grains')
    organs.append(grains)

    # create the roots
    roots = model.Roots(t=t, mstruct=0.504, sucrose=0, nitrates=250, amino_acids=0, name='Roots')
    organs.append(roots)

    # create the phloem
    phloem = model.Phloem(t=t, sucrose=0, amino_acids=0, name='Phloem')
    organs.append(phloem)

    # initialize the model
    cnwheat_ = simulation.CNWheat(organs=organs)

    # run the model
    start_time = 0
    stop_time = 48
    time_step = 4
     
    output_df_list = []
     
    for t in xrange(start_time, stop_time, time_step):
        # update the organs
        for organ in organs:
            organ.t = t
            if isinstance(organ, model.PhotosyntheticOrgan):
                organ.An = An_Tr_dict[organ.name]['An'][t]
                organ.Tr = An_Tr_dict[organ.name]['Tr'][t]
        output_df = cnwheat_.run(start_time=t, stop_time=t+time_step, number_of_output_steps=time_step+1)
        output_df_list.append(output_df)
             
    global_output_df = pd.concat(output_df_list, ignore_index=True)
    global_output_df.drop_duplicates(subset='t', inplace=True)

    compare_actual_to_desired(DATA_DIRPATH, global_output_df, True)


if __name__ == '__main__':
    test_run()


