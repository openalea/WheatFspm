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

PRECISION = 2
RELATIVE_TOLERANCE = 10**-PRECISION
ABSOLUTE_TOLERANCE = RELATIVE_TOLERANCE


def read_t_data(curr_data_dirpath, data_filename):
    data_filepath = os.path.join(curr_data_dirpath, data_filename)
    return pd.read_csv(data_filepath, sep=None, index_col='t', engine = 'python')


def compare_actual_to_desired(DATA_DIRPATH, actual_output_df, desired_output_filename, actual_output_filename, save_actual_output=False):
    # read desired output
    desired_output_filepath = os.path.join(DATA_DIRPATH, desired_output_filename)
    desired_output_df = pd.read_csv(desired_output_filepath)

    # keep only the rows to test
    actual_output_df = actual_output_df[actual_output_df['t'].isin(desired_output_df['t'])]

    # keep only the columns to test
    actual_output_df = actual_output_df[desired_output_df.columns]
    
    if save_actual_output:
        actual_output_filepath = os.path.join(DATA_DIRPATH, actual_output_filename)
        actual_output_df.to_csv(actual_output_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(PRECISION))
    
    # keep only numerical data
    for column in ('axis', 'organ', 'exposed'):
        if column in desired_output_df.columns:
            del desired_output_df[column]
            del actual_output_df[column]
    
    # compare to the desired output
    np.testing.assert_allclose(actual_output_df.values, desired_output_df.values, RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE)


def test_run():

    population = model.Population()
    
    plant = model.Plant()
    population.plants.append(plant)
    
    axis = model.Axis()
    plant.axes.append(axis)
    
    axis.grains = model.Grains(starch=0, structure=10850, proteins=0)

    axis.roots = model.Roots(mstruct=0.504, sucrose=0, nitrates=250, amino_acids=0)

    axis.phloem = model.Phloem(sucrose=0, amino_acids=0)
    
    phytomer1 = model.Phytomer(index=1)
    
    phytomer1.lamina = model.Lamina()
    
    lamina_element = model.LaminaElement(area=0.00346, mstruct=0.14, width= 0.018, 
                                    height=0.6, starch=0, sucrose=0, triosesP=0, 
                                    fructan=0, nitrates=0, amino_acids=0, proteins=0)
    
    phytomer1.lamina.elements.append(lamina_element)
    
    phytomer1.internode = model.Internode()
    
    internode_element = model.InternodeElement(area=0.0012, mstruct=0.148, width=0.042, 
                                        height=0.4, starch=0, sucrose=0, triosesP=0, 
                                        fructan=0, nitrates=0, amino_acids=0, proteins=0)
    
    phytomer1.internode.elements.append(internode_element)
    
    phytomer1.sheath = model.Sheath()
    
    sheath_element = model.SheathElement(area=0.0006, mstruct=0.103, width=0.042, 
                                    height=0.5, starch=0, sucrose=0, triosesP=0, 
                                    fructan=0, nitrates=0 , amino_acids=0, proteins=0)
    
    phytomer1.sheath.elements.append(sheath_element)
    
    axis.phytomers.append(phytomer1)
    
    phytomer2 = model.Phytomer(index=2)
    
    phytomer2.peduncle = model.Peduncle()
    
    peduncle_element = model.PeduncleElement(area=0.00155, mstruct=0.168, width= 0.031, 
                                        height=0.65, starch=0, sucrose=0, triosesP=0, 
                                        fructan=0, nitrates=0, amino_acids=0, proteins=0)
    
    phytomer2.peduncle.elements.append(peduncle_element)
    
    axis.phytomers.append(phytomer2)
    
    phytomer3 = model.Phytomer(index=3)
    
    phytomer3.chaff = model.Chaff()
    
    chaff_element = model.ChaffElement(area=0.00075, mstruct=0.21, width=0.02, 
                                       height= 0.7, starch=0, sucrose=0, triosesP=0, 
                                       fructan=0, nitrates=0, amino_acids=0, proteins=0)
    
    phytomer3.chaff.elements.append(chaff_element)
    
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
    all_elements_df_list = []
     
    for t in xrange(start_time, stop_time, time_step):
        # update the population
        for plant in population.plants:
            for axis in plant.axes:
                for phytomer in axis.phytomers:
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ is not None: 
                            for photosynthetic_organ_element in organ.elements:
                                photosynthetic_organ_element.An = An_Tr_dict[organ.__class__.__name__]['An'][t]
                                photosynthetic_organ_element.Tr = An_Tr_dict[organ.__class__.__name__]['Tr'][t]
        # run the model
        all_plants_df, all_axes_df, all_phytomers_df, all_organs_df, all_elements_df = cnwheat_.run(start_time=t, stop_time=t+time_step, number_of_output_steps=time_step+1)
        all_plants_df_list.append(all_plants_df)
        all_axes_df_list.append(all_axes_df)
        all_phytomers_df_list.append(all_phytomers_df)
        all_organs_df_list.append(all_organs_df)
        all_elements_df_list.append(all_elements_df)
    
    global_plants_df = pd.concat(all_plants_df_list, ignore_index=True)
    global_plants_df.drop_duplicates(subset=simulation.CNWheat.PLANTS_INDEXES, inplace=True)
    compare_actual_to_desired(DATA_DIRPATH, global_plants_df, DESIRED_PLANTS_OUTPUTS_FILENAME, ACTUAL_PLANTS_OUTPUTS_FILENAME, True)
    
    global_axes_df = pd.concat(all_axes_df_list, ignore_index=True)
    global_axes_df.drop_duplicates(subset=simulation.CNWheat.AXES_INDEXES, inplace=True)
    compare_actual_to_desired(DATA_DIRPATH, global_axes_df, DESIRED_AXES_OUTPUTS_FILENAME, ACTUAL_AXES_OUTPUTS_FILENAME, True)
    
    global_phytomers_df = pd.concat(all_phytomers_df_list, ignore_index=True)
    global_phytomers_df.drop_duplicates(subset=simulation.CNWheat.PHYTOMERS_INDEXES, inplace=True)
    compare_actual_to_desired(DATA_DIRPATH, global_phytomers_df, DESIRED_PHYTOMERS_OUTPUTS_FILENAME, ACTUAL_PHYTOMERS_OUTPUTS_FILENAME, True)
    
    global_organs_df = pd.concat(all_organs_df_list, ignore_index=True)
    global_organs_df.drop_duplicates(subset=simulation.CNWheat.ORGANS_INDEXES, inplace=True)
    compare_actual_to_desired(DATA_DIRPATH, global_organs_df, DESIRED_ORGANS_OUTPUTS_FILENAME, ACTUAL_ORGANS_OUTPUTS_FILENAME, True)
    
    global_elements_df = pd.concat(all_elements_df_list, ignore_index=True)
    global_elements_df.drop_duplicates(subset=simulation.CNWheat.ELEMENTS_INDEXES, inplace=True)
    compare_actual_to_desired(DATA_DIRPATH, global_elements_df, DESIRED_ELEMENTS_OUTPUTS_FILENAME, ACTUAL_ELEMENTS_OUTPUTS_FILENAME, True)


if __name__ == '__main__':
    test_run()


