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
from cnwheat import model as cnwheat_model

INPUTS_DIRPATH = 'inputs'

AN_TR_TS_GS_FILENAME = 'An_Tr_Ts_gs.csv'

OUTPUTS_DIRPATH = 'outputs'

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

PRECISION = 6
RELATIVE_TOLERANCE = 10**-PRECISION
ABSOLUTE_TOLERANCE = RELATIVE_TOLERANCE


def read_t_data(curr_data_dirpath, data_filename):
    data_filepath = os.path.join(curr_data_dirpath, data_filename)
    return pd.read_csv(data_filepath, sep=None, index_col='t', engine = 'python')


def compare_actual_to_desired(outputs_dirpath, actual_output_df, desired_output_filename, actual_output_filename, save_actual_output=False):
    # read desired output
    desired_output_filepath = os.path.join(outputs_dirpath, desired_output_filename)
    desired_output_df = pd.read_csv(desired_output_filepath)

    # keep only the rows to test
    actual_output_df = actual_output_df[actual_output_df['t'].isin(desired_output_df['t'])]

    # keep only the columns to test
    actual_output_df = actual_output_df[desired_output_df.columns]

    if save_actual_output:
        actual_output_filepath = os.path.join(outputs_dirpath, actual_output_filename)
        actual_output_df.to_csv(actual_output_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(PRECISION))

    # keep only numerical data
    for column in ('axis', 'organ', 'exposed'):
        if column in desired_output_df.columns:
            del desired_output_df[column]
            del actual_output_df[column]

    # compare to the desired output
    np.testing.assert_allclose(actual_output_df.values, desired_output_df.values, RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE)


def test_run():

    population = cnwheat_model.Population()

    plant = cnwheat_model.Plant(index=1)
    population.plants.append(plant)

    axis = cnwheat_model.Axis(axis_type='MS', index=0)
    plant.axes.append(axis)

    axis.grains = cnwheat_model.Grains(starch=0, structure=10850, proteins=170)

    axis.roots = cnwheat_model.Roots(mstruct=0.504, sucrose=0, nitrates=0, amino_acids=0)

    axis.phloem = cnwheat_model.Phloem(sucrose=0, amino_acids=0)

    # Phytomer 1
    phytomer1 = cnwheat_model.Phytomer(index=1)

    phytomer1.lamina = cnwheat_model.Lamina()
    lamina_element = cnwheat_model.LaminaElement(area=0.00346, mstruct=0.14, width= 0.018, height=0.6,
                                    starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0,
                                    amino_acids=0, proteins=380)
    phytomer1.lamina.elements.append(lamina_element)

    phytomer1.sheath = cnwheat_model.Sheath()
    sheath_element = cnwheat_model.SheathElement(area=0.0006, mstruct=0.103, width=0.0011, height=0.5,
                                    starch=0, sucrose=0, triosesP=0, fructan=0,
                                    nitrates=0 , amino_acids=0, proteins=130)
    phytomer1.sheath.elements.append(sheath_element)

    phytomer1.internode = cnwheat_model.Internode()
    internode_element1 = cnwheat_model.InternodeElement(area=0.0012, mstruct=0.148, width=0.00257, height=0.3,
                                          starch=0, sucrose=0, triosesP=0, fructan=0,
                                          nitrates=0, amino_acids=0, proteins=20, index=1, exposed=False)
    phytomer1.internode.elements.append(internode_element1)
    axis.phytomers.append(phytomer1)

    # Phytomer 4
    phytomer4 = cnwheat_model.Phytomer(index=4)
    phytomer4.peduncle = cnwheat_model.Peduncle()
    peduncle_element1 = cnwheat_model.PeduncleElement(area=0.00155, mstruct=0.168, width= 0.00349, height=0.65,
                                        starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0,
                                        amino_acids=0, proteins=30, index=1, exposed=False)
    phytomer4.peduncle.elements.append(peduncle_element1)

    # Exposed peduncle
    peduncle_element2 = cnwheat_model.PeduncleElement(area=0.00085, mstruct=0.089, width= 0.00349, height=0.5,
                                        starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0,
                                        amino_acids=0, proteins=180, index=2)
    phytomer4.peduncle.elements.append(peduncle_element2)
    axis.phytomers.append(phytomer4)

    # Phytomer 5
    phytomer5 = cnwheat_model.Phytomer(index=5)
    phytomer5.chaff = cnwheat_model.Chaff()
    chaff_element = cnwheat_model.ChaffElement(area=0.00075, mstruct=0.21, width=0.00265, height= 0.7, starch=0,
                                  sucrose=0, triosesP=0, fructan=0, nitrates=0, amino_acids=0,
                                  proteins=260)
    phytomer5.chaff.elements.append(chaff_element)
    axis.phytomers.append(phytomer5)

    # Get assimilation, transpiration, organ temperature and stomatal conductance data
    An_Tr_Ts_gs_filepath = os.path.join(INPUTS_DIRPATH, AN_TR_TS_GS_FILENAME)
    An_Tr_Ts_gs_df = pd.read_csv(An_Tr_Ts_gs_filepath)
    An_Tr_Ts_gs_grouped = An_Tr_Ts_gs_df.groupby(simulation.CNWheat.ELEMENTS_INDEXES)

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
        population.t = t
        for plant in population.plants:
            pid = plant.index
            for axis in plant.axes:
                axid = axis.id
                for phytomer in axis.phytomers:
                    phytoid = phytomer.index
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ is not None:
                            orgid = organ.__class__.__name__
                            for photosynthetic_organ_element in organ.elements:
                                eltid = photosynthetic_organ_element.index
                                exposed = photosynthetic_organ_element.exposed
                                group = An_Tr_Ts_gs_grouped.get_group((t, pid, axid, phytoid, orgid, eltid, exposed))
                                row_index = group.first_valid_index()
                                photosynthetic_organ_element.An = group.An[row_index]
                                photosynthetic_organ_element.Tr = group.Tr[row_index]
                                photosynthetic_organ_element.Ts = group.Ts[row_index]
                                photosynthetic_organ_element.gs = group.gs[row_index]

        # run the model of CN exchanges ; the population is internally updated by the model of CN exchanges
        all_plants_df, all_axes_df, all_phytomers_df, all_organs_df, all_elements_df = cnwheat_.run(start_time=t, stop_time=t+time_step, number_of_output_steps=time_step+1)
        all_plants_df_list.append(all_plants_df)
        all_axes_df_list.append(all_axes_df)
        all_phytomers_df_list.append(all_phytomers_df)
        all_organs_df_list.append(all_organs_df)
        all_elements_df_list.append(all_elements_df)

    global_plants_df = pd.concat(all_plants_df_list, ignore_index=True)
    global_plants_df.drop_duplicates(subset=simulation.CNWheat.PLANTS_INDEXES, inplace=True)
    compare_actual_to_desired(OUTPUTS_DIRPATH, global_plants_df, DESIRED_PLANTS_OUTPUTS_FILENAME, ACTUAL_PLANTS_OUTPUTS_FILENAME, True)

    global_axes_df = pd.concat(all_axes_df_list, ignore_index=True)
    global_axes_df.drop_duplicates(subset=simulation.CNWheat.AXES_INDEXES, inplace=True)
    compare_actual_to_desired(OUTPUTS_DIRPATH, global_axes_df, DESIRED_AXES_OUTPUTS_FILENAME, ACTUAL_AXES_OUTPUTS_FILENAME, True)

    global_phytomers_df = pd.concat(all_phytomers_df_list, ignore_index=True)
    global_phytomers_df.drop_duplicates(subset=simulation.CNWheat.PHYTOMERS_INDEXES, inplace=True)
    compare_actual_to_desired(OUTPUTS_DIRPATH, global_phytomers_df, DESIRED_PHYTOMERS_OUTPUTS_FILENAME, ACTUAL_PHYTOMERS_OUTPUTS_FILENAME, True)

    global_organs_df = pd.concat(all_organs_df_list, ignore_index=True)
    global_organs_df.drop_duplicates(subset=simulation.CNWheat.ORGANS_INDEXES, inplace=True)
    compare_actual_to_desired(OUTPUTS_DIRPATH, global_organs_df, DESIRED_ORGANS_OUTPUTS_FILENAME, ACTUAL_ORGANS_OUTPUTS_FILENAME, True)

    global_elements_df = pd.concat(all_elements_df_list, ignore_index=True)
    global_elements_df.drop_duplicates(subset=simulation.CNWheat.ELEMENTS_INDEXES, inplace=True)
    compare_actual_to_desired(OUTPUTS_DIRPATH, global_elements_df, DESIRED_ELEMENTS_OUTPUTS_FILENAME, ACTUAL_ELEMENTS_OUTPUTS_FILENAME, True)


if __name__ == '__main__':
    test_run()