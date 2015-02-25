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
import time, datetime

import logging
import logging.config
import json

import pandas as pd

from cnwheat import simulation
from cnwheat import model

t0 = time.time()

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

    population = model.Population()

    plant = model.Plant(index=1)
    population.plants.append(plant)

    axis = model.Axis(axis_type='MS', index=0)
    plant.axes.append(axis)

    axis.grains = model.Grains(starch=0, structure=10850, proteins=1180)

    axis.roots = model.Roots(mstruct=0.504, sucrose=0, nitrates=0, amino_acids=0)

    axis.phloem = model.Phloem(sucrose=0, amino_acids=0)

    # Phytomer 1
    phytomer1 = model.Phytomer(index=1)

    phytomer1.lamina = model.Lamina()
    lamina_element = model.LaminaElement(area=0.00346, mstruct=0.14, width= 0.018, height=0.6,
                                    starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0,
                                    amino_acids=0, proteins=2600)
    phytomer1.lamina.elements.append(lamina_element)

    phytomer1.sheath = model.Sheath()
    sheath_element = model.SheathElement(area=0.0006, mstruct=0.103, width=0.042, height=0.5,
                                    starch=0, sucrose=0, triosesP=0, fructan=0,
                                    nitrates=0 , amino_acids=0, proteins=325)
    phytomer1.sheath.elements.append(sheath_element)

    # Internode enclosed
    phytomer1.internode = model.Internode()
    internode_element1 = model.InternodeElement(area=0.0012, mstruct=0.148, width=0.042, height=0.4,
                                          starch=0, sucrose=0, triosesP=0, fructan=0,
                                          nitrates=0, amino_acids=0, proteins=605, index="Enclosed")
    phytomer1.internode.elements.append(internode_element1)

    # Internode exposed
    internode_element2 = model.InternodeElement(area=0.0003, mstruct=0.04, width=0.042, height=0.3,
                                          starch=0, sucrose=0, triosesP=0, fructan=0,
                                          nitrates=0, amino_acids=0, proteins=500, index="Exposed")
    phytomer1.internode.elements.append(internode_element2)

    axis.phytomers.append(phytomer1)
    print phytomer1.lamina.elements
    # Phytomer 2
    phytomer2 = model.Phytomer(index=2)

    phytomer2.lamina = model.Lamina()
    lamina_element = model.LaminaElement(area=0.0034, mstruct=0.09, width= 0.014, height=0.38,
                                    starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0,
                                    amino_acids=0, proteins=1420)
    print phytomer1.lamina.elements
    phytomer2.lamina.elements.append(lamina_element)
    print phytomer1.lamina.elements

    phytomer2.sheath = model.Sheath()
    sheath_element = model.SheathElement(area=0.0005, mstruct=0.069, width=0.043, height=0.3,
                                    starch=0, sucrose=0, triosesP=0, fructan=0,
                                    nitrates=0 , amino_acids=0, proteins=90)
    phytomer2.sheath.elements.append(sheath_element)

    phytomer2.internode = model.Internode()
    internode_element = model.InternodeElement(area=0.0004, mstruct=0.18, width=0.043, height=0.18,
                                          starch=0, sucrose=0, triosesP=0, fructan=0,
                                          nitrates=0, amino_acids=0, proteins=140)
    phytomer2.internode.elements.append(internode_element)

    axis.phytomers.append(phytomer2)

    # Phytomer 3
    phytomer3 = model.Phytomer(index=3)

    phytomer3.lamina = model.Lamina()
    lamina_element = model.LaminaElement(area=0.00228, mstruct=0.05, width= 0.0125, height=0.24,
                                    starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0,
                                    amino_acids=0, proteins=590)
    phytomer3.lamina.elements.append(lamina_element)

    phytomer3.sheath = model.Sheath()
    sheath_element = model.SheathElement(area=0.0004, mstruct=0.043, width=0.04, height=0.18,
                                    starch=0, sucrose=0, triosesP=0, fructan=0,
                                    nitrates=0 , amino_acids=0, proteins=10)
    phytomer3.sheath.elements.append(sheath_element)

    phytomer3.internode = model.Internode()
    internode_element = model.InternodeElement(area=0.00025, mstruct=0.154, width=0.04, height=0.08,
                                          starch=0, sucrose=0, triosesP=0, fructan=0,
                                          nitrates=0, amino_acids=0, proteins=0)
    phytomer3.internode.elements.append(internode_element)

    axis.phytomers.append(phytomer3)

    # Phytomer 4 (reproductive)
    phytomer4 = model.Phytomer(index=4)

    # Enclosed peduncle
    phytomer4.peduncle = model.Peduncle()
    peduncle_element1 = model.PeduncleElement(area=0.00155, mstruct=0.168, width= 0.031, height=0.65,
                                        starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0,
                                        amino_acids=0, proteins=1240, index="Enclosed")
    phytomer4.peduncle.elements.append(peduncle_element1)

    # Exposed peduncle
    peduncle_element2 = model.PeduncleElement(area=0.00085, mstruct=0.089, width= 0.031, height=0.5,
                                        starch=0, sucrose=0, triosesP=0, fructan=0, nitrates=0,
                                        amino_acids=0, proteins=1000, index="Exposed")
    phytomer4.peduncle.elements.append(peduncle_element2)

    phytomer4.chaff = model.Chaff()
    chaff_element = model.ChaffElement(area=0.00075, mstruct=0.21, width=0.02, height= 0.7, starch=0,
                                  sucrose=0, triosesP=0, fructan=0, nitrates=0, amino_acids=0,
                                  proteins=1800)
    phytomer4.chaff.elements.append(chaff_element)

    axis.phytomers.append(phytomer4)

    An_Tr_dict = {}
    for organ_name in ('Chaff', 'Internode', 'Lamina', 'Peduncle', 'Sheath'):
        An_Tr_dict[organ_name] = read_t_data(DATA_DIRPATH, 'An_Tr_%s.csv' % organ_name)

    # initialize the model
    cnwheat_ = simulation.CNWheat(population=population)

    start_time = 0
    stop_time = 500 # 960
    time_step = 125 # 241

    all_plants_df_list = []
    all_axes_df_list = []
    all_phytomers_df_list = []
    all_organs_df_list = []
    all_elements_df_list = []

    for t in xrange(start_time, stop_time, time_step):
        print 't=', t
        # update the population
        for plant in population.plants:
            for axis in plant.axes:
                for phytomer in axis.phytomers:
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ is not None:
                            for photosynthetic_organ_element in organ.elements:
                                print organ.__class__.__name__, phytomer.index
                                photosynthetic_organ_element.An = An_Tr_dict[organ.__class__.__name__]['An'][t]
                                photosynthetic_organ_element.Tr = An_Tr_dict[organ.__class__.__name__]['Tr'][t]
        # run the model
        print 'model Run'
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