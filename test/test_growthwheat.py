# -*- coding: latin-1 -*-
"""
    test_growthwheat
    ~~~~~~~~~~~~~~~~

    An example to show how to:

        * initialize and run the model Growth-Wheat,
        * format the outputs of Growth-Wheat.

    You must first install :mod:`growthwheat` and its dependencies
    before running this script with the command `python`.

    :copyright: Copyright 2014-2016 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.

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

import os

import numpy as np
import pandas as pd

from growthwheat import model, simulation, converter

INPUTS_DIRPATH = 'inputs'
HIDDENZONE_INPUTS_FILENAME = 'hiddenzones_inputs.csv'
ELEMENT_INPUTS_FILENAME = 'elements_inputs.csv'
ROOT_INPUTS_FILENAME = 'roots_inputs.csv'
SAM_INPUTS_FILENAME = 'SAM_inputs.csv'

OUTPUTS_DIRPATH = 'outputs'
DESIRED_HIDDENZONE_OUTPUTS_FILENAME = 'desired_hiddenzones_outputs.csv'
DESIRED_ELEMENT_OUTPUTS_FILENAME = 'desired_elements_outputs.csv'
DESIRED_ROOT_OUTPUTS_FILENAME = 'desired_roots_outputs.csv'
ACTUAL_HIDDENZONE_OUTPUTS_FILENAME = 'actual_hiddenzones_outputs.csv'
ACTUAL_ELEMENT_OUTPUTS_FILENAME = 'actual_elements_outputs.csv'
ACTUAL_ROOT_OUTPUTS_FILENAME = 'actual_roots_outputs.csv'

PRECISION = 6
RELATIVE_TOLERANCE = 10**-PRECISION
ABSOLUTE_TOLERANCE = RELATIVE_TOLERANCE

def compare_actual_to_desired(data_dirpath, actual_data_df, desired_data_filename, actual_data_filename=None):
    # read desired data
    desired_data_filepath = os.path.join(data_dirpath, desired_data_filename)
    desired_data_df = pd.read_csv(desired_data_filepath)

    if actual_data_filename is not None:
        actual_data_filepath = os.path.join(data_dirpath, actual_data_filename)
        actual_data_df.to_csv(actual_data_filepath, na_rep='NA', index=False)

    # keep only numerical data
    for column in ('axis', 'organ', 'element', 'leaf_is_emerged', 'internode_is_visible','is_growing'):
        if column in desired_data_df.columns:
            assert desired_data_df[column].equals(actual_data_df[column])
            del desired_data_df[column]
            del actual_data_df[column]

    # compare to the desired data
    np.testing.assert_allclose(actual_data_df.values, desired_data_df.values, RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE)


def test_run():
    # read growthwheat inputs at t0
    hiddenzones_inputs_t0 = pd.read_csv(os.path.join(INPUTS_DIRPATH, HIDDENZONE_INPUTS_FILENAME))
    element_inputs_t0 = pd.read_csv(os.path.join(INPUTS_DIRPATH, ELEMENT_INPUTS_FILENAME))
    root_inputs_t0 = pd.read_csv(os.path.join(INPUTS_DIRPATH, ROOT_INPUTS_FILENAME))
    SAM_inputs_t0 = pd.read_csv(os.path.join(INPUTS_DIRPATH, SAM_INPUTS_FILENAME))

    # Create population
    simulation_ = simulation.Simulation(delta_t=3600)
    # convert the dataframe to simulation inputs format
    all_inputs_dict = converter.from_dataframes(hiddenzones_inputs_t0, element_inputs_t0, root_inputs_t0, SAM_inputs_t0)
    hz_inputs_reconverted_df = all_inputs_dict['hiddenzone']
    element_inputs_reconverted_df = all_inputs_dict['elements']
    roots_inputs_reconverted_df = all_inputs_dict['roots']
    # compare inputs
##    compare_actual_to_desired('inputs', hz_inputs_reconverted_df, HIDDENZONE_INPUTS_FILENAME)
##    compare_actual_to_desired('inputs', organ_inputs_reconverted_df, ELEMENT_INPUTS_FILENAME)
##    compare_actual_to_desired('inputs', roots_inputs_reconverted_df, ROOT_INPUTS_FILENAME)
    # initialize the simulation with the inputs
    simulation_.initialize(all_inputs_dict)
    # run the simulation
    simulation_.run()
    # convert the outputs to Pandas dataframe
    hiddenzone_outputs_df, element_outputs_df, root_outputs_df = converter.to_dataframes(simulation_.outputs)
    # compare outputs
    compare_actual_to_desired('outputs', hiddenzone_outputs_df, DESIRED_HIDDENZONE_OUTPUTS_FILENAME, ACTUAL_HIDDENZONE_OUTPUTS_FILENAME)
    compare_actual_to_desired('outputs', element_outputs_df, DESIRED_ELEMENT_OUTPUTS_FILENAME, ACTUAL_ELEMENT_OUTPUTS_FILENAME)
    compare_actual_to_desired('outputs', root_outputs_df, DESIRED_ROOT_OUTPUTS_FILENAME, ACTUAL_ROOT_OUTPUTS_FILENAME)

if __name__ == '__main__':
    test_run()