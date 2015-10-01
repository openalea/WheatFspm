# -*- coding: latin-1 -*-
"""
    test_farquhar_wheat
    ~~~~~~~~~~~~~~~~~~~

    Test the Farquhar-Wheat model (standalone).

    You must first install :mod:`farquharwheat` and its dependencies
    before running this script with the command `python`.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.
    
    .. seealso:: Barillot et al. 2015.
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

from farquharwheat import model, simulation, converter

INPUTS_FILENAME = 'inputs.csv'
PAR_FILENAME = 'PAR.csv'

DESIRED_OUTPUTS_FILENAME = 'desired_outputs.csv'

ACTUAL_OUTPUTS_FILENAME = 'actual_outputs.csv'

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
    for column in ('axis', 'organ', 'element', 'label'):
        if column in desired_data_df.columns:
            del desired_data_df[column]
            del actual_data_df[column]

    # compare to the desired data
    np.testing.assert_allclose(actual_data_df.values, desired_data_df.values, RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE)


def test_run():

    # create a simulation
    simulation_ = simulation.Simulation()
    # read inputs from Pandas dataframe
    inputs_df = pd.read_csv(INPUTS_FILENAME)
    # convert the dataframe to simulation inputs format
    inputs = converter.from_dataframe(inputs_df)
    # initialize the simulation with the inputs
    simulation_.initialize(inputs)
    # convert the inputs to Pandas dataframe
    reconverted_inputs = converter.to_dataframe(simulation_.inputs)
    # compare inputs
    compare_actual_to_desired('.', reconverted_inputs, INPUTS_FILENAME)
    # get the PAR from dataframe
    PAR_df = pd.read_csv(PAR_FILENAME, index_col=converter.ELEMENTS_TOPOLOGY_COLUMNS)
    # compute incident PAR
    PARi_df = PAR_df
    PARi_df.PAR *= 0.9 * 0.95
    PARi_dict = PARi_df.to_dict()['PAR']
    # run the simulation
    simulation_.run(Ta=18.8, ambient_CO2=360, RH=0.68, Ur=3.171, PARi_dict=PARi_dict)
    # convert the outputs to Pandas dataframe
    outputs_df = converter.to_dataframe(simulation_.outputs)
    # compare outputs
    compare_actual_to_desired('.', outputs_df, DESIRED_OUTPUTS_FILENAME, ACTUAL_OUTPUTS_FILENAME)


if __name__ == '__main__':
    test_run()
