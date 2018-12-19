# -*- coding: latin-1 -*-
"""
    main
    ~~~~

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

from growthwheat import simulation, model, converter

INPUTS_DIRPATH = 'inputs'

# growthwheat inputs at t0
HIDDENZONE_INPUTS_FILEPATH = os.path.join(INPUTS_DIRPATH, 'hiddenzones_inputs.csv')
ELEMENT_INPUTS_FILEPATH = os.path.join(INPUTS_DIRPATH, 'elements_inputs.csv')
ROOT_INPUTS_FILEPATH = os.path.join(INPUTS_DIRPATH, 'roots_inputs.csv')
SAM_INPUTS_FILEPATH = os.path.join(INPUTS_DIRPATH, 'SAM_inputs.csv')

# growthwheat outputs
OUTPUTS_DIRPATH = 'outputs'
HGZ_OUTPUTS_FILENAME = 'hiddenzones_outputs.csv'
ELEMENT_OUTPUTS_FILENAME = 'elements_outputs.csv'
ROOT_OUTPUTS_FILENAME = 'roots_outputs.csv'

# read growthwheat inputs at t0
hiddenzones_inputs_t0 = pd.read_csv(HIDDENZONE_INPUTS_FILEPATH)
element_inputs_t0 = pd.read_csv(ELEMENT_INPUTS_FILEPATH)
root_inputs_t0 = pd.read_csv(ROOT_INPUTS_FILEPATH)
SAM_inputs_t0 = pd.read_csv(SAM_INPUTS_FILEPATH)

OUTPUTS_PRECISION = 6

if __name__ == '__main__':

    # Create population
    simulation_ = simulation.Simulation(delta_t=3600)
    # convert the dataframe to simulation inputs format
    inputs = converter.from_dataframes(hiddenzones_inputs_t0, element_inputs_t0, root_inputs_t0,SAM_inputs_t0)
    # initialize the simulation with the inputs
    simulation_.initialize(inputs)
    # run the simulation
    simulation_.run()
    # convert the outputs to Pandas dataframe
    hiddenzones_outputs, elements_outputs, roots_outputs = converter.to_dataframes(simulation_.outputs)
    # write the dataframe to CSV
    hiddenzones_outputs.to_csv(os.path.join(OUTPUTS_DIRPATH, HGZ_OUTPUTS_FILENAME), index=False, na_rep='NA', float_format='%.{}f'.format(OUTPUTS_PRECISION))
    elements_outputs.to_csv(os.path.join(OUTPUTS_DIRPATH, ELEMENT_OUTPUTS_FILENAME), index=False, na_rep='NA', float_format='%.{}f'.format(OUTPUTS_PRECISION))
    roots_outputs.to_csv(os.path.join(OUTPUTS_DIRPATH, ROOT_OUTPUTS_FILENAME), index=False, na_rep='NA', float_format='%.{}f'.format(OUTPUTS_PRECISION))