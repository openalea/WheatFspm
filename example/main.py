# -*- coding: latin-1 -*-

import pandas as pd

from farquharwheat import simulation, converter

'''
    main
    ~~~~

    An example to show how to:

        * initialize and run the model Farquhar-Wheat,
        * format the outputs of Farquhar-Wheat.

    You must first install :mod:`farquharwheat` and its dependencies
    before running this script with the command `python`.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: see LICENSE for details.

'''

'''
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
'''


INPUTS_ELEMENT_FILENAME = 'elements_inputs.csv'
INPUTS_SAM_FILENAME = 'SAMs_inputs.csv'
OUTPUTS_FILENAME = 'outputs.csv'

OUTPUTS_PRECISION = 6

if __name__ == '__main__':

    # create a simulation and a converter
    simulation_ = simulation.Simulation()
    # read inputs from Pandas dataframe
    elements_inputs_df = pd.read_csv(INPUTS_ELEMENT_FILENAME)
    SAMs_inputs_df = pd.read_csv(INPUTS_SAM_FILENAME)
    # convert the dataframe to simulation inputs format
    inputs = converter.from_dataframe(elements_inputs_df, SAMs_inputs_df)
    # initialize the simulation with the inputs
    simulation_.initialize(inputs)
    # run the simulation
    simulation_.run(Ta=18.8, ambient_CO2=360, RH=0.530000, Ur=2.200000)
    # convert the outputs to Pandas dataframe
    outputs_df = converter.to_dataframe(simulation_.outputs)
    # write the dataframe to CSV
    outputs_df.to_csv(OUTPUTS_FILENAME, index=False, na_rep='NA', float_format='%.{}f'.format(OUTPUTS_PRECISION))
