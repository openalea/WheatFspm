# -*- coding: latin-1 -*-

'''
    main
    ~~~~

    An example to show how to:
        
        * initialize and run the model Farquhar-Wheat,
        * format the outputs of Farquhar-Wheat.

    You must first install :mod:`farquharwheat` and its dependencies 
    before running this script with the command `python`.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2015.
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

import pandas as pd

from farquharwheat import model, simulation, converter 

INPUTS_FILENAME = 'inputs.csv'
PAR_FILENAME = 'PAR.csv'

OUTPUTS_FILENAME = 'outputs.csv'

OUTPUTS_PRECISION = 6

if __name__ == '__main__':
    
    # create a simulation
    simulation_ = simulation.Simulation()
    # read inputs from Pandas dataframe
    inputs_df = pd.read_csv(INPUTS_FILENAME)
    # convert the dataframe to simulation inputs format
    inputs = converter.from_dataframe(inputs_df)
    # initialize the simulation with the inputs
    simulation_.initialize(inputs)
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
    # write the dataframe to CSV
    outputs_df.to_csv(OUTPUTS_FILENAME, index=False, na_rep='NA', float_format='%.{}f'.format(OUTPUTS_PRECISION)) 
    
