# -*- coding: latin-1 -*-

'''
    main
    ~~~~

    An example to show how to:
        
        * initialize and run the model Farquhar-Wheat,
        * format the outputs of Farquhar-Wheat.

    You must first install :mod:`farquharwheat` and its dependencies 
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

import pandas as pd

from farquharwheat import model, simulation 

INPUTS_DIRPATH = 'inputs'

CSV_INPUTS_FILENAME = 'farquharwheat_in.csv'

OUTPUTS_DIRPATH = 'outputs'

CSV_OUTPUTS_FILENAME = 'farquharwheat_out.csv'

OUTPUTS_PRECISION = 6

if __name__ == '__main__':
    
    # create a simulation
    simulation_ = simulation.Simulation()
    # read inputs from Pandas dataframe
    inputs = pd.read_csv(os.path.join(INPUTS_DIRPATH, CSV_INPUTS_FILENAME))
    # initialize the simulation with the inputs
    simulation_.initialize(inputs)
    # run the simulation
    simulation_.run(Ta=18.8, ambient_CO2=360, RH=0.68, Ur=3.171, PARi=2262400)
    # format the outputs to Pandas dataframe
    outputs_df = simulation_.format_outputs()
    # write the dataframe to CSV
    outputs_df.to_csv(os.path.join(OUTPUTS_DIRPATH, CSV_OUTPUTS_FILENAME), index=False, na_rep='NA', float_format='%.{}f'.format(OUTPUTS_PRECISION)) 
    
