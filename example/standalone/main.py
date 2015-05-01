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

from farquharwheat import model, simulation 

INPUTS_DIRPATH = 'inputs'

OUTPUTS_DIRPATH = 'outputs'

CSV_OUTPUTS_FILENAME = 'farquharwheat_out.csv'

OUTPUTS_PRECISION = 6


def setup_inputs():
    """Setup the inputs.
    """
    return \
    {(1, 'MS', 9, 'blade', 'visible'): {'Na': 0.7912279999999999,
                                        'STAR': 1.91e-06,
                                        'organ_height': 32.12856529091365,
                                        'organ_type': 'blade',
                                        'organ_width': 0.018000000000000002},
     (1, 'MS', 9, 'sheath', 'visible'): {'Na': 1.0547030000000002,
                                         'STAR': 3.2000000000000006e-07,
                                         'organ_height': 18.58919516066019,
                                         'organ_type': 'sheath',
                                         'organ_width': 0.3609132887878119},
     (1, 'MS', 10, 'blade', 'visible'): {'Na': 1.15,
                                         'STAR': 6.210000000000001e-06,
                                         'organ_height': 42.96438200561979,
                                         'organ_type': 'blade',
                                         'organ_width': 0.018000000000000002},
     (1, 'MS', 10, 'internode', 'hidden'): {'Na': 2.26,
                                            'STAR': 3.200000000000001e-08,
                                            'organ_height': 18.58919516066019,
                                            'organ_type': 'internode',
                                            'organ_width': 0.22541436464088413},
     (1, 'MS', 10, 'sheath', 'visible'): {'Na': 1.96,
                                          'STAR': 1.1100000000000002e-06,
                                          'organ_height': 29.71785793205823,
                                          'organ_type': 'sheath',
                                          'organ_width': 0.3381029382567559},
     (1, 'MS', 11, 'blade', 'visible'): {'Na': 1.8971099999999999,
                                         'STAR': 2.61e-05,
                                         'organ_height': 59.86435547731891,
                                         'organ_type': 'blade',
                                         'organ_width': 0.018000000000000002},
     (1, 'MS', 11, 'internode', 'hidden'): {'Na': 1.596103,
                                            'STAR': 1.1100000000000003e-07,
                                            'organ_height': 29.71785793205823,
                                            'organ_type': 'internode',
                                            'organ_width': 0.34},
     (1, 'MS', 11, 'internode', 'visible'): {'Na': 1.584906,
                                             'STAR': 2.07e-06,
                                             'organ_height': 38.76237341840967,
                                             'organ_type': 'internode',
                                             'organ_width': 0.34},
     (1, 'MS', 11, 'sheath', 'visible'): {'Na': 4.446667,
                                          'STAR': 6.5199999999999986e-06,
                                          'organ_height': 48.15501051047851,
                                          'organ_type': 'sheath',
                                          'organ_width': 0.260878274577915},
     (1, 'MS', 12, 'peduncle', 'hidden'): {'Na': 1.7735849999999997,
                                           'STAR': 6.519999999999999e-07,
                                           'organ_height': 48.15501051047851,
                                           'organ_type': 'peduncle',
                                           'organ_width': 0.3},
     (1, 'MS', 12, 'peduncle', 'visible'): {'Na': 1.768919,
                                            'STAR': 1.321e-05,
                                            'organ_height': 58.0490108767635,
                                            'organ_type': 'peduncle',
                                            'organ_width': 0.3},
     (1, 'MS', 13, 'ear', 'visible'): {'Na': 6.7466669999999995,
                                       'STAR': 3.023e-05,
                                       'organ_height': 65.35037211630063,
                                       'organ_type': 'ear',
                                       'organ_width': 1.02}}


if __name__ == '__main__':
    
    # create a simulation
    simulation_ = simulation.Simulation()
    # setup the inputs 
    inputs = setup_inputs()
    # initialize the simulation with the inputs
    simulation_.initialize(inputs)
    # run the simulation
    simulation_.run(Ta=18.8, ambient_CO2=360, RH=0.68, Ur=3.171, PARi=2262400)
    # format the outputs to Pandas dataframe
    outputs_df = simulation_.format_outputs()
    # write the dataframe to CSV
    outputs_df.to_csv(os.path.join(OUTPUTS_DIRPATH, CSV_OUTPUTS_FILENAME), index=False, na_rep='NA', float_format='%.{}f'.format(OUTPUTS_PRECISION)) 
    
