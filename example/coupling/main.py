# -*- coding: latin-1 -*-

'''
    main
    ~~~~

    An example to show how to:
        
        * initialize and run the model Farquhar-Wheat from an MTG,
        * format the outputs of Farquhar-Wheat to both MTG and Pandas dataframes.

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

from farquharwheat import simulation, converter

INPUTS_FILENAME = 'inputs.csv'
OUTPUTS_FILENAME = 'outputs.csv'

OUTPUTS_PRECISION = 6


def setup_MTG():
    """Construct and fill a valid MTG.
    TODO: load the fulfilled MTG from file. 
    """
    import random
    import numpy as np
    import pandas as pd
    from alinea.adel.astk_interface import initialise_stand
    
    random.seed(1234)
    np.random.seed(1234)
    
    inputs_df = pd.read_csv(INPUTS_FILENAME)
    
    g, wheat, domain_area, domain, convUnit = initialise_stand(1500)
    
    # add the properties which do not exist yet
    for property_ in converter.FARQUHARWHEAT_INPUTS:
        if property_ not in g.properties():
            g.add_property(property_)
    
    plants_indexes_in_inputs_df = inputs_df.plant.unique()
    # traverse the MTG recursively from top
    for plant_vid in g.components_iter(g.root):
        plant_index = int(g.index(plant_vid))
        if plant_index not in plants_indexes_in_inputs_df:
            # TODO: remove plant_vid and all its components
            continue
        plant_inputs_df = inputs_df[inputs_df['plant'] == plant_index]
        for axis_vid in g.components_iter(plant_vid):
            axis_id = g.label(axis_vid)
            if axis_id not in plant_inputs_df.axis.unique():
                # TODO: remove axis_vid and all its components
                continue
            axis_inputs_df = plant_inputs_df[plant_inputs_df['axis'] == axis_id]
            for metamer_vid in g.components(axis_vid): 
                metamer_index = int(g.index(metamer_vid))
                if metamer_index not in axis_inputs_df['metamer'].unique():
                    # TODO: remove metamer_vid and all its components
                    continue
                metamer_inputs_df = axis_inputs_df[axis_inputs_df['metamer'] == metamer_index]
                for organ_vid in g.components_iter(metamer_vid):
                    organ_label = g.label(organ_vid)
                    if organ_label not in converter.FARQUHARWHEAT_ORGANS_NAMES:
                        continue
                    if organ_label not in metamer_inputs_df['organ'].unique():
                        # TODO: remove organ_vid and all its components
                        continue
                    organ_inputs_df = metamer_inputs_df[metamer_inputs_df['organ'] == organ_label]
                    for element_vid in g.components_iter(organ_vid):
                        element_label = g.label(element_vid)
                        if element_label not in organ_inputs_df['element'].unique():
                            # TODO: remove element_vid and all its components
                            continue
                        element_df = organ_inputs_df[organ_inputs_df['element'] == element_label]
                        element_series = element_df.loc[element_df.first_valid_index()]
                        for property_ in converter.FARQUHARWHEAT_ELEMENTS_INPUTS:
                            g.property(property_)[element_vid] = element_series[property_]
    return g

if __name__ == '__main__':
    
    # create a simulation
    simulation_ = simulation.Simulation()
    # setup a MTG 
    g = setup_MTG()
    # convert the MTG to Farquhar-Wheat inputs and initialize the simulation
    simulation_.initialize(converter.from_MTG(g))
    # run the simulation
    simulation_.run(Ta=18.8, ambient_CO2=360, RH=0.68, Ur=3.171, PARi=2262400)
    # update the MTG from Farquhar-Wheat outputs
    converter.update_MTG(simulation_.outputs, g)
    # format Farquhar-Wheat outputs to Pandas dataframe
    outputs_df = converter.to_dataframe(simulation_.outputs)
    # write the dataframe to CSV
    outputs_df.to_csv(OUTPUTS_FILENAME, index=False, na_rep='NA', float_format='%.{}f'.format(OUTPUTS_PRECISION)) 


