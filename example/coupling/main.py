# -*- coding: latin-1 -*-

'''
    main
    ~~~~

    An example to show to couple Farquhar-Wheat with others models, using MTG format 
    to share data between the models. 

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

from farquharwheat import simulation, converter

INPUTS_DIRPATH = 'inputs'

OUTPUTS_DIRPATH = 'outputs'

CSV_OUTPUTS_FILENAME = 'farquharwheat_out.csv'

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
    
    PROPERTIES_FILENAME = 'properties.csv'
    
    properties_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, PROPERTIES_FILENAME))
    
    ORGANS_MAPPING = {'internode': 'Internode', 'blade': 'Lamina', 'sheath': 'Sheath', 'peduncle': 'Peduncle', 'ear': 'Chaff'}
    
    g, wheat, domain_area, domain, convUnit = initialise_stand(1500)
    
    # add the properties which do not exist yet
    if not 'exposed_area' in g.properties():
        g.add_property('exposed_area')
    if not 'width' in g.properties():
        g.add_property('width')
    if not 'surfacic_nitrogen' in g.properties():
        g.add_property('surfacic_nitrogen')
    
    plants_indexes_in_properties_df = properties_df.plant.unique()
    # traverse the MTG recursively from top
    for plant_vid in g.components_iter(g.root):
        plant_index = int(g.index(plant_vid))
        if plant_index not in plants_indexes_in_properties_df:
            # TODO: remove plant_vid and all its components
            continue
        plant_properties_df = properties_df[properties_df['plant'] == plant_index]
        for axis_vid in g.components_iter(plant_vid):
            axis_id = g.label(axis_vid)
            if axis_id not in plant_properties_df.axis.unique():
                # TODO: remove axis_vid and all its components
                continue
            axis_properties_df = plant_properties_df[plant_properties_df['axis'] == axis_id]
            for metamer_vid in g.components(axis_vid): 
                metamer_index = int(g.index(metamer_vid))
                if metamer_index not in axis_properties_df['metamer'].unique():
                    # TODO: remove metamer_vid and all its components
                    continue
                metamer_properties_df = axis_properties_df[axis_properties_df['metamer'] == metamer_index]
                for organ_vid in g.components_iter(metamer_vid):
                    organ_label = g.label(organ_vid)
                    if organ_label not in ORGANS_MAPPING:
                        continue
                    organ_type = ORGANS_MAPPING[organ_label]
                    if organ_type not in metamer_properties_df['organ'].unique():
                        # TODO: remove organ_vid and all its components
                        continue
                    organ_properties_df = metamer_properties_df[metamer_properties_df['organ'] == organ_type]
                    # width
                    g.property('width')[organ_vid] = organ_properties_df['width'][organ_properties_df.first_valid_index()]
                    for element_vid in g.components_iter(organ_vid):
                        vertex_properties = g.get_vertex_property(element_vid)
                        element_label = g.label(element_vid)
                        if element_label.startswith('Hidden'):
                            element_type = 'enclosed'
                        else:
                            element_type = 'exposed'
                        if element_type not in organ_properties_df['element'].unique():
                            # TODO: remove element_vid and all its components
                            continue
                        element_df = organ_properties_df[organ_properties_df['element'] == element_type]
                        # exposed_area
                        element_area = vertex_properties['area'] / 10000.0 # conversion from cm2 to m2
                        element_STAR = element_df['STAR'][element_df.first_valid_index()]
                        g.property('exposed_area')[element_vid] = element_STAR * element_area
                        g.property('surfacic_nitrogen')[element_vid] = element_df['SLN'][element_df.first_valid_index()]
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
    outputs_df = simulation_.format_outputs()
    # write the dataframe to CSV
    outputs_df.to_csv(os.path.join(OUTPUTS_DIRPATH, CSV_OUTPUTS_FILENAME), index=False, na_rep='NA', float_format='%.{}f'.format(OUTPUTS_PRECISION)) 


