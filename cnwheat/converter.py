# -*- coding: latin-1 -*-

"""
    cnwheat.converter
    ~~~~~~~~~~~~~~~~~

    The module :mod:`cnwheat.converter` defines functions to:
    
        * convert :class:`dataframes <pandas.DataFrame>` to/from :class:`population <model.Population>`.

    :copyright: Copyright 2015 INRA-EcoSys, see AUTHORS.
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

import warnings

import numpy as np
import pandas as pd

from cnwheat import model, simulation


PLANTS_DATAFRAME_COLUMNS = simulation.Simulation.PLANTS_INPUTS_INDEXES + simulation.Simulation.PLANTS_STATE
AXES_DATAFRAME_COLUMNS = simulation.Simulation.AXES_INPUTS_INDEXES + simulation.Simulation.AXES_STATE
PHYTOMERS_DATAFRAME_COLUMNS = simulation.Simulation.PHYTOMERS_INPUTS_INDEXES + simulation.Simulation.PHYTOMERS_STATE
ORGANS_DATAFRAME_COLUMNS = simulation.Simulation.ORGANS_INPUTS_INDEXES + simulation.Simulation.ORGANS_STATE
ELEMENTS_DATAFRAME_COLUMNS = simulation.Simulation.ELEMENTS_INPUTS_INDEXES + simulation.Simulation.ELEMENTS_STATE


def from_dataframes(plants_dataframe, axes_dataframe, phytomers_dataframe, organs_dataframe, elements_dataframe):
    """
    Construct a :class:`model.Population` from `plants_dataframe`, `axes_dataframe`, `phytomers_dataframe`, 
    `organs_dataframe` and `elements_dataframe`.
    
    :Parameters:
        
        - `plants_dataframe` (:class:`pandas.DataFrame`) - Plants dataframe, with one line by plant.
        
        - `axes_dataframe` (:class:`pandas.DataFrame`) - Axes dataframe, with one line by axis.
        
        - `phytomers_dataframe` (:class:`pandas.DataFrame`) - Phytomers dataframe, with one line by phytomer.
        
        - `organs_dataframe` (:class:`pandas.DataFrame`) - Organs dataframe, with one line by organ.
        
        - `elements_dataframe` (:class:`pandas.DataFrame`) - Elements dataframe, with one line by element.
        
    :Returns:
        A new :class:`population <model.Population>`.
        
    :Returns Type:
        :class:`model.Population` 
        
    """
    
    population = model.Population()
    
    for plant_index in plants_dataframe.plant:
        # create a new plant
        plant = model.Plant()
        population.plants.append(plant)
        curr_axes_dataframe = axes_dataframe[axes_dataframe['plant'] == plant_index]
        for axis_id in curr_axes_dataframe.axis:
            # create a new axis
            axis = model.Axis()
            plant.axes.append(axis)
            
            curr_organs_dataframe = organs_dataframe[(organs_dataframe['plant'] == plant_index) & (organs_dataframe['axis'] == axis_id)]
            
            for axis_attribute_name, axis_attribute_class in (('grains', model.Grains), ('roots', model.Roots), ('soil', model.Soil), ('phloem', model.Phloem)):
                organ_class_name = axis_attribute_class.__name__
                organ_dataframe = curr_organs_dataframe[curr_organs_dataframe['organ'] == organ_class_name]
                # create a new organ
                organ = axis_attribute_class()
                organ_attributes_names = [state_var_name for state_var_name in simulation.Simulation.ORGANS_STATE if hasattr(organ, state_var_name)]
                organ_row = organ_dataframe.loc[organ_dataframe.first_valid_index()]
                organ_attributes_values = organ_row[organ_attributes_names].tolist()
                organ_attributes = dict(zip(organ_attributes_names, organ_attributes_values))
                organ.__dict__.update(organ_attributes)
                organ.initialize()
                setattr(axis, axis_attribute_name, organ)
                    
            curr_phytomers_dataframe = phytomers_dataframe[(phytomers_dataframe['plant'] == plant_index) & (phytomers_dataframe['axis'] == axis_id)]
            for phytomer_index in curr_phytomers_dataframe.phytomer:
                # create a new phytomer
                phytomer = model.Phytomer()
                axis.phytomers.append(phytomer)
                
                curr_organs_dataframe = organs_dataframe[(organs_dataframe['plant'] == plant_index) & (organs_dataframe['axis'] == axis_id) & (organs_dataframe['phytomer'] == phytomer_index)]
                
                for phytomer_attribute_name, phytomer_attribute_class, phytomer_attribute_element_class in \
                    (('chaff', model.Chaff, model.ChaffElement), 
                     ('lamina', model.Lamina, model.LaminaElement), 
                     ('internode', model.Internode, model.InternodeElement), 
                     ('peduncle', model.Peduncle, model.PeduncleElement),
                     ('sheath', model.Sheath, model.SheathElement)):
                    
                    if phytomer_attribute_name == 'peduncle':
                        pass
                    
                    organ_class_name = phytomer_attribute_class.__name__
                    curr_elements_dataframe = elements_dataframe[(elements_dataframe['plant'] == plant_index) & (elements_dataframe['axis'] == axis_id) & (elements_dataframe['phytomer'] == phytomer_index) & (elements_dataframe['organ'] == organ_class_name)]
                
                    if organ_class_name not in curr_organs_dataframe.organ.values and organ_class_name not in curr_elements_dataframe.organ.values:
                        continue
                    # create a new organ
                    organ = phytomer_attribute_class()
                    setattr(phytomer, phytomer_attribute_name, organ)
                    
                    for organ_attribute_name, organ_attribute_type in (('enclosed_element', 'enclosed'), ('exposed_element', 'exposed')):
                        element_dataframe = curr_elements_dataframe[curr_elements_dataframe['element'] == organ_attribute_type]
                        if len(element_dataframe) == 0:
                            continue
                        element_dataframe = element_dataframe.loc[:, simulation.Simulation.ELEMENTS_STATE]
                        # create a new element
                        element_dict = element_dataframe.loc[element_dataframe.first_valid_index()].to_dict()
                        element = phytomer_attribute_element_class(**element_dict)
                        setattr(organ, organ_attribute_name, element)
                        
    return population


def to_dataframes(population):
    """
    Convert `population` to Pandas dataframes.
    
    :Parameters:
        
        - `population` (:class:`model.Population`) - The population to convert.
    
    :Returns:
        :class:`dataframes <pandas.DataFrame>` describing internal state and compartments at each scale:

            * plant scale: plant index, state parameters and compartments of each plant (see :attr:`Simulation:PLANTS_DATAFRAME_COLUMNS`) 
            * axis scale: plant index, axis id, state parameters and compartments of each axis (see :attr:`Simulation:AXES_DATAFRAME_COLUMNS`)
            * phytomer scale: plant index, axis id, phytomer index, state parameters and compartments of phytomer plant (see :attr:`Simulation:PHYTOMERS_DATAFRAME_COLUMNS`)
            * organ scale: plant index, axis id, phytomer index, organ type, state parameters and compartments of each organ (see :attr:`Simulation:ORGANS_DATAFRAME_COLUMNS`)
            * and element scale: plant index, axis id, phytomer index, organ type, element type, state parameters and compartments of each element (see :attr:`Simulation:ELEMENTS_DATAFRAME_COLUMNS`)
    
    :Returns Type:
        :class:`tuple` of :class:`pandas.DataFrame`
    
    """
    
    # initialize the dataframes
    all_plants_df = pd.DataFrame(columns=PLANTS_DATAFRAME_COLUMNS)
    all_axes_df = pd.DataFrame(columns=AXES_DATAFRAME_COLUMNS)
    all_phytomers_df = pd.DataFrame(columns=PHYTOMERS_DATAFRAME_COLUMNS)
    all_organs_df = pd.DataFrame(columns=ORGANS_DATAFRAME_COLUMNS)
    all_elements_df = pd.DataFrame(columns=ELEMENTS_DATAFRAME_COLUMNS)
    
    def append_row(model_object, indexes, attributes_names, inputs_df):
        # function to append a row in a dataframe
        attributes_values = []
        for attribute_name in attributes_names:
            attributes_values.append(getattr(model_object, attribute_name, np.nan))
        if len(attributes_values) == 0 \
            or any(attribute_value is None for attribute_value in attributes_values) \
            or not np.isnan([attribute_value for attribute_value in attributes_values if attribute_value is not None]).all():
            inputs_df.loc[len(inputs_df),:] = indexes + attributes_values
    
    # run through the population tree and fill the dataframes
    plant_index = 1
    for plant in population.plants:
        append_row(plant, [plant_index], simulation.Simulation.PLANTS_STATE, all_plants_df)
        axis_index = 0
        for axis in plant.axes:
            axis_id = model.Axis.get_axis_id(axis_index)
            append_row(axis, [plant_index, axis_id], simulation.Simulation.AXES_STATE, all_axes_df)
            for organ in (axis.roots, axis.soil, axis.phloem, axis.grains):
                if organ is None:
                    continue
                append_row(organ, [plant_index, axis_id, np.nan, organ.__class__.__name__], simulation.Simulation.ORGANS_STATE, all_organs_df)
            phytomer_index = 1
            for phytomer in axis.phytomers:
                append_row(phytomer, [plant_index, axis_id, phytomer_index], simulation.Simulation.PHYTOMERS_STATE, all_phytomers_df)
                for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                    if organ is None:
                        continue
                    append_row(organ, [plant_index, axis_id, phytomer_index, organ.__class__.__name__], simulation.Simulation.ORGANS_STATE, all_organs_df)
                    for element, element_type in ((organ.exposed_element, 'exposed'), (organ.enclosed_element, 'enclosed')):
                        if element is None:
                            continue
                        append_row(element, [plant_index, axis_id, phytomer_index, organ.__class__.__name__, element_type], simulation.Simulation.ELEMENTS_STATE, all_elements_df)
                phytomer_index += 1
            axis_index += 1
        plant_index += 1
            
    # sort the rows of the dataframes by columns
    all_plants_df.sort_index(by=PLANTS_DATAFRAME_COLUMNS, inplace=True)
    all_axes_df.sort_index(by=AXES_DATAFRAME_COLUMNS, inplace=True)
    all_phytomers_df.sort_index(by=PHYTOMERS_DATAFRAME_COLUMNS, inplace=True)
    all_organs_df.sort_index(by=ORGANS_DATAFRAME_COLUMNS, inplace=True)
    all_elements_df.sort_index(by=ELEMENTS_DATAFRAME_COLUMNS, inplace=True)

    # infer the right types of the columns in the dataframes
    all_plants_df = all_plants_df.convert_objects(copy=False)
    all_axes_df = all_axes_df.convert_objects(copy=False)
    all_phytomers_df = all_phytomers_df.convert_objects(copy=False)
    all_organs_df = all_organs_df.convert_objects(copy=False)
    all_elements_df = all_elements_df.convert_objects(copy=False)

    # convert the indexes of plants, phytomers and elements to integers in the dataframes
    all_plants_df['plant'] = all_plants_df['plant'].astype(int)
    all_axes_df['plant'] = all_axes_df['plant'].astype(int)
    all_phytomers_df[['plant', 'phytomer']] = all_phytomers_df[['plant', 'phytomer']].astype(int)
    all_organs_df['plant'] = all_organs_df['plant'].astype(int)
    all_elements_df[['plant', 'phytomer']] = all_elements_df[['plant', 'phytomer']].astype(int)
    
    return all_plants_df, all_axes_df, all_phytomers_df, all_organs_df, all_elements_df

    