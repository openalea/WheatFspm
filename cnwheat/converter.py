# -*- coding: latin-1 -*-

"""
    cnwheat.converter
    ~~~~~~~~~~~~~~~~~

    The module :mod:`cnwheat.converter` defines functions to:
    
        * convert :class:`dataframes <pandas.DataFrame>` to/from a CNWheat :class:`population <model.Population>`.
        * convert a :class:`MTG <openalea.mtg.mtg.MTG>` to/from a CNWheat :class:`population <model.Population>`.
    
    Both dataframes and MTG follow AdelWheat naming convention.

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

import numpy as np
import pandas as pd
from openalea.mtg import io, fat_mtg

from cnwheat import model, simulation


#: the topology and state variables at PLANT scale
PLANTS_STATE_VARIABLES = simulation.Simulation.PLANTS_INPUTS_INDEXES + simulation.Simulation.PLANTS_STATE

#: the topology and state variables at AXIS scale
AXES_STATE_VARIABLES = simulation.Simulation.AXES_INPUTS_INDEXES + simulation.Simulation.AXES_STATE

#: the topology and state variables at PHYTOMER scale
PHYTOMERS_STATE_VARIABLES = simulation.Simulation.PHYTOMERS_INPUTS_INDEXES + simulation.Simulation.PHYTOMERS_STATE

#: the topology and state variables at ORGAN scale
ORGANS_STATE_VARIABLES = simulation.Simulation.ORGANS_INPUTS_INDEXES + simulation.Simulation.ORGANS_STATE

#: the topology and state variables at ELEMENT scale
ELEMENTS_STATE_VARIABLES = simulation.Simulation.ELEMENTS_INPUTS_INDEXES + simulation.Simulation.ELEMENTS_STATE

#: the topology and state variables of soils
SOILS_STATE_VARIABLES = simulation.Simulation.SOILS_INPUTS_INDEXES + simulation.Simulation.SOILS_STATE

#: the columns, in the input/output dataframes, which define the topology of the organs which belong to an axis
AXES_ORGANS_TOPOLOGY_COLUMNS = ['plant', 'axis', 'organ']

#: the columns, in the input/output dataframes, which define the topology of the soils
SOILS_TOPOLOGY_COLUMNS = ['plant', 'axis']

#: the mapping of CNWheat organ classes to the attributes in axis and phytomer which represent an organ
CNWHEAT_ATTRIBUTES_MAPPING = {model.Internode: 'internode', model.Lamina: 'lamina', 
                              model.Sheath: 'sheath', model.Peduncle: 'peduncle', model.Chaff: 'chaff',
                              model.Roots: 'roots', model.Grains: 'grains', model.Phloem: 'phloem'}

#: the mapping of the CNWheat organ classes to organ names in MTG
CNWHEAT_CLASSES_TO_MTG_ORGANS_MAPPING = {model.Internode: 'internode', model.Lamina: 'blade', 
                                         model.Sheath: 'sheath', model.Peduncle: 'peduncle', model.Chaff: 'ear',
                                         model.Roots: 'roots', model.Grains: 'grains', model.Phloem: 'phloem'}

#: the mapping of the name of each element, from MTG to CNWheat
MTG_TO_CNWHEAT_ELEMENTS_NAMES_MAPPING = {'HiddenElement': 'enclosed_element', 'StemElement': 'exposed_element', 'LeafElement1': 'exposed_element'}

#: the mapping of organs (which belong to an axis) labels in MTG to organ classes in CNWheat
MTG_TO_CNWHEAT_AXES_ORGANS_MAPPING = {'grains': model.Grains, 'phloem': model.Phloem, 'roots': model.Roots}

#: the mapping of organs (which belong to a phytomer) labels in MTG to organ classes in CNWheat
MTG_TO_CNWHEAT_PHYTOMERS_ORGANS_MAPPING = {'internode': model.Internode, 'blade': model.Lamina, 'sheath': model.Sheath, 'peduncle': model.Peduncle, 'ear': model.Chaff}

# the mapping of CNWheat photosynthetic organs to CNWheat photosynthetic organ elements 
CNWHEAT_ORGANS_TO_ELEMENTS_MAPPING = {model.Internode: model.InternodeElement, model.Lamina: model.LaminaElement, model.Sheath: model.SheathElement, model.Peduncle: model.PeduncleElement, model.Chaff: model.ChaffElement}

#: the parameters and variables which define the state of a CNWheat population
POPULATION_STATE_VARIABLE = set(simulation.Simulation.PLANTS_STATE + simulation.Simulation.AXES_STATE + 
                                simulation.Simulation.PHYTOMERS_STATE + simulation.Simulation.ORGANS_STATE +
                                simulation.Simulation.ELEMENTS_STATE)


def from_dataframes(plants_inputs=None, axes_inputs=None, metamers_inputs=None, organs_inputs=None, elements_inputs=None, soils_inputs=None):
    """
    If `plants_inputs`, `axes_inputs`, `metamers_inputs`, `organs_inputs` and `elements_inputs` are not `None`, convert `plants_inputs`, 
    `axes_inputs`, `metamers_inputs`, `organs_inputs` and  `elements_inputs` to a :class:`population <model.Population>`.
    
    If `soils_inputs` is not `None`, convert `soils_inputs` to a dictionary of :class:`soils <model.Soil>`.
    
    Data in `plants_inputs`, `axes_inputs`, `metamers_inputs`, `organs_inputs`, `elements_inputs` and `soils_inputs`
    respect the naming convention of AdelWheat.
    
    :Parameters:
        
        - `plants_inputs` (:class:`pandas.DataFrame`) - Plants inputs, with one line by plant.
        
        - `axes_inputs` (:class:`pandas.DataFrame`) - Axes inputs, with one line by axis.
        
        - `metamers_inputs` (:class:`pandas.DataFrame`) - Metamers inputs, with one line by metamer.
        
        - `organs_inputs` (:class:`pandas.DataFrame`) - Organs inputs, with one line by organ.
        
        - `elements_inputs` (:class:`pandas.DataFrame`) - Elements inputs, with one line by element.
        
        - `soils_inputs` (:class:`pandas.DataFrame`) - Soils inputs, with one line by soil.
        
    :Returns:
        If `plants_inputs`, `axes_inputs`, `metamers_inputs`, `organs_inputs` and `elements_inputs` are not `None`, return a :class:`population <model.Population>`,
        and/or 
        if `soils_inputs` is not `None`,  return a :class:`dict` of :class:`soils <model.Soil>`.
        
    :Returns Type:
        :class:`tuple` or :class:`population <model.Population>` or :class:`dict`.
        
    """
    
    convert_dataframes_to_population = plants_inputs is not None and axes_inputs is not None and metamers_inputs is not None and organs_inputs is not None and elements_inputs is not None
    convert_dataframe_to_soils_dict = soils_inputs is not None
    
    if convert_dataframes_to_population:
        population = model.Population()
        
        for plant_index in plants_inputs.plant:
            # create a new plant
            plant = model.Plant(plant_index)
            population.plants.append(plant)
            curr_axes_inputs = axes_inputs[axes_inputs['plant'] == plant_index]
            for axis_label in curr_axes_inputs.axis:
                # create a new axis
                axis = model.Axis(axis_label)
                curr_organs_inputs = organs_inputs[(organs_inputs['plant'] == plant_index) & (organs_inputs['axis'] == axis_label)]
                for axis_attribute_name, axis_attribute_class in (('grains', model.Grains), ('roots', model.Roots), ('phloem', model.Phloem)):
                    organ_label = CNWHEAT_CLASSES_TO_MTG_ORGANS_MAPPING[axis_attribute_class]
                    organ_inputs = curr_organs_inputs[curr_organs_inputs['organ'] == organ_label]
                    # create a new organ
                    organ = axis_attribute_class(organ_label)
                    organ_attributes_names = [state_var_name for state_var_name in simulation.Simulation.ORGANS_STATE if hasattr(organ, state_var_name)]
                    organ_row = organ_inputs.loc[organ_inputs.first_valid_index()]
                    organ_attributes_values = organ_row[organ_attributes_names].tolist()
                    organ_attributes = dict(zip(organ_attributes_names, organ_attributes_values))
                    organ.__dict__.update(organ_attributes)
                    organ.initialize()
                    setattr(axis, axis_attribute_name, organ)
                        
                curr_metamers_inputs = metamers_inputs[(metamers_inputs['plant'] == plant_index) & (metamers_inputs['axis'] == axis_label)]
                for metamer_index in curr_metamers_inputs.metamer:
                    # create a new phytomer
                    phytomer = model.Phytomer(metamer_index)
                    axis.phytomers.append(phytomer)
                    
                    curr_organs_inputs = organs_inputs[(organs_inputs['plant'] == plant_index) & (organs_inputs['axis'] == axis_label) & (organs_inputs['metamer'] == metamer_index)]
                    
                    for phytomer_attribute_name, phytomer_attribute_class, phytomer_attribute_element_class in \
                        (('chaff', model.Chaff, model.ChaffElement), 
                         ('lamina', model.Lamina, model.LaminaElement), 
                         ('internode', model.Internode, model.InternodeElement), 
                         ('peduncle', model.Peduncle, model.PeduncleElement),
                         ('sheath', model.Sheath, model.SheathElement)):
                        
                        organ_label = CNWHEAT_CLASSES_TO_MTG_ORGANS_MAPPING[phytomer_attribute_class]
                        curr_elements_inputs = elements_inputs[(elements_inputs['plant'] == plant_index) & (elements_inputs['axis'] == axis_label) & (elements_inputs['metamer'] == metamer_index) & (elements_inputs['organ'] == organ_label)]
                        if organ_label not in curr_organs_inputs.organ.values and organ_label not in curr_elements_inputs.organ.values:
                            continue
                        # create a new organ
                        organ = phytomer_attribute_class(organ_label)
                        organ.initialize()
                        setattr(phytomer, phytomer_attribute_name, organ)
                        
                        for mtg_element_label, cnwheat_element_name in MTG_TO_CNWHEAT_ELEMENTS_NAMES_MAPPING.iteritems():
                            element_inputs = curr_elements_inputs[curr_elements_inputs['element'] == mtg_element_label]
                            if len(element_inputs) == 0:
                                continue
                            element_inputs = element_inputs.loc[:, simulation.Simulation.ELEMENTS_STATE]
                            element_dict = element_inputs.loc[element_inputs.first_valid_index()].to_dict()
                            if element_dict['green_area'] == 0.0:
                                continue
                            # create a new element
                            element = phytomer_attribute_element_class(mtg_element_label, **element_dict)
                            setattr(organ, cnwheat_element_name, element)
                            
                plant.axes.append(axis)
                
    if convert_dataframe_to_soils_dict:
        soils = {}
        for soil_id, soil_group in soils_inputs.groupby(SOILS_TOPOLOGY_COLUMNS):
            # create a new soil
            soil_attributes = soil_group.loc[soil_group.first_valid_index(), simulation.Simulation.SOILS_STATE].to_dict()
            soil = model.Soil(**soil_attributes)
            soils[soil_id] = soil

    if convert_dataframes_to_population and convert_dataframe_to_soils_dict:
        return population, soils
    elif convert_dataframes_to_population:
        return population
    else:
        return soils
                        

def to_dataframes(population=None, soils=None):
    """
    If `population` is not None, convert `population` to Pandas dataframes.
    If `soils` is not None, convert `soils` to Pandas dataframe.
    Convert a CNWheat :class:`population <model.Population>` and a dictionary of :class:`soils <model.Soil>` to Pandas dataframes.
    
    Data in the resulting dataframe(s) respect(s) the naming convention of AdelWheat.
    
    :Parameters:
        
        - `population` (:class:`model.Population`) - The CNWheat population to convert.
        
        - `soils` (:class:`dict` of `model.Soil`) - The soils to convert.
    
    :Returns:
        If `population` is not None, return :class:`dataframes <pandas.DataFrame>` describing the internal state and compartments of the population at each scale:

            * plant scale: plant index, state parameters and compartments of each plant (see :attr:`Simulation:PLANTS_STATE_VARIABLES`) 
            * axis scale: plant index, axis id, state parameters and compartments of each axis (see :attr:`Simulation:AXES_STATE_VARIABLES`)
            * metamer scale: plant index, axis id, metamer index, state parameters and compartments of each metamer (see :attr:`Simulation:PHYTOMERS_STATE_VARIABLES`)
            * organ scale: plant index, axis id, metamer index, organ type, state parameters and compartments of each organ (see :attr:`Simulation:ORGANS_STATE_VARIABLES`)
            * and element scale: plant index, axis id, metamer index, organ type, element type, state parameters and compartments of each element (see :attr:`Simulation:ELEMENTS_STATE_VARIABLES`),
        
        and/or
        
        if `soils` is not None, return a :class:`dataframe <pandas.DataFrame>` describing internal state and compartments of the soils, with one line by soil.
        
    
    :Returns Type:
        :class:`tuple` of :class:`pandas.DataFrame` 
        or
        :class:`pandas.DataFrame` 
    
    """
    
    convert_population_to_dataframes = population is not None
    convert_soils_to_dataframe = soils is not None
    
    def append_row(model_object, indexes, attributes_names, inputs_df):
        # function to append a row to a dataframe
        attributes_values = []
        for attribute_name in attributes_names:
            attributes_values.append(getattr(model_object, attribute_name, np.nan))
        if len(attributes_values) == 0 \
            or any(attribute_value is None for attribute_value in attributes_values) \
            or not np.isnan([attribute_value for attribute_value in attributes_values if attribute_value is not None]).all():
            inputs_df.loc[len(inputs_df),:] = indexes + attributes_values
    
    if convert_population_to_dataframes:
        # initialize the dataframes
        all_plants_df = pd.DataFrame(columns=PLANTS_STATE_VARIABLES)
        all_axes_df = pd.DataFrame(columns=AXES_STATE_VARIABLES)
        all_metamers_df = pd.DataFrame(columns=PHYTOMERS_STATE_VARIABLES)
        all_organs_df = pd.DataFrame(columns=ORGANS_STATE_VARIABLES)
        all_elements_df = pd.DataFrame(columns=ELEMENTS_STATE_VARIABLES)
        
        # run through the population tree and fill the dataframes
        for plant in population.plants:
            append_row(plant, [plant.index], simulation.Simulation.PLANTS_STATE, all_plants_df)
            for axis in plant.axes:
                append_row(axis, [plant.index, axis.label], simulation.Simulation.AXES_STATE, all_axes_df)
                for organ in (axis.roots, axis.phloem, axis.grains):
                    append_row(organ, [plant.index, axis.label, np.nan, organ.label], simulation.Simulation.ORGANS_STATE, all_organs_df)
                for phytomer in axis.phytomers:
                    append_row(phytomer, [plant.index, axis.label, phytomer.index], simulation.Simulation.PHYTOMERS_STATE, all_metamers_df)
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ is None:
                            continue
                        append_row(organ, [plant.index, axis.label, phytomer.index, organ.label], simulation.Simulation.ORGANS_STATE, all_organs_df)
                        for element in (organ.exposed_element, organ.enclosed_element):
                            if element is None:
                                continue
                            append_row(element, [plant.index, axis.label, phytomer.index, organ.label, element.label], simulation.Simulation.ELEMENTS_STATE, all_elements_df)
                
        # sort the rows of the dataframes by columns
        all_plants_df.sort_index(by=PLANTS_STATE_VARIABLES, inplace=True)
        all_axes_df.sort_index(by=AXES_STATE_VARIABLES, inplace=True)
        all_metamers_df.sort_index(by=PHYTOMERS_STATE_VARIABLES, inplace=True)
        all_organs_df.sort_index(by=ORGANS_STATE_VARIABLES, inplace=True)
        all_elements_df.sort_index(by=ELEMENTS_STATE_VARIABLES, inplace=True)
    
        # infer the right types of the columns in the dataframes
        all_plants_df = all_plants_df.convert_objects(copy=False)
        all_axes_df = all_axes_df.convert_objects(copy=False)
        all_metamers_df = all_metamers_df.convert_objects(copy=False)
        all_organs_df = all_organs_df.convert_objects(copy=False)
        all_elements_df = all_elements_df.convert_objects(copy=False)
    
        # convert the indexes of plants, metamers and elements to integers in the dataframes
        all_plants_df['plant'] = all_plants_df['plant'].astype(int)
        all_axes_df['plant'] = all_axes_df['plant'].astype(int)
        all_metamers_df[['plant', 'metamer']] = all_metamers_df[['plant', 'metamer']].astype(int)
        all_organs_df['plant'] = all_organs_df['plant'].astype(int)
        all_elements_df[['plant', 'metamer']] = all_elements_df[['plant', 'metamer']].astype(int)
        
        all_plants_df.reset_index(drop=True, inplace=True)
        all_axes_df.reset_index(drop=True, inplace=True)
        all_metamers_df.reset_index(drop=True, inplace=True)
        all_organs_df.reset_index(drop=True, inplace=True)
        all_elements_df.reset_index(drop=True, inplace=True)
        
    if convert_soils_to_dataframe:
        all_soils_df = pd.DataFrame(columns=SOILS_STATE_VARIABLES)
        for soil_id, soil in soils.iteritems():
            append_row(soil, list(soil_id), simulation.Simulation.SOILS_STATE, all_soils_df)
        all_soils_df.sort_index(by=SOILS_STATE_VARIABLES, inplace=True)
        all_soils_df = all_soils_df.convert_objects(copy=False)
        all_soils_df['plant'] = all_soils_df['plant'].astype(int)
        all_soils_df.reset_index(drop=True, inplace=True)
    
    if convert_population_to_dataframes and convert_soils_to_dataframe:
        return all_plants_df, all_axes_df, all_metamers_df, all_organs_df, all_elements_df, all_soils_df
    elif convert_population_to_dataframes:
        return all_plants_df, all_axes_df, all_metamers_df, all_organs_df, all_elements_df
    else:
        return all_soils_df
    

def from_MTG(g, organs_inputs, elements_inputs):
    """
    Convert a MTG to a CNWheat :class:`population <model.Population>`. 
    Use data in `organs_inputs` and `elements_inputs` if `g` is incomplete.
    The property names in the MTG respect the naming convention of AdelWheat.
    
    :Parameters:
        
            - g (:class:`openalea.mtg.mtg.MTG`) - A MTG which contains the inputs
              of CN-Wheat. These inputs are: :mod:`CNWHEAT_INPUTS`.
              
            - `organs_inputs` (:class:`pandas.DataFrame`) - Organs dataframe, with one line by organ.
            
            - `elements_inputs` (:class:`pandas.DataFrame`) - Elements dataframe, with one line by element.
            
    :Returns:
        A CNWheat :class:`population <model.Population>`.
        
    :Returns Type:
        :class:`population <model.Population>`
        
    """
    
    population = model.Population()
    
    organs_inputs_grouped = organs_inputs.groupby(AXES_ORGANS_TOPOLOGY_COLUMNS)
    elements_inputs_grouped = elements_inputs.groupby(simulation.Simulation.ELEMENTS_INPUTS_INDEXES)
    
    # traverse the MTG recursively from top
    for plant_vid in g.components_iter(g.root):
        plant_index = int(g.index(plant_vid))
        # create a new plant
        plant = model.Plant(plant_index)
        is_valid_plant = False
        for axis_vid in g.components_iter(plant_vid):
            axis_label = g.label(axis_vid)
            # create a new axis
            axis = model.Axis(axis_label)
            axis_components_iter = g.components_iter(axis_vid)
            
            first_axis_component_vid = next(axis_components_iter)
            second_axis_component_vid = next(axis_components_iter)
            
            first_and_second_axis_components = {g.label(first_axis_component_vid): first_axis_component_vid, 
                                                g.label(second_axis_component_vid): second_axis_component_vid}
            
            is_valid_axis = True
            # create and initialize a roots and a phloem
            for organ_class in (model.Roots, model.Phloem):
                organ_label = CNWHEAT_CLASSES_TO_MTG_ORGANS_MAPPING[organ_class]
                organ_id = (plant_index, axis_label, organ_label)
                organ = MTG_TO_CNWHEAT_AXES_ORGANS_MAPPING[organ_label](organ_label)
                organ_inputs_names = [state_var_name for state_var_name in simulation.Simulation.ORGANS_STATE if hasattr(organ, state_var_name)]
                if organ_id in organs_inputs_grouped.groups:
                    organ_inputs_group = organs_inputs_grouped.get_group(organ_id)
                    organ_inputs_group_series = organ_inputs_group.loc[organ_inputs_group.first_valid_index(), organ_inputs_group.columns.intersection(organ_inputs_names)]
                else:
                    organ_inputs_group_series = pd.Series()
                is_valid_organ = True
                if organ_label in first_and_second_axis_components:
                    axis_component_vid = first_and_second_axis_components[organ_label]
                    organ_vid = axis_component_vid
                    vertex_properties = g.get_vertex_property(organ_vid)
                    organ_inputs_dict = {}
                    for organ_input_name in organ_inputs_names:
                        if organ_input_name in vertex_properties:
                            # use the property of the vertex
                            organ_inputs_dict[organ_input_name] = vertex_properties[organ_input_name]
                        else:
                            # use the value in organs_inputs
                            if organ_input_name in organ_inputs_group_series:
                                organ_inputs_dict[organ_input_name] = organ_inputs_group_series[organ_input_name]
                            else:
                                is_valid_organ = False
                                break
                else:
                    organ_inputs_dict = organ_inputs_group_series.to_dict()
                    if not set(organ_inputs_dict).issuperset(organ_inputs_names):
                        is_valid_organ = False
                if is_valid_organ:
                    organ.__dict__.update(organ_inputs_dict)
                    organ.initialize()
                    setattr(axis, CNWHEAT_ATTRIBUTES_MAPPING[organ_class], organ)
                else:
                    is_valid_axis = False
                    break
            
            if not is_valid_axis:
                continue
            
            mtg_has_grains = False
            grains_label = CNWHEAT_CLASSES_TO_MTG_ORGANS_MAPPING[model.Grains]
            grains_id = (plant_index, axis_label, grains_label)
            # create a new grains
            grains = MTG_TO_CNWHEAT_AXES_ORGANS_MAPPING[grains_label](grains_label)
            grains_inputs_names = [state_var_name for state_var_name in simulation.Simulation.ORGANS_STATE if hasattr(grains, state_var_name)]
            if grains_id in organs_inputs_grouped.groups:
                grains_inputs_group = organs_inputs_grouped.get_group(grains_id)
                grains_inputs_group_series = grains_inputs_group.loc[grains_inputs_group.first_valid_index(), grains_inputs_group.columns.intersection(grains_inputs_names)]
            else:
                grains_inputs_group_series = pd.Series()
            is_valid_grains = True
            
            has_valid_phytomer = False
            for axis_component_vid in g.components_iter(axis_vid):
                mtg_axis_component_label = g.label(axis_component_vid)
                if mtg_axis_component_label == grains_label:
                    mtg_has_grains = True
                    grains_vid = axis_component_vid
                    vertex_properties = g.get_vertex_property(grains_vid)
                    grains_inputs_dict = {}
                    for grains_input_name in grains_inputs_names:
                        if grains_input_name in vertex_properties:
                            # use the property of the vertex
                            grains_inputs_dict[grains_input_name] = vertex_properties[grains_input_name]
                        else:
                            # use the value in organs_inputs
                            if grains_input_name in grains_inputs_group_series:
                                grains_inputs_dict[grains_input_name] = grains_inputs_group_series[grains_input_name]
                            else:
                                is_valid_grains = False
                                break
                elif mtg_axis_component_label.startswith('metamer'):
                    metamer_vid = axis_component_vid
                    metamer_index = int(g.index(metamer_vid))
                    # create a new phytomer
                    phytomer = model.Phytomer(metamer_index)
                    has_valid_organ = False
                    for organ_vid in g.components_iter(metamer_vid):
                        organ_label = g.label(organ_vid)
                        if organ_label not in MTG_TO_CNWHEAT_PHYTOMERS_ORGANS_MAPPING: continue
                        # create a new organ
                        organ_class = MTG_TO_CNWHEAT_PHYTOMERS_ORGANS_MAPPING[organ_label]
                        organ = organ_class(organ_label)
                        organ.initialize()
                        has_valid_element = False
                        for element_vid in g.components_iter(organ_vid):
                            vertex_properties = g.get_vertex_property(element_vid)
                            element_label = g.label(element_vid)
                            element_id = (plant_index, axis_label, metamer_index, organ_label, element_label)
                            if element_id in elements_inputs_grouped.groups:
                                elements_inputs_group = elements_inputs_grouped.get_group(element_id)
                                elements_inputs_group_series = elements_inputs_group.loc[elements_inputs_group.first_valid_index(), elements_inputs_group.columns.intersection(simulation.Simulation.ELEMENTS_STATE)] 
                            else:
                                elements_inputs_group_series = pd.Series()
                            element_inputs = {}
                            is_valid_element = True
                            for element_input_name in simulation.Simulation.ELEMENTS_STATE:
                                if element_input_name in vertex_properties:
                                    # use the properties of the vertex
                                    element_inputs[element_input_name] = vertex_properties[element_input_name]
                                    if element_input_name == 'green_area':
                                        element_inputs[element_input_name] /= 10000.0 # convert from cm2 to m2 ; TODO: homogenize the units between the models 
                                else:
                                    # use the value in elements_inputs
                                    if element_input_name in elements_inputs_group_series:
                                        element_inputs[element_input_name] = elements_inputs_group_series[element_input_name]
                                    else:
                                        is_valid_element = False
                                        break
                                if element_input_name == 'green_area' and element_inputs[element_input_name] == 0.0:
                                    is_valid_element = False
                                    break
                            if is_valid_element:
                                has_valid_element = True
                                element = CNWHEAT_ORGANS_TO_ELEMENTS_MAPPING[organ_class](element_label, **element_inputs)
                                setattr(organ, MTG_TO_CNWHEAT_ELEMENTS_NAMES_MAPPING[element_label], element)
                        if has_valid_element:
                            has_valid_organ = True
                            setattr(phytomer, CNWHEAT_ATTRIBUTES_MAPPING[organ_class], organ)
                            
                    if has_valid_organ:  
                        axis.phytomers.append(phytomer)
                        has_valid_phytomer = True
                
            if not mtg_has_grains:
                grains_inputs_dict = grains_inputs_group_series.to_dict()
                if not set(grains_inputs_dict).issuperset(grains_inputs_names):
                    is_valid_grains = False
                
            if is_valid_grains:
                grains.__dict__.update(grains_inputs_dict)
                grains.initialize()
                setattr(axis, CNWHEAT_ATTRIBUTES_MAPPING[model.Grains], grains)
            else:
                is_valid_axis = False
                
            if not has_valid_phytomer:
                is_valid_axis = False
            
            if is_valid_axis:
                plant.axes.append(axis)
                is_valid_plant = True
                    
        if is_valid_plant:
            population.plants.append(plant)
        
    return population
                        
                    
def update_MTG(population, g):
    """
    Update a MTG from a CN-Wheat :class:`population <model.Population>`.
    The property names in the MTG respect the naming convention of AdelWheat. 
    
    :Parameters:
        
            - population (:class:`model.Population`) - a CN-Wheat :class:`population <model.Population>`.
        
            - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG to update from the CN-Wheat :class:`population <model.Population>`. 
            
    """
    
    # add the properties if needed
    property_names = g.property_names()
    for cnwheat_variable_name in POPULATION_STATE_VARIABLE:
        if cnwheat_variable_name not in property_names:
            g.add_property(cnwheat_variable_name)
    
    
    plants_iterator = g.components_iter(g.root)
    
    # traverse the population recursively from top
    for plant in population.plants:
        plant_index = plant.index
        while True:
            plant_vid = next(plants_iterator)
            if int(g.index(plant_vid)) == plant_index:
                break
        axes_iterator = g.components_iter(plant_vid)
        for axis in plant.axes:
            axis_label = axis.label
            while True:
                axis_vid = next(axes_iterator)
                if g.label(axis_vid) == axis_label:
                    break
            
            axis_components_iter = g.components_iter(axis_vid)
            
            first_axis_component_vid = next(axis_components_iter)
            second_axis_component_vid = next(axis_components_iter)
            
            first_and_second_axis_components = {g.label(first_axis_component_vid): first_axis_component_vid, 
                                                g.label(second_axis_component_vid): second_axis_component_vid}
            
            # insert a roots and a phloem if needed, and update them
            for organ_class in (model.Roots, model.Phloem):
                organ = getattr(axis, CNWHEAT_ATTRIBUTES_MAPPING[organ_class])
                organ_label = organ.label
                if organ_label in first_and_second_axis_components:
                    organ_vid = first_and_second_axis_components[organ_label]
                else:
                    first_axis_component_vid = organ_vid = insert_parent_at_all_scales(g, first_axis_component_vid, label=organ_label)[0]
                    g = fat_mtg(g)
                organ_property_names = [property_name for property_name in simulation.Simulation.ORGANS_STATE if hasattr(organ, property_name)]
                for organ_property_name in organ_property_names:
                    g.property(organ_property_name)[organ_vid] = getattr(organ, organ_property_name)
                
            # update the metamers
            for phytomer in axis.phytomers:
                phytomer_index = phytomer.index
                while True:
                    axis_component_vid = next(axis_components_iter)
                    if g.label(axis_component_vid).startswith('metamer') and int(g.index(axis_component_vid)) == phytomer_index:
                        break
                metamer_vid = axis_component_vid
                
                # update the organs
                for organ_vid in g.components_iter(metamer_vid):
                    organ_label = g.label(organ_vid)
                    if organ_label not in MTG_TO_CNWHEAT_PHYTOMERS_ORGANS_MAPPING: continue
                    organ = getattr(phytomer, CNWHEAT_ATTRIBUTES_MAPPING[MTG_TO_CNWHEAT_PHYTOMERS_ORGANS_MAPPING[organ_label]])
                    organ_property_names = [property_name for property_name in simulation.Simulation.ORGANS_STATE if hasattr(organ, property_name)]
                    for organ_property_name in organ_property_names:
                        g.property(organ_property_name)[organ_vid] = getattr(organ, organ_property_name)
                    
                    # update the elements
                    for element_vid in g.components_iter(organ_vid):
                        element_label = g.label(element_vid)
                        if element_label not in MTG_TO_CNWHEAT_ELEMENTS_NAMES_MAPPING: continue
                        element = getattr(organ, MTG_TO_CNWHEAT_ELEMENTS_NAMES_MAPPING[element_label])
                        element_property_names = [property_name for property_name in simulation.Simulation.ELEMENTS_STATE if hasattr(element, property_name)]
                        for element_property_name in element_property_names:
                            element_property_value = getattr(element, element_property_name)
                            if element_property_name == 'green_area':
                                element_property_value *= 10000.0 # convert from m2 to cm2 ; TODO: homogenize the units between the models 
                            g.property(element_property_name)[element_vid] = element_property_value
                    
            # insert a grains in the MTG if needed, and update it
            grains = axis.grains
            grains_label = grains.label
            last_axis_component_vid = g.components(axis_vid)[-1]
            if g.label(last_axis_component_vid) == grains_label:
                grains_vid = last_axis_component_vid
            else:
                grains_vid = add_child_at_all_scales(g, last_axis_component_vid, label=grains_label)[0]
            grains_property_names = [property_name for property_name in simulation.Simulation.ORGANS_STATE if hasattr(grains, property_name)]
            for grains_property_name in grains_property_names:
                g.property(grains_property_name)[grains_vid] = getattr(grains, grains_property_name)
                            


##########################################################################################################
####### TODO: move the following to openalea.mtg #########################################################
##########################################################################################################

def insert_parent_at_all_scales(g, parent_id, edge_type='+', label='roots'):
    added_vertices = []
    scale = g.scale(parent_id)
    max_scale = g.max_scale()
    edge_types = g.properties()['edge_type']
    
    # Add a parent
    if g.parent(parent_id) is None:
        vid = g.insert_parent(parent_id, label=label, edge_type='/')
        edge_types[parent_id] = edge_type
    else:
        return added_vertices
    
    # Update complex and components
    cid = g.complex(parent_id)
    croots = g._components[cid]
    if parent_id in croots:
        i = croots.index(parent_id)
        g._components[cid][i] = vid
        g._complex[vid] = cid
        print 'ADDED', vid
    else:
        print 'ERROR'

    added_vertices.append(vid)
    
    while scale+1 <= max_scale:
        pid = g.component_roots(parent_id)[0]
        component_id = g.add_component(vid)
        vid = g.insert_parent(pid, parent_id=component_id, label=label, edge_type='/')
        edge_types[pid] = edge_type
        scale += 1
        added_vertices.append(vid)
        parent_id = pid
        print 'ADDED', vid
    
    return added_vertices


def add_child_at_all_scales(g, parent_id, edge_type='<', label='grains'):
    added_vertices = []
    scale = g.scale(parent_id)
    max_scale = g.max_scale()
    edge_types = g.properties()['edge_type']
    
    # Add a parent
    parents = [parent_id]
    pid = parent_id
    while scale+1 <= max_scale:
        croots = g.component_roots(pid)
        if croots:
            pid = g.component_roots(pid)[0]
            parents.append(pid)
        else:
            break

    vid = g.add_child(parent_id, edge_type=edge_type, label=label)
    added_vertices.append(vid)
    print 'ADDED', vid
    
    for pid in parents[1:]:
        vid, cid = g.add_child_and_complex(pid, complex=vid, edge_type=edge_type, label=label)
        added_vertices.append(vid)
        print 'ADDED', vid
        
    return added_vertices
