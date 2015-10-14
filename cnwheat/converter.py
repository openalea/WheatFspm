# -*- coding: latin-1 -*-

"""
    cnwheat.converter
    ~~~~~~~~~~~~~~~~~

    The module :mod:`cnwheat.converter` defines functions to:
    
        * convert :class:`dataframes <pandas.DataFrame>` to/from a CNWheat :class:`population <model.Population>`.
        * convert a :class:`MTG <openalea.mtg.mtg.MTG>` to/from a CNWheat :class:`population <model.Population>`.

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

from cnwheat import model, simulation


class ConverterError(Exception): pass


class PropertyNotFoundError(ConverterError):
    '''Property not found in a vertex of a MTG.
    TODO: make this warning an error'''
    def __init__(self, property_, vertex_id):
        self.message = 'Property "{0}" not found in vertex {1}.'.format(property_, vertex_id)
    def __str__(self):
        return repr(self.message)
    
class NotModeledComponentError(ConverterError):
    '''MTG component not modeled by CNWheat.
    TODO: make this warning an error'''
    def __init__(self, component_label, component_vertex_id):
        self.message = 'MTG component "{0}" is not modeled by CNWheat (vertex {1}).'.format(component_label, component_vertex_id)
    def __str__(self):
        return repr(self.message)

class MismatchedTopologiesError(ConverterError):
    '''Topologies mismatched between CNWheat and MTG.
    TODO: make this warning an error'''
    def __init__(self, cnwheat_id, vertex_id):
        self.message = 'Topologies mismatched between CNWheat and MTG: no mapping between CNWheat object {0} and vertex {1}.'.format(cnwheat_id, vertex_id)
    def __str__(self):
        return repr(self.message)
    

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

#: the inputs needed by CNWheat
CNWHEAT_INPUTS = set(simulation.Simulation.PLANTS_STATE + simulation.Simulation.AXES_STATE + simulation.Simulation.PHYTOMERS_STATE + simulation.Simulation.ORGANS_STATE + simulation.Simulation.ELEMENTS_STATE)

#: the outputs computed by CNWheat
CNWHEAT_OUTPUTS = set([compartment_name for compartments_names in simulation.Simulation.MODEL_COMPARTMENTS_NAMES.values() for compartment_name in compartments_names])

#: the class corresponding to each organ modeled by CNWheat and which belongs to an axis
CNWHEAT_AXES_ORGANS_MAPPING = {'grains': model.Grains, 'phloem': model.Phloem, 'roots': model.Roots, 
                               'soil': model.Soil}

#: the class corresponding to each organ modeled by CNWheat and which belongs to a phytomer
CNWHEAT_PHYTOMERS_ORGANS_MAPPING = {'internode': model.Internode, 'lamina': model.Lamina, 
                                    'sheath': model.Sheath, 'peduncle': model.Peduncle, 'chaff': model.Chaff}

#: the mapping of the CNWheat organ classes to organ names in MTG, for organs which belong to a phytomer
CNWHEAT_CLASSES_TO_MTG_ORGANS_MAPPING = {model.Internode: 'internode', model.Lamina: 'blade', 
                                         model.Sheath: 'sheath', model.Peduncle: 'peduncle', model.Chaff: 'ear',
                                         model.Roots: 'roots', model.Grains: 'grains', model.Phloem: 'phloem', 
                                         model.Soil: 'soil'}

#: the mapping of the names of the attributes which designate the elements in the organ to the type of the elements
CNWHEAT_ELEMENTS_ATTRIBUTES_NAMES_TO_ELEMENTS_TYPES_MAPPING = {'enclosed_element': 'enclosed', 'exposed_element': 'exposed'}

#: the mapping of the name of each organ, from MTG to CNWheat, for organs which belong to an axis 
MTG_TO_CNWHEAT_AXES_ORGANS_NAMES_MAPPING = {'grains': 'grains', 'phloem': 'phloem', 'roots': 'roots', 
                                            'soil': 'soil'}

#: the mapping of the name of each organ, from MTG to CNWheat, for organs which belong to a phytomer 
MTG_TO_CNWHEAT_PHYTOMERS_ORGANS_NAMES_MAPPING = {'internode': 'internode', 'blade': 'lamina', 'sheath': 'sheath', 
                                                 'peduncle': 'peduncle', 'ear': 'chaff'}

#: the mapping of the name of each element, from MTG to CNWheat
MTG_TO_CNWHEAT_ELEMENTS_NAMES_MAPPING = {'HiddenElement': 'enclosed_element', 'StemElement': 'exposed_element', 'LeafElement1': 'exposed_element'}

#: the mapping of the name of each element, from CNWheat to MTG
CNWHEAT_TO_MTG_ELEMENTS_NAMES_MAPPING = {'exposed_element': {model.Lamina: 'LeafElement1', model.Internode: 'StemElement', model.Sheath: 'StemElement', 
                                                             model.Peduncle: 'StemElement', model.Chaff: 'StemElement'}, 
                                         'enclosed_element': {model.Lamina: 'HiddenElement', model.Internode: 'HiddenElement', model.Sheath: 'HiddenElement', 
                                                              model.Peduncle: 'HiddenElement', model.Chaff: 'HiddenElement'}}

#: the mapping of organ names to element classes in CNWheat 
CNWHEAT_ORGANS_TO_ELEMENTS_MAPPING = {'internode': model.InternodeElement, 'lamina': model.LaminaElement, 'sheath': model.SheathElement, 
                                      'peduncle': model.PeduncleElement, 'chaff': model.ChaffElement}


def from_dataframes(plants_dataframe, axes_dataframe, phytomers_dataframe, organs_dataframe, elements_dataframe):
    """
    Convert `plants_dataframe`, `axes_dataframe`, `phytomers_dataframe`, `organs_dataframe` and `elements_dataframe` 
    to a CNWheat :class:`population <model.Population>`.
    
    :Parameters:
        
        - `plants_dataframe` (:class:`pandas.DataFrame`) - Plants dataframe, with one line by plant.
        
        - `axes_dataframe` (:class:`pandas.DataFrame`) - Axes dataframe, with one line by axis.
        
        - `phytomers_dataframe` (:class:`pandas.DataFrame`) - Phytomers dataframe, with one line by phytomer.
        
        - `organs_dataframe` (:class:`pandas.DataFrame`) - Organs dataframe, with one line by organ.
        
        - `elements_dataframe` (:class:`pandas.DataFrame`) - Elements dataframe, with one line by element.
        
    :Returns:
        A CNWheat :class:`population <model.Population>`.
        
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
                    
                    organ_class_name = phytomer_attribute_class.__name__
                    curr_elements_dataframe = elements_dataframe[(elements_dataframe['plant'] == plant_index) & (elements_dataframe['axis'] == axis_id) & (elements_dataframe['phytomer'] == phytomer_index) & (elements_dataframe['organ'] == organ_class_name)]
                
                    if organ_class_name not in curr_organs_dataframe.organ.values and organ_class_name not in curr_elements_dataframe.organ.values:
                        continue
                    # create a new organ
                    organ = phytomer_attribute_class()
                    organ.initialize()
                    setattr(phytomer, phytomer_attribute_name, organ)
                    
                    for element_name, element_type in CNWHEAT_ELEMENTS_ATTRIBUTES_NAMES_TO_ELEMENTS_TYPES_MAPPING.iteritems():
                        element_dataframe = curr_elements_dataframe[curr_elements_dataframe['element'] == element_type]
                        if len(element_dataframe) == 0:
                            continue
                        element_dataframe = element_dataframe.loc[:, simulation.Simulation.ELEMENTS_STATE]
                        # create a new element
                        element_dict = element_dataframe.loc[element_dataframe.first_valid_index()].to_dict()
                        element = phytomer_attribute_element_class(**element_dict)
                        setattr(organ, element_name, element)
                        
    return population


def to_dataframes(population):
    """
    Convert a CNWheat :class:`population <model.Population>` to Pandas dataframes.
    
    :Parameters:
        
        - `population` (:class:`model.Population`) - The CNWheat population to convert.
    
    :Returns:
        :class:`dataframes <pandas.DataFrame>` describing internal state and compartments at each scale:

            * plant scale: plant index, state parameters and compartments of each plant (see :attr:`Simulation:PLANTS_STATE_VARIABLES`) 
            * axis scale: plant index, axis id, state parameters and compartments of each axis (see :attr:`Simulation:AXES_STATE_VARIABLES`)
            * phytomer scale: plant index, axis id, phytomer index, state parameters and compartments of phytomer plant (see :attr:`Simulation:PHYTOMERS_STATE_VARIABLES`)
            * organ scale: plant index, axis id, phytomer index, organ type, state parameters and compartments of each organ (see :attr:`Simulation:ORGANS_STATE_VARIABLES`)
            * and element scale: plant index, axis id, phytomer index, organ type, element type, state parameters and compartments of each element (see :attr:`Simulation:ELEMENTS_STATE_VARIABLES`)
    
    :Returns Type:
        :class:`tuple` of :class:`pandas.DataFrame`
    
    """
    
    # initialize the dataframes
    all_plants_df = pd.DataFrame(columns=PLANTS_STATE_VARIABLES)
    all_axes_df = pd.DataFrame(columns=AXES_STATE_VARIABLES)
    all_phytomers_df = pd.DataFrame(columns=PHYTOMERS_STATE_VARIABLES)
    all_organs_df = pd.DataFrame(columns=ORGANS_STATE_VARIABLES)
    all_elements_df = pd.DataFrame(columns=ELEMENTS_STATE_VARIABLES)
    
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
                append_row(organ, [plant_index, axis_id, np.nan, organ.__class__.__name__], simulation.Simulation.ORGANS_STATE, all_organs_df)
            phytomer_index = 9 # TODO: replace by 1
            for phytomer in axis.phytomers:
                append_row(phytomer, [plant_index, axis_id, phytomer_index], simulation.Simulation.PHYTOMERS_STATE, all_phytomers_df)
                for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                    if organ is None:
                        continue
                    append_row(organ, [plant_index, axis_id, phytomer_index, organ.__class__.__name__], simulation.Simulation.ORGANS_STATE, all_organs_df)
                    for element_name, element_type in CNWHEAT_ELEMENTS_ATTRIBUTES_NAMES_TO_ELEMENTS_TYPES_MAPPING.iteritems():
                        element = getattr(organ, element_name)
                        if element is None:
                            continue
                        append_row(element, [plant_index, axis_id, phytomer_index, organ.__class__.__name__, element_type], simulation.Simulation.ELEMENTS_STATE, all_elements_df)
                phytomer_index += 1
            axis_index += 1
        plant_index += 1
            
    # sort the rows of the dataframes by columns
    all_plants_df.sort_index(by=PLANTS_STATE_VARIABLES, inplace=True)
    all_axes_df.sort_index(by=AXES_STATE_VARIABLES, inplace=True)
    all_phytomers_df.sort_index(by=PHYTOMERS_STATE_VARIABLES, inplace=True)
    all_organs_df.sort_index(by=ORGANS_STATE_VARIABLES, inplace=True)
    all_elements_df.sort_index(by=ELEMENTS_STATE_VARIABLES, inplace=True)

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
    
    all_plants_df.reset_index(drop=True, inplace=True)
    all_axes_df.reset_index(drop=True, inplace=True)
    all_phytomers_df.reset_index(drop=True, inplace=True)
    all_organs_df.reset_index(drop=True, inplace=True)
    all_elements_df.reset_index(drop=True, inplace=True)
    
    return all_plants_df, all_axes_df, all_phytomers_df, all_organs_df, all_elements_df


def from_MTG(g, organs_dataframe, available_components):
    """
    Convert a MTG to a CNWheat :class:`population <model.Population>`. 
    
    :Parameters:
        
            - g (:class:`openalea.mtg.mtg.MTG`) - A MTG which contains the inputs
              of CN-Wheat. These inputs are: :mod:`CNWHEAT_INPUTS`.
              
            - `organs_dataframe` (:class:`pandas.DataFrame`) - Organs dataframe, with one line by organ.
            
            - `available_components` - TODO: remove this argument
              
    :Returns:
        A CNWheat :class:`population <model.Population>`.
        
    :Returns Type:
        :class:`model.Population`
        
    .. todo:: remove argument `organs_dataframe` as soon as :mod:`openalea.mtg` permits to add/remove components.
        
    """
    population = model.Population()
    
    # traverse the MTG recursively from top
    for plant_vertex_id in g.components_iter(g.root):
        if g.index(plant_vertex_id) not in available_components: continue
        # create a new plant
        plant = model.Plant()
        for axis_vertex_id in g.components_iter(plant_vertex_id):
            if (g.index(plant_vertex_id), g.label(axis_vertex_id)) not in available_components: continue
            # create a new axis
            axis = model.Axis()
            for axis_component_vertex_id in g.components_iter(axis_vertex_id):
                if (g.index(plant_vertex_id), g.label(axis_vertex_id), g.index(axis_component_vertex_id)) not in available_components: continue
                mtg_axis_component_label = g.label(axis_component_vertex_id)
                if mtg_axis_component_label in MTG_TO_CNWHEAT_AXES_ORGANS_NAMES_MAPPING: # grains, phloem, roots, soil
                    mtg_organ_label = mtg_axis_component_label
                    cnwheat_organ_name = MTG_TO_CNWHEAT_AXES_ORGANS_NAMES_MAPPING[mtg_organ_label]
                    organ_vertex_id = axis_component_vertex_id
                    # create a new organ
                    organ = CNWHEAT_AXES_ORGANS_MAPPING[cnwheat_organ_name]()
                    organ_attributes_names = [state_var_name for state_var_name in simulation.Simulation.ORGANS_STATE if hasattr(organ, state_var_name)]
                    attributes_found_in_vertex_properties = {}
                    vertex_properties = g.get_vertex_property(organ_vertex_id)
                    for organ_attribute_name in organ_attributes_names:
                        if organ_attribute_name not in vertex_properties:
                            raise PropertyNotFoundError(organ_attribute_name, organ_vertex_id)
                        attributes_found_in_vertex_properties[organ_attribute_name] = vertex_properties[organ_attribute_name]
                    # initialize the new organ and add it to the current axis
                    organ.__dict__.update(attributes_found_in_vertex_properties)
                    organ.initialize()
                    setattr(axis, cnwheat_organ_name, organ)
                    
                elif mtg_axis_component_label.startswith('metamer'): # phytomer
                    metamer_vertex_id = axis_component_vertex_id
                    # create a new phytomer
                    phytomer = model.Phytomer()
                    for organ_vertex_id in g.components_iter(metamer_vertex_id):
                        if (g.index(plant_vertex_id), g.label(axis_vertex_id), g.index(axis_component_vertex_id), g.label(organ_vertex_id)) not in available_components: continue
                        mtg_organ_label = g.label(organ_vertex_id)
                        if mtg_organ_label not in MTG_TO_CNWHEAT_PHYTOMERS_ORGANS_NAMES_MAPPING:
                            raise NotModeledComponentError(mtg_organ_label, organ_vertex_id)
                        cnwheat_organ_name = MTG_TO_CNWHEAT_PHYTOMERS_ORGANS_NAMES_MAPPING[mtg_organ_label]
                        # create a new organ
                        organ = CNWHEAT_PHYTOMERS_ORGANS_MAPPING[cnwheat_organ_name]()
                        organ_attributes_names = [state_var_name for state_var_name in simulation.Simulation.ORGANS_STATE if hasattr(organ, state_var_name)]
                        attributes_found_in_vertex_properties = {}
                        vertex_properties = g.get_vertex_property(organ_vertex_id)
                        for organ_attribute_name in organ_attributes_names:
                            if organ_attribute_name not in vertex_properties:
                                raise PropertyNotFoundError(organ_attribute_name, organ_vertex_id)
                            attributes_found_in_vertex_properties[organ_attribute_name] = vertex_properties[organ_attribute_name]
                        # initialize the new organ and add it to the current phytomer
                        organ.__dict__.update(attributes_found_in_vertex_properties)
                        organ.initialize()
                        for element_vertex_id in g.components_iter(organ_vertex_id):
                            if (g.index(plant_vertex_id), g.label(axis_vertex_id), g.index(axis_component_vertex_id), g.label(organ_vertex_id), g.label(element_vertex_id)) not in available_components: continue
                            mtg_element_label = g.label(element_vertex_id)
                            # create a new element
                            cnwheat_element_name = MTG_TO_CNWHEAT_ELEMENTS_NAMES_MAPPING[mtg_element_label]
                            attributes_found_in_vertex_properties = {}
                            vertex_properties = g.get_vertex_property(element_vertex_id)
                            for element_attribute_name in simulation.Simulation.ELEMENTS_STATE:
                                if element_attribute_name not in vertex_properties:
                                    raise PropertyNotFoundError(element_attribute_name, element_vertex_id)
                                attributes_found_in_vertex_properties[element_attribute_name] = vertex_properties[element_attribute_name]
                            element = CNWHEAT_ORGANS_TO_ELEMENTS_MAPPING[cnwheat_organ_name](**attributes_found_in_vertex_properties)
                            setattr(organ, cnwheat_element_name, element)
                        
                        setattr(phytomer, cnwheat_organ_name, organ)
                        
                    axis.phytomers.append(phytomer)
                else:
                    raise NotModeledComponentError(mtg_axis_component_label, axis_component_vertex_id)
            plant.axes.append(axis)
                
        population.plants.append(plant)
        
        #TODO: remove the following code as soon as :mod:`openalea.mtg` permits to add/remove components.
        plant_index = 1
        for plant in population.plants:
            axis_index = 0
            for axis in plant.axes:
                axis_id = model.Axis.get_axis_id(axis_index)
                for organ_name, organ_class in CNWHEAT_AXES_ORGANS_MAPPING.iteritems():
                    organ = organ_class()
                    organ_attributes_names = [state_var_name for state_var_name in simulation.Simulation.ORGANS_STATE if hasattr(organ, state_var_name)]
                    organ_dataframe = organs_dataframe[(organs_dataframe['plant'] == plant_index) & (organs_dataframe['axis'] == axis_id) & (organs_dataframe['organ'] == organ_class.__name__)]
                    organ_attributes_values = organ_dataframe.loc[organ_dataframe.first_valid_index(), organ_attributes_names].tolist()
                    organ_attributes = dict(zip(organ_attributes_names, organ_attributes_values))
                    organ.__dict__.update(organ_attributes)
                    organ.initialize()
                    setattr(axis, organ_name, organ)
                axis_index += 1
            plant_index += 1
                                
        return population         
                    
                    
def update_MTG(population, g, organs_dataframe, available_components):
    """
    Update a MTG from a CN-Wheat :class:`population <model.Population>`. 
    
    :Parameters:
        
            - population (:class:`model.Population`) - a CN-Wheat :class:`population <model.Population>`.
        
            - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG to update from the CN-Wheat :class:`population <model.Population>`. 
            
            - `organs_dataframe` (:class:`pandas.DataFrame`) - Organs dataframe, with one line by organ.
            
            - `available_components` - TODO: remove this argument
            
    .. todo:: remove argument `organs_dataframe` as soon as :mod:`openalea.mtg` permits to add/remove components.        
    
    .. todo:: update the MTG from CNWheat roots, phloem, grains and soil

    """
    # create a mapping of topology IDs to CNWheat objects
    cnwheat_objects_mapping = {}
    plant_index = 1
    for plant in population.plants:
        cnwheat_objects_mapping[(plant_index,)] = plant
        axis_index = 0
        for axis in plant.axes:
            axis_id = model.Axis.get_axis_id(axis_index)
            cnwheat_objects_mapping[(plant_index, axis_id)] = axis
            for organ in (axis.roots, axis.soil, axis.phloem, axis.grains):
                if organ is not None:
                    mtg_organ_label = CNWHEAT_CLASSES_TO_MTG_ORGANS_MAPPING[organ.__class__]
                    cnwheat_objects_mapping[(plant_index, axis_id, None, mtg_organ_label)] = organ
            phytomer_index = 9 # TODO: replace by phytomer_index = 1
            for phytomer in axis.phytomers:
                cnwheat_objects_mapping[(plant_index, axis_id, phytomer_index)] = phytomer
                for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                    if organ is None:
                        continue
                    mtg_organ_label = CNWHEAT_CLASSES_TO_MTG_ORGANS_MAPPING[organ.__class__]
                    cnwheat_objects_mapping[(plant_index, axis_id, phytomer_index, mtg_organ_label)] = organ
                    for cnwheat_element_name, cnwheat_classes_to_mtg_elements_names_mapping in CNWHEAT_TO_MTG_ELEMENTS_NAMES_MAPPING.iteritems():
                        element = getattr(organ, cnwheat_element_name)
                        if element is None:
                            continue
                        mtg_element_label = cnwheat_classes_to_mtg_elements_names_mapping[organ.__class__]
                        cnwheat_objects_mapping[(plant_index, axis_id, phytomer_index, mtg_organ_label, mtg_element_label)] = element
                        
                phytomer_index += 1
            axis_index += 1
        plant_index += 1
    
    # add the properties if needed
    property_names = g.property_names()
    for cnwheat_variable_name in CNWHEAT_OUTPUTS:
        if cnwheat_variable_name not in property_names:
            g.add_property(cnwheat_variable_name)
            
    # traverse the MTG recursively from top
    for plant_vertex_id in g.components_iter(g.root):
        plant_index = int(g.index(plant_vertex_id))
        if g.index(plant_vertex_id) not in available_components: continue
        plant_topology_id = (plant_index,)
        if plant_topology_id not in cnwheat_objects_mapping:
            raise MismatchedTopologiesError(plant_topology_id, plant_vertex_id)
        for axis_vertex_id in g.components_iter(plant_vertex_id):
            if (g.index(plant_vertex_id), g.label(axis_vertex_id)) not in available_components: continue
            axis_id = g.label(axis_vertex_id)
            axis_topology_id = (plant_index, axis_id)
            if axis_topology_id not in cnwheat_objects_mapping:
                raise MismatchedTopologiesError(axis_topology_id, axis_vertex_id)
            for axis_component_vertex_id in g.components_iter(axis_vertex_id):
                if (g.index(plant_vertex_id), g.label(axis_vertex_id), g.index(axis_component_vertex_id)) not in available_components: continue
                mtg_axis_component_label = g.label(axis_component_vertex_id)
                if mtg_axis_component_label in MTG_TO_CNWHEAT_AXES_ORGANS_NAMES_MAPPING: # grains, phloem, roots, soil
                    organ_vertex_id = axis_component_vertex_id
                    mtg_organ_label = mtg_axis_component_label
                    organ_topology_id = (plant_index, axis_id, None, mtg_organ_label)
                    if organ_topology_id not in cnwheat_objects_mapping:
                        raise MismatchedTopologiesError(organ_topology_id, organ_vertex_id)
                    organ = cnwheat_objects_mapping[(plant_index, axis_id, None, mtg_organ_label)]
                    organ_compartment_names = [compartment_name for compartment_name in simulation.Simulation.MODEL_COMPARTMENTS_NAMES[model.Organ] if hasattr(organ, compartment_name)]
                    #TODO: is it possible to set all properties in the MTG without iterating ?
                    for organ_compartment_name in organ_compartment_names:
                        g.property(organ_compartment_name)[organ_vertex_id] = getattr(organ, organ_compartment_name)
                elif mtg_axis_component_label.startswith('metamer'): # phytomer
                    metamer_vertex_id = axis_component_vertex_id
                    metamer_index = int(g.index(metamer_vertex_id))
                    metamer_topology_id = (plant_index, axis_id, metamer_index)
                    if metamer_topology_id not in cnwheat_objects_mapping:
                        raise MismatchedTopologiesError(metamer_topology_id, metamer_vertex_id)
                    for organ_vertex_id in g.components_iter(metamer_vertex_id):
                        if (g.index(plant_vertex_id), g.label(axis_vertex_id), g.index(axis_component_vertex_id), g.label(organ_vertex_id)) not in available_components: continue
                        mtg_organ_label = g.label(organ_vertex_id)
                        organ_topology_id = (plant_index, axis_id, metamer_index, mtg_organ_label)
                        if organ_topology_id not in cnwheat_objects_mapping:
                            raise MismatchedTopologiesError(organ_topology_id, organ_vertex_id)
                        for element_vertex_id in g.components_iter(organ_vertex_id):
                            if (g.index(plant_vertex_id), g.label(axis_vertex_id), g.index(axis_component_vertex_id), g.label(organ_vertex_id), g.label(element_vertex_id)) not in available_components: continue
                            mtg_element_label = g.label(element_vertex_id)
                            element_topology_id = (plant_index, axis_id, metamer_index, mtg_organ_label, mtg_element_label)
                            if element_topology_id not in cnwheat_objects_mapping:
                                raise MismatchedTopologiesError(element_topology_id, element_vertex_id)
                            element = cnwheat_objects_mapping[(plant_index, axis_id, metamer_index, mtg_organ_label, mtg_element_label)]
                            #TODO: is it possible to set all properties in the MTG without iterating ?
                            for element_compartment_name in simulation.Simulation.MODEL_COMPARTMENTS_NAMES[model.PhotosyntheticOrganElement]:
                                g.property(element_compartment_name)[element_vertex_id] = getattr(element, element_compartment_name)
    
    #TODO: remove the following code as soon as :mod:`openalea.mtg` permits to add/remove components.
    plant_index = 1
    for plant in population.plants:
        axis_index = 0
        for axis in plant.axes:
            axis_id = model.Axis.get_axis_id(axis_index)
            for organ in (axis.roots, axis.soil, axis.phloem, axis.grains):
                organ_attributes = dict([(state_var_name, getattr(organ, state_var_name)) for state_var_name in simulation.Simulation.ORGANS_STATE if hasattr(organ, state_var_name)])
                organs_dataframe.loc[(organs_dataframe['plant'] == plant_index) & (organs_dataframe['axis'] == axis_id) & (organs_dataframe['organ'] == organ.__class__.__name__), organ_attributes.keys()] = organ_attributes.values()
            axis_index += 1
        plant_index += 1
    
