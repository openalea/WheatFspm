
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
    
    for curr_plant_index in plants_dataframe.plant:
        # create a new plant
        curr_plant = model.Plant(index=curr_plant_index)
        population.plants.append(curr_plant)
        curr_axes_dataframe = axes_dataframe[axes_dataframe['plant'] == curr_plant_index]
        for curr_axis_id in curr_axes_dataframe.axis:
            # create a new axis
            curr_axis = model.Axis(axis_id=curr_axis_id)
            curr_plant.axes.append(curr_axis)
            
            curr_organs_dataframe = organs_dataframe[(organs_dataframe['plant'] == curr_plant_index) & (organs_dataframe['axis'] == curr_axis_id)]
            
            for curr_axis_attribute_name, curr_axis_attribute_class in (('grains', model.Grains), ('roots', model.Roots), ('soil', model.Soil), ('phloem', model.Phloem)):
                curr_organ_class_name = curr_axis_attribute_class.__name__
                curr_organ_dataframe = curr_organs_dataframe[curr_organs_dataframe['organ'] == curr_organ_class_name]
                if len(curr_organ_dataframe) != 0:
                    # create a new organ
                    curr_organ = curr_axis_attribute_class()
                    curr_organ_attributes_names = [state_var_name for state_var_name in simulation.Simulation.ORGANS_STATE if hasattr(curr_organ, state_var_name)]
                    curr_organ_row = curr_organ_dataframe.loc[curr_organ_dataframe.first_valid_index()]
                    curr_organ_attributes_values = curr_organ_row[curr_organ_attributes_names].tolist()
                    curr_organ_attributes = dict(zip(curr_organ_attributes_names, curr_organ_attributes_values))
                    curr_organ.__dict__.update(curr_organ_attributes)
                    curr_organ.initialize()
                    setattr(curr_axis, curr_axis_attribute_name, curr_organ)
                    
            curr_phytomers_dataframe = phytomers_dataframe[(phytomers_dataframe['plant'] == curr_plant_index) & (phytomers_dataframe['axis'] == curr_axis_id)]
            for curr_phytomer_index in curr_phytomers_dataframe.phytomer:
                # create a new phytomer
                curr_phytomer = model.Phytomer(index=curr_phytomer_index)
                curr_axis.phytomers.append(curr_phytomer)
                
                curr_organs_dataframe = organs_dataframe[(organs_dataframe['plant'] == curr_plant_index) & (organs_dataframe['axis'] == curr_axis_id) & (organs_dataframe['phytomer'] == curr_phytomer_index)]
                
                for curr_phytomer_attribute_name, curr_phytomer_attribute_class, curr_phytomer_attribute_element_class in \
                    (('chaff', model.Chaff, model.ChaffElement), 
                     ('lamina', model.Lamina, model.LaminaElement), 
                     ('internode', model.Internode, model.InternodeElement), 
                     ('peduncle', model.Peduncle, model.PeduncleElement),
                     ('sheath', model.Sheath, model.SheathElement)):
                    
                    curr_organ_class_name = curr_phytomer_attribute_class.__name__
                    curr_elements_dataframe = elements_dataframe[(elements_dataframe['plant'] == curr_plant_index) & (elements_dataframe['axis'] == curr_axis_id) & (elements_dataframe['phytomer'] == curr_phytomer_index) & (elements_dataframe['organ'] == curr_organ_class_name)]
                
                    if curr_organ_class_name not in curr_organs_dataframe.organ.values and curr_organ_class_name not in curr_elements_dataframe.organ.values:
                        continue
                    # create a new organ
                    curr_organ = curr_phytomer_attribute_class()
                    setattr(curr_phytomer, curr_phytomer_attribute_name, curr_organ)
                    
                    for curr_organ_attribute_name, curr_organ_attribute_type in (('enclosed_element', 'enclosed'), ('exposed_element', 'exposed')):
                        curr_element_dataframe = curr_elements_dataframe[curr_elements_dataframe['element'] == curr_organ_attribute_type]
                        if len(curr_element_dataframe) == 0:
                            continue
                        curr_element_dataframe = curr_element_dataframe.loc[:, simulation.Simulation.ELEMENTS_STATE]
                        # create a new element
                        curr_element_dict = curr_element_dataframe.loc[curr_element_dataframe.first_valid_index()].to_dict()
                        curr_element = curr_phytomer_attribute_element_class(**curr_element_dict)
                        setattr(curr_organ, curr_organ_attribute_name, curr_element)
                        
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
    for plant in population.plants:
        append_row(plant, [plant.index], simulation.Simulation.PLANTS_STATE, all_plants_df)
        for axis in plant.axes:
            append_row(axis, [plant.index, axis.id], simulation.Simulation.AXES_STATE, all_axes_df)
            for organ in (axis.roots, axis.soil, axis.phloem, axis.grains):
                if organ is None:
                    continue
                append_row(organ, [plant.index, axis.id, np.nan, organ.__class__.__name__], simulation.Simulation.ORGANS_STATE, all_organs_df)
            for phytomer in axis.phytomers:
                append_row(phytomer, [plant.index, axis.id, phytomer.index], simulation.Simulation.PHYTOMERS_STATE, all_phytomers_df)
                for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                    if organ is None:
                        continue
                    append_row(organ, [plant.index, axis.id, phytomer.index, organ.__class__.__name__], simulation.Simulation.ORGANS_STATE, all_organs_df)
                    for element, element_type in ((organ.exposed_element, 'exposed'), (organ.enclosed_element, 'enclosed')):
                        if element is None:
                            continue
                        append_row(element, [plant.index, axis.id, phytomer.index, organ.__class__.__name__, element_type], simulation.Simulation.ELEMENTS_STATE, all_elements_df)
            
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
    