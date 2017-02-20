# -*- coding: latin-1 -*-

"""
    cnwheat.converter
    ~~~~~~~~~~~~~~~~~

    The module :mod:`cnwheat.converter` defines functions to convert
    :class:`dataframes <pandas.DataFrame>` to/from CNWheat inputs or outputs format.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2016.
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


#: the topology and state variables at PLANT scale
PLANTS_STATE_VARIABLES = simulation.Simulation.PLANTS_INPUTS_INDEXES + simulation.Simulation.PLANTS_STATE

#: the topology and state variables at AXIS scale
AXES_STATE_VARIABLES = simulation.Simulation.AXES_INPUTS_INDEXES + simulation.Simulation.AXES_STATE

#: the topology and state variables at PHYTOMER scale
PHYTOMERS_STATE_VARIABLES = simulation.Simulation.PHYTOMERS_INPUTS_INDEXES + simulation.Simulation.PHYTOMERS_STATE

#: the topology and state variables at ORGAN scale
ORGANS_STATE_VARIABLES = simulation.Simulation.ORGANS_INPUTS_INDEXES + simulation.Simulation.ORGANS_STATE

#: the topology and state variables at ORGAN scale
HIDDENZONES_STATE_VARIABLES = simulation.Simulation.HIDDENZONE_INPUTS_INDEXES + simulation.Simulation.HIDDENZONE_STATE

#: the topology and state variables at ELEMENT scale
ELEMENTS_STATE_VARIABLES = simulation.Simulation.ELEMENTS_INPUTS_INDEXES + simulation.Simulation.ELEMENTS_STATE

#: the topology and state variables of soils
SOILS_STATE_VARIABLES = simulation.Simulation.SOILS_INPUTS_INDEXES + simulation.Simulation.SOILS_STATE

#: the columns, in the input/output dataframes, which define the topology of the soils
SOILS_TOPOLOGY_COLUMNS = ['plant', 'axis']

#: the mapping of the CNWheat organ classes to organ names in MTG
CNWHEAT_CLASSES_TO_MTG_ORGANS_MAPPING = {model.Internode: 'internode', model.Lamina: 'blade',
                                         model.Sheath: 'sheath', model.Peduncle: 'peduncle', model.Chaff: 'ear',
                                         model.Roots: 'roots', model.Grains: 'grains', model.Phloem: 'phloem', model.HiddenZone: 'hiddenzone'}

#: the mapping of the name of each element, from MTG to CNWheat
MTG_TO_CNWHEAT_ELEMENTS_NAMES_MAPPING = {'HiddenElement': 'enclosed_element', 'StemElement': 'exposed_element', 'LeafElement1': 'exposed_element'}


def from_dataframes(organs_inputs=None, hiddenzones_inputs=None, elements_inputs=None, soils_inputs=None):
    """
    If `organs_inputs`, `hiddenzones_inputs` and `elements_inputs` are not `None`,
    convert `organs_inputs`, `hiddenzones_inputs` and  `elements_inputs` to a :class:`population <model.Population>`.

    If `soils_inputs` is not `None`, convert `soils_inputs` to a dictionary of :class:`soils <model.Soil>`.

    :Parameters:

        - `organs_inputs` (:class:`pandas.DataFrame`) - Organs inputs, with one line by organ.

        - `hiddenzones_inputs` (:class:`pandas.DataFrame`) - Hidden zones inputs, with one line by hidden zone.

        - `elements_inputs` (:class:`pandas.DataFrame`) - Elements inputs, with one line by element.

        - `soils_inputs` (:class:`pandas.DataFrame`) - Soils inputs, with one line by soil.

    :Returns:
        If `organs_inputs`, `hiddenzones_inputs` and `elements_inputs` are not `None`, return a :class:`population <model.Population>`,
        and/or
        if `soils_inputs` is not `None`,  return a :class:`dict` of :class:`soils <model.Soil>`.

    :Returns Type:
        :class:`tuple` or :class:`dict`.

    """

    convert_dataframes_to_population = organs_inputs is not None and hiddenzones_inputs is not None and elements_inputs is not None
    convert_dataframe_to_soils_dict = soils_inputs is not None

    if convert_dataframes_to_population:
        population = model.Population()

        for plant_index in organs_inputs.plant.unique():
            # create a new plant
            plant = model.Plant(plant_index)
            population.plants.append(plant)
            curr_axes_labels = organs_inputs[organs_inputs['plant'] == plant_index].axis.unique()
            for axis_label in curr_axes_labels:
                # create a new axis
                axis = model.Axis(axis_label)
                curr_organs_inputs = organs_inputs[(organs_inputs['plant'] == plant_index) & (organs_inputs['axis'] == axis_label)]
                for axis_attribute_name, axis_attribute_class in (('roots', model.Roots), ('phloem', model.Phloem), ('grains', model.Grains)):
                    organ_label = CNWHEAT_CLASSES_TO_MTG_ORGANS_MAPPING[axis_attribute_class]
                    organ_inputs = curr_organs_inputs[curr_organs_inputs['organ'] == organ_label]
                    if not organ_inputs.empty:
                        # create a new organ
                        organ = axis_attribute_class(organ_label)
                        organ_attributes_names = [state_var_name for state_var_name in simulation.Simulation.ORGANS_STATE if hasattr(organ, state_var_name)]
                        organ_row = organ_inputs.loc[organ_inputs.first_valid_index()]
                        organ_attributes_values = organ_row[organ_attributes_names].tolist()
                        organ_attributes = dict(zip(organ_attributes_names, organ_attributes_values))
                        organ.__dict__.update(organ_attributes)
                        organ.initialize()
                        setattr(axis, axis_attribute_name, organ)

                curr_metamers_indexes_for_hiddenzones = hiddenzones_inputs[(hiddenzones_inputs['plant'] == plant_index) & (hiddenzones_inputs['axis'] == axis_label)].metamer.unique()
                curr_metamers_indexes_for_elements = elements_inputs[(elements_inputs['plant'] == plant_index) & (elements_inputs['axis'] == axis_label)].metamer.unique()
                curr_metamers_indexes = np.unique(np.concatenate((curr_metamers_indexes_for_hiddenzones,
                                                                  curr_metamers_indexes_for_elements)))
                for metamer_index in curr_metamers_indexes:
                    # create a new phytomer
                    phytomer = model.Phytomer(metamer_index)
                    axis.phytomers.append(phytomer)

                    for phytomer_attribute_name, phytomer_attribute_class, phytomer_attribute_element_class in \
                        (('chaff', model.Chaff, model.ChaffElement),
                         ('lamina', model.Lamina, model.LaminaElement),
                         ('internode', model.Internode, model.InternodeElement),
                         ('peduncle', model.Peduncle, model.PeduncleElement),
                         ('sheath', model.Sheath, model.SheathElement)):

                        organ_label = CNWHEAT_CLASSES_TO_MTG_ORGANS_MAPPING[phytomer_attribute_class]

                        if metamer_index in curr_metamers_indexes_for_elements:
                            curr_elements_inputs = elements_inputs[(elements_inputs['plant'] == plant_index) & (elements_inputs['axis'] == axis_label) & (elements_inputs['metamer'] == metamer_index) & (elements_inputs['organ'] == organ_label)]
                            if organ_label not in curr_elements_inputs.organ.values:
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
                                # create a new element
                                element = phytomer_attribute_element_class(mtg_element_label, **element_dict)
                                setattr(organ, cnwheat_element_name, element)

                    if metamer_index in curr_metamers_indexes_for_hiddenzones:
                        hiddenzone_inputs = hiddenzones_inputs[(hiddenzones_inputs['plant'] == plant_index) & (hiddenzones_inputs['axis'] == axis_label) & (hiddenzones_inputs['metamer'] == metamer_index)]
                        if len(hiddenzone_inputs) == 0:
                            continue
                        hiddenzone_inputs = hiddenzone_inputs.loc[:, simulation.Simulation.HIDDENZONE_STATE]
                        hiddenzone_dict = hiddenzone_inputs.loc[hiddenzone_inputs.first_valid_index()].to_dict()
                        # create a new hidden zone
                        hiddenzone = model.HiddenZone(CNWHEAT_CLASSES_TO_MTG_ORGANS_MAPPING[model.HiddenZone], **hiddenzone_dict)
                        hiddenzone.initialize()
                        phytomer.hiddenzone = hiddenzone

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

    :Parameters:

        - `population` (:class:`model.Population`) - The CNWheat population to convert.

        - `soils` (:class:`dict` of `model.Soil`) - The soils to convert.

    :Returns:
        If `population` is not None, return :class:`dataframes <pandas.DataFrame>` describing the internal state and compartments of the population at each scale:

            * axis scale: plant index, axis id, state parameters and compartments of each axis (see :attr:`Simulation:AXES_STATE_VARIABLES`)
            * organ scale:
                * hidden zones: plant index, axis id, metamer index, state parameters and compartments of each organ (see :attr:`Simulation:HIDDENZONE_STATE_VARIABLES`)
                * roots, phloem and grains: plant index, axis id, organ type, state parameters and compartments of each organ (see :attr:`Simulation:ORGANS_STATE_VARIABLES`)
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
        inputs_df.loc[len(inputs_df),:] = indexes + attributes_values

    if convert_population_to_dataframes:
        # initialize the dataframes
        all_axes_df = pd.DataFrame(columns=AXES_STATE_VARIABLES)
        all_organs_df = pd.DataFrame(columns=ORGANS_STATE_VARIABLES)
        all_hiddenzones_df = pd.DataFrame(columns=HIDDENZONES_STATE_VARIABLES)
        all_elements_df = pd.DataFrame(columns=ELEMENTS_STATE_VARIABLES)

        # run through the population tree and fill the dataframes
        for plant in population.plants:
            for axis in plant.axes:
                append_row(axis, [plant.index, axis.label], simulation.Simulation.AXES_STATE, all_axes_df)
                for organ in (axis.roots, axis.phloem, axis.grains):
                    if organ is not None:
                        append_row(organ, [plant.index, axis.label, organ.label], simulation.Simulation.ORGANS_STATE, all_organs_df)
                for phytomer in axis.phytomers:
                    if phytomer.hiddenzone is not None:
                        append_row(phytomer.hiddenzone, [plant.index, axis.label, phytomer.index], simulation.Simulation.HIDDENZONE_STATE, all_hiddenzones_df)
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ is None:
                            continue
                        for element in (organ.exposed_element, organ.enclosed_element):
                            if element is None:
                                continue
                            append_row(element, [plant.index, axis.label, phytomer.index, organ.label, element.label], simulation.Simulation.ELEMENTS_STATE, all_elements_df)

        # sort the rows of the dataframes by columns
        all_axes_df.sort_values(by=AXES_STATE_VARIABLES, inplace=True)
        all_organs_df.sort_values(by=ORGANS_STATE_VARIABLES, inplace=True)
        all_hiddenzones_df.sort_values(by=HIDDENZONES_STATE_VARIABLES, inplace=True)
        all_elements_df.sort_values(by=ELEMENTS_STATE_VARIABLES, inplace=True)

        # convert the indexes of plants, metamers and elements to integers in the dataframes
        all_axes_df['plant'] = all_axes_df['plant'].astype(int)
        all_organs_df['plant'] = all_organs_df['plant'].astype(int)
        all_hiddenzones_df[['plant', 'metamer']] = all_hiddenzones_df[['plant', 'metamer']].astype(int)
        all_elements_df[['plant', 'metamer']] = all_elements_df[['plant', 'metamer']].astype(int)

        all_axes_df.reset_index(drop=True, inplace=True)
        all_organs_df.reset_index(drop=True, inplace=True)
        all_hiddenzones_df.reset_index(drop=True, inplace=True)
        all_elements_df.reset_index(drop=True, inplace=True)

    if convert_soils_to_dataframe:
        all_soils_df = pd.DataFrame(columns=SOILS_STATE_VARIABLES)
        for soil_id, soil in soils.iteritems():
            append_row(soil, list(soil_id), simulation.Simulation.SOILS_STATE, all_soils_df)
        all_soils_df.sort_values(by=SOILS_STATE_VARIABLES, inplace=True)
        all_soils_df['plant'] = all_soils_df['plant'].astype(int)
        all_soils_df.reset_index(drop=True, inplace=True)

    if convert_population_to_dataframes and convert_soils_to_dataframe:
        return all_axes_df, all_organs_df, all_hiddenzones_df, all_elements_df, all_soils_df
    elif convert_population_to_dataframes:
        return all_axes_df, all_organs_df, all_hiddenzones_df, all_elements_df
    else:
        return all_soils_df

