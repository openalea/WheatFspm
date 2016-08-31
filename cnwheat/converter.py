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

#: the topology and state variables at ORGAN scale
HGZS_STATE_VARIABLES = simulation.Simulation.HIDDENGROWINGZONE_INPUTS_INDEXES + simulation.Simulation.HIDDENGROWINGZONE_STATE

#: the topology and state variables at ELEMENT scale
ELEMENTS_STATE_VARIABLES = simulation.Simulation.ELEMENTS_INPUTS_INDEXES + simulation.Simulation.ELEMENTS_STATE

#: the topology and state variables of soils
SOILS_STATE_VARIABLES = simulation.Simulation.SOILS_INPUTS_INDEXES + simulation.Simulation.SOILS_STATE

#: the columns, in the input/output dataframes, which define the topology of the soils
SOILS_TOPOLOGY_COLUMNS = ['plant', 'axis']

#: the mapping of CNWheat organ classes to the attributes in axis and phytomer which represent an organ
CNWHEAT_ATTRIBUTES_MAPPING = {model.Internode: 'internode', model.Lamina: 'lamina',
                              model.Sheath: 'sheath', model.Peduncle: 'peduncle', model.Chaff: 'chaff',
                              model.Roots: 'roots', model.Grains: 'grains', model.Phloem: 'phloem',
                              model.HiddenGrowingZone: 'hgz'}

#: the mapping of the CNWheat organ classes to organ names in MTG
CNWHEAT_CLASSES_TO_MTG_ORGANS_MAPPING = {model.Internode: 'internode', model.Lamina: 'blade',
                                         model.Sheath: 'sheath', model.Peduncle: 'peduncle', model.Chaff: 'ear',
                                         model.Roots: 'roots', model.Grains: 'grains', model.Phloem: 'phloem', model.HiddenGrowingZone: 'hgz'}

#: the mapping of the name of each element, from MTG to CNWheat
MTG_TO_CNWHEAT_ELEMENTS_NAMES_MAPPING = {'HiddenElement': 'enclosed_element', 'StemElement': 'exposed_element', 'LeafElement1': 'exposed_element'}

#: the mapping of organs (which belong to an axis) labels in MTG to organ classes in CNWheat
MTG_TO_CNWHEAT_AXES_ORGANS_MAPPING = {'grains': model.Grains, 'phloem': model.Phloem, 'roots': model.Roots}

#: the mapping of organs (which belong to a phytomer) labels in MTG to organ classes in CNWheat
MTG_TO_CNWHEAT_PHYTOMERS_ORGANS_MAPPING = {'internode': model.Internode, 'blade': model.Lamina, 'sheath': model.Sheath, 'peduncle': model.Peduncle, 'ear': model.Chaff, 'hgz': model.HiddenGrowingZone}

# the mapping of CNWheat photosynthetic organs to CNWheat photosynthetic organ elements
CNWHEAT_ORGANS_TO_ELEMENTS_MAPPING = {model.Internode: model.InternodeElement, model.Lamina: model.LaminaElement, model.Sheath: model.SheathElement, model.Peduncle: model.PeduncleElement, model.Chaff: model.ChaffElement}

#: the parameters and variables which define the state of a CNWheat population
POPULATION_STATE_VARIABLE = set(simulation.Simulation.PLANTS_STATE + simulation.Simulation.AXES_STATE +
                                simulation.Simulation.PHYTOMERS_STATE + simulation.Simulation.ORGANS_STATE +
                                simulation.Simulation.HIDDENGROWINGZONE_STATE + simulation.Simulation.ELEMENTS_STATE)


def from_dataframes(plants_inputs=None, axes_inputs=None, metamers_inputs=None, organs_inputs=None, hgzs_inputs=None, elements_inputs=None, soils_inputs=None):
    """
    If `plants_inputs`, `axes_inputs`, `metamers_inputs`, `organs_inputs`, `hgzs_inputs` and `elements_inputs` are not `None`, convert `plants_inputs`,
    `axes_inputs`, `metamers_inputs`, `organs_inputs`, `hgzs_inputs` and  `elements_inputs` to a :class:`population <model.Population>`.

    If `soils_inputs` is not `None`, convert `soils_inputs` to a dictionary of :class:`soils <model.Soil>`.

    :Parameters:

        - `plants_inputs` (:class:`pandas.DataFrame`) - Plants inputs, with one line by plant.

        - `axes_inputs` (:class:`pandas.DataFrame`) - Axes inputs, with one line by axis.

        - `metamers_inputs` (:class:`pandas.DataFrame`) - Metamers inputs, with one line by metamer.

        - `organs_inputs` (:class:`pandas.DataFrame`) - Organs inputs, with one line by organ.

        - `hgzs_inputs` (:class:`pandas.DataFrame`) - Hidden growing zones inputs, with one line by hidden growing zone.

        - `elements_inputs` (:class:`pandas.DataFrame`) - Elements inputs, with one line by element.

        - `soils_inputs` (:class:`pandas.DataFrame`) - Soils inputs, with one line by soil.

    :Returns:
        If `plants_inputs`, `axes_inputs`, `metamers_inputs`, `organs_inputs`, `hgzs_inputs` and `elements_inputs` are not `None`, return a :class:`population <model.Population>`,
        and/or
        if `soils_inputs` is not `None`,  return a :class:`dict` of :class:`soils <model.Soil>`.

    :Returns Type:
        :class:`tuple` or :class:`dict`.

    """

    convert_dataframes_to_population = plants_inputs is not None and axes_inputs is not None and metamers_inputs is not None and organs_inputs is not None and hgzs_inputs is not None and elements_inputs is not None
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

                    for phytomer_attribute_name, phytomer_attribute_class, phytomer_attribute_element_class in \
                        (('chaff', model.Chaff, model.ChaffElement),
                         ('lamina', model.Lamina, model.LaminaElement),
                         ('internode', model.Internode, model.InternodeElement),
                         ('peduncle', model.Peduncle, model.PeduncleElement),
                         ('sheath', model.Sheath, model.SheathElement)):

                        organ_label = CNWHEAT_CLASSES_TO_MTG_ORGANS_MAPPING[phytomer_attribute_class]

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

                    hgz_inputs = hgzs_inputs[(hgzs_inputs['plant'] == plant_index) & (hgzs_inputs['axis'] == axis_label) & (hgzs_inputs['metamer'] == metamer_index)]
                    if len(hgz_inputs) == 0:
                        continue
                    hgz_inputs = hgz_inputs.loc[:, simulation.Simulation.HIDDENGROWINGZONE_STATE]
                    hgz_dict = hgz_inputs.loc[hgz_inputs.first_valid_index()].to_dict()
                    # create a new hidden growing zone
                    hgz = model.HiddenGrowingZone(CNWHEAT_CLASSES_TO_MTG_ORGANS_MAPPING[model.HiddenGrowingZone], **hgz_dict)
                    hgz.initialize()
                    phytomer.hiddengrowingzone = hgz

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

            * plant scale: plant index, state parameters and compartments of each plant (see :attr:`Simulation:PLANTS_STATE_VARIABLES`)
            * axis scale: plant index, axis id, state parameters and compartments of each axis (see :attr:`Simulation:AXES_STATE_VARIABLES`)
            * metamer scale: plant index, axis id, metamer index, state parameters and compartments of each metamer (see :attr:`Simulation:PHYTOMERS_STATE_VARIABLES`)
            * organ scale:
                * hidden growing zones: plant index, axis id, metamer index, state parameters and compartments of each organ (see :attr:`Simulation:HIDDENGROWINGZONE_STATE_VARIABLES`)
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
        all_plants_df = pd.DataFrame(columns=PLANTS_STATE_VARIABLES)
        all_axes_df = pd.DataFrame(columns=AXES_STATE_VARIABLES)
        all_metamers_df = pd.DataFrame(columns=PHYTOMERS_STATE_VARIABLES)
        all_organs_df = pd.DataFrame(columns=ORGANS_STATE_VARIABLES)
        all_hgzs_df = pd.DataFrame(columns=HGZS_STATE_VARIABLES)
        all_elements_df = pd.DataFrame(columns=ELEMENTS_STATE_VARIABLES)

        # run through the population tree and fill the dataframes
        for plant in population.plants:
            append_row(plant, [plant.index], simulation.Simulation.PLANTS_STATE, all_plants_df)
            for axis in plant.axes:
                append_row(axis, [plant.index, axis.label], simulation.Simulation.AXES_STATE, all_axes_df)
                for organ in (axis.roots, axis.phloem, axis.grains):
                    append_row(organ, [plant.index, axis.label, organ.label], simulation.Simulation.ORGANS_STATE, all_organs_df)
                for phytomer in axis.phytomers:
                    append_row(phytomer, [plant.index, axis.label, phytomer.index], simulation.Simulation.PHYTOMERS_STATE, all_metamers_df)
                    if phytomer.hiddengrowingzone is not None:
                        append_row(phytomer.hiddengrowingzone, [plant.index, axis.label, phytomer.index], simulation.Simulation.HIDDENGROWINGZONE_STATE, all_hgzs_df)
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ is None:
                            continue
                        for element in (organ.exposed_element, organ.enclosed_element):
                            if element is None:
                                continue
                            append_row(element, [plant.index, axis.label, phytomer.index, organ.label, element.label], simulation.Simulation.ELEMENTS_STATE, all_elements_df)

        # sort the rows of the dataframes by columns
        all_plants_df.sort_values(by=PLANTS_STATE_VARIABLES, inplace=True)
        all_axes_df.sort_values(by=AXES_STATE_VARIABLES, inplace=True)
        all_metamers_df.sort_values(by=PHYTOMERS_STATE_VARIABLES, inplace=True)
        all_organs_df.sort_values(by=ORGANS_STATE_VARIABLES, inplace=True)
        all_hgzs_df.sort_values(by=HGZS_STATE_VARIABLES, inplace=True)
        all_elements_df.sort_values(by=ELEMENTS_STATE_VARIABLES, inplace=True)

        # infer the right types of the columns in the dataframes
        all_plants_df = all_plants_df.convert_objects(copy=False)
        all_axes_df = all_axes_df.convert_objects(copy=False)
        all_metamers_df = all_metamers_df.convert_objects(copy=False)
        all_organs_df = all_organs_df.convert_objects(copy=False)
        all_hgzs_df = all_hgzs_df.convert_objects(copy=False)
        all_elements_df = all_elements_df.convert_objects(copy=False)

        # convert the indexes of plants, metamers and elements to integers in the dataframes
        all_plants_df['plant'] = all_plants_df['plant'].astype(int)
        all_axes_df['plant'] = all_axes_df['plant'].astype(int)
        all_metamers_df[['plant', 'metamer']] = all_metamers_df[['plant', 'metamer']].astype(int)
        all_organs_df['plant'] = all_organs_df['plant'].astype(int)
        all_hgzs_df[['plant', 'metamer']] = all_hgzs_df[['plant', 'metamer']].astype(int)
        all_elements_df[['plant', 'metamer']] = all_elements_df[['plant', 'metamer']].astype(int)

        all_plants_df.reset_index(drop=True, inplace=True)
        all_axes_df.reset_index(drop=True, inplace=True)
        all_metamers_df.reset_index(drop=True, inplace=True)
        all_organs_df.reset_index(drop=True, inplace=True)
        all_hgzs_df.reset_index(drop=True, inplace=True)
        all_elements_df.reset_index(drop=True, inplace=True)

    if convert_soils_to_dataframe:
        all_soils_df = pd.DataFrame(columns=SOILS_STATE_VARIABLES)
        for soil_id, soil in soils.iteritems():
            append_row(soil, list(soil_id), simulation.Simulation.SOILS_STATE, all_soils_df)
        all_soils_df.sort_values(by=SOILS_STATE_VARIABLES, inplace=True)
        all_soils_df = all_soils_df.convert_objects(copy=False)
        all_soils_df['plant'] = all_soils_df['plant'].astype(int)
        all_soils_df.reset_index(drop=True, inplace=True)

    if convert_population_to_dataframes and convert_soils_to_dataframe:
        return all_plants_df, all_axes_df, all_metamers_df, all_organs_df, all_hgzs_df, all_elements_df, all_soils_df
    elif convert_population_to_dataframes:
        return all_plants_df, all_axes_df, all_metamers_df, all_organs_df, all_hgzs_df, all_elements_df
    else:
        return all_soils_df


def from_MTG(g, organs_inputs, hgzs_inputs, elements_inputs): #TODO: update doc
    """
    Convert a MTG to a CN-Wheat :class:`population <model.Population>`.
    Use data in `organs_inputs` and `elements_inputs` if `g` is incomplete.

    :Parameters:

            - g (:class:`openalea.mtg.mtg.MTG`) - A MTG which contains the inputs
              of CN-Wheat. These inputs are: :mod:`POPULATION_STATE_VARIABLE`.

            - `organs_inputs` (:class:`pandas.DataFrame`) - Organs dataframe, with one line by organ.

            - `hgzs_inputs` (:class:`pandas.DataFrame`) - Hidden growing zones inputs, with one line by hidden growing zone.

            - `elements_inputs` (:class:`pandas.DataFrame`) - Elements dataframe, with one line by element.

    :Returns:
        A CN-Wheat :class:`population <model.Population>`.

    :Returns Type:
        :class:`population <model.Population>`

    """

    population = model.Population()

    organs_inputs_grouped = organs_inputs.groupby(simulation.Simulation.ORGANS_INPUTS_INDEXES)
    hgzs_inputs_grouped = hgzs_inputs.groupby(simulation.Simulation.HIDDENGROWINGZONE_INPUTS_INDEXES)
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
            is_valid_axis = True
            for organ_class in (model.Roots, model.Phloem, model.Grains):
                organ_label = CNWHEAT_CLASSES_TO_MTG_ORGANS_MAPPING[organ_class]
                organ_id = (plant_index, axis_label, organ_label)
                # create a new organ
                organ = organ_class(organ_label)
                organ_input_names = set(simulation.Simulation.ORGANS_STATE).intersection(organ.__dict__)
                if organ_id in organs_inputs_grouped.groups:
                    organ_inputs_group = organs_inputs_grouped.get_group(organ_id)
                    organ_inputs_group_series = organ_inputs_group.loc[organ_inputs_group.first_valid_index(), organ_inputs_group.columns.intersection(organ_input_names)]
                else:
                    organ_inputs_group_series = pd.Series()
                axis_properties = g.get_vertex_property(axis_vid)
                if organ_label in axis_properties:
                    organ_properties = axis_properties[organ_label]
                    is_valid_organ = True
                    organ_inputs_dict = {}
                    for organ_input_name in organ_input_names:
                        if organ_input_name in organ_properties:
                            # use the input from the MTG
                            organ_inputs_dict[organ_input_name] = organ_properties[organ_input_name]
                        else:
                            # use the input from the series
                            if organ_input_name in organ_inputs_group_series:
                                organ_inputs_dict[organ_input_name] = organ_inputs_group_series[organ_input_name]
                            else:
                                is_valid_organ = False
                                break
                else:
                    organ_inputs_dict = organ_inputs_group_series.to_dict()
                    if not set(organ_inputs_dict).issuperset(organ_input_names):
                        is_valid_organ = False
                if is_valid_organ:
                    organ.__dict__.update(organ_inputs_dict)
                    organ.initialize()
                    # add the new organ to current axis
                    setattr(axis, organ_label, organ)
                elif organ_class is not model.Grains:
                    is_valid_axis = False
                    break

            if not is_valid_axis:
                continue

            has_valid_phytomer = False
            for metamer_vid in g.components_iter(axis_vid):
                metamer_index = int(g.index(metamer_vid))
                # create a new phytomer
                phytomer = model.Phytomer(metamer_index)

                hgz_id = (plant_index, axis_label, metamer_index)
                if hgz_id in hgzs_inputs_grouped.groups:
                    hgz_inputs_group = hgzs_inputs_grouped.get_group(hgz_id)
                    hgz_inputs_group_series = hgz_inputs_group.loc[hgz_inputs_group.first_valid_index(), hgz_inputs_group.columns.intersection(simulation.Simulation.HIDDENGROWINGZONE_STATE)]
                else:
                    hgz_inputs_group_series = pd.Series()
                hgz_label = CNWHEAT_CLASSES_TO_MTG_ORGANS_MAPPING[model.HiddenGrowingZone]
                metamer_properties = g.get_vertex_property(metamer_vid)
                if hgz_label in metamer_properties:
                    hgz_properties = metamer_properties[hgz_label]
                    is_valid_hgz = True
                    hgz_inputs_dict = {}
                    for hgz_input_name in simulation.Simulation.HIDDENGROWINGZONE_STATE:
                        if hgz_input_name in hgz_properties:
                            # use the input from the MTG
                            hgz_inputs_dict[hgz_input_name] = hgz_properties[hgz_input_name]
                        else:
                            # use the input from the series
                            if hgz_input_name in hgz_inputs_group_series:
                                hgz_inputs_dict[hgz_input_name] = hgz_inputs_group_series[hgz_input_name]
                            else:
                                is_valid_hgz = False
                                break
                else:
                    hgz_inputs_dict = hgz_inputs_group_series.to_dict()
                    if not set(hgz_inputs_dict).issuperset(simulation.Simulation.HIDDENGROWINGZONE_STATE):
                        is_valid_hgz = False
                if is_valid_hgz:
                    # create a new hgz
                    hgz = model.HiddenGrowingZone(hgz_label, **hgz_inputs_dict)
                    hgz.initialize()
                    # add the new hgz to current phytomer
                    setattr(phytomer, hgz_label, hgz)

                has_valid_organ = False

                for organ_vid in g.components_iter(metamer_vid):
                    organ_label = g.label(organ_vid)
                    if organ_label not in MTG_TO_CNWHEAT_PHYTOMERS_ORGANS_MAPPING or g.get_vertex_property(organ_vid)['length']==0: continue
                    # create a new organ
                    organ_class = MTG_TO_CNWHEAT_PHYTOMERS_ORGANS_MAPPING[organ_label]
                    organ = organ_class(organ_label)
                    organ.initialize()
                    has_valid_element = False
                    for element_vid in g.components_iter(organ_vid):
                        vertex_properties = g.get_vertex_property(element_vid)
                        element_label = g.label(element_vid)
                        if element_label not in MTG_TO_CNWHEAT_ELEMENTS_NAMES_MAPPING: continue
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
                                # use the input from the MTG
                                element_inputs[element_input_name] = vertex_properties[element_input_name]
                                if element_input_name == 'green_area':
                                    element_inputs[element_input_name] /= 10000.0 # convert from cm2 to m2 ; TODO: homogenize the units between the models
                            else:
                                # use the input from the series
                                if element_input_name in elements_inputs_group_series:
                                    element_inputs[element_input_name] = elements_inputs_group_series[element_input_name]
                                else:
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

    :Parameters:

            - population (:class:`model.Population`) - a CN-Wheat :class:`population <model.Population>`.

            - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG to update from the CN-Wheat :class:`population <model.Population>`.

    .. seealso:: see :mod:`model` for the structure of the population.

    """

    # add the missing properties
    property_names = g.property_names()
    for cnwheat_data_name in POPULATION_STATE_VARIABLE:
        if cnwheat_data_name not in property_names:
            g.add_property(cnwheat_data_name)
    for organ_label in MTG_TO_CNWHEAT_AXES_ORGANS_MAPPING.keys() + [CNWHEAT_CLASSES_TO_MTG_ORGANS_MAPPING[model.HiddenGrowingZone]]:
        if organ_label not in property_names:
            g.add_property(organ_label)

    plants_iterator = g.components_iter(g.root)
    # traverse CN-Wheat population from top
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
            for organ_label in MTG_TO_CNWHEAT_AXES_ORGANS_MAPPING.iterkeys():
                if organ_label not in g.get_vertex_property(axis_vid):
                    # Add a property describing the organ to the current axis of the MTG
                    g.property(organ_label)[axis_vid] = {}
                # Update the property describing the organ of the current axis in the MTG
                organ = getattr(axis, organ_label)
                organ_properties = g.get_vertex_property(axis_vid)[organ_label]
                for property_name in simulation.Simulation.ORGANS_STATE:
                    if hasattr(organ, property_name):
                        organ_properties[property_name] = getattr(organ, property_name)
            metamers_iterator = g.components_iter(axis_vid)
            for phytomer in axis.phytomers:
                phytomer_index = phytomer.index
                while True:
                    metamer_vid = next(metamers_iterator)
                    if int(g.index(metamer_vid)) == phytomer_index:
                        break
                if phytomer.hiddengrowingzone is not None:
                    hgz_label = CNWHEAT_CLASSES_TO_MTG_ORGANS_MAPPING[model.HiddenGrowingZone]
                    if hgz_label not in g.get_vertex_property(metamer_vid):
                        # Add a property describing the hgz to the current metamer of the MTG
                        g.property(hgz_label)[metamer_vid] = {}
                    # Update the property describing the hgz of the current metamer in the MTG
                    hgz_properties = g.get_vertex_property(metamer_vid)[hgz_label]
                    hgz_properties.update(phytomer.hiddengrowingzone.__dict__)
                for organ_vid in g.components_iter(metamer_vid):
                    organ_label = g.label(organ_vid)
                    if organ_label not in MTG_TO_CNWHEAT_PHYTOMERS_ORGANS_MAPPING: continue
                    organ = getattr(phytomer, CNWHEAT_ATTRIBUTES_MAPPING[MTG_TO_CNWHEAT_PHYTOMERS_ORGANS_MAPPING[organ_label]])
                    if organ is None: continue
                    organ_property_names = [property_name for property_name in simulation.Simulation.ORGANS_STATE if hasattr(organ, property_name)]
                    for organ_property_name in organ_property_names:
                        g.property(organ_property_name)[organ_vid] = getattr(organ, organ_property_name)
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

