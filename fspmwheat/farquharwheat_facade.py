# -*- coding: latin-1 -*-

"""
    fspmwheat.farquharwheat_facade
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`fspmwheat.farquharwheat_facade` is a facade of the model FarquharWheat.

    This module permits to initialize and run the model FarquharWheat from a :class:`MTG <openalea.mtg.mtg.MTG>`
    in a convenient and transparent way, wrapping all the internal complexity of the model, and dealing
    with all the tedious initialization and conversion processes.

    :copyright: Copyright 2014-2016 INRA-ECOSYS, see AUTHORS.
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

import pandas as pd

from farquharwheat import converter, simulation

import tools

#: the name of the organs modeled by FarquharWheat
FARQUHARWHEAT_ORGANS_NAMES = set(['internode', 'blade', 'sheath', 'peduncle', 'ear'])

#: the columns which define the topology in the elements scale dataframe shared between all models
SHARED_ELEMENTS_INPUTS_OUTPUTS_INDEXES = ['plant', 'axis', 'metamer', 'organ', 'element']


class FarquharWheatFacade(object):
    """
    The FarquharWheatFacade class permits to initialize, run the model FarquharWheat
    from a :class:`MTG <openalea.mtg.mtg.MTG>`, and update the MTG and the dataframes
    shared between all models.

    Use :meth:`run` to run the model.

    :Parameters:

        - `shared_mtg` (:class:`openalea.mtg.mtg.MTG`) - The MTG shared between all models.
        - `model_elements_inputs_df` (:class:`pandas.DataFrame`) - the inputs of the model at elements scale.
        - `shared_elements_inputs_outputs_df` (:class:`pandas.DataFrame`) - the dataframe of inputs and outputs at elements scale shared between all models.
    """

    def __init__(self, shared_mtg,
                 model_elements_inputs_df,
                 shared_elements_inputs_outputs_df):

        self._shared_mtg = shared_mtg #: the MTG shared between all models

        self._simulation = simulation.Simulation() #: the simulator to use to run the model

        all_farquharwheat_inputs_dict = converter.from_dataframe(model_elements_inputs_df)
        self._update_shared_MTG(all_farquharwheat_inputs_dict)

        self._shared_elements_inputs_outputs_df = shared_elements_inputs_outputs_df #: the dataframe at elements scale shared between all models
        self._update_shared_dataframes(model_elements_inputs_df)


    def run(self, Ta, ambient_CO2, RH, Ur):
        """
        Run the model and update the MTG and the dataframes shared between all models.

        :Parameters:

            - `Ta` (:class:`float`) - air temperature at t (degree Celsius)

            - `ambient_CO2` (:class:`float`) - air CO2 at t (µmol mol-1)

            - `RH` (:class:`float`) - relative humidity at t (decimal fraction)

            - `Ur` (:class:`float`) - wind speed at the top of the canopy at t (m s-1)

        """
        self._initialize_model()
        self._simulation.run(Ta, ambient_CO2, RH, Ur)
        self._update_shared_MTG(self._simulation.outputs)
        farquharwheat_elements_outputs_df = converter.to_dataframe(self._simulation.outputs)
        self._update_shared_dataframes(farquharwheat_elements_outputs_df)


    def _initialize_model(self):
        """
        Initialize the inputs of the model from the MTG shared between all models.
        """
        all_farquharwheat_elements_inputs_dict = {}

        # traverse the MTG recursively from top ...
        for mtg_plant_vid in self._shared_mtg.components_iter(self._shared_mtg.root):
            mtg_plant_index = int(self._shared_mtg.index(mtg_plant_vid))
            for mtg_axis_vid in self._shared_mtg.components_iter(mtg_plant_vid):
                mtg_axis_label = self._shared_mtg.label(mtg_axis_vid)
                for mtg_metamer_vid in self._shared_mtg.components_iter(mtg_axis_vid):
                    mtg_metamer_index = int(self._shared_mtg.index(mtg_metamer_vid))
                    for mtg_organ_vid in self._shared_mtg.components_iter(mtg_metamer_vid):
                        mtg_organ_label = self._shared_mtg.label(mtg_organ_vid)
                        if mtg_organ_label not in FARQUHARWHEAT_ORGANS_NAMES: continue
                        for mtg_element_vid in self._shared_mtg.components_iter(mtg_organ_vid):
                            mtg_element_properties = self._shared_mtg.get_vertex_property(mtg_element_vid)
                            mtg_element_label = self._shared_mtg.label(mtg_element_vid)
                            element_id = (mtg_plant_index, mtg_axis_label, mtg_metamer_index, mtg_organ_label, mtg_element_label)
                            if set(mtg_element_properties).issuperset(converter.FARQUHARWHEAT_INPUTS):
                                farquharwheat_element_inputs_dict = {}
                                for farquharwheat_element_input_name in converter.FARQUHARWHEAT_INPUTS:
                                    farquharwheat_element_inputs_dict[farquharwheat_element_input_name] = mtg_element_properties[farquharwheat_element_input_name]
                                all_farquharwheat_elements_inputs_dict[element_id] = farquharwheat_element_inputs_dict
                                # TODO: temporary ; replace 'FARQUHARWHEAT_ELEMENT_PROPERTIES_TEMP' by default values
                                FARQUHARWHEAT_ELEMENT_PROPERTIES_TEMP = {'width': 0, 'height': 0.6, 'PARa': 0}
                                farquharwheat_element_inputs_dict = {}
                                for farquharwheat_element_input_name in converter.FARQUHARWHEAT_INPUTS:
                                    if farquharwheat_element_input_name in mtg_element_properties:
                                        farquharwheat_element_inputs_dict[farquharwheat_element_input_name] = mtg_element_properties[farquharwheat_element_input_name]
                                    else:
                                        farquharwheat_element_inputs_dict[farquharwheat_element_input_name] = FARQUHARWHEAT_ELEMENT_PROPERTIES_TEMP[farquharwheat_element_input_name]
                                all_farquharwheat_elements_inputs_dict[element_id] = farquharwheat_element_inputs_dict

        self._simulation.initialize(all_farquharwheat_elements_inputs_dict)


    def _update_shared_MTG(self, farquharwheat_elements_data_dict):
        """
        Update the MTG shared between all models from the inputs or the outputs of the model.
        """
        # add the properties if needed
        mtg_property_names = self._shared_mtg.property_names()
        for farquharwheat_elements_data_name in converter.FARQUHARWHEAT_INPUTS_OUTPUTS:
            if farquharwheat_elements_data_name not in mtg_property_names:
                self._shared_mtg.add_property(farquharwheat_elements_data_name)

        # traverse the MTG recursively from top ...
        for mtg_plant_vid in self._shared_mtg.components_iter(self._shared_mtg.root):
            mtg_plant_index = int(self._shared_mtg.index(mtg_plant_vid))
            for mtg_axis_vid in self._shared_mtg.components_iter(mtg_plant_vid):
                mtg_axis_label = self._shared_mtg.label(mtg_axis_vid)
                for mtg_metamer_vid in self._shared_mtg.components_iter(mtg_axis_vid):
                    mtg_metamer_index = int(self._shared_mtg.index(mtg_metamer_vid))
                    for mtg_organ_vid in self._shared_mtg.components_iter(mtg_metamer_vid):
                        mtg_organ_label = self._shared_mtg.label(mtg_organ_vid)
                        if mtg_organ_label not in FARQUHARWHEAT_ORGANS_NAMES: continue
                        for mtg_element_vid in self._shared_mtg.components_iter(mtg_organ_vid):
                            mtg_element_label = self._shared_mtg.label(mtg_element_vid)
                            element_id = (mtg_plant_index, mtg_axis_label, mtg_metamer_index, mtg_organ_label, mtg_element_label)
                            if element_id not in farquharwheat_elements_data_dict: continue
                            # update the element in the MTG
                            farquharwheat_element_data_dict = farquharwheat_elements_data_dict[element_id]
                            for farquharwheat_element_data_name, farquharwheat_element_data_value in farquharwheat_element_data_dict.iteritems():
                                self._shared_mtg.property(farquharwheat_element_data_name)[mtg_element_vid] = farquharwheat_element_data_value


    def _update_shared_dataframes(self, farquharwheat_elements_data_df):
        """
        Update the dataframes shared between all models from the inputs dataframes or the outputs dataframes of the model.
        """
        tools.combine_dataframes_inplace(farquharwheat_elements_data_df, SHARED_ELEMENTS_INPUTS_OUTPUTS_INDEXES, self._shared_elements_inputs_outputs_df)

