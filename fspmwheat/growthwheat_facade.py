# -*- coding: latin-1 -*-

"""
    fspmwheat.growthwheat_facade
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`fspmwheat.growthwheat_facade` is a facade of the model GrowthWheat.

    This module permits to initialize and run the model GrowthWheat from a :class:`MTG <openalea.mtg.mtg.MTG>`
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

from growthwheat import converter, simulation
import tools

EMERGED_GROWING_ORGAN_LABELS = ['StemElement', 'LeafElement1']

SHARED_ORGANS_INPUTS_OUTPUTS_INDEXES = ['plant', 'axis', 'organ']

SHARED_HIDDENZONES_INPUTS_OUTPUTS_INDEXES = ['plant', 'axis', 'metamer']

SHARED_ELEMENTS_INPUTS_OUTPUTS_INDEXES = ['plant', 'axis', 'metamer', 'organ', 'element']

class GrowthWheatFacade(object):
    """
    The GrowthWheatFacade class permits to initialize, run the model GrowthWheat
    from a :class:`MTG <openalea.mtg.mtg.MTG>`, and update the MTG and the dataframes
    shared between all models.

    Use :meth:`run` to run the model.

    :Parameters:

        - `shared_mtg` (:class:`openalea.mtg.mtg.MTG`) - The MTG shared between all models.
        - `delta_t` (:class:`int`) - The delta between two runs, in seconds.
        - `model_hiddenzones_inputs_df` (:class:`pandas.DataFrame`) - the inputs of the model at hiddenzones scale.
        - `model_organs_inputs_df` (:class:`pandas.DataFrame`) - the inputs of the model at organs scale.
        - `model_roots_inputs_df` (:class:`pandas.DataFrame`) - the inputs of the model at roots scale.
        - `shared_organs_inputs_outputs_df` (:class:`pandas.DataFrame`) - the dataframe of inputs and outputs at organs scale shared between all models.
        - `shared_hiddenzones_inputs_outputs_df` (:class:`pandas.DataFrame`) - the dataframe of inputs and outputs at hiddenzones scale shared between all models.
        - `shared_elements_inputs_outputs_df` (:class:`pandas.DataFrame`) - the dataframe of inputs and outputs at elements scale shared between all models.
    """

    def __init__(self, shared_mtg, delta_t,
                 model_hiddenzones_inputs_df,
                 model_organs_inputs_df,
                 model_roots_inputs_df,
                 shared_organs_inputs_outputs_df,
                 shared_hiddenzones_inputs_outputs_df,
                 shared_elements_inputs_outputs_df):

        self._shared_mtg = shared_mtg #: the MTG shared between all models

        self._simulation = simulation.Simulation(delta_t=delta_t) #: the simulator to use to run the model

        all_growthwheat_inputs_dict = converter.from_dataframes(model_hiddenzones_inputs_df, model_organs_inputs_df, model_roots_inputs_df)

        self._update_shared_MTG(all_growthwheat_inputs_dict['hiddenzone'], all_growthwheat_inputs_dict['organs'], all_growthwheat_inputs_dict['roots'])

        self._shared_organs_inputs_outputs_df = shared_organs_inputs_outputs_df #: the dataframe at organs scale shared between all models
        self._shared_hiddenzones_inputs_outputs_df = shared_hiddenzones_inputs_outputs_df #: the dataframe at hiddenzones scale shared between all models
        self._shared_elements_inputs_outputs_df = shared_elements_inputs_outputs_df #: the dataframe at elements scale shared between all models
        self._update_shared_dataframes(model_hiddenzones_inputs_df, model_organs_inputs_df, model_roots_inputs_df)


    def run(self):
        """
        Run the model and update the MTG and the dataframes shared between all models.
        """
        self._initialize_model()
        self._simulation.run()
        self._update_shared_MTG(self._simulation.outputs['hiddenzone'], self._simulation.outputs['organs'], self._simulation.outputs['roots'])
        growthwheat_hiddenzones_outputs_df, growthwheat_organs_outputs_df, growthwheat_roots_outputs_df = converter.to_dataframes(self._simulation.outputs)
        self._update_shared_dataframes(growthwheat_hiddenzones_outputs_df, growthwheat_organs_outputs_df, growthwheat_roots_outputs_df)


    def _initialize_model(self):
        """
        Initialize the inputs of the model from the MTG shared between all models.
        """

        all_growthwheat_hiddenzones_inputs_dict = {}
        all_growthwheat_organs_inputs_dict = {}
        all_growthwheat_roots_inputs_dict = {}

        for mtg_plant_vid in self._shared_mtg.components_iter(self._shared_mtg.root):
            mtg_plant_index = int(self._shared_mtg.index(mtg_plant_vid))
            for mtg_axis_vid in self._shared_mtg.components_iter(mtg_plant_vid):
                mtg_axis_label = self._shared_mtg.label(mtg_axis_vid)

                # Roots
                mtg_axis_properties = self._shared_mtg.get_vertex_property(mtg_axis_vid)
                if 'roots' in mtg_axis_properties:
                    roots_id = (mtg_plant_index, mtg_axis_label, 'roots')
                    mtg_roots_properties = mtg_axis_properties['roots']

                    if set(mtg_roots_properties).issuperset(simulation.ROOT_INPUTS):
                        growthwheat_roots_inputs_dict = {}
                        for growthwheat_roots_input_name in simulation.ROOT_INPUTS:
                            growthwheat_roots_inputs_dict[growthwheat_roots_input_name] = mtg_roots_properties[growthwheat_roots_input_name]
                        all_growthwheat_roots_inputs_dict[roots_id] = growthwheat_roots_inputs_dict

                for mtg_metamer_vid in self._shared_mtg.components_iter(mtg_axis_vid):

                    mtg_metamer_index = int(self._shared_mtg.index(mtg_metamer_vid))
                    for mtg_organ_vid in self._shared_mtg.components_iter(mtg_metamer_vid):
                        mtg_organ_label = self._shared_mtg.label(mtg_organ_vid)
                        mtg_organ_properties = self._shared_mtg.get_vertex_property(mtg_organ_vid)
                        organ_id = (mtg_plant_index, mtg_axis_label, mtg_metamer_index, mtg_organ_label)
                        #: TODO: a changer car ici je recopie les prop de l'élement dans le dict des organes
                        for mtg_element_vid in self._shared_mtg.components_iter(mtg_organ_vid):
                            mtg_element_properties = self._shared_mtg.get_vertex_property(mtg_element_vid)
                            if set(mtg_element_properties).issuperset(simulation.ORGAN_INPUTS) and mtg_element_properties['label'] in EMERGED_GROWING_ORGAN_LABELS:
                                growthwheat_organ_inputs_dict = {}
                                for growthwheat_organ_input_name in simulation.ORGAN_INPUTS:
                                    growthwheat_organ_inputs_dict[growthwheat_organ_input_name] = mtg_element_properties[growthwheat_organ_input_name]
                                all_growthwheat_organs_inputs_dict[organ_id] = growthwheat_organ_inputs_dict

                    mtg_metamer_properties = self._shared_mtg.get_vertex_property(mtg_metamer_vid)
                    if 'hiddenzone' in mtg_metamer_properties:
                        hiddenzone_id = (mtg_plant_index, mtg_axis_label, mtg_metamer_index)
                        mtg_hiddenzone_properties = mtg_metamer_properties['hiddenzone']

                        if set(mtg_hiddenzone_properties).issuperset(simulation.HIDDENZONE_INPUTS):
                            growthwheat_hiddenzone_inputs_dict = {}
                            for growthwheat_hiddenzone_input_name in simulation.HIDDENZONE_INPUTS:
                                growthwheat_hiddenzone_inputs_dict[growthwheat_hiddenzone_input_name] = mtg_hiddenzone_properties[growthwheat_hiddenzone_input_name]
                            all_growthwheat_hiddenzones_inputs_dict[hiddenzone_id] = growthwheat_hiddenzone_inputs_dict

        self._simulation.initialize({'hiddenzone': all_growthwheat_hiddenzones_inputs_dict, 'organs': all_growthwheat_organs_inputs_dict, 'roots': all_growthwheat_roots_inputs_dict})


    def _update_shared_MTG(self, all_growthwheat_hiddenzones_data_dict, all_growthwheat_organs_data_dict, all_growthwheat_roots_data_dict):
        """
        Update the MTG shared between all models from the inputs or the outputs of the model.
        """
        # add the properties if needed
        mtg_property_names = self._shared_mtg.property_names()
        if 'roots' not in mtg_property_names:
            self._shared_mtg.add_property('roots')
        for growthwheat_data_name in set(simulation.HIDDENZONE_INPUTS_OUTPUTS + simulation.ORGAN_INPUTS_OUTPUTS):
            if growthwheat_data_name not in mtg_property_names:
                self._shared_mtg.add_property(growthwheat_data_name)

        # update the properties of the MTG
        for mtg_plant_vid in self._shared_mtg.components_iter(self._shared_mtg.root):
            mtg_plant_index = int(self._shared_mtg.index(mtg_plant_vid))
            for mtg_axis_vid in self._shared_mtg.components_iter(mtg_plant_vid):
                mtg_axis_label = self._shared_mtg.label(mtg_axis_vid)
                roots_id = (mtg_plant_index, mtg_axis_label, 'roots')
                if roots_id in all_growthwheat_roots_data_dict:
                    growthwheat_roots_data_dict = all_growthwheat_roots_data_dict[roots_id]
                    mtg_axis_properties = self._shared_mtg.get_vertex_property(mtg_axis_vid)
                    if 'roots' not in mtg_axis_properties:
                        self._shared_mtg.property('roots')[mtg_axis_vid] = {}
                    for roots_data_name, roots_data_value in growthwheat_roots_data_dict.iteritems():
                        self._shared_mtg.property('roots')[mtg_axis_vid][roots_data_name] = roots_data_value

                for mtg_metamer_vid in self._shared_mtg.components_iter(mtg_axis_vid):
                    mtg_metamer_index = int(self._shared_mtg.index(mtg_metamer_vid))
                    hiddenzone_id = (mtg_plant_index, mtg_axis_label, mtg_metamer_index)
                    if hiddenzone_id in all_growthwheat_hiddenzones_data_dict:
                        growthwheat_hiddenzone_data_dict = all_growthwheat_hiddenzones_data_dict[hiddenzone_id]
                        mtg_metamer_properties = self._shared_mtg.get_vertex_property(mtg_metamer_vid)
                        if 'hiddenzone' not in mtg_metamer_properties:
                            self._shared_mtg.property('hiddenzone')[mtg_metamer_vid] = {}
                        for hiddenzone_data_name, hiddenzone_data_value in growthwheat_hiddenzone_data_dict.iteritems():
                            self._shared_mtg.property('hiddenzone')[mtg_metamer_vid][hiddenzone_data_name] = hiddenzone_data_value

                    for mtg_organ_vid in self._shared_mtg.components_iter(mtg_metamer_vid):
                        mtg_organ_label = self._shared_mtg.label(mtg_organ_vid)
                        organ_id = (mtg_plant_index, mtg_axis_label, mtg_metamer_index, mtg_organ_label)

                        if organ_id in all_growthwheat_organs_data_dict:
                            growthwheat_organ_data_dict = all_growthwheat_organs_data_dict[organ_id]

                            for organ_data_name, organ_data_value in growthwheat_organ_data_dict.iteritems():
                                self._shared_mtg.property(organ_data_name)[mtg_organ_vid] = organ_data_value

                                # Write properties in element scale #: TODO: voir ce qui ce passe avec +sieurs élements
                                for mtg_element_vid in self._shared_mtg.components_iter(mtg_organ_vid):
                                    if self._shared_mtg.get_vertex_property(mtg_element_vid)['label'] in ('StemElement', 'LeafElement1'):
                                        self._shared_mtg.property(organ_data_name)[mtg_element_vid] = growthwheat_organ_data_dict.get(organ_data_name)


    def _update_shared_dataframes(self, growthwheat_hiddenzones_data_df, growthwheat_organs_data_df, growthwheat_roots_data_df):
        """
        Update the dataframes shared between all models from the inputs dataframes or the outputs dataframes of the model.
        """

        for growthwheat_data_df, \
            shared_inputs_outputs_indexes, \
            shared_inputs_outputs_df in ((growthwheat_hiddenzones_data_df, SHARED_HIDDENZONES_INPUTS_OUTPUTS_INDEXES, self._shared_hiddenzones_inputs_outputs_df),
                                         (growthwheat_organs_data_df, SHARED_ELEMENTS_INPUTS_OUTPUTS_INDEXES, self._shared_elements_inputs_outputs_df),
                                         (growthwheat_roots_data_df, SHARED_ORGANS_INPUTS_OUTPUTS_INDEXES, self._shared_organs_inputs_outputs_df)):

            if growthwheat_data_df is growthwheat_roots_data_df:
                growthwheat_data_df = growthwheat_data_df.copy()
                growthwheat_data_df.loc[:, 'organ'] = 'roots'
            elif growthwheat_data_df is growthwheat_organs_data_df:
                growthwheat_data_df = growthwheat_data_df.copy()
                growthwheat_data_df['element'] = "" # TODO: check with Camille
                growthwheat_data_df.loc[growthwheat_data_df.organ == 'blade', 'element'] = 'LeafElement1'
                growthwheat_data_df.loc[growthwheat_data_df.organ != 'blade', 'element'] = 'StemElement'

            tools.combine_dataframes_inplace(growthwheat_data_df, shared_inputs_outputs_indexes, shared_inputs_outputs_df)
