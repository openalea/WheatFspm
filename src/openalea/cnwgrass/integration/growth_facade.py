# -*- coding: latin-1 -*-

from openalea.cnwgrass.growth import converter, simulation, parameters
from openalea.cnwgrass.integration import tools

"""
    integration.growth_facade
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`integration.growth_facade` is a facade of the model Growth.

    This module permits to initialize and run the model Growth from a :class:`MTG <openalea.mtg.mtg.MTG>`
    in a convenient and transparent way, wrapping all the internal complexity of the model, and dealing
    with all the tedious initialization and conversion processes.

"""

LEAF_LABELS = ['blade', 'sheath']

EMERGED_GROWING_ORGAN_LABELS = ['StemElement', 'LeafElement1']
ELEMENT_LABELS = ['StemElement', 'LeafElement1', 'HiddenElement']

SHARED_AXES_INPUTS_OUTPUTS_INDEXES = ['plant', 'axis']

SHARED_ORGANS_INPUTS_OUTPUTS_INDEXES = ['plant', 'axis', 'organ']

SHARED_HIDDENZONES_INPUTS_OUTPUTS_INDEXES = ['plant', 'axis', 'metamer']

SHARED_ELEMENTS_INPUTS_OUTPUTS_INDEXES = ['plant', 'axis', 'metamer', 'organ', 'element']


class GrowthFacade(object):
    """
    The GrowthFacade class permits to initialize, run the model Growth
    from a :class:`MTG <openalea.mtg.mtg.MTG>`, and update the MTG and the dataframes
    shared between all models.

    Use :meth:`run` to run the model.

    :Parameters:

"""

    #: hiddenzone-scale inputs declared as required by the Growth model but not unconditionally needed:
    #: 'leaf_Lmax'/'internode_Lmax' are never actually read by any computation in growth/simulation.py;
    #: 'is_over' is only ever written by Growth, never read as an input; 'delta_internode_pseudo_age' /
    #: 'delta_leaf_pseudo_age' are declared in morphogenesis's own `hiddenzone_outputs` (written) but NOT in its
    #: `hiddenzone_inputs` (read from the MTG) -- they only exist in the MTG once morphogenesis's real,
    #: MS-only computation has freshly added them that same time step (internode past ligulation / leaf past
    #: emergence), which is also exactly the only case Growth itself ever reads them (see the
    #: `internode_pseudo_age >= internode_rapid_growth_t` / `leaf_is_growing` gates in growth/simulation.py) --
    #: so requiring them here before that state is ever reached would permanently exclude the hiddenzone.
    #: See the note where this is used in `_initialize_model`.
    HIDDENZONE_INPUTS_NOT_REQUIRED = {'leaf_Lmax', 'internode_Lmax', 'is_over', 'delta_internode_pseudo_age', 'delta_leaf_pseudo_age'}

    def __init__(self, shared_mtg, delta_t,
                 model_hiddenzones_inputs_df,
                 model_elements_inputs_df,
                 model_roots_inputs_df,
                 model_axes_inputs_df,
                 shared_organs_inputs_outputs_df,
                 shared_hiddenzones_inputs_outputs_df,
                 shared_elements_inputs_outputs_df,
                 shared_axes_inputs_outputs_df,
                 hydraulics=False,
                 explicit_tillers=False,
                 update_parameters=None,
                 update_shared_df=True,
                 cnwgrass_roots=True):

        """
        :param openalea.mtg.mtg.MTG shared_mtg: The MTG shared between all models.
        :param int delta_t: The delta between two runs, in seconds.
        :param pandas.DataFrame model_hiddenzones_inputs_df: the inputs of the model at hiddenzones scale.
        :param pandas.DataFrame model_elements_inputs_df: the inputs of the model at elements scale.
        :param pandas.DataFrame model_roots_inputs_df: the inputs of the model at roots scale.
        :param pandas.DataFrame model_axes_inputs_df: the inputs of the model at axes scale.
        :param pandas.DataFrame shared_organs_inputs_outputs_df: the dataframe of inputs and outputs at organs scale shared between all models.
        :param pandas.DataFrame shared_hiddenzones_inputs_outputs_df: the dataframe of inputs and outputs at hiddenzones scale shared between all models.
        :param pandas.DataFrame shared_elements_inputs_outputs_df: the dataframe of inputs and outputs at elements scale shared between all models.
        :param pandas.DataFrame shared_axes_inputs_outputs_df: the dataframe of inputs and outputs at axis scale shared between all models.
        :param bool hydraulics: if True the model will assume the coupling to the turgor-driven growth model
        :param bool explicit_tillers: if True, compute real dimensional/mass growth for every tiller's hiddenzones/elements found in the MTG, instead of only the main stem's.
        :param dict or None update_parameters: A dictionary with the parameters to update, should have the form {'param1': value1, 'param2': value2, ...}.
        :param bool update_shared_df: If `True`  update the shared dataframes at init and at each run (unless stated otherwise)
        """

        self._shared_mtg = shared_mtg  #: the MTG shared between all models
        self.explicit_tillers = explicit_tillers

        self._simulation = simulation.Simulation(delta_t=delta_t, hydraulics=hydraulics, update_parameters=update_parameters, cnwgrass_roots=cnwgrass_roots)  #: the simulator to use to run the model

        all_growth_inputs_dict = converter.from_dataframes(model_hiddenzones_inputs_df, model_elements_inputs_df, model_roots_inputs_df, model_axes_inputs_df)

        self._update_shared_MTG(all_growth_inputs_dict['hiddenzone'], all_growth_inputs_dict['elements'], all_growth_inputs_dict['roots'], all_growth_inputs_dict['axes'])

        self._shared_organs_inputs_outputs_df = shared_organs_inputs_outputs_df  #: the dataframe at organs scale shared between all models
        self._shared_hiddenzones_inputs_outputs_df = shared_hiddenzones_inputs_outputs_df  #: the dataframe at hiddenzones scale shared between all models
        self._shared_elements_inputs_outputs_df = shared_elements_inputs_outputs_df  #: the dataframe at elements scale shared between all models
        self._shared_axes_inputs_outputs_df = shared_axes_inputs_outputs_df  #: the dataframe at axis scale shared between all models
        self._update_shared_df = update_shared_df
        if self._update_shared_df:
            self._update_shared_dataframes(model_hiddenzones_inputs_df, model_elements_inputs_df, model_roots_inputs_df, model_axes_inputs_df)

        self.seed_is_moistened = True

    def run(self, postflowering_stages=False, update_shared_df=None):
        """
        Run the model and update the MTG and the dataframes shared between all models.
        :param bool postflowering_stages: if True the model will calculate root growth with the parameters calibrated for post flowering stages
        :param bool update_shared_df: if 'True', update the shared dataframes at this time step.
        """
        self._initialize_model()
        if self.seed_is_moistened:
            self._simulation.run(postflowering_stages)
            self._update_shared_MTG(self._simulation.outputs['hiddenzone'], self._simulation.outputs['elements'], self._simulation.outputs['roots'], self._simulation.outputs['axes'])

            if update_shared_df or (update_shared_df is None and self._update_shared_df):
                growth_hiddenzones_outputs_df, growth_elements_outputs_df, growth_roots_outputs_df, growth_axes_outputs_df = converter.to_dataframes(self._simulation.outputs,
                                                                                                                                                     self._simulation.axis_outputs,
                                                                                                                                                     self._simulation.hiddenzone_outputs,
                                                                                                                                                     self._simulation.element_outputs,
                                                                                                                                                     self._simulation.root_outputs)
                self._update_shared_dataframes(growth_hiddenzones_outputs_df, growth_elements_outputs_df, growth_roots_outputs_df, growth_axes_outputs_df)

    def _initialize_model(self):
        """
        Initialize the inputs of the model from the MTG shared between all models.
        """

        all_growth_hiddenzones_inputs_dict = {}
        all_growth_elements_inputs_dict = {}
        all_growth_roots_inputs_dict = {}
        all_growth_axes_inputs_dict = {}

        for mtg_plant_vid in self._shared_mtg.components_iter(self._shared_mtg.root):
            mtg_plant_index = int(self._shared_mtg.index(mtg_plant_vid))
            for mtg_axis_vid in self._shared_mtg.components_iter(mtg_plant_vid):
                mtg_axis_label = self._shared_mtg.label(mtg_axis_vid)
                if isinstance(mtg_axis_label, bytes):
                    #: Some MTGs carry byte-string axis labels (numpy.bytes_ included) -- normalize to str so
                    #: 'MS' comparisons behave consistently, and so hiddenzone_id/element_id tuples built from it
                    #: don't mix str and bytes axis labels within the same `sorted(...)` call in growth/simulation.py
                    #: (str/bytes are not orderable against each other and sorted() raises TypeError otherwise).
                    mtg_axis_label = mtg_axis_label.decode('UTF-8')
                if mtg_axis_label != 'MS' and not self.explicit_tillers:
                    continue

                mtg_axis_properties = self._shared_mtg.get_vertex_property(mtg_axis_vid)

                if 'endosperm' in mtg_axis_properties and mtg_axis_properties['endosperm']['moistening'] < 1:
                    self.seed_is_moistened = False
                    continue
                elif mtg_axis_label =='MS':
                    self.seed_is_moistened = True

                axis_id = (mtg_plant_index, mtg_axis_label)
                if set(mtg_axis_properties).issuperset(self._simulation.axis_inputs):
                    growth_axis_inputs_dict = {}
                    for growth_axis_input_name in self._simulation.axis_inputs:
                        growth_axis_inputs_dict[growth_axis_input_name] = mtg_axis_properties[growth_axis_input_name]
                    all_growth_axes_inputs_dict[axis_id] = growth_axis_inputs_dict

                # Roots
                if 'roots' in mtg_axis_properties:
                    roots_id = (mtg_plant_index, mtg_axis_label, 'roots')
                    mtg_roots_properties = mtg_axis_properties['roots']

                    if set(mtg_roots_properties).issuperset(self._simulation.root_inputs):
                        growth_roots_inputs_dict = {}
                        for growth_roots_input_name in self._simulation.root_inputs:
                            growth_roots_inputs_dict[growth_roots_input_name] = mtg_roots_properties[growth_roots_input_name]
                        all_growth_roots_inputs_dict[roots_id] = growth_roots_inputs_dict

                for mtg_metamer_vid in self._shared_mtg.components_iter(mtg_axis_vid):

                    mtg_metamer_index = int(self._shared_mtg.index(mtg_metamer_vid))

                    mtg_metamer_properties = self._shared_mtg.get_vertex_property(mtg_metamer_vid)
                    if 'hiddenzone' in mtg_metamer_properties:
                        hiddenzone_id = (mtg_plant_index, mtg_axis_label, mtg_metamer_index)
                        mtg_hiddenzone_properties = mtg_metamer_properties['hiddenzone']

                        #: 'leaf_Lmax' / 'internode_Lmax' are declared as required hiddenzone inputs but never
                        #: actually read by any computation in growth/simulation.py, and 'is_over' is Growth's
                        #: own output (only ever written, never read as an input) -- so none of the three need to
                        #: already be present in the MTG. This matters for a brand-new tiller hiddenzone: under
                        #: hydraulics=True, morphogenesis itself stops tracking 'leaf_Lmax' (see
                        #: morphogenesis/model.py's `variables_to_delete`), so a tiller can never inherit it from
                        #: MS via the copy mechanism, and 'is_over' is never set until Growth itself first
                        #: processes the hiddenzone -- requiring them here would permanently exclude every tiller.
                        required_hiddenzone_inputs = set(self._simulation.hiddenzone_inputs) - GrowthFacade.HIDDENZONE_INPUTS_NOT_REQUIRED
                        if set(mtg_hiddenzone_properties).issuperset(required_hiddenzone_inputs):  # Initial values are set by morphogenesis
                            growth_hiddenzone_inputs_dict = {}
                            for growth_hiddenzone_input_name in self._simulation.hiddenzone_inputs:
                                #: default to 0.0 (not None) when absent: 'delta_internode_pseudo_age' / 'delta_leaf_pseudo_age'
                                #: can genuinely be read by a numeric computation once the corresponding hiddenzone
                                #: reaches the relevant developmental stage (see HIDDENZONE_INPUTS_NOT_REQUIRED above);
                                #: 0.0 ("no progress this step") is a safe fallback there, whereas None would raise.
                                growth_hiddenzone_inputs_dict[growth_hiddenzone_input_name] = mtg_hiddenzone_properties.get(growth_hiddenzone_input_name, 0.0)
                            all_growth_hiddenzones_inputs_dict[hiddenzone_id] = growth_hiddenzone_inputs_dict

                        # We take only the elements of growing metamers i.e. the ones with hiddenzones
                        for mtg_organ_vid in self._shared_mtg.components_iter(mtg_metamer_vid):
                            mtg_organ_label = self._shared_mtg.label(mtg_organ_vid)

                            for mtg_element_vid in self._shared_mtg.components_iter(mtg_organ_vid):
                                mtg_element_label = self._shared_mtg.label(mtg_element_vid)
                                element_id = (mtg_plant_index, mtg_axis_label, mtg_metamer_index, mtg_organ_label, mtg_element_label)
                                mtg_element_properties = self._shared_mtg.get_vertex_property(mtg_element_vid)

                                if mtg_element_label in ELEMENT_LABELS and \
                                        mtg_element_properties.get('length', 0) > 0:  # Note : ADEL puts length to positive value after updates even for HiddenElement.

                                    growth_element_inputs_dict = {}

                                    # Exclude the HiddenElement apart from remobilization cases
                                    remobilisation = False
                                    if mtg_element_label == 'HiddenElement':
                                        if growth_hiddenzone_inputs_dict['leaf_is_remobilizing'] or growth_hiddenzone_inputs_dict['internode_is_remobilizing']:
                                            remobilisation = True
                                        else:
                                            continue

                                    for growth_element_input_name in self._simulation.element_inputs:
                                        mtg_element_input = mtg_element_properties.get(growth_element_input_name)
                                        if mtg_element_input is None:
                                            mtg_element_input = parameters.OrganInit().__dict__[growth_element_input_name]
                                        growth_element_inputs_dict[growth_element_input_name] = mtg_element_input
                                        if remobilisation:
                                            # Needed later on for CN-Metabolism calculation. TODO: Should it be in morphogenesis_facade instead ? (MG)
                                            growth_element_inputs_dict['green_area'] = mtg_element_properties.get('area')
                                    all_growth_elements_inputs_dict[element_id] = growth_element_inputs_dict

        self._simulation.initialize({'hiddenzone': all_growth_hiddenzones_inputs_dict, 'elements': all_growth_elements_inputs_dict,
                                     'roots': all_growth_roots_inputs_dict, 'axes': all_growth_axes_inputs_dict})

    def _update_shared_MTG(self, all_growth_hiddenzones_data_dict, all_growth_elements_data_dict, all_growth_roots_data_dict, all_growth_axes_data_dict):
        """
        Update the MTG shared between all models from the inputs or the outputs of the model.

        :param dict all_growth_hiddenzones_data_dict: Growth outputs at hidden zone scale
        :param dict all_growth_elements_data_dict: Growth outputs at element scale
        :param dict all_growth_roots_data_dict: Growth outputs at root scale
        :param dict all_growth_axes_data_dict: Growth outputs at axis scale
        """

        # add the properties if needed
        mtg_property_names = self._shared_mtg.property_names()
        if 'roots' not in mtg_property_names:
            self._shared_mtg.add_property('roots')
        for growth_data_name in set(self._simulation.hiddenzone_inputs_outputs + self._simulation.element_inputs_outputs):
            if growth_data_name not in mtg_property_names:
                self._shared_mtg.add_property(growth_data_name)

        # update the properties of the MTG
        for mtg_plant_vid in self._shared_mtg.components_iter(self._shared_mtg.root):
            mtg_plant_index = int(self._shared_mtg.index(mtg_plant_vid))
            for mtg_axis_vid in self._shared_mtg.components_iter(mtg_plant_vid):
                mtg_axis_label = self._shared_mtg.label(mtg_axis_vid)
                if isinstance(mtg_axis_label, bytes):
                    mtg_axis_label = mtg_axis_label.decode('UTF-8')
                axis_id = (mtg_plant_index, mtg_axis_label)

                if mtg_axis_label != 'MS' and not self.explicit_tillers:
                    continue

                if axis_id not in all_growth_axes_data_dict:
                    #: A tiller freshly created by morphogenesis this same time step, or otherwise not yet
                    #: carrying the axis-scale inputs Growth requires, was excluded from _initialize_model's
                    #: inputs (see the `set(...).issuperset(...)` check above) and so has no output either --
                    #: skip it for this time step rather than crashing on a missing dict entry.
                    continue

                growth_axis_data_dict = all_growth_axes_data_dict[axis_id]
                for axis_data_name, axis_data_value in growth_axis_data_dict.items():
                    self._shared_mtg.property(axis_data_name)[mtg_axis_vid] = axis_data_value

                roots_id = (mtg_plant_index, mtg_axis_label, 'roots')
                if roots_id in all_growth_roots_data_dict:
                    growth_roots_data_dict = all_growth_roots_data_dict[roots_id]
                    mtg_axis_properties = self._shared_mtg.get_vertex_property(mtg_axis_vid)
                    if 'roots' not in mtg_axis_properties:
                        self._shared_mtg.property('roots')[mtg_axis_vid] = {}
                    for roots_data_name, roots_data_value in growth_roots_data_dict.items():
                        self._shared_mtg.property('roots')[mtg_axis_vid][roots_data_name] = roots_data_value

                #: Metamer scale
                for mtg_metamer_vid in self._shared_mtg.components_iter(mtg_axis_vid):
                    mtg_metamer_index = int(self._shared_mtg.index(mtg_metamer_vid))
                    hiddenzone_id = (mtg_plant_index, mtg_axis_label, mtg_metamer_index)
                    if hiddenzone_id in all_growth_hiddenzones_data_dict:
                        growth_hiddenzone_data_dict = all_growth_hiddenzones_data_dict[hiddenzone_id]
                        mtg_metamer_properties = self._shared_mtg.get_vertex_property(mtg_metamer_vid)
                        if 'hiddenzone' not in mtg_metamer_properties:  # MG : when is it used ?
                            self._shared_mtg.property('hiddenzone')[mtg_metamer_vid] = {}
                        for hiddenzone_data_name, hiddenzone_data_value in growth_hiddenzone_data_dict.items():
                            self._shared_mtg.property('hiddenzone')[mtg_metamer_vid][hiddenzone_data_name] = hiddenzone_data_value

                    elif mtg_axis_label == 'MS' and 'hiddenzone' in self._shared_mtg.get_vertex_property(mtg_metamer_vid):
                        #: remove the 'hiddenzone' property from this metamer -- only valid for the main stem,
                        #: whose hiddenzones are self-sustained by Growth's own writes every time step (see the
                        #: `issuperset` check in _initialize_model): there, hiddenzone_id being absent from
                        #: Growth's outputs really does mean this leaf has fully emerged and its hiddenzone is
                        #: legitimately gone. A tiller's hiddenzone is never self-sustained that way (it is
                        #: seeded from morphogenesis's copy-from-MS mechanism, not from Growth), so its
                        #: hiddenzone_id can be absent from Growth's outputs simply because Growth doesn't have
                        #: all the inputs it wants for it yet -- deleting it in that case would silently wipe
                        #: the tiller's own geometry, not reflect a real developmental event.
                        del self._shared_mtg.property('hiddenzone')[mtg_metamer_vid]

                    #: Organ scale
                    for mtg_organ_vid in self._shared_mtg.components_iter(mtg_metamer_vid):
                        mtg_organ_label = self._shared_mtg.label(mtg_organ_vid)

                        #: Element scale
                        for mtg_element_vid in self._shared_mtg.components_iter(mtg_organ_vid):
                            mtg_element_label = self._shared_mtg.label(mtg_element_vid)
                            element_id = (mtg_plant_index, mtg_axis_label, mtg_metamer_index, mtg_organ_label, mtg_element_label)

                            if element_id in all_growth_elements_data_dict:
                                growth_element_data_dict = all_growth_elements_data_dict[element_id]

                                for element_data_name, element_data_value in growth_element_data_dict.items():
                                    self._shared_mtg.property(element_data_name)[mtg_element_vid] = element_data_value

    def _update_shared_dataframes(self, growth_hiddenzones_data_df, growth_elements_data_df, growth_roots_data_df, growth_axes_data_df):
        """
        Update the dataframes shared between all models from the inputs dataframes or the outputs dataframes of the model.

        :param pandas.DataFrame growth_hiddenzones_data_df: Growth shared dataframe at hidden zone scale
        :param pandas.DataFrame growth_elements_data_df: Growth shared dataframe at element scale
        :param pandas.DataFrame growth_roots_data_df: Growth shared dataframe at roots scale
        :param pandas.DataFrame growth_axes_data_df: Growth shared dataframe at axis scale
        """

        for growth_data_df, \
            shared_inputs_outputs_indexes, \
            shared_inputs_outputs_df in ((growth_hiddenzones_data_df, SHARED_HIDDENZONES_INPUTS_OUTPUTS_INDEXES, self._shared_hiddenzones_inputs_outputs_df),
                                         (growth_elements_data_df, SHARED_ELEMENTS_INPUTS_OUTPUTS_INDEXES, self._shared_elements_inputs_outputs_df),
                                         (growth_roots_data_df, SHARED_ORGANS_INPUTS_OUTPUTS_INDEXES, self._shared_organs_inputs_outputs_df),
                                         (growth_axes_data_df, SHARED_AXES_INPUTS_OUTPUTS_INDEXES, self._shared_axes_inputs_outputs_df)):

            if growth_data_df is growth_roots_data_df:
                growth_data_df = growth_data_df.copy()
                growth_data_df.loc[:, 'organ'] = 'roots'

            tools.combine_dataframes_inplace(growth_data_df, shared_inputs_outputs_indexes, shared_inputs_outputs_df)
