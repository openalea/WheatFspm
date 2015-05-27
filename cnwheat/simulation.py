# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    cnwheat.simulation
    ~~~~~~~~~~~~~~~~~~

    The module :mod:`cnwheat.simulation` is the front-end to run the CN-Wheat model.

    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2014.
"""

"""
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
"""

import sys
import logging

import numpy as np
import pandas as pd
from scipy.integrate import odeint
from respiwheat.model import RespirationModel

import model

class SimulationError(Exception): pass
class SimulationRunError(SimulationError): pass
class SimulationInputsError(SimulationError): pass

class Simulation(object):
    """
    The Simulation class permits to initialize and run the model.

    Use :meth:`run` to run the model.

    :Parameters:

        - `population` (:class:`list`) - the :class:`population <cnwheat.model.Population` in the system.

    """

    #: the name of the compartments attributes in the model.
    MODEL_COMPARTMENTS_NAMES = {model.Plant: [],
                                model.Axis: [],
                                model.Phytomer: [],
                                model.Organ: ['sucrose', 'amino_acids', 'nitrates', 'structure', 'starch', 'proteins', 'mstruct', 'Nstruct'],
                                model.PhotosyntheticOrganElement: ['nitrates', 'starch', 'amino_acids', 'proteins',
                                                                   'sucrose', 'triosesP', 'fructan', 'mstruct', 'Nstruct']}

    T_INDEX = 't'
    
    PLANTS_INPUTS_INDEXES = ['plant']
    PLANTS_OUTPUTS_INDEXES = [T_INDEX] + PLANTS_INPUTS_INDEXES
    PLANTS_INTERNAL_STATE_PARAMETERS = []
    PLANTS_EXTERNAL_STATE_PARAMETERS = []
    PLANTS_STATE_PARAMETERS = PLANTS_INTERNAL_STATE_PARAMETERS + PLANTS_EXTERNAL_STATE_PARAMETERS
    PLANTS_STATE = PLANTS_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.Plant, [])
    PLANTS_INTERMEDIATE_VARIABLES = []
    PLANTS_FLUXES = []
    PLANTS_INTEGRATIVE_VARIABLES = []
    PLANTS_RUN_VARIABLES = PLANTS_OUTPUTS_INDEXES + PLANTS_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.Plant, []) + PLANTS_INTERMEDIATE_VARIABLES + PLANTS_FLUXES + PLANTS_INTEGRATIVE_VARIABLES
    PLANTS_POSTPROCESSING_VARIABLES = []
    PLANTS_FORMAT_INPUTS = PLANTS_INPUTS_INDEXES + PLANTS_STATE
    PLANTS_FORMAT_OUTPUTS = PLANTS_RUN_VARIABLES + PLANTS_POSTPROCESSING_VARIABLES

    AXES_INPUTS_INDEXES = ['plant', 'axis']
    AXES_OUTPUTS_INDEXES = [T_INDEX] + AXES_INPUTS_INDEXES
    AXES_INTERNAL_STATE_PARAMETERS = []
    AXES_EXTERNAL_STATE_PARAMETERS = []
    AXES_STATE_PARAMETERS = AXES_INTERNAL_STATE_PARAMETERS + AXES_EXTERNAL_STATE_PARAMETERS
    AXES_STATE = AXES_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.Axis, [])
    AXES_INTERMEDIATE_VARIABLES = []
    AXES_FLUXES = []
    AXES_INTEGRATIVE_VARIABLES = ['Total_transpiration']
    AXES_RUN_VARIABLES = AXES_OUTPUTS_INDEXES + AXES_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.Axis, []) + AXES_INTERMEDIATE_VARIABLES + AXES_FLUXES + AXES_INTEGRATIVE_VARIABLES
    AXES_POSTPROCESSING_VARIABLES = []
    AXES_FORMAT_INPUTS = AXES_INPUTS_INDEXES + AXES_STATE
    AXES_FORMAT_OUTPUTS = AXES_RUN_VARIABLES + AXES_POSTPROCESSING_VARIABLES
    
    PHYTOMERS_INPUTS_INDEXES = ['plant', 'axis', 'phytomer']
    PHYTOMERS_OUTPUTS_INDEXES = [T_INDEX] + PHYTOMERS_INPUTS_INDEXES
    PHYTOMERS_INTERNAL_STATE_PARAMETERS = []
    PHYTOMERS_EXTERNAL_STATE_PARAMETERS = []
    PHYTOMERS_STATE_PARAMETERS = PHYTOMERS_INTERNAL_STATE_PARAMETERS + PHYTOMERS_EXTERNAL_STATE_PARAMETERS
    PHYTOMERS_STATE = PHYTOMERS_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.Phytomer, [])
    PHYTOMERS_INTERMEDIATE_VARIABLES = []
    PHYTOMERS_FLUXES = []
    PHYTOMERS_INTEGRATIVE_VARIABLES = []
    PHYTOMERS_RUN_VARIABLES = PHYTOMERS_OUTPUTS_INDEXES + PHYTOMERS_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.Phytomer, []) + PHYTOMERS_INTERMEDIATE_VARIABLES + PHYTOMERS_FLUXES + PHYTOMERS_INTEGRATIVE_VARIABLES
    PHYTOMERS_POSTPROCESSING_VARIABLES = []
    PHYTOMERS_FORMAT_INPUTS = PHYTOMERS_INPUTS_INDEXES + PHYTOMERS_STATE
    PHYTOMERS_FORMAT_OUTPUTS = PHYTOMERS_RUN_VARIABLES + PHYTOMERS_POSTPROCESSING_VARIABLES

    ORGANS_INPUTS_INDEXES = ['plant', 'axis', 'phytomer', 'organ']
    ORGANS_OUTPUTS_INDEXES = [T_INDEX] + ORGANS_INPUTS_INDEXES
    ORGANS_INTERNAL_STATE_PARAMETERS = ['volume']
    ORGANS_EXTERNAL_STATE_PARAMETERS = ['mstruct_growth', 'Nstruct_N_growth', 'mstruct_senescence']
    ORGANS_STATE_PARAMETERS = ORGANS_INTERNAL_STATE_PARAMETERS + ORGANS_EXTERNAL_STATE_PARAMETERS
    ORGANS_STATE = ORGANS_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.Organ, [])
    ORGANS_INTERMEDIATE_VARIABLES = ['Conc_Nitrates_Soil', 'RGR_Structure', 'R_Nnit_upt', 'R_Nnit_red', 'R_residual', 'R_grain_growth_struct', 'R_grain_growth_starch', 'R_growth', 
                                     'C_exudation', 'N_exudation']
    ORGANS_FLUXES = ['Unloading_Sucrose', 'Export_Amino_Acids', 'Potential_Uptake_Nitrates', 'S_Proteins', 'S_Amino_Acids', 'Unloading_Amino_Acids', 
                     'S_grain_starch', 'Uptake_Nitrates', 'S_grain_structure']
    ORGANS_INTEGRATIVE_VARIABLES = ['total_nitrogen']
    ORGANS_RUN_VARIABLES = ORGANS_OUTPUTS_INDEXES + ORGANS_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.Organ, []) + ORGANS_INTERMEDIATE_VARIABLES + ORGANS_FLUXES + ORGANS_INTEGRATIVE_VARIABLES
    ORGANS_POSTPROCESSING_VARIABLES = ['Conc_Nitrates', 'Conc_Amino_Acids', 'Dry_Mass', 'Conc_Sucrose', 'Proteins_N_Mass']
    ORGANS_FORMAT_INPUTS = ORGANS_INPUTS_INDEXES + ORGANS_STATE
    ORGANS_FORMAT_OUTPUTS = ORGANS_RUN_VARIABLES + ORGANS_POSTPROCESSING_VARIABLES
    
    ELEMENTS_INPUTS_INDEXES = ['plant', 'axis', 'phytomer', 'organ', 'element']
    ELEMENTS_OUTPUTS_INDEXES = [T_INDEX] + ELEMENTS_INPUTS_INDEXES
    ELEMENTS_INTERNAL_STATE_PARAMETERS = ['area', 'green_area', 'width', 'height']
    ELEMENTS_EXTERNAL_STATE_PARAMETERS = ['An', 'Tr', 'Ag', 'Rd', 'Ts', 'gs', 'surfacic_nitrogen', 'relative_delta_green_area']
    ELEMENTS_STATE_PARAMETERS = ELEMENTS_INTERNAL_STATE_PARAMETERS + ELEMENTS_EXTERNAL_STATE_PARAMETERS
    ELEMENTS_STATE = ELEMENTS_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.PhotosyntheticOrganElement, [])
    ELEMENTS_INTERMEDIATE_VARIABLES = ['Transpiration', 'R_phloem_loading', 'R_Nnit_red', 'R_residual', 'Photosynthesis']
    ELEMENTS_FLUXES = ['Loading_Sucrose', 'Regul_S_Fructan', 'S_Starch', 'D_Starch', 'S_Sucrose', 'S_Fructan', 'D_Fructan', 'remob_starch_senescence', 
                       'remob_fructan_senescence', 'remob_proteins_senescence']
    ELEMENTS_INTEGRATIVE_VARIABLES = ['total_nitrogen']
    ELEMENTS_RUN_VARIABLES = ELEMENTS_OUTPUTS_INDEXES + ELEMENTS_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.PhotosyntheticOrganElement, []) + ELEMENTS_INTERMEDIATE_VARIABLES + ELEMENTS_FLUXES + ELEMENTS_INTEGRATIVE_VARIABLES
    ELEMENTS_POSTPROCESSING_VARIABLES = ['Conc_TriosesP', 'Conc_Starch', 'Conc_Sucrose', 'Conc_Fructan', 'Conc_Nitrates', 'Conc_Amino_Acids', 'Conc_Proteins', 
                                         'SLN', 'Nitrates_import', 'Amino_Acids_import', 'S_Amino_Acids', 'S_Proteins', 'D_Proteins', 'Loading_Amino_Acids']
    ELEMENTS_FORMAT_INPUTS = ELEMENTS_INPUTS_INDEXES + ELEMENTS_STATE
    ELEMENTS_FORMAT_OUTPUTS = ELEMENTS_RUN_VARIABLES + ELEMENTS_POSTPROCESSING_VARIABLES
    
    LOGGERS_NAMES = {'compartments': {model.Plant: 'cnwheat.compartments.plants',
                                      model.Axis: 'cnwheat.compartments.axes',
                                      model.Phytomer: 'cnwheat.compartments.phytomers',
                                      model.Organ: 'cnwheat.compartments.organs',
                                      model.PhotosyntheticOrganElement: 'cnwheat.compartments.elements'},
                     'derivatives': {model.Plant: 'cnwheat.derivatives.plants',
                                     model.Axis: 'cnwheat.derivatives.axes',
                                     model.Phytomer: 'cnwheat.derivatives.phytomers',
                                     model.Organ: 'cnwheat.derivatives.organs',
                                     model.PhotosyntheticOrganElement: 'cnwheat.derivatives.elements'}}

    def __init__(self):

        self.population = model.Population() #: the population to simulate on
        
        self.initial_conditions = [] #: the initial conditions of the compartments in the population
        self.initial_conditions_mapping = {} #: dictionary to map the compartments to their indexes in :attr:`initial_conditions`
        
        self._time_grid = np.array([]) #: the time grid of the simulation
        self._solver_output = np.array([]) #: the value of the compartments for each time step, with the initial conditions in the first row

        self.progressbar = ProgressBar(title='Solver progress') #: progress bar to show the progress of the solver
        self.show_progressbar = False #: True: show the progress bar ; False: DO NOT show the progress bar

    def initialize(self, population=None, plants_inputs=None, axes_inputs=None, phytomers_inputs=None, organs_inputs=None, elements_inputs=None):
        """
        Initialize :attr:`population` from `population` (if given), and complete it with `plants_inputs`, `axes_inputs`, `phytomers_inputs`, 
        `organs_inputs` and `elements_inputs` (if given).
        
        :Parameters:
            
            - `population` (:class:`model.Population`) - a population of plants.
            
            - `plants_inputs` (:class:`pandas.DataFrame`) - Plants inputs, with one line by plant.
            
            - `axes_inputs` (:class:`pandas.DataFrame`) - Axes inputs, with one line by axis.
            
            - `phytomers_inputs` (:class:`pandas.DataFrame`) - Phytomers inputs, with one line by phytomer.
            
            - `organs_inputs` (:class:`pandas.DataFrame`) - Organs inputs, with one line by organ.
            
            - `elements_inputs` (:class:`pandas.DataFrame`) - Elements inputs, with one line by element.
            
            .. seealso:: See :mod:`cnwheat.model` for more information about the inputs.  
            
        """
        
        logger = logging.getLogger(__name__)

        logger.info('Initialization of the simulation...')
        
        if (plants_inputs is None or axes_inputs is None or phytomers_inputs is None or organs_inputs is None or elements_inputs is None) and population is None:
            raise SimulationInputsError('You must pass either {} or {} or both.'.format('population', ('plants_inputs', 'axes_inputs', 'phytomers_inputs', 'organs_inputs', 'elements_inputs')))
            
        # clean the attributes of the simulation
        del self.population.plants[:]
        del self.initial_conditions[:]
        self.initial_conditions_mapping.clear()
        
        # build the population
        if population is not None:
            self.population.plants.extend(population.plants)
            
        if not (plants_inputs is None or axes_inputs is None or phytomers_inputs is None or organs_inputs is None or elements_inputs is None):
            
            elements_attributes_names = [state_var_name for state_var_name in Simulation.ELEMENTS_STATE if state_var_name in elements_inputs]
            
            for curr_plant_index in plants_inputs.plant:
                # create a new plant
                curr_plant = model.Plant(index=curr_plant_index)
                self.population.plants.append(curr_plant)
                curr_axes_inputs = axes_inputs[axes_inputs['plant'] == curr_plant_index]
                for curr_axis_id in curr_axes_inputs.axis:
                    # create a new axis
                    curr_axis = model.Axis(axis_id=curr_axis_id)
                    curr_plant.axes.append(curr_axis)
                    
                    curr_organs_inputs = organs_inputs[(organs_inputs['plant'] == curr_plant_index) & (organs_inputs['axis'] == curr_axis_id)]
                    
                    for curr_axis_attribute_name, curr_axis_attribute_class in (('grains', model.Grains), ('roots', model.Roots), ('soil', model.Soil), ('phloem', model.Phloem)):
                        curr_organ_class_name = curr_axis_attribute_class.__name__
                        curr_organ_inputs = curr_organs_inputs[curr_organs_inputs['organ'] == curr_organ_class_name]
                        if len(curr_organ_inputs) != 0:
                            # create a new organ
                            curr_organ = curr_axis_attribute_class()
                            curr_organ_attributes_names = [state_var_name for state_var_name in Simulation.ORGANS_STATE if hasattr(curr_organ, state_var_name)]
                            curr_organ_inputs_row = curr_organ_inputs.loc[curr_organ_inputs.first_valid_index()]
                            curr_organ_attributes_values = curr_organ_inputs_row[curr_organ_attributes_names].tolist()
                            curr_organ_attributes = dict(zip(curr_organ_attributes_names, curr_organ_attributes_values))
                            curr_organ.__dict__.update(curr_organ_attributes)
                            setattr(curr_axis, curr_axis_attribute_name, curr_organ)
                            
                    curr_phytomers_inputs = phytomers_inputs[(phytomers_inputs['plant'] == curr_plant_index) & (phytomers_inputs['axis'] == curr_axis_id)]
                    for curr_phytomer_index in curr_phytomers_inputs.phytomer:
                        # create a new phytomer
                        curr_phytomer = model.Phytomer(index=curr_phytomer_index)
                        curr_axis.phytomers.append(curr_phytomer)
                        
                        curr_organs_inputs = organs_inputs[(organs_inputs['plant'] == curr_plant_index) & (organs_inputs['axis'] == curr_axis_id) & (organs_inputs['phytomer'] == curr_phytomer_index)]
                        
                        for curr_phytomer_attribute_name, curr_phytomer_attribute_class, curr_phytomer_attribute_element_class in \
                            (('chaff', model.Chaff, model.ChaffElement), 
                             ('lamina', model.Lamina, model.LaminaElement), 
                             ('internode', model.Internode, model.InternodeElement), 
                             ('peduncle', model.Peduncle, model.PeduncleElement),
                             ('sheath', model.Sheath, model.SheathElement)):
                            
                            curr_organ_class_name = curr_phytomer_attribute_class.__name__
                            curr_elements_inputs = elements_inputs[(elements_inputs['plant'] == curr_plant_index) & (elements_inputs['axis'] == curr_axis_id) & (elements_inputs['phytomer'] == curr_phytomer_index) & (elements_inputs['organ'] == curr_organ_class_name)]
                        
                            if curr_organ_class_name not in curr_organs_inputs.organ.values and curr_organ_class_name not in curr_elements_inputs.organ.values:
                                continue
                            # create a new organ
                            curr_organ = curr_phytomer_attribute_class()
                            setattr(curr_phytomer, curr_phytomer_attribute_name, curr_organ)
                            
                            for curr_organ_attribute_name, curr_organ_attribute_type in (('enclosed_element', 'enclosed'), ('exposed_element', 'exposed')):
                                curr_element_inputs = curr_elements_inputs[curr_elements_inputs['element'] == curr_organ_attribute_type][elements_attributes_names]
                                if len(curr_element_inputs) == 0:
                                    continue
                                # create a new element
                                curr_element_inputs_dict = curr_element_inputs.loc[curr_element_inputs.first_valid_index()].to_dict()
                                curr_element = curr_phytomer_attribute_element_class(**curr_element_inputs_dict)
                                setattr(curr_organ, curr_organ_attribute_name, curr_element)
                                
                        
        # initialize initial conditions
        def _init_initial_conditions(model_object, i):
            class_ = model_object.__class__
            if issubclass(class_, model.Organ):
                class_ = model.Organ
            elif issubclass(class_, model.PhotosyntheticOrganElement):
                class_ = model.PhotosyntheticOrganElement
            compartments_names = Simulation.MODEL_COMPARTMENTS_NAMES[class_]
            self.initial_conditions_mapping[model_object] = {}
            for compartment_name in compartments_names:
                if hasattr(model_object, compartment_name):
                    self.initial_conditions_mapping[model_object][compartment_name] = i
                    self.initial_conditions.append(0)
                    i += 1
            return i

        i = 0

        for plant in self.population.plants:
            i = _init_initial_conditions(plant, i)
            for axis in plant.axes:
                i = _init_initial_conditions(axis, i)
                for organ in (axis.roots, axis.soil, axis.phloem, axis.grains):
                    if organ is None:
                        continue
                    i = _init_initial_conditions(organ, i)
                for phytomer in axis.phytomers:
                    i = _init_initial_conditions(phytomer, i)
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ is None:
                            continue
                        i = _init_initial_conditions(organ, i)
                        for element in (organ.exposed_element, organ.enclosed_element):
                            if element is None:
                                continue
                            i = _init_initial_conditions(element, i)
        
        #TODO: check the consistency of the population
        logger.info('Initialization of the simulation DONE')

    def run(self, start_time, stop_time, number_of_output_steps, odeint_mxstep=0, show_progressbar=False):
        """
        Compute CN exchanges in :attr:`population` from `start_time` to `stop_time`, for `number_of_output_steps` steps.

        :Parameters:

            - `start_time` (:class:`int`) - The starting of the time grid.

            - `stop_time` (:class:`int`) - The end of the time grid.

            - `number_of_output_steps` (:class:`int`) - Number of time points for which to compute the CN exchanges in :attr:`population`.

            - `odeint_mxstep` (:class:`int`) - Maximum number of (internally defined) steps allowed for each integration point in time grid.
              `odeint_mxstep` is passed to :func:`scipy.integrate.odeint` as `mxstep`. If `odeint_mxstep` = 0 (the default), then `mxstep` is determined by the solver.
              Normally, the `mxstep` determined by the solver permits to solve the current model. User can try to increase this value if a more complex model is defined
              and if the integration failed. However, take care that the origin of an integration failure could be a discontinuity in the RHS function used
              by :func:`scipy.integrate.odeint`, and that this discontinuity could be due to a bug in your model. To summary: if the integration failed, first
              check the logs.
              
            - `show_progressbar` (:class:`bool`) - True: show the progress bar ; False: do not show the progress bar.
            
        :Returns:
            Dictionary containing output information from the solver. 
            This is the dictionary returned by :func:`scipy.integrate.odeint` as second output. 
            See the documentation of :func:`scipy.integrate.odeint` for more information. 

        :Returns Type:
            :class:`dict`

        """
        logger = logging.getLogger(__name__)
        logger.info('Run of CN-Wheat from {} to {}...'.format(start_time, stop_time))

        t = np.linspace(start_time, stop_time, number_of_output_steps)
        
        self.show_progressbar = show_progressbar
        if self.show_progressbar:
            self.progressbar.set_t_max(stop_time)

        compartments_logger = logging.getLogger('cnwheat.compartments')
        derivatives_logger = logging.getLogger('cnwheat.derivatives')
        if compartments_logger.isEnabledFor(logging.DEBUG) or derivatives_logger.isEnabledFor(logging.DEBUG):
            sep = ','
            if compartments_logger.isEnabledFor(logging.DEBUG):
                plants_compartments_logger = logging.getLogger('cnwheat.compartments.plants')
                plants_compartments_logger.debug(sep.join(Simulation.PLANTS_OUTPUTS_INDEXES + Simulation.MODEL_COMPARTMENTS_NAMES[model.Plant]))
                axes_compartments_logger = logging.getLogger('cnwheat.compartments.axes')
                axes_compartments_logger.debug(sep.join(Simulation.AXES_OUTPUTS_INDEXES + Simulation.MODEL_COMPARTMENTS_NAMES[model.Axis]))
                phytomers_compartments_logger = logging.getLogger('cnwheat.compartments.phytomers')
                phytomers_compartments_logger.debug(sep.join(Simulation.PHYTOMERS_OUTPUTS_INDEXES + Simulation.MODEL_COMPARTMENTS_NAMES[model.Phytomer]))
                organs_compartments_logger = logging.getLogger('cnwheat.compartments.organs')
                organs_compartments_logger.debug(sep.join(Simulation.ORGANS_OUTPUTS_INDEXES + Simulation.MODEL_COMPARTMENTS_NAMES[model.Organ]))
                elements_compartments_logger = logging.getLogger('cnwheat.compartments.elements')
                elements_compartments_logger.debug(sep.join(Simulation.ELEMENTS_OUTPUTS_INDEXES + Simulation.MODEL_COMPARTMENTS_NAMES[model.PhotosyntheticOrganElement]))
            if derivatives_logger.isEnabledFor(logging.DEBUG):
                plants_derivatives_logger = logging.getLogger('cnwheat.derivatives.plants')
                plants_derivatives_logger.debug(sep.join(Simulation.PLANTS_OUTPUTS_INDEXES + Simulation.MODEL_COMPARTMENTS_NAMES[model.Plant]))
                axes_derivatives_logger = logging.getLogger('cnwheat.derivatives.axes')
                axes_derivatives_logger.debug(sep.join(Simulation.AXES_OUTPUTS_INDEXES + Simulation.MODEL_COMPARTMENTS_NAMES[model.Axis]))
                phytomers_derivatives_logger = logging.getLogger('cnwheat.derivatives.phytomers')
                phytomers_derivatives_logger.debug(sep.join(Simulation.PHYTOMERS_OUTPUTS_INDEXES + Simulation.MODEL_COMPARTMENTS_NAMES[model.Phytomer]))
                organs_derivatives_logger = logging.getLogger('cnwheat.derivatives.organs')
                organs_derivatives_logger.debug(sep.join(Simulation.ORGANS_OUTPUTS_INDEXES + Simulation.MODEL_COMPARTMENTS_NAMES[model.Organ]))
                elements_derivatives_logger = logging.getLogger('cnwheat.derivatives.elements')
                elements_derivatives_logger.debug(sep.join(Simulation.ELEMENTS_OUTPUTS_INDEXES + Simulation.MODEL_COMPARTMENTS_NAMES[model.PhotosyntheticOrganElement]))

        self._update_initial_conditions()

        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(
                """Run the solver with:
                    - time grid = %s,
                    - odeint mxstep = %s""",
                t, odeint_mxstep)

        soln, infodict = odeint(self._calculate_all_derivatives, self.initial_conditions, t, full_output=True, mxstep=odeint_mxstep)

        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(
                """Run of the solver DONE: infodict = %s""",
                infodict)
        
        # update self._time_grid
        self._time_grid.resize(t.shape)
        np.copyto(self._time_grid, t)
        
        # update self._solver_output
        self._solver_output.resize(soln.shape)
        np.copyto(self._solver_output, soln)

        if not set(infodict['mused']).issubset([1,2]):
            message = "Integration failed. See the logs of lsoda or try to increase the value of 'mxstep'."
            logger.exception(message)
            raise SimulationRunError(message)
        
        last_compartments_values = self._solver_output[-1]
        self._update_population(last_compartments_values)

        self.population.calculate_integrative_variables()

        logger.info('Run of CN-Wheat from {} to {} DONE'.format(start_time, stop_time))

        return infodict


    def _update_initial_conditions(self):
        """Update the compartments values in :attr:`initial_conditions` from the compartments values of :attr:`population`.
        """
        for model_object, compartments in self.initial_conditions_mapping.iteritems():
            for compartment_name, compartment_index in compartments.iteritems():
                self.initial_conditions[compartment_index] = getattr(model_object, compartment_name)
                

    def _log_compartments(self, t, y, loggers_names):
        """Log the values in `y` to the loggers in `loggers_names`.
        """
        def update_rows(model_object, indexes, rows, i):
            row = []
            class_ = model_object.__class__
            if issubclass(class_, model.Organ):
                class_ = model.Organ
            elif issubclass(class_, model.PhotosyntheticOrganElement):
                class_ = model.PhotosyntheticOrganElement
            compartments_names = Simulation.MODEL_COMPARTMENTS_NAMES[class_]
            for compartment_name in compartments_names:
                if hasattr(model_object, compartment_name):
                    row.append(str(y[i]))
                    i += 1
                else:
                    row.append('NA')
            if row.count('NA') < len(row):
                rows.append([str(index) for index in indexes] + row)
            return i

        i = 0
        all_rows = dict([(class_, []) for class_ in loggers_names])
        for plant in self.population.plants:
            i = update_rows(plant, [t, plant.index], all_rows[model.Plant], i)
            for axis in plant.axes:
                i = update_rows(axis, [t, plant.index, axis.id], all_rows[model.Axis], i)
                for organ in (axis.roots, axis.soil, axis.phloem, axis.grains):
                    if organ is None:
                        continue
                    i = update_rows(organ, [t, plant.index, axis.id, 'NA', organ.__class__.__name__], all_rows[model.Organ], i)
                for phytomer in axis.phytomers:
                    i = update_rows(phytomer, [t, plant.index, axis.id, phytomer.index], all_rows[model.Phytomer], i)
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ is None:
                            continue
                        i = update_rows(organ, [t, plant.index, axis.id, phytomer.index, organ.__class__.__name__], all_rows[model.Organ], i)
                        for element, element_type in ((organ.exposed_element, 'exposed'), (organ.enclosed_element, 'enclosed')):
                            if element is None:
                                continue
                            i = update_rows(element, [t, plant.index, axis.id, phytomer.index, organ.__class__.__name__, element_type], all_rows[model.PhotosyntheticOrganElement], i)

        row_sep = '\n'
        column_sep = ','
        for class_, logger_name in loggers_names.iteritems():
            compartments_logger = logging.getLogger(logger_name)
            formatted_initial_conditions = row_sep.join([column_sep.join(row) for row in all_rows[class_]])
            compartments_logger.debug(formatted_initial_conditions)


    def _calculate_all_derivatives(self, y, t):
        """Compute the derivative of `y` at `t`.

        :meth:`_calculate_all_derivatives` is passed as **func** argument to
        :func:`odeint(func, y0, t, args=(),...) <scipy.integrate.odeint>`.
        :meth:`_calculate_all_derivatives` is called automatically by
        :func:`scipy.integrate.odeint <scipy.integrate.odeint>`.

        First call to :meth:`_calculate_all_derivatives` uses `y` = **y0** and
        `t` = **t** [0], where **y0** and **t** are arguments passed to :func:`odeint(func, y0, t, args=(),...) <scipy.integrate.odeint>`.

        Following calls to :meth:`_calculate_all_derivatives` use `t` in [min( **t** ), max( **t** ) + 1] where
        **t** is an argument passed to :func:`odeint(func, y0, t, args=(),...) <scipy.integrate.odeint>`. `y` is
        computed automatically by the solver.

        :Parameters:

            - `y` (:class:`list`) - The current y values. `y` is automatically set by
              :func:`scipy.integrate.odeint`. User does not have control over `y`.
              At first call to :meth:`_calculate_all_derivatives` by :func:`scipy.integrate.odeint`, `y` = **y0**
              where **y0** is one of the arguments passed to :func:`odeint(func, y0, t, args=(),...) <scipy.integrate.odeint>`.
              For each following call to :meth:`_calculate_all_derivatives`, `y` is
              computed automatically by the solver.

            - `t` (:class:`float`) - The current t at which we want to compute the derivatives.
              `t` is automatically set by :func:`scipy.integrate.odeint`.
              User does not have control over `t`.
              At first call to :meth:`_calculate_all_derivatives` :func:`scipy.integrate.odeint`,
              `t` = **t** [0], where **t** is one of the arguments passed to :func:`odeint(func, y0, t, args=(),...) <scipy.integrate.odeint>`.
              For each following call to :meth:`_calculate_all_derivatives`, `t` belongs
              to the interval [min( **t** ), max( **t** ) + 1], where **t** is an
              argument passed to :func:`odeint(func, y0, t, args=(),...) <scipy.integrate.odeint>`.

        :Returns:
            The derivatives of `y` at `t`.

        :Returns Type:
            :class:`list`


        """
        logger = logging.getLogger(__name__)

        logger.debug('t = {}'.format(t))

        compartments_logger = logging.getLogger('cnwheat.compartments')
        if compartments_logger.isEnabledFor(logging.DEBUG):
            self._log_compartments(t, y, Simulation.LOGGERS_NAMES['compartments'])

        # check that the solver is not crashed
        y_isnan = np.isnan(y)
        if y_isnan.any():
            message = 'The solver did not manage to compute a compartment. See the logs.'
            logger.exception(message)
            raise SimulationRunError(message)

        y_derivatives = np.zeros_like(y)

        for plant in self.population.plants:
            for axis in plant.axes:

                phloem_contributors = []

                axis.phloem.sucrose = y[self.initial_conditions_mapping[axis.phloem]['sucrose']]

                axis.phloem.amino_acids = y[self.initial_conditions_mapping[axis.phloem]['amino_acids']]

                axis.roots.nitrates = y[self.initial_conditions_mapping[axis.roots]['nitrates']]
                axis.roots.amino_acids = y[self.initial_conditions_mapping[axis.roots]['amino_acids']]
                axis.soil.nitrates = y[self.initial_conditions_mapping[axis.soil]['nitrates']]
                phloem_contributors.append(axis.roots)

                # compute the total transpiration at t_inf
                total_transpiration = 0.0
                transpiration_mapping = {}
                for phytomer in axis.phytomers:
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ is not None:
                            for element in (organ.exposed_element, organ.enclosed_element):
                                if element is not None:
                                    transpiration_mapping[element] = element.calculate_transpiration(element.Tr, element.green_area)
                                    total_transpiration += transpiration_mapping[element]
                                    phloem_contributors.append(element)

                # compute the flows from/to the roots to/from photosynthetic organs
                conc_nitrates_soil = axis.soil.calculate_conc_nitrates(axis.soil.nitrates)
                roots_uptake_nitrate, potential_roots_uptake_nitrates = axis.roots.calculate_uptake_nitrates(conc_nitrates_soil, axis.roots.nitrates, total_transpiration)
                roots_export_amino_acids = axis.roots.calculate_export_amino_acids(axis.roots.amino_acids, total_transpiration)

                # compute the derivative of each photosynthetic organ element compartment
                for phytomer in axis.phytomers:
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):

                        if organ is None:
                            continue

                        for element in (organ.exposed_element, organ.enclosed_element):
                            if element is None:
                                continue

                            element.starch = y[self.initial_conditions_mapping[element]['starch']]
                            element.sucrose = y[self.initial_conditions_mapping[element]['sucrose']]
                            element.triosesP = y[self.initial_conditions_mapping[element]['triosesP']]
                            element.fructan = y[self.initial_conditions_mapping[element]['fructan']]
                            element.nitrates = y[self.initial_conditions_mapping[element]['nitrates']]
                            element.amino_acids = y[self.initial_conditions_mapping[element]['amino_acids']]
                            element.proteins = y[self.initial_conditions_mapping[element]['proteins']]

                            # intermediate variables
                            photosynthesis = element.calculate_photosynthesis(element.Ag, element.green_area)
                            element_transpiration = transpiration_mapping[element]

                            # flows
                            s_starch = element.calculate_s_starch(element.triosesP)
                            d_starch = element.calculate_d_starch(element.starch)
                            remob_starch_senescence = element.calculate_remob_starch_senescence(element.starch, element.relative_delta_green_area)
                            s_sucrose = element.calculate_s_sucrose(element.triosesP)
                            element.loading_sucrose = element.calculate_loading_sucrose(element.sucrose, axis.phloem.sucrose)
                            regul_s_fructan = element.calculate_regul_s_fructan(element.loading_sucrose)
                            d_fructan = element.calculate_d_fructan(element.sucrose, element.fructan)
                            remob_fructan_senescence = element.calculate_remob_fructan_senescence(element.fructan, element.relative_delta_green_area)
                            s_fructan = element.calculate_s_fructan(element.sucrose, regul_s_fructan)
                            nitrates_import = element.calculate_nitrates_import(roots_uptake_nitrate, element_transpiration, total_transpiration)
                            amino_acids_import = element.calculate_amino_acids_import(roots_export_amino_acids, element_transpiration, total_transpiration)
                            s_amino_acids = element.calculate_s_amino_acids(element.nitrates, element.triosesP)
                            s_proteins = element.calculate_s_proteins(element.amino_acids)
                            d_proteins = element.calculate_d_proteins(element.proteins)
                            remob_proteins_senescence = element.calculate_remob_proteins_senescence(element.proteins, element.relative_delta_green_area)
                            element.loading_amino_acids = element.calculate_loading_amino_acids(element.amino_acids, axis.phloem.amino_acids)

                            # compartments derivatives
                            starch_derivative = element.calculate_starch_derivative(s_starch, d_starch, remob_starch_senescence)
                            R_phloem_loading = RespirationModel.R_phloem(element.loading_sucrose, element.sucrose, element.mstruct*element.__class__.PARAMETERS.ALPHA)
                            R_Nnit_red = RespirationModel.R_Nnit_red(s_amino_acids, element.sucrose, element.mstruct*element.__class__.PARAMETERS.ALPHA)
                            element.total_nitrogen = element.calculate_total_nitrogen(element.nitrates, element.amino_acids, element.proteins, element.Nstruct)
                            R_residual = RespirationModel.R_residual(element.sucrose, element.mstruct*element.__class__.PARAMETERS.ALPHA, element.total_nitrogen, organ.__class__.PARAMETERS.DELTA_T)
                            sucrose_derivative = element.calculate_sucrose_derivative(s_sucrose, d_starch, remob_starch_senescence, element.loading_sucrose, s_fructan, d_fructan, remob_fructan_senescence, R_phloem_loading, R_Nnit_red, R_residual)
                            triosesP_derivative = element.calculate_triosesP_derivative(photosynthesis, s_sucrose, s_starch, s_amino_acids)
                            fructan_derivative = element.calculate_fructan_derivative(s_fructan, d_fructan, remob_fructan_senescence)
                            nitrates_derivative = element.calculate_nitrates_derivative(nitrates_import, s_amino_acids)
                            amino_acids_derivative = element.calculate_amino_acids_derivative(amino_acids_import, s_amino_acids, s_proteins, d_proteins, remob_proteins_senescence, element.loading_amino_acids)
                            proteins_derivative = element.calculate_proteins_derivative(s_proteins, d_proteins, remob_proteins_senescence)

                            y_derivatives[self.initial_conditions_mapping[element]['starch']] = starch_derivative
                            y_derivatives[self.initial_conditions_mapping[element]['sucrose']] = sucrose_derivative
                            y_derivatives[self.initial_conditions_mapping[element]['triosesP']] = triosesP_derivative
                            y_derivatives[self.initial_conditions_mapping[element]['fructan']] = fructan_derivative
                            y_derivatives[self.initial_conditions_mapping[element]['nitrates']] = nitrates_derivative
                            y_derivatives[self.initial_conditions_mapping[element]['amino_acids']] = amino_acids_derivative
                            y_derivatives[self.initial_conditions_mapping[element]['proteins']] = proteins_derivative

                if axis.grains is not None:
                    phloem_contributors.append(axis.grains)
                    # compute the derivative of each compartment of grains
                    axis.grains.structure = y[self.initial_conditions_mapping[axis.grains]['structure']]
                    axis.grains.starch = y[self.initial_conditions_mapping[axis.grains]['starch']]
                    axis.grains.proteins = y[self.initial_conditions_mapping[axis.grains]['proteins']]

                    # intermediate variables
                    RGR_structure = axis.grains.calculate_RGR_structure(axis.phloem.sucrose)
                    structural_dry_mass = axis.grains.calculate_structural_dry_mass(axis.grains.structure)

                    # flows
                    axis.grains.s_grain_structure = axis.grains.calculate_s_grain_structure(t, axis.grains.structure, RGR_structure)
                    axis.grains.s_grain_starch = axis.grains.calculate_s_grain_starch(t, axis.phloem.sucrose)
                    axis.grains.s_proteins = axis.grains.calculate_s_proteins(axis.grains.s_grain_structure, axis.grains.s_grain_starch, axis.phloem.amino_acids, axis.phloem.sucrose, structural_dry_mass)

                    # compartments derivatives
                    R_grain_growth_struct, R_grain_growth_starch = RespirationModel.R_grain_growth(axis.grains.s_grain_structure, axis.grains.s_grain_starch, structural_dry_mass)
                    structure_derivative = axis.grains.calculate_structure_derivative(axis.grains.s_grain_structure, R_grain_growth_struct)
                    starch_derivative = axis.grains.calculate_starch_derivative(axis.grains.s_grain_starch, structural_dry_mass, R_grain_growth_starch)
                    proteins_derivative = axis.grains.calculate_proteins_derivative(axis.grains.s_proteins)
                    y_derivatives[self.initial_conditions_mapping[axis.grains]['structure']] = structure_derivative
                    y_derivatives[self.initial_conditions_mapping[axis.grains]['starch']] = starch_derivative
                    y_derivatives[self.initial_conditions_mapping[axis.grains]['proteins']] = proteins_derivative

                # compute the derivative of each compartment of roots and soil
                y_derivatives[self.initial_conditions_mapping[axis.soil]['nitrates']] = axis.soil.calculate_nitrates_derivative(roots_uptake_nitrate)

                axis.roots.sucrose = y[self.initial_conditions_mapping[axis.roots]['sucrose']]
                mstruct_C_growth = axis.roots.mstruct_C_growth
                Nstruct_N_growth = axis.roots.Nstruct_N_growth

                # flows
                axis.roots.unloading_sucrose = axis.roots.calculate_unloading_sucrose(axis.phloem.sucrose)
                axis.roots.unloading_amino_acids = axis.roots.calculate_unloading_amino_acids(axis.roots.unloading_sucrose, axis.phloem.sucrose, axis.phloem.amino_acids)
                axis.roots.s_amino_acids = axis.roots.calculate_s_amino_acids(axis.roots.nitrates, axis.roots.sucrose)
                C_exudated, N_exudated = axis.roots.calculate_exudation(axis.roots.unloading_sucrose, axis.phloem.sucrose, axis.phloem.amino_acids)

                # compartments derivatives
                axis.roots.total_nitrogen = axis.roots.calculate_total_nitrogen(axis.roots.nitrates, axis.roots.amino_acids, axis.roots.Nstruct)
                R_Nnit_upt = RespirationModel.R_Nnit_upt(roots_uptake_nitrate, axis.roots.sucrose)
                R_Nnit_red = RespirationModel.R_Nnit_red(axis.roots.s_amino_acids, axis.roots.sucrose, axis.roots.mstruct*model.Roots.PARAMETERS.ALPHA, root=True)
                R_residual = RespirationModel.R_residual(axis.roots.sucrose, axis.roots.mstruct*model.Roots.PARAMETERS.ALPHA, axis.roots.total_nitrogen, axis.roots.PARAMETERS.DELTA_T)
                R_roots_growth = RespirationModel.R_growth(mstruct_C_growth, axis.roots.mstruct)
                sucrose_derivative = axis.roots.calculate_sucrose_derivative(axis.roots.unloading_sucrose, axis.roots.s_amino_acids, mstruct_C_growth, C_exudated, R_Nnit_upt, R_Nnit_red, R_residual, R_roots_growth)
                nitrates_derivative = axis.roots.calculate_nitrates_derivative(roots_uptake_nitrate, axis.roots.s_amino_acids)
                amino_acids_derivative = axis.roots.calculate_amino_acids_derivative(axis.roots.unloading_amino_acids, axis.roots.s_amino_acids, roots_export_amino_acids, Nstruct_N_growth, N_exudated)
                y_derivatives[self.initial_conditions_mapping[axis.roots]['sucrose']] = sucrose_derivative
                y_derivatives[self.initial_conditions_mapping[axis.roots]['nitrates']] = nitrates_derivative
                y_derivatives[self.initial_conditions_mapping[axis.roots]['amino_acids']] = amino_acids_derivative

                # compute the derivative of each compartment of phloem
                sucrose_phloem_derivative = axis.phloem.calculate_sucrose_derivative(phloem_contributors)
                amino_acids_phloem_derivative = axis.phloem.calculate_amino_acids_derivative(phloem_contributors)
                y_derivatives[self.initial_conditions_mapping[axis.phloem]['sucrose']] = sucrose_phloem_derivative
                y_derivatives[self.initial_conditions_mapping[axis.phloem]['amino_acids']] = amino_acids_phloem_derivative

        if self.show_progressbar:
            self.progressbar.update(t)

        derivatives_logger = logging.getLogger('cnwheat.derivatives')
        if derivatives_logger.isEnabledFor(logging.DEBUG):
            self._log_compartments(t, y_derivatives, Simulation.LOGGERS_NAMES['derivatives'])

        return y_derivatives


    def _update_population(self, compartments_values):
        """Update the state of :attr:`population` from the values in `compartments_values`.
        """
        logger = logging.getLogger(__name__)
        logger.debug('Updating the state of the population...')
        for model_object, compartments in self.initial_conditions_mapping.iteritems():
            for compartment_name, compartment_index in compartments.iteritems():
                setattr(model_object, compartment_name, compartments_values[compartment_index])
        logger.debug('Updating the state of the population DONE')
        
    
    def format_inputs(self):
        """
        Format internal state parameters and compartments of :attr:`population` to Pandas dataframes.
        
        :Returns:
            :class:`dataframes <pandas.DataFrame>` of internal state parameters and compartments at different scales:

                * plant: plant index, state parameters and compartments of each plant (see :attr:`Simulation:PLANTS_FORMAT_INPUTS`) 
                * axis: plant index, axis id, state parameters and compartments of each axis (see :attr:`Simulation:AXES_FORMAT_INPUTS`)
                * phytomer: plant index, axis id, phytomer index, state parameters and compartments of phytomer plant (see :attr:`Simulation:PHYTOMERS_FORMAT_INPUTS`)
                * organ: plant index, axis id, phytomer index, organ type, state parameters and compartments of each organ (see :attr:`Simulation:ORGANS_FORMAT_INPUTS`)
                * and element: plant index, axis id, phytomer index, organ type, element type, state parameters and compartments of each element (see :attr:`Simulation:ELEMENTS_FORMAT_INPUTS`)
        
        :Returns Type:
            :class:`tuple` of :class:`pandas.DataFrame`
        
        """
        
        logger = logging.getLogger(__name__)
        logger.debug('Formatting of inputs...')
        
        all_plants_df = pd.DataFrame(columns=Simulation.PLANTS_FORMAT_INPUTS)
        all_axes_df = pd.DataFrame(columns=Simulation.AXES_FORMAT_INPUTS)
        all_phytomers_df = pd.DataFrame(columns=Simulation.PHYTOMERS_FORMAT_INPUTS)
        all_organs_df = pd.DataFrame(columns=Simulation.ORGANS_FORMAT_INPUTS)
        all_elements_df = pd.DataFrame(columns=Simulation.ELEMENTS_FORMAT_INPUTS)
        
        def append_row(model_object, indexes, attributes_names, inputs_df):
            attributes_values = []
            for attribute_name in attributes_names:
                attributes_values.append(getattr(model_object, attribute_name, np.nan))
            if len(attributes_values) == 0 \
                or any(attribute_value is None for attribute_value in attributes_values) \
                or not np.isnan([attribute_value for attribute_value in attributes_values if attribute_value is not None]).all():
                inputs_df.loc[len(inputs_df),:] = indexes + attributes_values
                    
        for plant in self.population.plants:
            append_row(plant, [plant.index], Simulation.PLANTS_STATE, all_plants_df)
            for axis in plant.axes:
                append_row(axis, [plant.index, axis.id], Simulation.AXES_STATE, all_axes_df)
                for organ in (axis.roots, axis.soil, axis.phloem, axis.grains):
                    if organ is None:
                        continue
                    append_row(organ, [plant.index, axis.id, np.nan, organ.__class__.__name__], Simulation.ORGANS_STATE, all_organs_df)
                for phytomer in axis.phytomers:
                    append_row(phytomer, [plant.index, axis.id, phytomer.index], Simulation.PHYTOMERS_STATE, all_phytomers_df)
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ is None:
                            continue
                        append_row(organ, [plant.index, axis.id, phytomer.index, organ.__class__.__name__], Simulation.ORGANS_STATE, all_organs_df)
                        for element, element_type in ((organ.exposed_element, 'exposed'), (organ.enclosed_element, 'enclosed')):
                            if element is None:
                                continue
                            append_row(element, [plant.index, axis.id, phytomer.index, organ.__class__.__name__, element_type], Simulation.ELEMENTS_STATE, all_elements_df)
                
        # sort the rows by the columns
        all_plants_df.sort_index(by=Simulation.PLANTS_FORMAT_INPUTS, inplace=True)
        all_axes_df.sort_index(by=Simulation.AXES_FORMAT_INPUTS, inplace=True)
        all_phytomers_df.sort_index(by=Simulation.PHYTOMERS_FORMAT_INPUTS, inplace=True)
        all_organs_df.sort_index(by=Simulation.ORGANS_FORMAT_INPUTS, inplace=True)
        all_elements_df.sort_index(by=Simulation.ELEMENTS_FORMAT_INPUTS, inplace=True)

        # infer the right types of the columns
        all_plants_df = all_plants_df.convert_objects(copy=False)
        all_axes_df = all_axes_df.convert_objects(copy=False)
        all_phytomers_df = all_phytomers_df.convert_objects(copy=False)
        all_organs_df = all_organs_df.convert_objects(copy=False)
        all_elements_df = all_elements_df.convert_objects(copy=False)

        # convert the indexes of plants, phytomers and elements to integers
        all_plants_df['plant'] = all_plants_df['plant'].astype(int)
        all_axes_df['plant'] = all_axes_df['plant'].astype(int)
        all_phytomers_df[['plant', 'phytomer']] = all_phytomers_df[['plant', 'phytomer']].astype(int)
        all_organs_df['plant'] = all_organs_df['plant'].astype(int)
        all_elements_df[['plant', 'phytomer']] = all_elements_df[['plant', 'phytomer']].astype(int)
        
        logger.debug('Formatting of inputs DONE')
        
        return all_plants_df, all_axes_df, all_phytomers_df, all_organs_df, all_elements_df


    def format_outputs(self):
        """
        Create :class:`dataframes <pandas.DataFrame>` of outputs at different scales.
        
        :Returns:
            :class:`dataframes <pandas.DataFrame>` of outputs at different scales:

                * plant: plant index, state parameters and compartments of each plant (see :attr:`Simulation:PLANTS_FORMAT_INPUTS`) 
                * axis: plant index, axis id, state parameters and compartments of each axis (see :attr:`Simulation:AXES_FORMAT_INPUTS`)
                * phytomer: plant index, axis id, phytomer index, state parameters and compartments of phytomer plant (see :attr:`Simulation:PHYTOMERS_FORMAT_INPUTS`)
                * organ: plant index, axis id, phytomer index, organ type, state parameters and compartments of each organ (see :attr:`Simulation:ORGANS_FORMAT_INPUTS`)
                * and element: plant index, axis id, phytomer index, organ type, element type, state parameters and compartments of each element (see :attr:`Simulation:ELEMENTS_FORMAT_INPUTS`)
        
        :Returns Type:
            :class:`tuple` of :class:`pandas.DataFrame`
        
        """
        logger = logging.getLogger(__name__)
        logger.debug('Formatting of outputs...')

        solver_output_transposed = self._solver_output.T

        all_plants_df = pd.DataFrame(columns=Simulation.PLANTS_FORMAT_OUTPUTS)
        all_axes_df = pd.DataFrame(columns=Simulation.AXES_FORMAT_OUTPUTS)
        all_phytomers_df = pd.DataFrame(columns=Simulation.PHYTOMERS_FORMAT_OUTPUTS)
        all_organs_df = pd.DataFrame(columns=Simulation.ORGANS_FORMAT_OUTPUTS)
        all_elements_df = pd.DataFrame(columns=Simulation.ELEMENTS_FORMAT_OUTPUTS)

        for plant in self.population.plants:

            plants_df = pd.DataFrame(columns=all_plants_df.columns)
            plants_df['t'] = self._time_grid
            plants_df['plant'] = plant.index

            for axis in plant.axes:

                axes_df = pd.DataFrame(columns=all_axes_df.columns)
                axes_df['t'] = self._time_grid
                axes_df['plant'] = plant.index
                axes_df['axis'] = axis.id

                # compute the total transpiration
                total_transpiration = np.zeros_like(self._time_grid)

                transpiration_mapping = {}
                for phytomer in axis.phytomers:
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ is not None:
                            for element in (organ.exposed_element, organ.enclosed_element):
                                if element is not None:
                                    transpiration_mapping[element] = map(element.calculate_transpiration, [element.Tr] * len(self._time_grid), [element.green_area] * len(self._time_grid))
                                    total_transpiration += transpiration_mapping[element]

                axes_df['Total_transpiration'] = total_transpiration

                # format phloem outputs
                organs_df = pd.DataFrame(columns=all_organs_df.columns)
                organs_df['t'] = self._time_grid
                organs_df['plant'] = plant.index
                organs_df['axis'] = axis.id
                organs_df['phytomer'] = np.nan
                organs_df['organ'] = axis.phloem.__class__.__name__
                phloem_sucrose = solver_output_transposed[self.initial_conditions_mapping[axis.phloem]['sucrose']]
                organs_df['sucrose'] = phloem_sucrose
                phloem_amino_acids = solver_output_transposed[self.initial_conditions_mapping[axis.phloem]['amino_acids']]
                organs_df['amino_acids'] = phloem_amino_acids
                organs_df['Conc_Sucrose'] = axis.phloem.calculate_conc_sucrose(organs_df['sucrose'])
                organs_df['Conc_Amino_Acids'] = axis.phloem.calculate_conc_amino_acids(organs_df['amino_acids'])
                all_organs_df = all_organs_df.append(organs_df, ignore_index=True)

                # format soil output
                organs_df = pd.DataFrame(columns=all_organs_df.columns)
                organs_df['t'] = self._time_grid
                organs_df['plant'] = plant.index
                organs_df['axis'] = axis.id
                organs_df['organ'] = axis.soil.__class__.__name__
                organs_df['Conc_Nitrates'] = axis.soil.calculate_conc_nitrates(solver_output_transposed[self.initial_conditions_mapping[axis.soil]['nitrates']])
                all_organs_df = all_organs_df.append(organs_df, ignore_index=True)

                # format roots outputs
                organs_df = pd.DataFrame(columns=all_organs_df.columns)
                organs_df['t'] = self._time_grid
                organs_df['plant'] = plant.index
                organs_df['axis'] = axis.id
                organs_df['phytomer'] = np.nan
                organs_df['organ'] = axis.roots.__class__.__name__
                organs_df['sucrose'] = solver_output_transposed[self.initial_conditions_mapping[axis.roots]['sucrose']]
                organs_df['nitrates'] = solver_output_transposed[self.initial_conditions_mapping[axis.roots]['nitrates']]
                organs_df['amino_acids'] = solver_output_transposed[self.initial_conditions_mapping[axis.roots]['amino_acids']]
                organs_df['mstruct'] = axis.roots.mstruct
                organs_df['Nstruct'] = axis.roots.Nstruct
                organs_df['Conc_Sucrose'] = axis.roots.calculate_conc_sucrose(organs_df['sucrose'])
                organs_df['Conc_Nitrates'] = axis.roots.calculate_conc_nitrates(organs_df['nitrates'])
                organs_df['Conc_Amino_Acids'] = axis.roots.calculate_conc_amino_acids(organs_df['amino_acids'])
                organs_df['Unloading_Sucrose'] = map(axis.roots.calculate_unloading_sucrose, phloem_sucrose)
                organs_df['Unloading_Amino_Acids'] = map(axis.roots.calculate_unloading_amino_acids, organs_df['Unloading_Sucrose'], phloem_sucrose, phloem_amino_acids)
                roots_uptake_nitrates, roots_potential_uptake_nitrates = axis.roots.calculate_uptake_nitrates(organs_df['Conc_Nitrates_Soil'], organs_df['nitrates'], total_transpiration)
                organs_df['Uptake_Nitrates'] = roots_uptake_nitrates
                organs_df['Potential_Uptake_Nitrates'] = map(axis.roots.calculate_export_amino_acids, organs_df['amino_acids'], total_transpiration)
                roots_export_amino_acids = map(axis.roots.calculate_export_amino_acids, organs_df['amino_acids'], total_transpiration)
                organs_df['Export_Amino_Acids'] = roots_export_amino_acids
                organs_df['S_Amino_Acids'] = map(axis.roots.calculate_s_amino_acids, organs_df['nitrates'], organs_df['sucrose'])
                organs_df['R_Nnit_upt'] = map(RespirationModel.R_Nnit_upt, organs_df['Uptake_Nitrates'], organs_df['sucrose'])
                organs_df['R_Nnit_red'] = map(RespirationModel.R_Nnit_red, organs_df['S_Amino_Acids'], organs_df['sucrose'], [axis.roots.mstruct*axis.roots.PARAMETERS.ALPHA] * len(self._time_grid), [True] * len(self._time_grid))
                total_nitrogen =  map(axis.roots.calculate_total_nitrogen, organs_df['nitrates'], organs_df['amino_acids'], [axis.roots.Nstruct] * len(self._time_grid))
                organs_df['total_nitrogen'] = total_nitrogen
                organs_df['R_residual'] = map(RespirationModel.R_residual, organs_df['sucrose'], [axis.roots.mstruct*axis.roots.PARAMETERS.ALPHA] * len(self._time_grid), total_nitrogen, [axis.roots.PARAMETERS.DELTA_T]*len(self._time_grid))
                organs_df['mstruct_growth'] = axis.roots.mstruct_C_growth
                organs_df['R_growth'] = map(RespirationModel.R_growth, organs_df['mstruct_growth'], [axis.roots.mstruct*axis.roots.PARAMETERS.ALPHA] * len(self._time_grid))
                organs_df['mstruct_senescence'] = axis.roots.mstruct_senescence
                organs_df['Nstruct_N_growth'] = axis.roots.Nstruct_N_growth
                exudation = np.array(map(axis.roots.calculate_exudation, organs_df['Unloading_Sucrose'], phloem_sucrose, phloem_amino_acids))
                organs_df['C_exudation'], organs_df['N_exudation'] = exudation[:,0], exudation[:,1]
                all_organs_df = all_organs_df.append(organs_df, ignore_index=True)

                # format photosynthetic organs elements outputs
                for phytomer in axis.phytomers:
                    phytomers_df = pd.DataFrame(columns=all_phytomers_df.columns)
                    phytomers_df['t'] = self._time_grid
                    phytomers_df['plant'] = plant.index
                    phytomers_df['axis'] = axis.id
                    phytomers_df['phytomer'] = phytomer.index

                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ is None:
                            continue
                        for element, element_type in ((organ.exposed_element, 'exposed'), (organ.enclosed_element, 'enclosed')):
                            if element is None:
                                continue

                            elements_df = pd.DataFrame(columns=all_elements_df.columns)
                            elements_df['t'] = self._time_grid
                            elements_df['plant'] = plant.index
                            elements_df['axis'] = axis.id
                            elements_df['phytomer'] = phytomer.index
                            elements_df['organ'] = organ.__class__.__name__
                            elements_df['element'] = element_type
                            elements_df['area'] = element.area
                            elements_df['green_area'] = element.green_area
                            elements_df['height'] = element.height
                            elements_df['width'] = element.width
                            elements_df['surfacic_nitrogen'] = element.surfacic_nitrogen
                            elements_df['relative_delta_green_area'] = element.relative_delta_green_area
                            elements_df['mstruct'] = element.mstruct
                            elements_df['Nstruct'] = element.Nstruct
                            elements_df['triosesP'] = solver_output_transposed[self.initial_conditions_mapping[element]['triosesP']]
                            elements_df['starch'] = solver_output_transposed[self.initial_conditions_mapping[element]['starch']]
                            elements_df['sucrose'] = solver_output_transposed[self.initial_conditions_mapping[element]['sucrose']]
                            elements_df['fructan'] = solver_output_transposed[self.initial_conditions_mapping[element]['fructan']]
                            elements_df['nitrates'] = solver_output_transposed[self.initial_conditions_mapping[element]['nitrates']]
                            elements_df['amino_acids'] = solver_output_transposed[self.initial_conditions_mapping[element]['amino_acids']]
                            elements_df['proteins'] = solver_output_transposed[self.initial_conditions_mapping[element]['proteins']]
                            elements_df['Loading_Sucrose'] = map(element.calculate_loading_sucrose, elements_df['sucrose'], phloem_sucrose)
                            elements_df['Regul_S_Fructan'] = map(element.calculate_regul_s_fructan, elements_df['Loading_Sucrose'])
                            elements_df['Ag'] = element.Ag
                            elements_df['An'] = element.An
                            elements_df['Rd'] = map(element.calculate_total_Rd, [element.Rd]*len(self._time_grid), [element.green_area] * len(self._time_grid))
                            elements_df['R_phloem_loading'] = map(RespirationModel.R_phloem, elements_df['Loading_Sucrose'], elements_df['sucrose'], [element.mstruct*element.__class__.PARAMETERS.ALPHA] * len(self._time_grid))
                            elements_df['S_Amino_Acids'] = map(element.calculate_s_amino_acids, elements_df['nitrates'], elements_df['triosesP'])
                            elements_df['R_Nnit_red'] = map(RespirationModel.R_Nnit_red, elements_df['S_Amino_Acids'], elements_df['sucrose'], [element.mstruct*element.__class__.PARAMETERS.ALPHA] * len(self._time_grid))
                            total_nitrogen = map(element.calculate_total_nitrogen, elements_df['nitrates'], elements_df['amino_acids'], elements_df['proteins'], [element.Nstruct] * len(self._time_grid))
                            elements_df['total_nitrogen'] = total_nitrogen
                            elements_df['R_residual'] = map(RespirationModel.R_residual, elements_df['sucrose'], [element.mstruct*element.__class__.PARAMETERS.ALPHA] * len(self._time_grid), total_nitrogen, [organ.__class__.PARAMETERS.DELTA_T] * len(self._time_grid))
                            elements_df['Tr'] = element.Tr
                            elements_df['Ts'] = element.Ts
                            elements_df['gs'] = element.gs
                            elements_df['Photosynthesis'] = map(element.calculate_photosynthesis, [element.Ag] * len(self._time_grid), [element.green_area] * len(self._time_grid))
                            elements_df['Transpiration'] = transpiration_mapping[element]
                            elements_df['Conc_TriosesP'] = element.calculate_conc_triosesP(elements_df['triosesP'])
                            elements_df['Conc_Starch'] = element.calculate_conc_starch(elements_df['starch'])
                            elements_df['Conc_Sucrose'] = element.calculate_conc_sucrose(elements_df['sucrose'])
                            elements_df['Conc_Fructan'] = element.calculate_conc_fructan(elements_df['fructan'])
                            elements_df['Conc_Nitrates'] = element.calculate_conc_nitrates(elements_df['nitrates'])
                            elements_df['Conc_Amino_Acids'] = element.calculate_conc_amino_acids(elements_df['amino_acids'])
                            elements_df['Conc_Proteins'] = element.calculate_conc_proteins(elements_df['proteins'])
                            elements_df['SLN'] = map(element.calculate_surfacic_nitrogen, elements_df['nitrates'], elements_df['amino_acids'], elements_df['proteins'], [element.Nstruct] * len(self._time_grid), [element.green_area] * len(self._time_grid))
                            elements_df['S_Starch'] = map(element.calculate_s_starch, elements_df['triosesP'])
                            elements_df['D_Starch'] = map(element.calculate_d_starch, elements_df['starch'])
                            elements_df['S_Sucrose'] = map(element.calculate_s_sucrose, elements_df['triosesP'])
                            elements_df['S_Fructan'] = map(element.calculate_s_fructan, elements_df['sucrose'], elements_df['Regul_S_Fructan'])
                            elements_df['D_Fructan'] = map(element.calculate_d_fructan, elements_df['sucrose'], elements_df['fructan'])
                            elements_df['Nitrates_import'] = map(element.calculate_nitrates_import, roots_uptake_nitrates, transpiration_mapping[element], total_transpiration)
                            elements_df['Amino_Acids_import'] = map(element.calculate_amino_acids_import, roots_export_amino_acids, transpiration_mapping[element], total_transpiration)
                            elements_df['S_Proteins'] = map(element.calculate_s_proteins, elements_df['amino_acids'])
                            elements_df['D_Proteins'] = map(element.calculate_d_proteins, elements_df['proteins'])
                            elements_df['Loading_Amino_Acids'] = map(element.calculate_loading_amino_acids, elements_df['amino_acids'], phloem_amino_acids)
                            elements_df['remob_starch_senescence'] = map(element.calculate_remob_starch_senescence, elements_df['starch'], [element.relative_delta_green_area] * len(self._time_grid))
                            elements_df['remob_fructan_senescence'] = map(element.calculate_remob_fructan_senescence, elements_df['fructan'], [element.relative_delta_green_area] * len(self._time_grid))
                            elements_df['remob_proteins_senescence'] = map(element.calculate_remob_proteins_senescence, elements_df['proteins'], [element.relative_delta_green_area] * len(self._time_grid))

                            all_elements_df = all_elements_df.append(elements_df, ignore_index=True)

                    all_phytomers_df = all_phytomers_df.append(phytomers_df, ignore_index=True)

                # format grains outputs
                if axis.grains is None:
                    continue

                organs_df = pd.DataFrame(columns=all_organs_df.columns)
                organs_df['t'] = self._time_grid
                organs_df['plant'] = plant.index
                organs_df['axis'] = axis.id
                organs_df['phytomer'] = np.nan
                organs_df['organ'] = axis.grains.__class__.__name__
                organs_df['structure'] = solver_output_transposed[self.initial_conditions_mapping[axis.grains]['structure']]
                organs_df['starch'] = solver_output_transposed[self.initial_conditions_mapping[axis.grains]['starch']]
                organs_df['proteins'] = solver_output_transposed[self.initial_conditions_mapping[axis.grains]['proteins']]
                organs_df['RGR_Structure'] = map(axis.grains.calculate_RGR_structure, phloem_sucrose)
                organs_df['S_grain_structure'] = map(axis.grains.calculate_s_grain_structure, self._time_grid, organs_df['structure'], organs_df['RGR_Structure'])
                structural_dry_mass = map(axis.grains.calculate_structural_dry_mass, organs_df['structure'])
                organs_df['S_grain_starch'] = map(axis.grains.calculate_s_grain_starch, self._time_grid, phloem_sucrose)
                organs_df['Dry_Mass'] = axis.grains.calculate_dry_mass(organs_df['structure'], organs_df['starch'], organs_df['proteins'])
                organs_df['Proteins_N_Mass'] = axis.grains.calculate_protein_mass(organs_df['proteins'])
                organs_df['Unloading_Sucrose'] = map(axis.grains.calculate_unloading_sucrose, organs_df['S_grain_structure'], organs_df['S_grain_starch'], structural_dry_mass)
                organs_df['S_Proteins'] = map(axis.grains.calculate_s_proteins, organs_df['S_grain_structure'], organs_df['S_grain_starch'], phloem_amino_acids, phloem_sucrose, structural_dry_mass)
                R_grain_growth =np.array(map(RespirationModel.R_grain_growth, organs_df['S_grain_structure'], organs_df['S_grain_starch'], structural_dry_mass))
                organs_df['R_grain_growth_struct'], organs_df['R_grain_growth_starch']  = R_grain_growth[:,0], R_grain_growth[:,1]

                all_organs_df = all_organs_df.append(organs_df, ignore_index=True)
                all_axes_df = all_axes_df.append(axes_df, ignore_index=True)

            all_plants_df = all_plants_df.append(plants_df, ignore_index=True)

        # set the order of the columns
        all_plants_df = all_plants_df.reindex_axis(Simulation.PLANTS_FORMAT_OUTPUTS, axis=1, copy=False)
        all_axes_df = all_axes_df.reindex_axis(Simulation.AXES_FORMAT_OUTPUTS, axis=1, copy=False)
        all_phytomers_df = all_phytomers_df.reindex_axis(Simulation.PHYTOMERS_FORMAT_OUTPUTS, axis=1, copy=False)
        all_organs_df = all_organs_df.reindex_axis(Simulation.ORGANS_FORMAT_OUTPUTS, axis=1, copy=False)
        all_elements_df = all_elements_df.reindex_axis(Simulation.ELEMENTS_FORMAT_OUTPUTS, axis=1, copy=False)

        # sort the rows by the columns
        all_plants_df.sort_index(by=Simulation.PLANTS_OUTPUTS_INDEXES, inplace=True)
        all_axes_df.sort_index(by=Simulation.AXES_OUTPUTS_INDEXES, inplace=True)
        all_phytomers_df.sort_index(by=Simulation.PHYTOMERS_OUTPUTS_INDEXES, inplace=True)
        all_organs_df.sort_index(by=Simulation.ORGANS_OUTPUTS_INDEXES, inplace=True)
        all_elements_df.sort_index(by=Simulation.ELEMENTS_OUTPUTS_INDEXES, inplace=True)

        # infer the right types of the columns
        all_plants_df = all_plants_df.convert_objects(copy=False)
        all_axes_df = all_axes_df.convert_objects(copy=False)
        all_phytomers_df = all_phytomers_df.convert_objects(copy=False)
        all_organs_df = all_organs_df.convert_objects(copy=False)
        all_elements_df = all_elements_df.convert_objects(copy=False)

        # convert the indexes of plants, phytomers and elements to integers
        all_plants_df['plant'] = all_plants_df['plant'].astype(int)
        all_axes_df['plant'] = all_axes_df['plant'].astype(int)
        all_phytomers_df[['plant', 'phytomer']] = all_phytomers_df[['plant', 'phytomer']].astype(int)
        all_organs_df['plant'] = all_organs_df['plant'].astype(int)
        all_elements_df[['plant', 'phytomer']] = all_elements_df[['plant', 'phytomer']].astype(int)

        logger.debug('Formatting of outputs DONE')

        return all_plants_df, all_axes_df, all_phytomers_df, all_organs_df, all_elements_df


class ProgressBarError(Exception): pass

class ProgressBar(object):
    """
    Display a console progress bar.
    """

    def __init__(self, bar_length=20, title='', block_character='#', uncomplete_character='-'):
        if bar_length <= 0:
            raise ProgressBarError('bar_length <= 0')
        self.bar_length = bar_length #: the number of blocks in the progress bar. MUST BE GREATER THAN ZERO !
        self.t_max = 1 #: the maximum t that the progress bar can display. MUST BE GREATER THAN ZERO !
        self.block_interval = 1 #: the time interval of each block. MUST BE GREATER THAN ZERO !
        self.last_upper_t = 0 #: the last upper t displayed by the progress bar
        self.progress_mapping = {} #: a mapping to optimize the refresh rate
        self.title = title #: the title to write on the left side of the progress bar
        self.block_character = block_character #: the character to represent a block
        self.uncomplete_character = uncomplete_character #: the character to represent the uncompleted part of the progress bar

    def set_t_max(self, t_max):
        """"Set :attr:`t_max` and update other attributes accordingly.
        """
        if t_max <= 0:
            raise ProgressBarError('t_max <= 0')
        self.t_max = t_max
        self.block_interval = self.t_max / self.bar_length
        self.last_upper_t = 0
        self.progress_mapping.clear()

    def update(self, t):
        """Update the progress bar if needed.
        """
        t = min(t, self.t_max)
        if t < self.last_upper_t:
            return
        else:
            self.last_upper_t = t
        t_inf = t // self.block_interval * self.block_interval
        if t_inf not in self.progress_mapping:
            progress = t / self.t_max
            block = int(round(self.bar_length * progress))
            text = "\r{0}: [{1}] {2:>5d}% ".format(self.title,
                                                  self.block_character * block + self.uncomplete_character * (self.bar_length - block),
                                                  int(progress*100))
            self.progress_mapping[t_inf] = text
            sys.stdout.write(self.progress_mapping[t_inf])
            sys.stdout.flush()

