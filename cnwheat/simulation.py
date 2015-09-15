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

import logging

import numpy as np
import pandas as pd
from scipy.integrate import odeint
from respiwheat.model import RespirationModel

import model, tools

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
                                model.Organ: ['sucrose', 'amino_acids', 'nitrates', 'structure', 'starch', 'proteins', 'mstruct', 'Nstruct', 'cytokinines'],
                                model.PhotosyntheticOrganElement: ['nitrates', 'starch', 'amino_acids', 'proteins',
                                                                   'sucrose', 'triosesP', 'fructan', 'mstruct', 'Nstruct', 'cytokinines']}

    T_INDEX = 't'
    
    PLANTS_INPUTS_INDEXES = ['plant']
    PLANTS_OUTPUTS_INDEXES = [T_INDEX] + PLANTS_INPUTS_INDEXES
    PLANTS_STATE_PARAMETERS = []
    PLANTS_STATE = PLANTS_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.Plant, [])
    PLANTS_INTERMEDIATE_VARIABLES = []
    PLANTS_FLUXES = []
    PLANTS_INTEGRATIVE_VARIABLES = []
    PLANTS_RUN_VARIABLES = PLANTS_OUTPUTS_INDEXES + PLANTS_STATE + PLANTS_INTERMEDIATE_VARIABLES + PLANTS_FLUXES + PLANTS_INTEGRATIVE_VARIABLES
    PLANTS_POSTPROCESSING_VARIABLES = []
    PLANTS_ALL_VARIABLES = PLANTS_RUN_VARIABLES + PLANTS_POSTPROCESSING_VARIABLES

    AXES_INPUTS_INDEXES = ['plant', 'axis']
    AXES_OUTPUTS_INDEXES = [T_INDEX] + AXES_INPUTS_INDEXES
    AXES_STATE_PARAMETERS = []
    AXES_STATE = AXES_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.Axis, [])
    AXES_INTERMEDIATE_VARIABLES = []
    AXES_FLUXES = []
    AXES_INTEGRATIVE_VARIABLES = ['Total_transpiration']
    AXES_RUN_VARIABLES = AXES_OUTPUTS_INDEXES + AXES_STATE + AXES_INTERMEDIATE_VARIABLES + AXES_FLUXES + AXES_INTEGRATIVE_VARIABLES
    AXES_POSTPROCESSING_VARIABLES = []
    AXES_ALL_VARIABLES = AXES_RUN_VARIABLES + AXES_POSTPROCESSING_VARIABLES

    PHYTOMERS_INPUTS_INDEXES = ['plant', 'axis', 'phytomer']
    PHYTOMERS_OUTPUTS_INDEXES = [T_INDEX] + PHYTOMERS_INPUTS_INDEXES
    PHYTOMERS_STATE_PARAMETERS = []
    PHYTOMERS_STATE = PHYTOMERS_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.Phytomer, [])
    PHYTOMERS_INTERMEDIATE_VARIABLES = []
    PHYTOMERS_FLUXES = []
    PHYTOMERS_INTEGRATIVE_VARIABLES = []
    PHYTOMERS_RUN_VARIABLES = PHYTOMERS_OUTPUTS_INDEXES + PHYTOMERS_STATE + PHYTOMERS_INTERMEDIATE_VARIABLES + PHYTOMERS_FLUXES + PHYTOMERS_INTEGRATIVE_VARIABLES
    PHYTOMERS_POSTPROCESSING_VARIABLES = []
    PHYTOMERS_ALL_VARIABLES = PHYTOMERS_RUN_VARIABLES + PHYTOMERS_POSTPROCESSING_VARIABLES

    ORGANS_INPUTS_INDEXES = ['plant', 'axis', 'phytomer', 'organ']
    ORGANS_OUTPUTS_INDEXES = [T_INDEX] + ORGANS_INPUTS_INDEXES
    ORGANS_STATE_PARAMETERS = ['volume', 'Tsoil', 'mstruct_C_growth', 'Nstruct_N_growth']
    ORGANS_STATE = ORGANS_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.Organ, [])
    ORGANS_INTERMEDIATE_VARIABLES = ['Conc_Nitrates_Soil', 'RGR_Structure', 'R_Nnit_upt', 'R_Nnit_red', 'R_residual', 'R_maintenance', 'R_grain_growth_struct', 'R_grain_growth_starch', 'R_growth',
                                     'C_exudation', 'N_exudation', 'regul_transpiration', 'HATS_LATS']
    ORGANS_FLUXES = ['Unloading_Sucrose', 'Export_Nitrates', 'Export_Amino_Acids', 'S_Proteins', 'S_Amino_Acids', 'Unloading_Amino_Acids',
                     'S_grain_starch', 'Uptake_Nitrates', 'S_grain_structure', 'S_cytokinines', 'D_cytokinines', 'Export_cytokinines']
    ORGANS_INTEGRATIVE_VARIABLES = ['total_organic_nitrogen']
    ORGANS_RUN_VARIABLES = ORGANS_OUTPUTS_INDEXES + ORGANS_STATE + ORGANS_INTERMEDIATE_VARIABLES + ORGANS_FLUXES + ORGANS_INTEGRATIVE_VARIABLES
    ORGANS_POSTPROCESSING_VARIABLES = ['Conc_Nitrates', 'Conc_Amino_Acids', 'Dry_Mass', 'Conc_Sucrose', 'Proteins_N_Mass', 'Conc_cytokinines']
    ORGANS_ALL_VARIABLES = ORGANS_RUN_VARIABLES + ORGANS_POSTPROCESSING_VARIABLES

    ELEMENTS_INPUTS_INDEXES = ['plant', 'axis', 'phytomer', 'organ', 'element']
    ELEMENTS_OUTPUTS_INDEXES = [T_INDEX] + ELEMENTS_INPUTS_INDEXES
    ELEMENTS_STATE_PARAMETERS = ['green_area', 'Ag', 'Rd', 'Tr', 'Ts']
    ELEMENTS_STATE = ELEMENTS_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.PhotosyntheticOrganElement, [])
    ELEMENTS_INTERMEDIATE_VARIABLES = ['Transpiration', 'R_phloem_loading', 'R_Nnit_red', 'R_residual', 'R_maintenance', 'Photosynthesis']
    ELEMENTS_FLUXES = ['Loading_Sucrose', 'Regul_S_Fructan', 'S_Starch', 'D_Starch', 'S_Sucrose', 'S_Fructan', 'D_Fructan']
    ELEMENTS_INTEGRATIVE_VARIABLES = ['total_organic_nitrogen']
    ELEMENTS_RUN_VARIABLES = ELEMENTS_OUTPUTS_INDEXES + ELEMENTS_STATE + ELEMENTS_INTERMEDIATE_VARIABLES + ELEMENTS_FLUXES + ELEMENTS_INTEGRATIVE_VARIABLES
    ELEMENTS_POSTPROCESSING_VARIABLES = ['Conc_TriosesP', 'Conc_Starch', 'Conc_Sucrose', 'Conc_Fructan', 'Conc_Nitrates', 'Conc_Amino_Acids', 'Conc_Proteins',
                                         'SLN', 'Nitrates_import', 'Amino_Acids_import', 'S_Amino_Acids', 'S_Proteins', 'D_Proteins', 'Loading_Amino_Acids', 
                                         'Conc_cytokinines', 'D_cytokinines', 'cytokinines_import', 'k_proteins', 'Rd_total']
    ELEMENTS_ALL_VARIABLES = ELEMENTS_RUN_VARIABLES + ELEMENTS_POSTPROCESSING_VARIABLES

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

        self.progressbar = tools.ProgressBar(title='Solver progress') #: progress bar to show the progress of the solver
        self.show_progressbar = False #: True: show the progress bar ; False: DO NOT show the progress bar


    def initialize(self, population):
        """
        Initialize:
            
            * :attr:`population`, 
            * :attr:`initial_conditions_mapping`, 
            * and :attr:`initial_conditions`
        
        from `population`.
        
        :Parameters:
            
            - `population` (:class:`model.Population`) - a population of plants.
            
        """
        
        logger = logging.getLogger(__name__)

        logger.info('Initialization of the simulation...')
        
        # clean the attributes of the simulation
        del self.population.plants[:]
        del self.initial_conditions[:]
        self.initial_conditions_mapping.clear()
        
        self.population.plants.extend(population.plants)
            
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
                axis.roots.sucrose = y[self.initial_conditions_mapping[axis.roots]['sucrose']]

                # test
                axis.roots.cytokinines = y[self.initial_conditions_mapping[axis.roots]['cytokinines']]

                axis.soil.nitrates = y[self.initial_conditions_mapping[axis.soil]['nitrates']]
                phloem_contributors.append(axis.roots)

                # compute the total transpiration at t_inf
                total_transpiration = 0.0 # mmol s-1
                total_green_area = 0.0 # m2
                transpiration_mapping = {}
                for phytomer in axis.phytomers:
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ is not None:
                            for element in (organ.exposed_element, organ.enclosed_element):
                                if element is not None:
                                    transpiration_mapping[element] = element.calculate_total_transpiration(element.Tr, element.green_area)
                                    total_transpiration += transpiration_mapping[element]
                                    total_green_area += element.green_area
                                    phloem_contributors.append(element)
                total_surfacic_transpiration = total_transpiration/total_green_area #: total transpiration rate of plant per unit area (mmol m-2 s-1)

                # compute the flows from/to the roots to/from photosynthetic organs
                conc_nitrates_soil = axis.soil.calculate_conc_nitrates(axis.soil.nitrates)
                roots_uptake_nitrate, _, _ = axis.roots.calculate_uptake_nitrates(conc_nitrates_soil, axis.roots.nitrates, total_surfacic_transpiration, axis.roots.sucrose)
                R_Nnit_upt = RespirationModel.R_Nnit_upt(roots_uptake_nitrate, axis.roots.sucrose)
                roots_export_nitrates = axis.roots.calculate_export_nitrates(axis.roots.nitrates, total_surfacic_transpiration)
                roots_export_amino_acids = axis.roots.calculate_export_amino_acids(axis.roots.amino_acids, total_surfacic_transpiration)

                # test
                roots_export_cytokinines = axis.roots.calculate_export_cytokinines(axis.roots.cytokinines, total_surfacic_transpiration)

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

                            # test
                            element.cytokinines = y[self.initial_conditions_mapping[element]['cytokinines']]

                            # intermediate variables
                            photosynthesis = element.calculate_photosynthesis(element.Ag, element.green_area)
                            element_transpiration = transpiration_mapping[element] # mmol s-1

                            # flows
                            s_starch = element.calculate_s_starch(element.triosesP)
                            d_starch = element.calculate_d_starch(element.starch)
                            s_sucrose = element.calculate_s_sucrose(element.triosesP)
                            element.loading_sucrose = element.calculate_loading_sucrose(element.sucrose, axis.phloem.sucrose)
                            R_phloem_loading, element.loading_sucrose = RespirationModel.R_phloem(element.loading_sucrose, element.sucrose, element.mstruct*element.__class__.PARAMETERS.ALPHA)
                            regul_s_fructan = element.calculate_regul_s_fructan(element.loading_sucrose)
                            d_fructan = element.calculate_d_fructan(element.sucrose, element.fructan)
                            s_fructan = element.calculate_s_fructan(element.sucrose, regul_s_fructan)
                            nitrates_import = element.calculate_nitrates_import(roots_uptake_nitrate, roots_export_nitrates, element_transpiration, total_transpiration)
                            amino_acids_import = element.calculate_amino_acids_import(roots_export_amino_acids, element_transpiration, total_transpiration)
                            s_amino_acids = element.calculate_s_amino_acids(element.nitrates, element.triosesP)
                            R_Nnit_red, s_amino_acids = RespirationModel.R_Nnit_red(s_amino_acids, element.sucrose, element.mstruct*element.__class__.PARAMETERS.ALPHA)
                            s_proteins = element.calculate_s_proteins(element.amino_acids)
                            k, d_proteins = element.calculate_d_proteins(element.proteins, element.cytokinines)
                            element.loading_amino_acids = element.calculate_loading_amino_acids(element.amino_acids, axis.phloem.amino_acids)

                            # test
                            cytokinines_import = element.calculate_cytokinines_import(roots_export_cytokinines, element_transpiration, total_transpiration)
                            d_cytokinines = element.calculate_d_cytokinines(element.cytokinines)
                            cytokinines_derivative = element.calculate_cytokinines_derivative(cytokinines_import, d_cytokinines)


                            # compartments derivatives
                            starch_derivative = element.calculate_starch_derivative(s_starch, d_starch)
                            element.total_organic_nitrogen = element.calculate_total_organic_nitrogen(element.amino_acids, element.proteins, element.Nstruct)
                            R_residual,_ = RespirationModel.R_residual(element.sucrose, element.mstruct*element.__class__.PARAMETERS.ALPHA, element.total_organic_nitrogen, organ.__class__.PARAMETERS.DELTA_T, element.Ts)
                            sum_respi = R_phloem_loading + R_Nnit_red + R_residual
                            sucrose_derivative = element.calculate_sucrose_derivative(s_sucrose, d_starch, element.loading_sucrose, s_fructan, d_fructan, sum_respi)
                            triosesP_derivative = element.calculate_triosesP_derivative(photosynthesis, s_sucrose, s_starch, s_amino_acids)
                            fructan_derivative = element.calculate_fructan_derivative(s_fructan, d_fructan)
                            nitrates_derivative = element.calculate_nitrates_derivative(nitrates_import, s_amino_acids)
                            amino_acids_derivative = element.calculate_amino_acids_derivative(amino_acids_import, s_amino_acids, s_proteins, d_proteins, element.loading_amino_acids)
                            proteins_derivative = element.calculate_proteins_derivative(s_proteins, d_proteins)

                            y_derivatives[self.initial_conditions_mapping[element]['starch']] = starch_derivative
                            y_derivatives[self.initial_conditions_mapping[element]['sucrose']] = sucrose_derivative
                            y_derivatives[self.initial_conditions_mapping[element]['triosesP']] = triosesP_derivative
                            y_derivatives[self.initial_conditions_mapping[element]['fructan']] = fructan_derivative
                            y_derivatives[self.initial_conditions_mapping[element]['nitrates']] = nitrates_derivative
                            y_derivatives[self.initial_conditions_mapping[element]['amino_acids']] = amino_acids_derivative
                            y_derivatives[self.initial_conditions_mapping[element]['proteins']] = proteins_derivative
                            # test
                            y_derivatives[self.initial_conditions_mapping[element]['cytokinines']] = cytokinines_derivative

                if axis.grains is not None:
                    phloem_contributors.append(axis.grains)
                    # compute the derivative of each compartment of grains
                    axis.grains.structure = y[self.initial_conditions_mapping[axis.grains]['structure']]
                    axis.grains.starch = y[self.initial_conditions_mapping[axis.grains]['starch']]
                    axis.grains.proteins = y[self.initial_conditions_mapping[axis.grains]['proteins']]

                    # intermediate variables
                    RGR_structure = axis.grains.calculate_RGR_structure(axis.phloem.sucrose)
                    axis.grains.structural_dry_mass = axis.grains.calculate_structural_dry_mass(axis.grains.structure)

                    # flows
                    axis.grains.s_grain_structure = axis.grains.calculate_s_grain_structure(t, axis.grains.structure, RGR_structure)
                    axis.grains.s_grain_starch = axis.grains.calculate_s_grain_starch(t, axis.phloem.sucrose)
                    axis.grains.s_proteins = axis.grains.calculate_s_proteins(axis.grains.s_grain_structure, axis.grains.s_grain_starch, axis.phloem.amino_acids, axis.phloem.sucrose, axis.grains.structural_dry_mass)
                    # compartments derivatives
                    R_grain_growth_struct, R_grain_growth_starch = RespirationModel.R_grain_growth(axis.grains.s_grain_structure, axis.grains.s_grain_starch, axis.grains.structural_dry_mass)
                    structure_derivative = axis.grains.calculate_structure_derivative(axis.grains.s_grain_structure, R_grain_growth_struct)
                    starch_derivative = axis.grains.calculate_starch_derivative(axis.grains.s_grain_starch, axis.grains.structural_dry_mass, R_grain_growth_starch)
                    proteins_derivative = axis.grains.calculate_proteins_derivative(axis.grains.s_proteins)
                    y_derivatives[self.initial_conditions_mapping[axis.grains]['structure']] = structure_derivative
                    y_derivatives[self.initial_conditions_mapping[axis.grains]['starch']] = starch_derivative
                    y_derivatives[self.initial_conditions_mapping[axis.grains]['proteins']] = proteins_derivative

                # compute the derivative of each compartment of roots and soil
                y_derivatives[self.initial_conditions_mapping[axis.soil]['nitrates']] = axis.soil.calculate_nitrates_derivative(roots_uptake_nitrate)

                mstruct_C_growth = axis.roots.mstruct_C_growth
                Nstruct_N_growth = axis.roots.Nstruct_N_growth

                # flows
                axis.roots.unloading_sucrose = axis.roots.calculate_unloading_sucrose(axis.phloem.sucrose)
                axis.roots.unloading_amino_acids = axis.roots.calculate_unloading_amino_acids(axis.roots.unloading_sucrose, axis.phloem.sucrose, axis.phloem.amino_acids)
                axis.roots.s_amino_acids = axis.roots.calculate_s_amino_acids(axis.roots.nitrates, axis.roots.sucrose)
                R_Nnit_red, axis.roots.s_amino_acids = RespirationModel.R_Nnit_red(axis.roots.s_amino_acids, axis.roots.sucrose, axis.roots.mstruct*model.Roots.PARAMETERS.ALPHA, root=True)
                C_exudated, N_exudated = axis.roots.calculate_exudation(axis.roots.unloading_sucrose, axis.roots.sucrose, axis.phloem.sucrose, axis.roots.amino_acids, axis.phloem.amino_acids)

                # test
                s_cytokinines = axis.roots.calculate_s_cytokinines(axis.roots.sucrose)
                d_cytokinines =  axis.roots.calculate_d_cytokinines(axis.roots.cytokinines)
                cytokinines_derivative = axis.roots.calculate_cytokinines_derivative(s_cytokinines, d_cytokinines, roots_export_cytokinines)

                # compartments derivatives
                axis.roots.total_organic_nitrogen = axis.roots.calculate_total_organic_nitrogen(axis.roots.amino_acids, axis.roots.Nstruct)
                R_residual,_ = RespirationModel.R_residual(axis.roots.sucrose, axis.roots.mstruct*model.Roots.PARAMETERS.ALPHA, axis.roots.total_organic_nitrogen, axis.roots.PARAMETERS.DELTA_T, axis.soil.Tsoil)
                R_roots_growth = RespirationModel.R_growth(mstruct_C_growth, axis.roots.mstruct)
                sum_respi = R_Nnit_upt + R_Nnit_red + R_residual + R_roots_growth
                sucrose_derivative = axis.roots.calculate_sucrose_derivative(axis.roots.unloading_sucrose, axis.roots.s_amino_acids, mstruct_C_growth, C_exudated, sum_respi)
                nitrates_derivative = axis.roots.calculate_nitrates_derivative(roots_uptake_nitrate, roots_export_nitrates, axis.roots.s_amino_acids)
                amino_acids_derivative = axis.roots.calculate_amino_acids_derivative(axis.roots.unloading_amino_acids, axis.roots.s_amino_acids, roots_export_amino_acids, Nstruct_N_growth, N_exudated)
                y_derivatives[self.initial_conditions_mapping[axis.roots]['sucrose']] = sucrose_derivative
                y_derivatives[self.initial_conditions_mapping[axis.roots]['nitrates']] = nitrates_derivative
                y_derivatives[self.initial_conditions_mapping[axis.roots]['amino_acids']] = amino_acids_derivative
                #test
                y_derivatives[self.initial_conditions_mapping[axis.roots]['cytokinines']] = cytokinines_derivative

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
        
    
    def postprocessings(self):
        """
        Compute:
        
            * intermediate variables (see :attr:`Simulation:PLANTS_INTERMEDIATE_VARIABLES`, :attr:`Simulation:AXES_INTERMEDIATE_VARIABLES`, :attr:`Simulation:PHYTOMERS_INTERMEDIATE_VARIABLES`, :attr:`Simulation:ORGANS_INTERMEDIATE_VARIABLES` and :attr:`Simulation:ELEMENTS_INTERMEDIATE_VARIABLES`),
            * fluxes (see :attr:`Simulation:PLANTS_FLUXES`, :attr:`Simulation:AXES_FLUXES`, :attr:`Simulation:PHYTOMERS_FLUXES`, :attr:`Simulation:ORGANS_FLUXES` and :attr:`Simulation:ELEMENTS_FLUXES`),
            * and integrative variables (see :attr:`Simulation:PLANTS_INTEGRATIVE_VARIABLES`, :attr:`Simulation:AXES_INTEGRATIVE_VARIABLES`, :attr:`Simulation:PHYTOMERS_INTEGRATIVE_VARIABLES`, :attr:`Simulation:ORGANS_INTEGRATIVE_VARIABLES` and :attr:`Simulation:ELEMENTS_INTEGRATIVE_VARIABLES`),
            
        from :attr:`_solver_output` and format them to :class:`dataframes <pandas.DataFrame>`.
        
        :Returns:
            :class:`dataframes <pandas.DataFrame>` of post-processing outputs at each scale: 

                * plant (see :attr:`Simulation:PLANTS_ALL_VARIABLES`) 
                * axis (see :attr:`Simulation:AXES_ALL_VARIABLES`)
                * phytomer (see :attr:`Simulation:PHYTOMERS_ALL_VARIABLES`)
                * organ (see :attr:`Simulation:ORGANS_ALL_VARIABLES`)
                * and element (see :attr:`Simulation:ELEMENTS_ALL_VARIABLES`)
                
        :Returns Type:
            :class:`tuple` of :class:`pandas.DataFrame`
        
        """
        logger = logging.getLogger(__name__)
        logger.debug('Formatting of outputs...')

        solver_output_transposed = self._solver_output.T

        all_plants_df = pd.DataFrame(columns=Simulation.PLANTS_ALL_VARIABLES)
        all_axes_df = pd.DataFrame(columns=Simulation.AXES_ALL_VARIABLES)
        all_phytomers_df = pd.DataFrame(columns=Simulation.PHYTOMERS_ALL_VARIABLES)
        all_organs_df = pd.DataFrame(columns=Simulation.ORGANS_ALL_VARIABLES)
        all_elements_df = pd.DataFrame(columns=Simulation.ELEMENTS_ALL_VARIABLES)

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
                total_green_area = np.zeros_like(self._time_grid)
                transpiration_mapping = {}
                for phytomer in axis.phytomers:
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ is not None:
                            for element in (organ.exposed_element, organ.enclosed_element):
                                if element is not None:
                                    transpiration_mapping[element] = map(element.calculate_total_transpiration, [element.Tr] * len(self._time_grid), [element.green_area] * len(self._time_grid))
                                    total_transpiration += transpiration_mapping[element]
                                    total_green_area += [element.green_area] * len(self._time_grid)
                total_surfacic_transpiration = total_transpiration / total_green_area

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
                conc_nitrates_soil = axis.soil.calculate_conc_nitrates(solver_output_transposed[self.initial_conditions_mapping[axis.soil]['nitrates']])
                organs_df['Conc_Nitrates'] = conc_nitrates_soil
                Tsoil = axis.soil.Tsoil
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
                uptake = np.array(map(axis.roots.calculate_uptake_nitrates, conc_nitrates_soil, organs_df['nitrates'], total_surfacic_transpiration, organs_df['sucrose']))
                roots_uptake_nitrates, organs_df['HATS_LATS'], organs_df['regul_transpiration'] = uptake[:,0], uptake[:,1], uptake[:,2]
                organs_df['Uptake_Nitrates'] = roots_uptake_nitrates
                roots_export_nitrates = map(axis.roots.calculate_export_nitrates, organs_df['nitrates'], total_surfacic_transpiration)
                organs_df['Export_Nitrates'] = roots_export_nitrates
                roots_export_amino_acids = map(axis.roots.calculate_export_amino_acids, organs_df['amino_acids'], total_surfacic_transpiration)
                organs_df['Export_Amino_Acids'] = roots_export_amino_acids
                organs_df['S_Amino_Acids'] = map(axis.roots.calculate_s_amino_acids, organs_df['nitrates'], organs_df['sucrose'])
                organs_df['R_Nnit_upt'] = map(RespirationModel.R_Nnit_upt, organs_df['Uptake_Nitrates'], organs_df['sucrose'])
                output_respi_model = np.array(map(RespirationModel.R_Nnit_red, organs_df['S_Amino_Acids'], organs_df['sucrose'], [axis.roots.mstruct*axis.roots.PARAMETERS.ALPHA] * len(self._time_grid), [True] * len(self._time_grid)))
                organs_df['R_Nnit_red'], organs_df['S_Amino_Acids'] = output_respi_model[:,0], output_respi_model[:,1]
                total_organic_nitrogen =  map(axis.roots.calculate_total_organic_nitrogen, organs_df['amino_acids'], [axis.roots.Nstruct] * len(self._time_grid))
                organs_df['total_organic_nitrogen'] = total_organic_nitrogen
                R_residual = np.array(map(RespirationModel.R_residual, organs_df['sucrose'], [axis.roots.mstruct*axis.roots.PARAMETERS.ALPHA] * len(self._time_grid), total_organic_nitrogen, [axis.roots.PARAMETERS.DELTA_T]*len(self._time_grid), [Tsoil] * len(self._time_grid)))
                organs_df['R_residual'], organs_df['R_maintenance'] = R_residual[:,0], R_residual[:,1]
                organs_df['mstruct_C_growth'] = axis.roots.mstruct_C_growth
                organs_df['R_growth'] = map(RespirationModel.R_growth, organs_df['mstruct_C_growth'], [axis.roots.mstruct*axis.roots.PARAMETERS.ALPHA] * len(self._time_grid))
                organs_df['Nstruct_N_growth'] = axis.roots.Nstruct_N_growth
                exudation = np.array(map(axis.roots.calculate_exudation, organs_df['Unloading_Sucrose'], organs_df['sucrose'], phloem_sucrose, organs_df['amino_acids'], phloem_amino_acids))
                organs_df['C_exudation'], organs_df['N_exudation'] = exudation[:,0], exudation[:,1]

                # test
                organs_df['cytokinines'] = solver_output_transposed[self.initial_conditions_mapping[axis.roots]['cytokinines']]
                organs_df['Conc_cytokinines'] = axis.roots.calculate_conc_cytokinines(organs_df['cytokinines'])
                organs_df['S_cytokinines'] = map(axis.roots.calculate_s_cytokinines, organs_df['sucrose'])
                organs_df['D_cytokinines'] = map(axis.roots.calculate_d_cytokinines, organs_df['cytokinines'])
                root_export_cytokinines = map(axis.roots.calculate_export_cytokinines, organs_df['cytokinines'], total_surfacic_transpiration)
                organs_df['Export_cytokinines'] = root_export_cytokinines

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
                            elements_df['green_area'] = element.green_area
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
                            output_respi_model =np.array(map(RespirationModel.R_phloem, elements_df['Loading_Sucrose'], elements_df['sucrose'], [element.mstruct*element.__class__.PARAMETERS.ALPHA] * len(self._time_grid)))
                            elements_df['R_phloem_loading'], elements_df['Loading_Sucrose'] = output_respi_model[:,0], output_respi_model[:,1]
                            elements_df['Regul_S_Fructan'] = map(element.calculate_regul_s_fructan, elements_df['Loading_Sucrose'])
                            elements_df['Ag'] = element.Ag
                            elements_df['Rd'] = element.Rd
                            elements_df['Rd_total'] = map(element.calculate_total_Rd, [element.Rd]*len(self._time_grid), [element.green_area] * len(self._time_grid))
                            elements_df['S_Amino_Acids'] = map(element.calculate_s_amino_acids, elements_df['nitrates'], elements_df['triosesP'])
                            output_respi_model = np.array(map(RespirationModel.R_Nnit_red, elements_df['S_Amino_Acids'], elements_df['sucrose'], [element.mstruct*element.__class__.PARAMETERS.ALPHA] * len(self._time_grid)))
                            elements_df['R_Nnit_red'], elements_df['S_Amino_Acids'] = output_respi_model[:,0], output_respi_model[:,1]
                            total_organic_nitrogen = map(element.calculate_total_organic_nitrogen, elements_df['amino_acids'], elements_df['proteins'], [element.Nstruct] * len(self._time_grid))
                            elements_df['total_organic_nitrogen'] = total_organic_nitrogen
                            R_residual = np.array(map(RespirationModel.R_residual, elements_df['sucrose'], [element.mstruct*element.__class__.PARAMETERS.ALPHA] * len(self._time_grid), total_organic_nitrogen, [organ.__class__.PARAMETERS.DELTA_T] * len(self._time_grid), [element.Ts] * len(self._time_grid)))
                            elements_df['R_residual'],  elements_df['R_maintenance'] = R_residual[:,0], R_residual[:,1]
                            elements_df['Tr'] = element.Tr
                            elements_df['Ts'] = element.Ts
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
                            elements_df['Nitrates_import'] = map(element.calculate_nitrates_import, roots_uptake_nitrates, roots_export_nitrates, transpiration_mapping[element], total_transpiration)
                            elements_df['Amino_Acids_import'] = map(element.calculate_amino_acids_import, roots_export_amino_acids, transpiration_mapping[element], total_transpiration)
                            elements_df['S_Proteins'] = map(element.calculate_s_proteins, elements_df['amino_acids'])
                            d_proteins = np.array(map(element.calculate_d_proteins, elements_df['proteins'], solver_output_transposed[self.initial_conditions_mapping[element]['cytokinines']]))
                            elements_df['k_proteins'] = d_proteins[:,0]
                            elements_df['D_Proteins'] = d_proteins[:,1]
                            elements_df['Loading_Amino_Acids'] = map(element.calculate_loading_amino_acids, elements_df['amino_acids'], phloem_amino_acids)

                            # test
                            elements_df['cytokinines'] = solver_output_transposed[self.initial_conditions_mapping[element]['cytokinines']]
                            elements_df['cytokinines_import'] = map(element.calculate_cytokinines_import, root_export_cytokinines, transpiration_mapping[element], total_transpiration)
                            elements_df['D_cytokinines'] = map(element.calculate_d_cytokinines, elements_df['cytokinines'])
                            elements_df['Conc_cytokinines'] = element.calculate_conc_cytokinines(elements_df['cytokinines'])



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
                organs_df['S_Proteins'] = map(axis.grains.calculate_s_proteins, organs_df['S_grain_structure'], organs_df['S_grain_starch'], phloem_amino_acids, phloem_sucrose, structural_dry_mass)
                R_grain_growth =np.array(map(RespirationModel.R_grain_growth, organs_df['S_grain_structure'], organs_df['S_grain_starch'], structural_dry_mass))
                organs_df['R_grain_growth_struct'], organs_df['R_grain_growth_starch']  = R_grain_growth[:,0], R_grain_growth[:,1]

                all_organs_df = all_organs_df.append(organs_df, ignore_index=True)
                all_axes_df = all_axes_df.append(axes_df, ignore_index=True)

            all_plants_df = all_plants_df.append(plants_df, ignore_index=True)

        # set the order of the columns
        all_plants_df = all_plants_df.reindex_axis(Simulation.PLANTS_ALL_VARIABLES, axis=1, copy=False)
        all_axes_df = all_axes_df.reindex_axis(Simulation.AXES_ALL_VARIABLES, axis=1, copy=False)
        all_phytomers_df = all_phytomers_df.reindex_axis(Simulation.PHYTOMERS_ALL_VARIABLES, axis=1, copy=False)
        all_organs_df = all_organs_df.reindex_axis(Simulation.ORGANS_ALL_VARIABLES, axis=1, copy=False)
        all_elements_df = all_elements_df.reindex_axis(Simulation.ELEMENTS_ALL_VARIABLES, axis=1, copy=False)

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

