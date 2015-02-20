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

import numpy as np
import pandas as pd
from scipy.integrate import odeint
import logging

import model

class CNWheatError(Exception): pass
class CNWheatInitError(Exception): pass
class CNWheatInputError(CNWheatError): pass
class CNWheatRunError(CNWheatError): pass

class CNWheat(object):
    """
    The CNWheat class permits to initialize and run the model.

    Use :meth:`run` to run the model.

    :Parameters:

        - `population` (:class:`list`) - the :class:`population <cnwheat.model.Population` in the system.

    """

    #: the name of the compartments attributes in the model. 
    MODEL_COMPARTMENTS_NAMES = {model.Plant: [],
                                model.Axis: [],
                                model.Phytomer: [],
                                model.Organ: ['nitrates', 'starch', 'amino_acids', 'proteins',
                                              'sucrose', 'triosesP', 'fructan', 'structure']}
    
    PLANTS_INDEXES = ['t', 'plant']
    PLANTS_OUTPUTS = PLANTS_INDEXES + MODEL_COMPARTMENTS_NAMES.get(model.Plant, [])
    
    AXES_INDEXES = ['t', 'plant', 'axis']
    AXES_OUTPUTS = AXES_INDEXES + MODEL_COMPARTMENTS_NAMES.get(model.Axis, []) + ['Total_transpiration']
    
    PHYTOMERS_INDEXES = ['t', 'plant', 'axis', 'phytomer']
    PHYTOMERS_OUTPUTS = PHYTOMERS_INDEXES + MODEL_COMPARTMENTS_NAMES.get(model.Phytomer, [])
    
    ORGANS_INDEXES = ['t', 'plant', 'axis', 'phytomer', 'organ']
    ORGANS_OUTPUTS = ORGANS_INDEXES + MODEL_COMPARTMENTS_NAMES.get(model.Organ, []) + ['Amino_Acids_import', 'Conc_C_Sucrose', 'Unloading_Sucrose', 'D_Starch', 'Tr', 
                      'Potential_Uptake_Nitrates', 'S_Proteins', 'S_Amino_Acids', 'S_Sucrose', 'Uptake_Nitrates', 'S_grain_structure', 'S_Fructan', 'Conc_TriosesP', 
                      'Conc_Nitrates', 'Conc_Proteins', 'Loading_Sucrose', 'D_Proteins', 'S_Starch', 'Conc_Starch', 'Nitrates_import', 'Unloading_Amino_Acids', 
                      'Conc_Nitrates_Soil', 'S_grain_starch', 'Loading_Amino_Acids', 'Photosynthesis', 'D_Fructan', 'Transpiration', 'Conc_Fructan', 
                      'Export_Amino_Acids', 'Conc_Amino_Acids', 'An', 'RGR_Structure', 'Conc_Sucrose', 'Regul_S_Fructan', 'Proteins_N_Mass', 'Dry_Mass']
    
    LOGGERS_NAMES = {'compartments': {model.Plant: 'cnwheat.compartments.plants',
                                      model.Axis: 'cnwheat.compartments.axes',
                                      model.Phytomer: 'cnwheat.compartments.phytomers',
                                      model.Organ: 'cnwheat.compartments.organs'},
                     'derivatives': {model.Plant: 'cnwheat.derivatives.plants',
                                     model.Axis: 'cnwheat.derivatives.axes',
                                     model.Phytomer: 'cnwheat.derivatives.phytomers',
                                     model.Organ: 'cnwheat.derivatives.organs'}}
    
    def __init__(self, population):

        logger = logging.getLogger(__name__)

        logger.info('Initialization of the simulation...')

        self.population = population #: the population used in the model

        # construct the list of initial conditions
        self.initial_conditions = [] #: the initial conditions of the compartments in the population
        self.initial_conditions_mapping = {} #: dictionary to map the compartments to their indexes in :attr:`initial_conditions`
        
        def _init_initial_conditions(model_object, i):
            if model_object is not None:
                class_ = model_object.__class__
                if issubclass(model_object.__class__, model.Organ):
                    class_ = model.Organ
                compartments_names = CNWheat.MODEL_COMPARTMENTS_NAMES[class_]
                self.initial_conditions_mapping[model_object] = {}
                for compartment_name in compartments_names:
                    if hasattr(model_object, compartment_name):
                        self.initial_conditions_mapping[model_object][compartment_name] = i
                        self.initial_conditions.append(0)
                        i += 1
            return i
        
        i = 0
        
        for plant in population.plants:
            i = _init_initial_conditions(plant, i)
            for axis in plant.axes:
                i = _init_initial_conditions(axis, i)
                for organ in (axis.roots, axis.phloem, axis.grains):
                    i = _init_initial_conditions(organ, i)
                for phytomer in axis.phytomers:
                    i = _init_initial_conditions(phytomer, i)
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                        i = _init_initial_conditions(organ, i)
        
        #TODO: check the validity of the population         

        self.progressbar = ProgressBar(title='Solver progress') #: progress bar to show the progress of the solver
        self.show_progressbar = False #: True: show the progress bar ; False: DO NOT show the progress bar
        logger.info('Initialization of the simulation DONE')


    def run(self, start_time, stop_time, number_of_output_steps, odeint_mxstep=5000, show_progressbar=False):
        """
        Compute CN exchanges between organs in :attr:`population`.
        The computation is done between `start_time` and `stop_time`, for `number_of_output_steps` steps.

        :Parameters:

            - `start_time` (:class:`int`) - The starting of the time grid.

            - `stop_time` (:class:`int`) - The end of the time grid.

            - `number_of_output_steps` (:class:`int`) - Number of time points for which to compute the CN exchanges in :attr:`population`.

            - `odeint_mxstep` (:class:`int`) - Maximum number of (internally defined) steps allowed for each integration point in time grid.
              `odeint_mxstep` is passed to :func:`scipy.integrate.odeint` as `mxstep`. If `odeint_mxstep` = 0, then `mxstep` is determined by the solver.
              The default value ( `5000` ) normally permits to solve the current model. User should increased this value if a more complex model is defined
              and if this model make the integration failed.

            - `show_progressbar` (:class:`bool`) - True: show the progress bar ; False: do not show the progress bar.

        :Returns:
            The :class:`dataframes <pandas.DataFrame>` of CN exchanges for each desired time step at different scales:
        
                * plant: t, plant index, outputs of plant,
                * axis: t, plant index, axis index, outputs of axis,
                * phytomer: t, plant index, axis index, phytomer index, outputs of phytomer,
                * and organ: t, plant index, axis index, phytomer index, organ name, outputs of organ.
        
        :Returns Type:
            :class:`tuple` of :class:`pandas.DataFrame`

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
                plants_compartments_logger.debug(sep.join(CNWheat.PLANTS_INDEXES + CNWheat.MODEL_COMPARTMENTS_NAMES[model.Plant]))
                axes_compartments_logger = logging.getLogger('cnwheat.compartments.axes')
                axes_compartments_logger.debug(sep.join(CNWheat.AXES_INDEXES + CNWheat.MODEL_COMPARTMENTS_NAMES[model.Axis]))
                phytomers_compartments_logger = logging.getLogger('cnwheat.compartments.phytomers')
                phytomers_compartments_logger.debug(sep.join(CNWheat.PHYTOMERS_INDEXES + CNWheat.MODEL_COMPARTMENTS_NAMES[model.Phytomer]))
                organs_compartments_logger = logging.getLogger('cnwheat.compartments.organs')
                organs_compartments_logger.debug(sep.join(CNWheat.ORGANS_INDEXES + CNWheat.MODEL_COMPARTMENTS_NAMES[model.Organ]))
            if derivatives_logger.isEnabledFor(logging.DEBUG):
                plants_derivatives_logger = logging.getLogger('cnwheat.derivatives.plants')
                plants_derivatives_logger.debug(sep.join(CNWheat.PLANTS_INDEXES + CNWheat.MODEL_COMPARTMENTS_NAMES[model.Plant]))
                axes_derivatives_logger = logging.getLogger('cnwheat.derivatives.axes')
                axes_derivatives_logger.debug(sep.join(CNWheat.AXES_INDEXES + CNWheat.MODEL_COMPARTMENTS_NAMES[model.Axis]))
                phytomers_derivatives_logger = logging.getLogger('cnwheat.derivatives.phytomers')
                phytomers_derivatives_logger.debug(sep.join(CNWheat.PHYTOMERS_INDEXES + CNWheat.MODEL_COMPARTMENTS_NAMES[model.Phytomer]))
                organs_derivatives_logger = logging.getLogger('cnwheat.derivatives.organs')
                organs_derivatives_logger.debug(sep.join(CNWheat.ORGANS_INDEXES + CNWheat.MODEL_COMPARTMENTS_NAMES[model.Organ]))
        
        self._update_initial_conditions() 
            
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(
                """Run the solver with:
                    - time grid = %s,
                    - odeint mxstep = %s""",
                t, odeint_mxstep)

        soln, infodict = odeint(self._calculate_all_derivatives, self.initial_conditions, t, full_output=True, mxstep=odeint_mxstep)
        
        if not set(infodict['mused']).issubset([1,2]):
            message = "Integration failed. See the logs of lsoda or try to increase the value of 'mxstep'."
            logger.exception(message)
            raise CNWheatRunError(message)
        
        last_compartments_values = soln[-1]
        self._update_population(last_compartments_values)
        
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(
                """Run of the solver DONE: infodict = %s""",
                infodict)
        
        all_plants_df, all_axes_df, all_phytomers_df, all_organs_df = self._format_solver_output(t, soln)
        
        logger.info('Run of CN-Wheat from {} to {} DONE'.format(start_time, stop_time))
        
        return all_plants_df, all_axes_df, all_phytomers_df, all_organs_df
    
    
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
            if model_object is not None:
                row = []
                class_ = model_object.__class__
                if issubclass(model_object.__class__, model.Organ):
                    class_ = model.Organ
                compartments_names = CNWheat.MODEL_COMPARTMENTS_NAMES[class_]
                for compartment_name in compartments_names:
                    if hasattr(model_object, compartment_name):
                        row.append(str(y[i]))
                        i += 1
                    else:
                        row.append('NA')
                rows.append([str(index) for index in indexes] + row) 
            return i
        
        i = 0
        all_rows = dict([(class_, []) for class_ in loggers_names])
        for plant in self.population.plants:
            i = update_rows(plant, [t, plant.index], all_rows[model.Plant], i)
            for axis in plant.axes:
                i = update_rows(axis, [t, plant.index, axis.index], all_rows[model.Axis], i)
                for organ in (axis.roots, axis.phloem, axis.grains):
                    i = update_rows(organ, [t, plant.index, axis.index, None, organ.__class__.__name__], all_rows[model.Organ], i)
                for phytomer in axis.phytomers:
                    i = update_rows(phytomer, [t, plant.index, axis.index, phytomer.index], all_rows[model.Phytomer], i)
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                        i = update_rows(organ, [t, plant.index, axis.index, phytomer.index, organ.__class__.__name__], all_rows[model.Organ], i)
        
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
            self._log_compartments(t, y, CNWheat.LOGGERS_NAMES['compartments'])

        # check that the solver is not crashed
        y_isnan = np.isnan(y)
        if y_isnan.any():
            message = 'The solver did not manage to compute a compartment. See the logs.'
            logger.exception(message)
            raise CNWheatRunError(message)

        y_derivatives = np.zeros_like(y)
        
        for plant in self.population.plants:
            for axis in plant.axes:
                
                axis_photosynthetic_organs = []
                
                axis.phloem.sucrose = y[self.initial_conditions_mapping[axis.phloem]['sucrose']]
                
                axis.phloem.amino_acids = y[self.initial_conditions_mapping[axis.phloem]['amino_acids']]
                
                axis.roots.nitrates = y[self.initial_conditions_mapping[axis.roots]['nitrates']]
                axis.roots.amino_acids = y[self.initial_conditions_mapping[axis.roots]['amino_acids']]
                axis_photosynthetic_organs.append(axis.roots)
                
                # compute the total transpiration at t_inf
                total_transpiration = 0.0
                transpiration_mapping = {}
                for phytomer in axis.phytomers:
                    phytomer_photosynthetic_organs = (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath)
                    for organ in phytomer_photosynthetic_organs:
                        if organ is not None:
                            transpiration_mapping[organ] = organ.calculate_transpiration(t, organ.Tr, phytomer.index)
                            total_transpiration += transpiration_mapping[organ]
                            axis_photosynthetic_organs.append(organ)

                # compute the flows from/to the roots to/from photosynthetic organs
                conc_nitrates_soil = axis.roots.calculate_conc_nitrates_soil(t)
                roots_uptake_nitrate, potential_roots_uptake_nitrates = axis.roots.calculate_uptake_nitrates(conc_nitrates_soil, axis.roots.nitrates, total_transpiration)
                roots_export_amino_acids = axis.roots.calculate_export_amino_acids(axis.roots.amino_acids, total_transpiration)

                # compute the derivative of each photosynthetic organ compartment
                for phytomer in axis.phytomers:
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                        
                        if organ is None:
                            continue
                        
                        organ.starch = y[self.initial_conditions_mapping[organ]['starch']]
                        organ.sucrose = y[self.initial_conditions_mapping[organ]['sucrose']]
                        organ.triosesP = y[self.initial_conditions_mapping[organ]['triosesP']]
                        organ.fructan = y[self.initial_conditions_mapping[organ]['fructan']]
                        organ.nitrates = y[self.initial_conditions_mapping[organ]['nitrates']]
                        organ.amino_acids = y[self.initial_conditions_mapping[organ]['amino_acids']]
                        organ.proteins = y[self.initial_conditions_mapping[organ]['proteins']]
            
                        # intermediate variables
                        photosynthesis = organ.calculate_photosynthesis(t, organ.An, phytomer.index)
                        organ_transpiration = transpiration_mapping[organ]
            
                        # flows
                        s_starch = organ.calculate_s_starch(organ.triosesP)
                        d_starch = organ.calculate_d_starch(organ.starch)
                        s_sucrose = organ.calculate_s_sucrose(organ.triosesP)
                        organ.loading_sucrose = organ.calculate_loading_sucrose(organ.sucrose, axis.phloem.sucrose)
                        regul_s_fructan = organ.calculate_regul_s_fructan(organ.loading_sucrose)
                        d_fructan = organ.calculate_d_fructan(organ.sucrose, organ.fructan)
                        s_fructan = organ.calculate_s_fructan(organ.sucrose, regul_s_fructan)
                        nitrates_import = organ.calculate_nitrates_import(roots_uptake_nitrate, organ_transpiration, total_transpiration)
                        amino_acids_import = organ.calculate_amino_acids_import(roots_export_amino_acids, organ_transpiration, total_transpiration)
                        s_amino_acids = organ.calculate_s_amino_acids(organ.nitrates, organ.triosesP)
                        s_proteins = organ.calculate_s_proteins(organ.amino_acids)
                        d_proteins = organ.calculate_d_proteins(organ.proteins)
                        organ.loading_amino_acids = organ.calculate_loading_amino_acids(organ.amino_acids, axis.phloem.amino_acids)
            
                        # compartments derivatives
                        starch_derivative = organ.calculate_starch_derivative(s_starch, d_starch)
                        sucrose_derivative = organ.calculate_sucrose_derivative(s_sucrose, d_starch, organ.loading_sucrose, s_fructan, d_fructan)
                        triosesP_derivative = organ.calculate_triosesP_derivative(photosynthesis, s_sucrose, s_starch, s_amino_acids)
                        fructan_derivative = organ.calculate_fructan_derivative(s_fructan, d_fructan)
                        nitrates_derivative = organ.calculate_nitrates_derivative(nitrates_import, s_amino_acids)
                        amino_acids_derivative = organ.calculate_amino_acids_derivative(amino_acids_import, s_amino_acids, s_proteins, d_proteins, organ.loading_amino_acids)
                        proteins_derivative = organ.calculate_proteins_derivative(s_proteins, d_proteins)
            
                        y_derivatives[self.initial_conditions_mapping[organ]['starch']] = starch_derivative
                        y_derivatives[self.initial_conditions_mapping[organ]['sucrose']] = sucrose_derivative
                        y_derivatives[self.initial_conditions_mapping[organ]['triosesP']] = triosesP_derivative
                        y_derivatives[self.initial_conditions_mapping[organ]['fructan']] = fructan_derivative
                        y_derivatives[self.initial_conditions_mapping[organ]['nitrates']] = nitrates_derivative
                        y_derivatives[self.initial_conditions_mapping[organ]['amino_acids']] = amino_acids_derivative
                        y_derivatives[self.initial_conditions_mapping[organ]['proteins']] = proteins_derivative

                if axis.grains is not None:
                    axis_photosynthetic_organs.append(axis.grains)
                    # compute the derivative of each compartment of grains
                    axis.grains.structure = y[self.initial_conditions_mapping[axis.grains]['structure']]
                    axis.grains.starch = y[self.initial_conditions_mapping[axis.grains]['starch']]
                    axis.grains.proteins = y[self.initial_conditions_mapping[axis.grains]['proteins']]
    
                    # intermediate variables
                    RGR_structure = axis.grains.calculate_RGR_structure(axis.phloem.sucrose)
    
                    # flows
                    axis.grains.s_grain_structure = axis.grains.calculate_s_grain_structure(t, axis.grains.structure, RGR_structure)
                    axis.grains.s_grain_starch = axis.grains.calculate_s_grain_starch(t, axis.phloem.sucrose)
                    axis.grains.s_proteins = axis.grains.calculate_s_proteins(axis.grains.s_grain_structure, axis.grains.s_grain_starch, axis.phloem.amino_acids, axis.phloem.sucrose, axis.grains.structure)
            
                    # compartments derivatives
                    structure_derivative = axis.grains.calculate_structure_derivative(axis.grains.s_grain_structure)
                    starch_derivative = axis.grains.calculate_starch_derivative(axis.grains.s_grain_starch, axis.grains.structure)
                    proteins_derivative = axis.grains.calculate_proteins_derivative(axis.grains.s_proteins)
                    y_derivatives[self.initial_conditions_mapping[axis.grains]['structure']] = structure_derivative
                    y_derivatives[self.initial_conditions_mapping[axis.grains]['starch']] = starch_derivative
                    y_derivatives[self.initial_conditions_mapping[axis.grains]['proteins']] = proteins_derivative
        
                # compute the derivative of each compartment of roots
                axis.roots.sucrose = y[self.initial_conditions_mapping[axis.roots]['sucrose']]
                # flows
                axis.roots.unloading_sucrose = axis.roots.calculate_unloading_sucrose(axis.phloem.sucrose)
                axis.roots.unloading_amino_acids = axis.roots.calculate_unloading_amino_acids(axis.phloem.amino_acids)
                axis.roots.s_amino_acids = axis.roots.calculate_s_amino_acids(axis.roots.nitrates, axis.roots.sucrose)
        
                # compartments derivatives
                sucrose_derivative = axis.roots.calculate_sucrose_derivative(axis.roots.unloading_sucrose, axis.roots.s_amino_acids)
                nitrates_derivative = axis.roots.calculate_nitrates_derivative(roots_uptake_nitrate, axis.roots.s_amino_acids)
                amino_acids_derivative = axis.roots.calculate_amino_acids_derivative(axis.roots.unloading_amino_acids, axis.roots.s_amino_acids, roots_export_amino_acids)
                y_derivatives[self.initial_conditions_mapping[axis.roots]['sucrose']] = sucrose_derivative
                y_derivatives[self.initial_conditions_mapping[axis.roots]['nitrates']] = nitrates_derivative
                y_derivatives[self.initial_conditions_mapping[axis.roots]['amino_acids']] = amino_acids_derivative
        
                # compute the derivative of each compartment of phloem
                sucrose_phloem_derivative = axis.phloem.calculate_sucrose_derivative(axis_photosynthetic_organs)
                amino_acids_phloem_derivative = axis.phloem.calculate_amino_acids_derivative(axis_photosynthetic_organs)
                y_derivatives[self.initial_conditions_mapping[axis.phloem]['sucrose']] = sucrose_phloem_derivative
                y_derivatives[self.initial_conditions_mapping[axis.phloem]['amino_acids']] = amino_acids_phloem_derivative

        if self.show_progressbar:
            self.progressbar.update(t)

        derivatives_logger = logging.getLogger('cnwheat.derivatives')
        if derivatives_logger.isEnabledFor(logging.DEBUG):
            self._log_compartments(t, y_derivatives, CNWheat.LOGGERS_NAMES['derivatives'])

        return y_derivatives

    
    def _update_population(self, compartments_values):
        """Update the state of :attr:`population` from the values in `compartments_values`.
        """
        logger = logging.getLogger(__name__)
        logger.debug('Updating the state of the population...')
        for organ, compartments in self.initial_conditions_mapping.iteritems():
            for compartment_name, compartment_index in compartments.iteritems():
                setattr(organ, compartment_name, compartments_values[compartment_index])
        logger.debug('Updating the state of the population DONE')
        

    def _format_solver_output(self, t, solver_output):
        """
        Create :class:`dataframes <pandas.DataFrame>` of outputs at different scales:
        
            * plant: t, plant index, outputs of plant,
            * axis: t, plant index, axis index, outputs of axis,
            * phytomer: t, plant index, axis index, phytomer index, outputs of phytomer,
            * and organ: t, plant index, axis index, phytomer index, organ name, outputs of organ.
        """
        logger = logging.getLogger(__name__)
        logger.debug('Formatting of solver output...')

        solver_output = solver_output.T

        all_plants_df = pd.DataFrame(columns=CNWheat.PLANTS_OUTPUTS)
        all_axes_df = pd.DataFrame(columns=CNWheat.AXES_OUTPUTS)
        all_phytomers_df = pd.DataFrame(columns=CNWheat.PHYTOMERS_OUTPUTS)
        all_organs_df = pd.DataFrame(columns=CNWheat.ORGANS_OUTPUTS)
        
        for plant in self.population.plants:
            
            plants_df = pd.DataFrame(columns=all_plants_df.columns)
            plants_df['t'] = t
            plants_df['plant'] = plant.index
            
            for axis in plant.axes:
                
                axes_df = pd.DataFrame(columns=all_axes_df.columns)
                axes_df['t'] = t
                axes_df['plant'] = plant.index
                axes_df['axis'] = axis.id
                
                # compute the total transpiration
                total_transpiration = np.zeros_like(t)
                
                transpiration_mapping = {}
                for phytomer in axis.phytomers:
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ is not None:
                            transpiration_mapping[organ] = map(organ.calculate_transpiration, t, [organ.Tr] * len(t), [phytomer.index] * len(t))
                            total_transpiration += transpiration_mapping[organ]
                
                axes_df['Total_transpiration'] = total_transpiration
          
                # format phloem outputs
                organs_df = pd.DataFrame(columns=all_organs_df.columns)
                organs_df['t'] = t
                organs_df['plant'] = plant.index
                organs_df['axis'] = axis.id
                organs_df['organ'] = axis.phloem.__class__.__name__
                phloem_sucrose = solver_output[self.initial_conditions_mapping[axis.phloem]['sucrose']]
                organs_df['sucrose'] = phloem_sucrose
                phloem_amino_acids = solver_output[self.initial_conditions_mapping[axis.phloem]['amino_acids']]
                organs_df['amino_acids'] = phloem_amino_acids
                organs_df['Conc_Sucrose'] = axis.phloem.calculate_conc_sucrose(organs_df['sucrose'])
                organs_df['Conc_C_Sucrose'] = axis.phloem.calculate_conc_c_sucrose(organs_df['sucrose'])
                organs_df['Conc_Amino_Acids'] = axis.phloem.calculate_conc_amino_acids(organs_df['amino_acids'])
                all_organs_df = all_organs_df.append(organs_df, ignore_index=True)
        
                # format roots outputs
                organs_df = pd.DataFrame(columns=all_organs_df.columns)
                organs_df['t'] = t
                organs_df['plant'] = plant.index
                organs_df['axis'] = axis.id
                organs_df['organ'] = axis.roots.__class__.__name__
                organs_df['sucrose'] = solver_output[self.initial_conditions_mapping[axis.roots]['sucrose']]
                organs_df['nitrates'] = solver_output[self.initial_conditions_mapping[axis.roots]['nitrates']]
                organs_df['amino_acids'] = solver_output[self.initial_conditions_mapping[axis.roots]['amino_acids']]
                organs_df['Conc_Sucrose'] = axis.roots.calculate_conc_sucrose(organs_df['sucrose'])
                organs_df['Conc_Nitrates'] = axis.roots.calculate_conc_nitrates(organs_df['nitrates'])
                organs_df['Conc_Amino_Acids'] = axis.roots.calculate_conc_amino_acids(organs_df['amino_acids'])
                organs_df['Conc_Nitrates_Soil'] = map(axis.roots.calculate_conc_nitrates_soil,t)
                organs_df['Unloading_Sucrose'] = map(axis.roots.calculate_unloading_sucrose, phloem_sucrose)
                organs_df['Unloading_Amino_Acids'] = map(axis.roots.calculate_unloading_amino_acids, phloem_amino_acids)
                roots_uptake_nitrates, roots_potential_uptake_nitrates = axis.roots.calculate_uptake_nitrates(organs_df['Conc_Nitrates_Soil'], organs_df['nitrates'], total_transpiration) 
                organs_df['Uptake_Nitrates'] = roots_uptake_nitrates
                organs_df['Potential_Uptake_Nitrates'] = map(axis.roots.calculate_export_amino_acids, organs_df['amino_acids'], total_transpiration)
                roots_export_amino_acids = map(axis.roots.calculate_export_amino_acids, organs_df['amino_acids'], total_transpiration)
                organs_df['Export_Amino_Acids'] = roots_export_amino_acids
                organs_df['S_Amino_Acids'] = map(axis.roots.calculate_s_amino_acids, organs_df['nitrates'], organs_df['sucrose'])
                all_organs_df = all_organs_df.append(organs_df, ignore_index=True)
        
                # format photosynthetic organs outputs
                for phytomer in axis.phytomers:
                    phytomers_df = pd.DataFrame(columns=all_phytomers_df.columns)
                    phytomers_df['t'] = t
                    phytomers_df['plant'] = plant.index
                    phytomers_df['axis'] = axis.id
                    phytomers_df['phytomer'] = phytomer.index
                    
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ is None:
                            continue
                        organs_df = pd.DataFrame(columns=all_organs_df.columns)
                        organs_df['t'] = t
                        organs_df['plant'] = plant.index
                        organs_df['axis'] = axis.id
                        organs_df['phytomer'] = phytomer.index
                        organs_df['organ'] = organ.__class__.__name__
                        organs_df['triosesP'] = solver_output[self.initial_conditions_mapping[organ]['triosesP']]
                        organs_df['starch'] = solver_output[self.initial_conditions_mapping[organ]['starch']]
                        organs_df['sucrose'] = solver_output[self.initial_conditions_mapping[organ]['sucrose']]
                        organs_df['fructan'] = solver_output[self.initial_conditions_mapping[organ]['fructan']]
                        organs_df['nitrates'] = solver_output[self.initial_conditions_mapping[organ]['nitrates']]
                        organs_df['amino_acids'] = solver_output[self.initial_conditions_mapping[organ]['amino_acids']]
                        organs_df['proteins'] = solver_output[self.initial_conditions_mapping[organ]['proteins']]
                        organs_df['Loading_Sucrose'] = map(organ.calculate_loading_sucrose, organs_df['sucrose'], phloem_sucrose)
                        organs_df['Regul_S_Fructan'] = map(organ.calculate_regul_s_fructan, organs_df['Loading_Sucrose'])
                        organs_df['An'] = organ.An
                        organs_df['Tr'] = organ.Tr
                        organs_df['Photosynthesis'] = map(organ.calculate_photosynthesis, t, [organ.An] * len(t), [phytomer.index] * len(t))
                        organs_df['Transpiration'] = transpiration_mapping[organ]
                        organs_df['Conc_TriosesP'] = organ.calculate_conc_triosesP(organs_df['triosesP'])
                        organs_df['Conc_Starch'] = organ.calculate_conc_starch(organs_df['starch'])
                        organs_df['Conc_Sucrose'] = organ.calculate_conc_sucrose(organs_df['sucrose'])
                        organs_df['Conc_Fructan'] = organ.calculate_conc_fructan(organs_df['fructan'])
                        organs_df['Conc_Nitrates'] = organ.calculate_conc_nitrates(organs_df['nitrates'])
                        organs_df['Conc_Amino_Acids'] = organ.calculate_conc_amino_acids(organs_df['amino_acids'])
                        organs_df['Conc_Proteins'] = organ.calculate_conc_proteins(organs_df['proteins'])
                        organs_df['S_Starch'] = map(organ.calculate_s_starch, organs_df['triosesP'])
                        organs_df['D_Starch'] = map(organ.calculate_d_starch, organs_df['starch'])
                        organs_df['S_Sucrose'] = map(organ.calculate_s_sucrose, organs_df['triosesP'])
                        organs_df['S_Fructan'] = map(organ.calculate_s_fructan, organs_df['sucrose'], organs_df['Regul_S_Fructan'])
                        organs_df['D_Fructan'] = map(organ.calculate_d_fructan, organs_df['sucrose'], organs_df['fructan'])
                        organs_df['Nitrates_import'] = map(organ.calculate_nitrates_import, roots_uptake_nitrates, transpiration_mapping[organ], total_transpiration)
                        organs_df['Amino_Acids_import'] = map(organ.calculate_amino_acids_import, roots_export_amino_acids, transpiration_mapping[organ], total_transpiration)
                        organs_df['S_Amino_Acids'] = map(organ.calculate_s_amino_acids, organs_df['nitrates'], organs_df['triosesP'])
                        organs_df['S_Proteins'] = map(organ.calculate_s_proteins, organs_df['amino_acids'])
                        organs_df['D_Proteins'] = map(organ.calculate_d_proteins, organs_df['proteins'])
                        organs_df['Loading_Amino_Acids_'] = map(organ.calculate_loading_amino_acids, organs_df['amino_acids'], phloem_amino_acids)
                        all_organs_df = all_organs_df.append(organs_df, ignore_index=True)
                        
                    all_phytomers_df = all_phytomers_df.append(phytomers_df, ignore_index=True)
        
                # format grains outputs
                if axis.grains is None:
                    continue
                
                organs_df = pd.DataFrame(columns=all_organs_df.columns)
                organs_df['t'] = t
                organs_df['plant'] = plant.index
                organs_df['axis'] = axis.id
                organs_df['organ'] = axis.grains.__class__.__name__
                organs_df['structure'] = solver_output[self.initial_conditions_mapping[axis.grains]['structure']]
                organs_df['starch'] = solver_output[self.initial_conditions_mapping[axis.grains]['starch']]
                organs_df['proteins'] = solver_output[self.initial_conditions_mapping[axis.grains]['proteins']]
                organs_df['RGR_Structure'] = map(axis.grains.calculate_RGR_structure, phloem_sucrose)
                organs_df['S_grain_structure'] = map(axis.grains.calculate_s_grain_structure, t, organs_df['structure'], organs_df['RGR_Structure'])
                organs_df['S_grain_starch'] = map(axis.grains.calculate_s_grain_starch, t, phloem_sucrose)
                organs_df['Dry_Mass'] = axis.grains.calculate_dry_mass(organs_df['structure'], organs_df['starch'])
                organs_df['Proteins_N_Mass'] = axis.grains.calculate_protein_mass(organs_df['proteins'])
                organs_df['Unloading_Sucrose_'] = map(axis.grains.calculate_unloading_sucrose, organs_df['S_grain_structure'], organs_df['S_grain_starch'], organs_df['structure'])
                organs_df['S_Proteins'] = map(axis.grains.calculate_s_proteins, organs_df['S_grain_structure'], organs_df['S_grain_starch'], phloem_amino_acids, phloem_sucrose, organs_df['structure'])
                all_organs_df = all_organs_df.append(organs_df, ignore_index=True)
                
                all_axes_df = all_axes_df.append(axes_df, ignore_index=True)
            
            all_plants_df = all_plants_df.append(plants_df, ignore_index=True)
            
        # set the order of the columns
        all_plants_df = all_plants_df.reindex_axis(CNWheat.PLANTS_OUTPUTS, axis=1, copy=False)
        all_axes_df = all_axes_df.reindex_axis(CNWheat.AXES_OUTPUTS, axis=1, copy=False)
        all_phytomers_df = all_phytomers_df.reindex_axis(CNWheat.PHYTOMERS_OUTPUTS, axis=1, copy=False)
        all_organs_df = all_organs_df.reindex_axis(CNWheat.ORGANS_OUTPUTS, axis=1, copy=False)
        
        # sort the rows by the columns
        all_plants_df.sort_index(by=CNWheat.PLANTS_INDEXES, inplace=True)
        all_axes_df.sort_index(by=CNWheat.AXES_INDEXES, inplace=True)
        all_phytomers_df.sort_index(by=CNWheat.PHYTOMERS_INDEXES, inplace=True)
        all_organs_df.sort_index(by=CNWheat.ORGANS_INDEXES, inplace=True)
        
        all_plants_df = all_plants_df.convert_objects(copy=False)
        all_axes_df = all_axes_df.convert_objects(copy=False)
        all_phytomers_df = all_phytomers_df.convert_objects(copy=False)
        all_organs_df = all_organs_df.convert_objects(copy=False)
        
        logger.debug('Formatting of solver output DONE')
        
        return all_plants_df, all_axes_df, all_phytomers_df, all_organs_df


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

