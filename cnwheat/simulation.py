# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    cnwheat.simulation
    ~~~~~~~~~~~~~~~~~~

    Front-end to run the model.

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

        - `organs` (:class:`list`) - List of :class:`cnwheat.model.Organ`.
          Must contain at least :

              * a :obj:`photosynthetic organ <cnwheat.model.PhotosyntheticOrgan>`,
              * a :obj:`phloem <cnwheat.model.Phloem>`,
              * a :obj:`roots <cnwheat.model.Roots>`,
              * a :obj:`grains <cnwheat.model.Grains>`.

    """

    #: the name of the compartments attributes in the organs. 
    MODEL_COMPARTMENTS_NAMES = {model.PhotosyntheticOrgan: ['triosesP', 'starch', 'sucrose', 'fructan', 'nitrates', 'amino_acids', 'proteins'], 
                                model.Phloem: ['sucrose', 'amino_acids'], 
                                model.Grains: ['starch', 'structure', 'proteins'], 
                                model.Roots: ['sucrose', 'nitrates', 'amino_acids']}

    def __init__(self, organs):

        logger = logging.getLogger(__name__)

        logger.info('Initialization of the simulation...')

        self.organs = organs #: the organs used in the model

        # construct the list of initial conditions
        self.initial_conditions = [] #: the initial conditions of the compartments in the organs
        self.initial_conditions_mapping = {} #: dictionary to map the compartments of each organ to their indexes in :attr:`initial_conditions`
        self.compartments_absolute_names = [] #: list of the compartments absolute names. The absolute name is the concatenation of the organ name and the compartment name. 
        self.photosynthetic_organs = [] #: the photosynthetic_organs
        
        i = 0
        for organ in organs:
            if isinstance(organ, model.PhotosyntheticOrgan):
                self.photosynthetic_organs.append(organ)
                compartments_names = CNWheat.MODEL_COMPARTMENTS_NAMES[model.PhotosyntheticOrgan]
            elif isinstance(organ, model.Phloem):
                self.phloem = organ # the phloem
                compartments_names = CNWheat.MODEL_COMPARTMENTS_NAMES[model.Phloem]
            elif isinstance(organ, model.Roots):
                self.roots = organ # the roots
                compartments_names = CNWheat.MODEL_COMPARTMENTS_NAMES[model.Roots]
            elif isinstance(organ, model.Grains):
                self.grains = organ # the grains
                compartments_names = CNWheat.MODEL_COMPARTMENTS_NAMES[model.Grains]
            self.initial_conditions_mapping[organ] = {}
            for compartment_name in compartments_names:
                self.initial_conditions_mapping[organ][compartment_name] = i
                self.compartments_absolute_names.append('.'.join([organ.name, compartment_name]))
                self.initial_conditions.append(0)
                i += 1
                
        try:
            if len(self.photosynthetic_organs) == 0:
                raise CNWheatInitError("No photosynthetic organ in 'organs'.")
            if self.phloem is None:
                raise CNWheatInitError("No phloem in 'organs'.")
            if self.roots is None:
                raise CNWheatInitError("No roots in 'organs'.")
            if self.grains is None:
                raise CNWheatInitError("No grains in 'organs'.")
        except CNWheatInitError, e:
            e.message += " 'organs' must contain at least: a photosynthetic organ, a phloem, a roots and a grains."
            logger.exception(e.message)
            raise e

        self.progressbar = ProgressBar(title='Solver progress') #: progress bar to show the progress of the solver
        self.show_progressbar = False #: True: show the progress bar ; False: DO NOT show the progress bar
        logger.info('Initialization of the simulation DONE')


    def run(self, start_time, stop_time, number_of_output_steps, odeint_mxstep=5000, show_progressbar=False):
        """
        Compute CN exchanges between organs :attr:`organs`.
        The computation is done between `start_time` and `stop_time`, for `number_of_output_steps` steps.

        :Parameters:

            - `start_time` (:class:`int`) - The starting of the time grid.

            - `stop_time` (:class:`int`) - The end of the time grid.

            - `number_of_output_steps` (:class:`int`) - Number of time points for which to compute the CN exchanges in the system.

            - `odeint_mxstep` (:class:`int`) - Maximum number of (internally defined) steps allowed for each integration point in time grid.
              `odeint_mxstep` is passed to :func:`scipy.integrate.odeint` as `mxstep`. If `odeint_mxstep` = 0, then `mxstep` is determined by the solver.
              The default value ( `5000` ) normally permits to solve the current model. User should increased this value if a more complex model is defined
              and if this model make the integration failed.

            - `show_progressbar` (:class:`bool`) - True: show the progress bar ; False: do not show the progress bar.

        :Returns:
            Dataframe containing the CN exchanges between organs for each desired time step.

        :Returns Type:
            :class:`pandas.DataFrame`

        """
        logger = logging.getLogger(__name__)
        logger.info('Run of the simulation...')

        t = np.linspace(start_time, stop_time, number_of_output_steps)

        self.show_progressbar = show_progressbar
        if self.show_progressbar:
            self.progressbar.set_t_max(stop_time)

        compartments_logger = logging.getLogger('cnwheat.compartments')
        derivatives_logger = logging.getLogger('cnwheat.derivatives')
        if compartments_logger.isEnabledFor(logging.DEBUG) or derivatives_logger.isEnabledFor(logging.DEBUG):
            formatted_compartment_names = ','.join(['t'] + self.compartments_absolute_names)
            if compartments_logger.isEnabledFor(logging.DEBUG):
                compartments_logger.debug(formatted_compartment_names)
            if derivatives_logger.isEnabledFor(logging.DEBUG):
                derivatives_logger.debug(formatted_compartment_names)
        
        self._update_initial_conditions() 
            
        if logger.isEnabledFor(logging.DEBUG):
            initial_conditions = dict(zip(self.compartments_absolute_names, self.initial_conditions))
            logger.debug(
                """Run the solver with:
                    - initial conditions = %s,
                    - time grid = %s,
                    - odeint mxstep = %s""",
                initial_conditions, t, odeint_mxstep)

        soln, infodict = odeint(self._calculate_all_derivatives, self.initial_conditions, t, full_output=True, mxstep=odeint_mxstep)

        if not set(infodict['mused']).issubset([1,2]): # I'm not sure if this test is robust or not... Beware especially when scipy is updated.
            message = "Integration failed. See the logs of lsoda or try to increase the value of 'mxstep'."
            logger.exception(message)
            raise CNWheatRunError(message)
        
        last_compartments_values = soln[len(soln)-1]
        self._update_organs(last_compartments_values)
        
        if logger.isEnabledFor(logging.DEBUG):
            compartments = dict(zip(self.compartments_absolute_names, soln[len(soln)-1]))
            logger.debug(
                """Run of the solver DONE:
                    - compartments(t=%s) = %s,
                    - infodict = %s""",
                stop_time, compartments, infodict)
        
        cnwheat_output_df = self._format_solver_output(t, soln)
        
        logger.info('Run of the simulation DONE')

        return cnwheat_output_df
    
    
    def _update_initial_conditions(self):
        """Update the compartments values in :attr:`initial_conditions` from the compartments values of :attr:`organs`.
        """
        for organ, compartments in self.initial_conditions_mapping.iteritems():
            for compartment_name, compartment_index in compartments.iteritems():
                self.initial_conditions[compartment_index] = getattr(organ, compartment_name)


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

            - `y` (:class:`list`) - The current y values. The y values of :attr:`phloem`
              must appear first. `y` is automatically set by
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
            formatted_initial_conditions = ','.join(map(str, [t] + y.tolist()))
            compartments_logger.debug(formatted_initial_conditions)

        # check that the solver is not crashed
        y_isnan = np.isnan(y)
        if y_isnan.any():
            y_isnan_indices = np.where(y_isnan)
            nan_compartments = []
            for organ, compartments in self.initial_conditions_mapping.iteritems():
                for compartment_name, compartment_index in compartments.iteritems():
                    if np.in1d(y_isnan_indices, compartment_index).any():
                        nan_compartments.append(organ.name + '.' + compartment_name)
            message = 'The solver did not manage to compute the compartments {}'.format(nan_compartments)
            logger.exception(message)
            raise CNWheatRunError(message)

        y_derivatives = np.zeros_like(y)

        self.phloem.sucrose = y[self.initial_conditions_mapping[self.phloem]['sucrose']]
        self.phloem.amino_acids = y[self.initial_conditions_mapping[self.phloem]['amino_acids']]

        self.roots.nitrates = y[self.initial_conditions_mapping[self.roots]['nitrates']]
        self.roots.amino_acids = y[self.initial_conditions_mapping[self.roots]['amino_acids']]

        # compute the total transpiration at t_inf
        total_transpiration = 0.0
        transpiration_mapping = dict.fromkeys(self.photosynthetic_organs)
        for organ in self.photosynthetic_organs:
            transpiration_mapping[organ] = organ.calculate_transpiration(t, organ.Tr)
            total_transpiration += transpiration_mapping[organ]
        #print 'total_transpiration', total_transpiration

        # compute the flows from/to the roots to/from photosynthetic organs
        conc_nitrates_soil = self.roots.calculate_conc_nitrates_soil(t)
        roots_uptake_nitrate, potential_roots_uptake_nitrates = self.roots.calculate_uptake_nitrates(conc_nitrates_soil, self.roots.nitrates, total_transpiration)
        roots_export_amino_acids = self.roots.calculate_export_amino_acids(self.roots.amino_acids, total_transpiration)
        #print 'roots_export_amino_acids',roots_export_amino_acids

        # compute the derivative of each photosynthetic organ compartment
        for organ in self.photosynthetic_organs:
            organ.starch = y[self.initial_conditions_mapping[organ]['starch']]
            organ.sucrose = y[self.initial_conditions_mapping[organ]['sucrose']]
            organ.triosesP = y[self.initial_conditions_mapping[organ]['triosesP']]
            organ.fructan = y[self.initial_conditions_mapping[organ]['fructan']]
            organ.nitrates = y[self.initial_conditions_mapping[organ]['nitrates']]
            organ.amino_acids = y[self.initial_conditions_mapping[organ]['amino_acids']]
            organ.proteins = y[self.initial_conditions_mapping[organ]['proteins']]

            # intermediate variables
            photosynthesis = organ.calculate_photosynthesis(t, organ.An)
            organ_transpiration = transpiration_mapping[organ]

            # flows
            s_starch = organ.calculate_s_starch(organ.triosesP)
            d_starch = organ.calculate_d_starch(organ.starch)
            s_sucrose = organ.calculate_s_sucrose(organ.triosesP)
            organ.loading_sucrose = organ.calculate_loading_sucrose(organ.sucrose, self.phloem.sucrose)
            regul_s_fructan = organ.calculate_regul_s_fructan(organ.loading_sucrose)
            d_fructan = organ.calculate_d_fructan(organ.sucrose, organ.fructan)
            s_fructan = organ.calculate_s_fructan(organ.sucrose, regul_s_fructan)
            nitrates_import = organ.calculate_nitrates_import(roots_uptake_nitrate, organ_transpiration, total_transpiration)
            amino_acids_import = organ.calculate_amino_acids_import(roots_export_amino_acids, organ_transpiration, total_transpiration)
            s_amino_acids = organ.calculate_s_amino_acids(organ.nitrates, organ.triosesP)
            s_proteins = organ.calculate_s_proteins(organ.amino_acids)
            d_proteins = organ.calculate_d_proteins(organ.proteins)
            organ.loading_amino_acids = organ.calculate_loading_amino_acids(organ.amino_acids, self.phloem.amino_acids)

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


        # compute the derivative of each compartment of grains
        self.grains.structure = y[self.initial_conditions_mapping[self.grains]['structure']]
        self.grains.starch = y[self.initial_conditions_mapping[self.grains]['starch']]
        self.grains.proteins = y[self.initial_conditions_mapping[self.grains]['proteins']]

        # intermediate variables
        RGR_structure = self.grains.calculate_RGR_structure(self.phloem.sucrose)

        # flows
        self.grains.s_grain_structure = self.grains.calculate_s_grain_structure(t, self.grains.structure, RGR_structure)
        self.grains.s_grain_starch = self.grains.calculate_s_grain_starch(t, self.phloem.sucrose)
        self.grains.s_proteins = self.grains.calculate_s_proteins(self.grains.s_grain_structure, self.grains.s_grain_starch, self.phloem.amino_acids, self.phloem.sucrose, self.grains.structure)

        # compartments derivatives
        structure_derivative = self.grains.calculate_structure_derivative(self.grains.s_grain_structure)
        starch_derivative = self.grains.calculate_starch_derivative(self.grains.s_grain_starch, self.grains.structure)
        proteins_derivative = self.grains.calculate_proteins_derivative(self.grains.s_proteins)
        y_derivatives[self.initial_conditions_mapping[self.grains]['structure']] = structure_derivative
        y_derivatives[self.initial_conditions_mapping[self.grains]['starch']] = starch_derivative
        y_derivatives[self.initial_conditions_mapping[self.grains]['proteins']] = proteins_derivative

        # compute the derivative of each compartment of roots
        self.roots.sucrose = y[self.initial_conditions_mapping[self.roots]['sucrose']]
        # flows
        self.roots.unloading_sucrose = self.roots.calculate_unloading_sucrose(self.phloem.sucrose)
        self.roots.unloading_amino_acids = self.roots.calculate_unloading_amino_acids(self.phloem.amino_acids)
        self.roots.s_amino_acids = self.roots.calculate_s_amino_acids(self.roots.nitrates, self.roots.sucrose)

        # compartments derivatives
        sucrose_derivative = self.roots.calculate_sucrose_derivative(self.roots.unloading_sucrose, self.roots.s_amino_acids)
        nitrates_derivative = self.roots.calculate_nitrates_derivative(roots_uptake_nitrate, self.roots.s_amino_acids)
        amino_acids_derivative = self.roots.calculate_amino_acids_derivative(self.roots.unloading_amino_acids, self.roots.s_amino_acids, roots_export_amino_acids)
        y_derivatives[self.initial_conditions_mapping[self.roots]['sucrose']] = sucrose_derivative
        y_derivatives[self.initial_conditions_mapping[self.roots]['nitrates']] = nitrates_derivative
        y_derivatives[self.initial_conditions_mapping[self.roots]['amino_acids']] = amino_acids_derivative

        # compute the derivative of each compartment of phloem
        sucrose_phloem_derivative = self.phloem.calculate_sucrose_derivative(self.organs)
        amino_acids_phloem_derivative = self.phloem.calculate_amino_acids_derivative(self.organs)
        y_derivatives[self.initial_conditions_mapping[self.phloem]['sucrose']] = sucrose_phloem_derivative
        y_derivatives[self.initial_conditions_mapping[self.phloem]['amino_acids']] = amino_acids_phloem_derivative

        if self.show_progressbar:
            self.progressbar.update(t)

        derivatives_logger = logging.getLogger('cnwheat.derivatives')
        if derivatives_logger.isEnabledFor(logging.DEBUG):
            formatted_derivatives = ','.join(map(str, [t] + y_derivatives.tolist()))
            derivatives_logger.debug(formatted_derivatives)

        return y_derivatives

    
    def _update_organs(self, compartments_values):
        """Update the organs in :attr:`organs` from the values in `compartments`.
        """
        logger = logging.getLogger(__name__)
        logger.debug('Updating of organs compartments...')
        for organ, compartments in self.initial_conditions_mapping.iteritems():
            for compartment_name, compartment_index in compartments.iteritems():
                setattr(organ, compartment_name, compartments_values[compartment_index])
        logger.debug('Updating of organs compartments DONE')
        

    def _format_solver_output(self, t, solver_output):
        """
        Create a :class:`pandas.DataFrame` wich columns are:

            * the time grid `t`,
            * the output of the solver `solver_output`,
            * and intermediate and post-processed variables useful for debug and validation.
            
        This is mainly use for debugging.
        """
        logger = logging.getLogger(__name__)
        logger.debug('Formatting of solver output...')

        solver_output = solver_output.T

        result_items = [('t', t)]
        
        # compute the total transpiration
        total_transpiration = np.zeros_like(t)
        transpiration_mapping = dict.fromkeys(self.photosynthetic_organs)
        for organ in self.photosynthetic_organs:
            transpiration_mapping[organ] = map(organ.calculate_transpiration, t, [organ.Tr] * len(t))
            total_transpiration += transpiration_mapping[organ]
        result_items.append(('Total_transpiration', total_transpiration))

        # format phloem outputs
        sucrose_phloem = solver_output[self.initial_conditions_mapping[self.phloem]['sucrose']]
        amino_acids_phloem = solver_output[self.initial_conditions_mapping[self.phloem]['amino_acids']]
        variables = [('Conc_Sucrose_{}'.format(self.phloem.name), self.phloem.calculate_conc_sucrose(sucrose_phloem)),
                     ('Conc_C_Sucrose_{}'.format(self.phloem.name), self.phloem.calculate_conc_c_sucrose(sucrose_phloem)),
                     ('Conc_Amino_Acids_{}'.format(self.phloem.name), self.phloem.calculate_conc_amino_acids(amino_acids_phloem))]
        compartments = [('Sucrose_{}'.format(self.phloem.name), sucrose_phloem),
                        ('Amino_Acids_{}'.format(self.phloem.name), amino_acids_phloem)]
        result_items.extend(variables + compartments)

        # format roots outputs
        sucrose_roots = solver_output[self.initial_conditions_mapping[self.roots]['sucrose']]
        unloading_sucrose = map(self.roots.calculate_unloading_sucrose, sucrose_phloem)
        nitrates_roots = solver_output[self.initial_conditions_mapping[self.roots]['nitrates']]
        conc_nitrates_soil = map(self.roots.calculate_conc_nitrates_soil,t)
        uptake_nitrates, potential_roots_uptake_nitrates = self.roots.calculate_uptake_nitrates(conc_nitrates_soil, nitrates_roots, total_transpiration)
        amino_acids_roots = solver_output[self.initial_conditions_mapping[self.roots]['amino_acids']]
        roots_export_amino_acids =  map(self.roots.calculate_export_amino_acids, amino_acids_roots, total_transpiration)
        variables = [('Conc_Sucrose_{}'.format(self.roots.name), self.roots.calculate_conc_sucrose(sucrose_roots)),
                     ('Conc_Nitrates_{}'.format(self.roots.name), self.roots.calculate_conc_nitrates(nitrates_roots)),
                     ('Conc_Amino_Acids_{}'.format(self.roots.name), self.roots.calculate_conc_amino_acids(amino_acids_roots)),
                     ('Conc_Nitrates_Soil_{}'.format(self.roots.name), conc_nitrates_soil)]
        flows = [('Unloading_Sucrose_{}'.format(self.roots.name), unloading_sucrose),
                 ('Unloading_Amino_Acids_{}'.format(self.roots.name), map(self.roots.calculate_unloading_amino_acids, amino_acids_phloem)),
                 ('Uptake_Nitrates_{}'.format(self.roots.name), uptake_nitrates),
                 ('Potential_Uptake_Nitrates_{}'.format(self.roots.name), potential_roots_uptake_nitrates),
                 ('Export_Amino_Acids_{}'.format(self.roots.name),roots_export_amino_acids),
                 ('S_Amino_Acids_{}'.format(self.roots.name), map(self.roots.calculate_s_amino_acids, nitrates_roots, sucrose_roots))]
        compartments = [('Sucrose_{}'.format(self.roots.name), sucrose_roots),
                        ('Nitrates_{}'.format(self.roots.name), nitrates_roots),
                        ('Amino_Acids_{}'.format(self.roots.name), amino_acids_roots)]
        result_items.extend(variables + flows + compartments)

        # format photosynthetic organs outputs
        for organ in self.photosynthetic_organs:
            triosesP = solver_output[self.initial_conditions_mapping[organ]['triosesP']]
            starch = solver_output[self.initial_conditions_mapping[organ]['starch']]
            sucrose = solver_output[self.initial_conditions_mapping[organ]['sucrose']]
            fructan = solver_output[self.initial_conditions_mapping[organ]['fructan']]
            loading_sucrose = map(organ.calculate_loading_sucrose, sucrose, sucrose_phloem)
            regul_s_fructan = map(organ.calculate_regul_s_fructan, loading_sucrose)

            nitrates = solver_output[self.initial_conditions_mapping[organ]['nitrates']]
            amino_acids = solver_output[self.initial_conditions_mapping[organ]['amino_acids']]
            proteins = solver_output[self.initial_conditions_mapping[organ]['proteins']]

            variables = [('An_{}'.format(organ.name), organ.An),
                         ('Tr_{}'.format(organ.name), organ.Tr),
                         ('Photosynthesis_{}'.format(organ.name), map(organ.calculate_photosynthesis, t, [organ.An] * len(t))),
                         ('Transpiration_{}'.format(organ.name), transpiration_mapping[organ]),
                         ('Conc_TriosesP_{}'.format(organ.name), organ.calculate_conc_triosesP(triosesP)),
                         ('Conc_Starch_{}'.format(organ.name), organ.calculate_conc_starch(starch)),
                         ('Conc_Sucrose_{}'.format(organ.name), organ.calculate_conc_sucrose(sucrose)),
                         ('Conc_Fructan_{}'.format(organ.name), organ.calculate_conc_fructan(fructan)),
                         ('Regul_S_Fructan_{}'.format(organ.name), regul_s_fructan),
                         ('Conc_Nitrates_{}'.format(organ.name), organ.calculate_conc_nitrates(nitrates)),
                         ('Conc_Amino_Acids_{}'.format(organ.name), organ.calculate_conc_amino_acids(amino_acids)),
                         ('Conc_Proteins_{}'.format(organ.name), organ.calculate_conc_proteins(proteins))]


            flows = [('S_Starch_{}'.format(organ.name), map(organ.calculate_s_starch, triosesP)),
                     ('D_Starch_{}'.format(organ.name), map(organ.calculate_d_starch, starch)),
                     ('S_Sucrose_{}'.format(organ.name), map(organ.calculate_s_sucrose, triosesP)),
                     ('Loading_Sucrose_{}'.format(organ.name), loading_sucrose),
                     ('S_Fructan_{}'.format(organ.name), map(organ.calculate_s_fructan, sucrose, regul_s_fructan)),
                     ('D_Fructan_{}'.format(organ.name), map(organ.calculate_d_fructan, sucrose, fructan)),
                     ('Nitrates_import_{}'.format(organ.name), map(organ.calculate_nitrates_import, uptake_nitrates, transpiration_mapping[organ], total_transpiration)),
                     ('Amino_Acids_import_{}'.format(organ.name), map(organ.calculate_amino_acids_import, roots_export_amino_acids, transpiration_mapping[organ], total_transpiration)),
                     ('S_Amino_Acids_{}'.format(organ.name), map(organ.calculate_s_amino_acids, nitrates, triosesP)),
                     ('S_Proteins_{}'.format(organ.name), map(organ.calculate_s_proteins, amino_acids)),
                     ('D_Proteins_{}'.format(organ.name), map(organ.calculate_d_proteins, proteins)),
                     ('Loading_Amino_Acids_{}'.format(organ.name), map(organ.calculate_loading_amino_acids, amino_acids, amino_acids_phloem))]

            compartments = [('TriosesP_{}'.format(organ.name), triosesP),
                            ('Starch_{}'.format(organ.name), starch),
                            ('Sucrose_{}'.format(organ.name), sucrose),
                            ('Fructan_{}'.format(organ.name), fructan),
                            ('Nitrates_{}'.format(organ.name), nitrates),
                            ('Amino_Acids_{}'.format(organ.name), amino_acids),
                            ('Proteins_{}'.format(organ.name), proteins)]

            result_items.extend(variables + flows + compartments)

        # format grains outputs
        structure_grains = solver_output[self.initial_conditions_mapping[self.grains]['structure']]
        starch_grains = solver_output[self.initial_conditions_mapping[self.grains]['starch']]
        proteins_grains = solver_output[self.initial_conditions_mapping[self.grains]['proteins']]

        RGR_structure_grains = map(self.grains.calculate_RGR_structure, sucrose_phloem)

        s_grain_structure = map(self.grains.calculate_s_grain_structure, t, structure_grains, RGR_structure_grains)
        s_grain_starch = map(self.grains.calculate_s_grain_starch, t, sucrose_phloem)

        variables = [('Dry_Mass_{}'.format(self.grains.name), self.grains.calculate_dry_mass(structure_grains, starch_grains)),
                     ('RGR_Structure_{}'.format(self.grains.name), RGR_structure_grains),
                     ('Proteins_N_Mass_{}'.format(self.grains.name), self.grains.calculate_protein_mass(proteins_grains)),
                     ('Unloading_Sucrose_{}'.format(self.grains.name), map(self.grains.calculate_unloading_sucrose, s_grain_structure, s_grain_starch, structure_grains))]

        flows = [('S_grain_structure{}'.format(self.grains.name), s_grain_structure),
                 ('S_grain_starch_{}'.format(self.grains.name), s_grain_starch),
                 ('S_Proteins_{}'.format(self.grains.name), map(self.grains.calculate_s_proteins, s_grain_structure, s_grain_starch, amino_acids_phloem, sucrose_phloem, structure_grains))]


        compartments = [('Structure_{}'.format(self.grains.name), structure_grains),
                        ('Starch_{}'.format(self.grains.name), starch_grains),
                        ('Proteins_{}'.format(self.grains.name), proteins_grains)]

        result_items.extend(variables + flows + compartments)

        logger.debug('Formatting of solver output DONE')

        return pd.DataFrame.from_items(result_items)


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

