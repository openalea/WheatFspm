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
from scipy.interpolate import interp1d

import model
from farquharwheat import model as photosynthesis_model

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

        - `meteo` (:class:`pandas.DataFrame`) - a :class:`pandas.DataFrame` which index
          is time in hours and which columns are "air_temperature", "humidity", "ambient_CO2" and "Wind".

    """

    def __init__(self, organs, meteo):

        self.organs = organs #: the organs used in the model
        self.meteo = meteo #: the meteo data used in the model

        # interpolate meteo data

        #: A dictionary which contains, for each data of meteo data, a function which
        #: permits to find the value of new points using interpolation. Keys
        #: are the name of meteo data (:class:`str`). Values are :class:`scipy.interpolate.interp1d`.
        self.meteo_interpolations = {}

        for column in meteo.columns:
            self.meteo_interpolations[column] = interp1d(meteo.index, meteo[column])

        # interpolate the PAR and construct the list of initial conditions
        self.initial_conditions = [] #: the initial conditions of the compartments in the organs
        self.initial_conditions_mapping = {} #: dictionary to map the compartments of each organ to their indexes in :attr:`initial_conditions`
        self.PAR_linear_interpolation = {} # the linear interpolation of PAR for each photosynthetic organ
        self.photosynthesis_mapping = {} # mapping to store the computed photosynthesis and transpiration for each photosynthetic organ

        i = 0
        for organ in organs:
            if isinstance(organ, model.PhotosyntheticOrgan):
                self.PAR_linear_interpolation[organ] = interp1d(organ.PAR.index, organ.PAR)
                self.photosynthesis_mapping[organ] = {}
            elif isinstance(organ, model.Phloem):
                self.phloem = organ # the phloem
            elif isinstance(organ, model.Roots):
                self.roots = organ # the roots
            elif isinstance(organ, model.Grains):
                self.grains = organ # the grains
            self.initial_conditions_mapping[organ] = {}
            for compartment_name, compartment_initial_value in organ.initial_conditions.iteritems():
                self.initial_conditions_mapping[organ][compartment_name] = i
                self.initial_conditions.append(compartment_initial_value)
                i += 1
        
        try:
            if len(self.photosynthesis_mapping) == 0:
                raise CNWheatInitError("No photosynthetic organ in 'organs'.")
            if self.phloem is None:
                raise CNWheatInitError("No phloem in 'organs'.")
            if self.roots is None:
                raise CNWheatInitError("No roots in 'organs'.")
            if self.grains is None:
                raise CNWheatInitError("No grains in 'organs'.")
        except CNWheatInitError, e:
            e.message += " 'organs' must contain at least: a photosynthetic organ, a phloem, a roots and a grains."
            raise e
        
        self.progressbar = ProgressBar(title='Solver progress') #: progress bar to show the progress of the solver
        self.show_progressbar = False #: True: show the progress bar ; False: DO NOT show the progress bar


    def run(self, start_time, stop_time, number_of_output_steps, photosynthesis_computation_interval=0, odeint_mxstep=5000, show_progressbar=False):
        """
        Compute CN exchanges between organs :attr:`organs`, using interpolated meteo data :attr:`meteo_interpolations`
        The computation is done between `start_time` and `stop_time`, for `number_of_output_steps` steps.

        :Parameters:

            - `start_time` (:class:`int`) - The starting of the time grid.

            - `stop_time` (:class:`int`) - The end of the time grid.

            - `number_of_output_steps` (:class:`int`) - Number of time points for which to compute the CN exchanges in the system.

            - `photosynthesis_computation_interval` (:class:`int`) - The interval which defined the time steps
              at which photosynthesis is computed. For example, if `photosynthesis_computation_interval` = 4,
              then photosynthesis is computed at t=0, t=4, t=8, ...
              This permits to save computation time for large simulation.
              If `photosynthesis_computation_interval` = 0 (the default), then photosynthesis
              is computed for each time step demanded by the solver.

            - `odeint_mxstep` (:class:`int`) - Maximum number of (internally defined) steps allowed for each integration point in time grid.
              `odeint_mxstep` is passed to :func:`scipy.integrate.odeint` as `mxstep`. If `odeint_mxstep` = 0, then `mxstep` is determined by the solver.
              The default value ( `5000` ) normally permits to solve the current model. User should increased this value if a more complex model is defined
              and if this model make the integration failed.

            - `show_progressbar` (:class:`bool`) - True: show the progress bar ; False: do not show the progress bar.

        :Returns:
            Dataframe containing the CN exchanges between organs for each desired time step.

        :Returns Type:
            :class:`pandas.DataFrame`

        .. warning:: due to the solver of :func:`scipy.integrate.odeint`, :attr:`meteo` must provide data for t = `stop_time` + 1.
                     For the same reason, the :attr:`PAR` of each organ of :attr:`organs` must also provide data for t = `stop_time` + 1.
                     This is automatically checked by the current function.

        """

        self._check_inputs_consistency(start_time, stop_time)

        t = np.linspace(start_time, stop_time, number_of_output_steps)

        self.show_progressbar = show_progressbar
        if self.show_progressbar:
            self.progressbar.set_t_max(stop_time)

        self._clear_organs()

        soln, infodict = odeint(self._calculate_all_derivatives, self.initial_conditions, t, (photosynthesis_computation_interval,), full_output=True, mxstep=odeint_mxstep)

        if not set(infodict['mused']).issubset([1,2]): # I'm not sure if this test is robust or not... Beware especially when scipy is updated.
            raise CNWheatRunError("Integration failed. See the logs of lsoda or try to increase the value of 'mxstep'.")

        cnwheat_output_df = self._format_solver_output(t, soln)

        return cnwheat_output_df


    def _check_inputs_consistency(self, start_time, stop_time):
        """
        Check the consistency between:
        
            * `start_time`, `stop_time` and :attr:`meteo`,
            * `start_time`, `stop_time` and the :attr:`PAR` of each organ of :attr:`organs`.
            
        Raise an exception if :attr:`meteo` or :attr:`PAR` is not valid.
        """

        # check the consistency of meteo
        lowest_t = self.meteo.first_valid_index()
        if start_time < lowest_t:
            raise CNWheatInputError('the lowest t ({}) in meteo data is greater than start_time ({}).'.format(lowest_t, start_time))

        solver_upper_boundary = stop_time + 1
        highest_t = self.meteo.last_valid_index()
        if highest_t < solver_upper_boundary:
            raise CNWheatInputError("""the highest t ({}) in meteo data is lower than stop_time + 1 = {}.
                            scipy.integrate.odeint requires the highest t to be equal or
                            greater than stop_time + 1""".format(highest_t, solver_upper_boundary))

        # check the consistency of the PAR
        for organ in self.PAR_linear_interpolation.iterkeys():
            lowest_t = organ.PAR.first_valid_index()
            if start_time < lowest_t:
                raise CNWheatInputError('the lowest t ({}) in the PAR of {} is greater than start_time ({}).'.format(lowest_t, organ.name, start_time))
            highest_t = organ.PAR.last_valid_index()
            if highest_t < solver_upper_boundary:
                raise CNWheatInputError("""the highest t ({}) in the PAR of {} is lower than stop_time + 1 = {}.
                                scipy.integrate.odeint requires the highest t to be equal or
                                greater than stop_time + 1""".format(highest_t, organ.name, solver_upper_boundary))


    def _clear_organs(self):
        """
        Clear the computed photosynthesis and transpiration for each photosynthetic organ in :attr:`photosynthesis_mapping`.
        """
        for organ_photosynthesis_mapping in self.photosynthesis_mapping.itervalues():
            organ_photosynthesis_mapping.clear()


    def _calculate_all_derivatives(self, y, t, photosynthesis_computation_interval):
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

            - `photosynthesis_computation_interval` (:class:`int`) - The interval which defined the time steps
              at which photosynthesis is computed. For example, if `photosynthesis_computation_interval` = 4,
              then photosynthesis is computed at t=0, t=4, t=8, ...
              This permits to save computation time for large simulation.
              If `photosynthesis_computation_interval` = 0 (the default), then photosynthesis is computed
              for each time step demanded by the solver.

        :Returns:
            The derivatives of `y` at `t`.

        :Returns Type:
            :class:`list`


        """
        #print 't = ', t
        y_isnan = np.isnan(y)
        if y_isnan.any():
            y_isnan_indices = np.where(y_isnan)
            nan_compartments = []
            for organ, compartments in self.initial_conditions_mapping.iteritems():
                for compartment_name, compartment_index in compartments.iteritems():
                    if np.in1d(y_isnan_indices, compartment_index).any():
                        nan_compartments.append(organ.name + '.' + compartment_name)
            raise CNWheatRunError('nan_compartments({}) = {}'.format(t, (nan_compartments)))

        y_derivatives = np.zeros_like(y)

        sucrose_phloem = y[self.initial_conditions_mapping[self.phloem]['sucrose']]
        amino_acids_phloem = y[self.initial_conditions_mapping[self.phloem]['amino_acids']]

        nitrates_roots = y[self.initial_conditions_mapping[self.roots]['nitrates']]
        amino_acids_roots = y[self.initial_conditions_mapping[self.roots]['amino_acids']]
        
        # calculate the time step at which we want to compute the photosynthesis
        if photosynthesis_computation_interval == 0: # i.e. we want to compute the photosynthesis for each time step required by the solver
            t_inf = t
        else:
            t_inf  = t // photosynthesis_computation_interval * photosynthesis_computation_interval

        # call the photosynthesis model
        for organ, PAR_linear_interpolation in self.PAR_linear_interpolation.iteritems():
            organ_photosynthesis_mapping = self.photosynthesis_mapping[organ]
            # calculate the photosynthesis of organ only if it has not been already calculated at t_inf
            if t_inf not in organ_photosynthesis_mapping:
                organ_photosynthesis_mapping[t_inf] = {}
                air_temperature_t_inf = self.meteo_interpolations['air_temperature'](t_inf)
                humidity_t_inf = self.meteo_interpolations['humidity'](t_inf)
                ambient_CO2_t_inf = self.meteo_interpolations['ambient_CO2'](t_inf)
                Wind_top_canopy_t_inf = self.meteo_interpolations['Wind'](t_inf)
                An, Tr = photosynthesis_model.PhotosynthesisModel.calculate_An(t_inf, organ.width, organ.height, PAR_linear_interpolation, 
                                                                               air_temperature_t_inf, ambient_CO2_t_inf, 
                                                                               humidity_t_inf, Wind_top_canopy_t_inf) # TODO: add dependancy to nitrogen
                #print 't=', t, 'PAR=', PAR_linear_interpolation(t_inf), 'Tr_{} = '.format(organ.name), Tr
                organ_photosynthesis_mapping[t_inf]['An'] = An
                organ_photosynthesis_mapping[t_inf]['Tr'] = Tr

        # compute the total transpiration at t_inf
        total_transpiration = 0.0
        transpiration_mapping = dict.fromkeys(self.photosynthesis_mapping)
        for organ, mapping in self.photosynthesis_mapping.iteritems():
            transpiration_mapping[organ] = organ.calculate_transpiration(t_inf, mapping[t_inf]['Tr'])
            total_transpiration += transpiration_mapping[organ]
        #print 'total_transpiration', total_transpiration

        # compute the flows from/to the roots to/from photosynthetic organs  
        conc_nitrates_soil = self.roots.calculate_conc_nitrates_soil(t)
        roots_uptake_nitrate, potential_roots_uptake_nitrates = self.roots.calculate_uptake_nitrates(conc_nitrates_soil, nitrates_roots, total_transpiration)
        roots_export_amino_acids = self.roots.calculate_export_amino_acids(amino_acids_roots, total_transpiration)
        #print 'roots_export_amino_acids',roots_export_amino_acids

        # compute the derivative of each photosynthetic organ compartment 
        for organ in self.PAR_linear_interpolation.iterkeys():
            starch = y[self.initial_conditions_mapping[organ]['starch']]
            sucrose = y[self.initial_conditions_mapping[organ]['sucrose']]
            triosesP = y[self.initial_conditions_mapping[organ]['triosesP']]
            fructan = y[self.initial_conditions_mapping[organ]['fructan']]
            nitrates = y[self.initial_conditions_mapping[organ]['nitrates']]
            amino_acids = y[self.initial_conditions_mapping[organ]['amino_acids']]
            proteins = y[self.initial_conditions_mapping[organ]['proteins']]

            # intermediate variables
            photosynthesis = organ.calculate_photosynthesis(t_inf, self.photosynthesis_mapping[organ][t_inf]['An'])
            organ_transpiration = transpiration_mapping[organ]

            # flows
            s_starch = organ.calculate_s_starch(triosesP)
            d_starch = organ.calculate_d_starch(starch)
            s_sucrose = organ.calculate_s_sucrose(triosesP)
            organ.loading_sucrose = organ.calculate_loading_sucrose(sucrose, sucrose_phloem)
            regul_s_fructan = organ.calculate_regul_s_fructan(organ.loading_sucrose)
            d_fructan = organ.calculate_d_fructan(sucrose, fructan)
            s_fructan = organ.calculate_s_fructan(sucrose, regul_s_fructan)
            nitrates_import = organ.calculate_nitrates_import(roots_uptake_nitrate, organ_transpiration, total_transpiration)
            amino_acids_import = organ.calculate_amino_acids_import(roots_export_amino_acids, organ_transpiration, total_transpiration)
            s_amino_acids = organ.calculate_s_amino_acids(nitrates, triosesP)
            #print 'Organ SAA=', organ.name , s_amino_acids, nitrates
            s_proteins = organ.calculate_s_proteins(amino_acids)
            d_proteins = organ.calculate_d_proteins(proteins)
            organ.loading_amino_acids = organ.calculate_loading_amino_acids(amino_acids, amino_acids_phloem)

            # compartments derivatives
            starch_derivative = organ.calculate_starch_derivative(s_starch, d_starch)
            sucrose_derivative = organ.calculate_sucrose_derivative(s_sucrose, d_starch, organ.loading_sucrose, s_fructan, d_fructan)
            triosesP_derivative = organ.calculate_triosesP_derivative(photosynthesis, s_sucrose, s_starch, s_amino_acids)
            fructan_derivative = organ.calculate_fructan_derivative(s_fructan, d_fructan)
            nitrates_derivative = organ.calculate_nitrates_derivative(nitrates_import, s_amino_acids)
            #print 'Organ:', organ.name, 'nitrates_derivative:',nitrates_derivative
            amino_acids_derivative = organ.calculate_amino_acids_derivative(amino_acids_import, s_amino_acids, s_proteins, d_proteins, organ.loading_amino_acids)
            #print 'Organ:', organ.name, 'amino_acids_derivative:',amino_acids_derivative
            proteins_derivative = organ.calculate_proteins_derivative(s_proteins, d_proteins)
            y_derivatives[self.initial_conditions_mapping[organ]['starch']] = starch_derivative
            y_derivatives[self.initial_conditions_mapping[organ]['sucrose']] = sucrose_derivative
            y_derivatives[self.initial_conditions_mapping[organ]['triosesP']] = triosesP_derivative
            y_derivatives[self.initial_conditions_mapping[organ]['fructan']] = fructan_derivative
            y_derivatives[self.initial_conditions_mapping[organ]['nitrates']] = nitrates_derivative
            y_derivatives[self.initial_conditions_mapping[organ]['amino_acids']] = amino_acids_derivative
            y_derivatives[self.initial_conditions_mapping[organ]['proteins']] = proteins_derivative


        # compute the derivative of each compartment of grains
        grains_structure = y[self.initial_conditions_mapping[self.grains]['structure']]
        grains_starch = y[self.initial_conditions_mapping[self.grains]['starch']]
        grains_proteins = y[self.initial_conditions_mapping[self.grains]['proteins']]

        # intermediate variables
        RGR_structure = self.grains.calculate_RGR_structure(sucrose_phloem)

        # flows
        self.grains.s_grain_structure = self.grains.calculate_s_grain_structure(t, grains_structure, RGR_structure)
        self.grains.s_grain_starch = self.grains.calculate_s_grain_starch(t, sucrose_phloem)
        self.grains.s_proteins = self.grains.calculate_s_proteins(self.grains.s_grain_structure, self.grains.s_grain_starch, amino_acids_phloem, sucrose_phloem, grains_structure)

        # compartments derivatives
        structure_derivative = self.grains.calculate_structure_derivative(self.grains.s_grain_structure)
        starch_derivative = self.grains.calculate_starch_derivative(self.grains.s_grain_starch, grains_structure)
        proteins_derivative = self.grains.calculate_proteins_derivative(self.grains.s_proteins)
        y_derivatives[self.initial_conditions_mapping[self.grains]['structure']] = structure_derivative
        y_derivatives[self.initial_conditions_mapping[self.grains]['starch']] = starch_derivative
        y_derivatives[self.initial_conditions_mapping[self.grains]['proteins']] = proteins_derivative

        # compute the derivative of each compartment of roots
        sucrose_roots = y[self.initial_conditions_mapping[self.roots]['sucrose']]
        # flows
        self.roots.unloading_sucrose = self.roots.calculate_unloading_sucrose(sucrose_phloem)
        self.roots.unloading_amino_acids = self.roots.calculate_unloading_amino_acids(amino_acids_phloem)
        self.roots.s_amino_acids = self.roots.calculate_s_amino_acids(nitrates_roots, sucrose_roots)
        #print 'Roots SAA', self.roots.s_amino_acids,nitrates_roots, sucrose_roots

        # compartments derivatives
        sucrose_derivative = self.roots.calculate_sucrose_derivative(self.roots.unloading_sucrose, self.roots.s_amino_acids)
        nitrates_derivative = self.roots.calculate_nitrates_derivative(roots_uptake_nitrate, self.roots.s_amino_acids)
        amino_acids_derivative = self.roots.calculate_amino_acids_derivative(self.roots.unloading_amino_acids, self.roots.s_amino_acids, roots_export_amino_acids)
        #print 'NO3-, AA derivative roots',nitrates_derivative,amino_acids_derivative
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

        return y_derivatives


    def _format_solver_output(self, t, solver_output):
        """
        Create a :class:`pandas.DataFrame` wich columns are:

            * the time grid `t`,
            * the output of the solver `solver_output`,
            * and intermediate and post-processed variables useful for debug and validation.
        """
        solver_output = solver_output.T

        result_items = [('t', t)]
        
        # call the photosynthesis model
        organs_An = dict.fromkeys(self.PAR_linear_interpolation)
        organs_Tr = dict.fromkeys(self.PAR_linear_interpolation)
        for organ, PAR_linear_interpolation in self.PAR_linear_interpolation.iteritems():
            An_Tr_list = []
            for t_ in t:
                An_Tr_list.append(photosynthesis_model.PhotosynthesisModel.calculate_An(
                                  t_,
                                  organ.width, 
                                  organ.height,
                                  PAR_linear_interpolation,
                                  self.meteo_interpolations['air_temperature'](t_),
                                  self.meteo_interpolations['ambient_CO2'](t_),
                                  self.meteo_interpolations['humidity'](t_),
                                  self.meteo_interpolations['Wind'](t_)))
            An_Tr = np.array(An_Tr_list)
            organs_An[organ] = An_Tr[:,0]
            organs_Tr[organ] = An_Tr[:,1]
            
        # compute the total transpiration
        total_transpiration = np.zeros_like(t)
        transpiration_mapping = dict.fromkeys(self.photosynthesis_mapping)
        for organ in self.photosynthesis_mapping.iterkeys():
            transpiration_mapping[organ] = map(organ.calculate_transpiration, t, organs_Tr[organ])
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
        variables = [('Dry_Mass_{}'.format(self.roots.name), self.roots.calculate_dry_mass(sucrose_roots)),
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
        for organ in self.PAR_linear_interpolation.iterkeys():
            triosesP = solver_output[self.initial_conditions_mapping[organ]['triosesP']]
            starch = solver_output[self.initial_conditions_mapping[organ]['starch']]
            sucrose = solver_output[self.initial_conditions_mapping[organ]['sucrose']]
            fructan = solver_output[self.initial_conditions_mapping[organ]['fructan']]
            loading_sucrose = map(organ.calculate_loading_sucrose, sucrose, sucrose_phloem)
            regul_s_fructan = map(organ.calculate_regul_s_fructan, loading_sucrose)

            nitrates = solver_output[self.initial_conditions_mapping[organ]['nitrates']]
            amino_acids = solver_output[self.initial_conditions_mapping[organ]['amino_acids']]
            proteins = solver_output[self.initial_conditions_mapping[organ]['proteins']]
            
            variables = [('Photosynthesis_Surfacic_Rate_{}'.format(organ.name), organs_An[organ]),
                         ('Tr_{}'.format(organ.name), organs_Tr[organ]),
                         ('Photosynthesis_{}'.format(organ.name), map(organ.calculate_photosynthesis, t, organs_An[organ])),
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

        variables = [('Dry_Mass_{}'.format(self.grains.name), self.grains.calculate_dry_mass(structure_grains, starch_grains)),
                     ('RGR_Structure_{}'.format(self.grains.name), RGR_structure_grains),
                     ('Proteins_Mass_{}'.format(self.grains.name), self.grains.calculate_protein_mass(proteins_grains, structure_grains))]

        s_grain_structure = map(self.grains.calculate_s_grain_structure, t, structure_grains, RGR_structure_grains)
        s_grain_starch = map(self.grains.calculate_s_grain_starch, t, sucrose_phloem)

        flows = [('S_grain_structure{}'.format(self.grains.name), s_grain_structure),
                 ('S_grain_starch_{}'.format(self.grains.name), s_grain_starch),
                 ('S_Proteins_{}'.format(self.grains.name), map(self.grains.calculate_s_proteins, s_grain_structure, s_grain_starch, amino_acids_phloem, sucrose_phloem, structure_grains))]

        compartments = [('Structure_{}'.format(self.grains.name), structure_grains),
                        ('Starch_{}'.format(self.grains.name), starch_grains),
                        ('Proteins_{}'.format(self.grains.name), proteins_grains)]

        result_items.extend(variables + flows + compartments)

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

