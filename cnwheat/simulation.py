# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    cnwheat.cnwheat
    ~~~~~~~~~~~~~~~

    Front-end to run the model.

    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.
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
class CNWheatInputError(CNWheatError): pass
class CNWheatRunError(CNWheatError): pass

class CNWheat(object):
    """
    Compute CN exchanges in wheat architecture defined by lamina, sheaths, internodes, peduncles, a chaff, a phloem, roots and grains.
    This class permits to initialize and run the model.

    :Parameters:

        - `organs` (:class:`list`) - List of :class:`cnwheat.model.Organ`.

        - `meteo` (:class:`pandas.DataFrame`) - a :class:`pandas.DataFrame` which index
          is time in hours and which columns are "air_temperature", "humidity", "ambient_CO2" and "Wind". Due to the
          solver of :func:`scipy.integrate.odeint`, `meteo` must provide data for t = `stop_time` + 1.

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
        for organ_ in organs:
            if isinstance(organ_, model.PhotosyntheticOrgan):
                self.PAR_linear_interpolation[organ_] = interp1d(organ_.PAR.index, organ_.PAR)
                self.photosynthesis_mapping[organ_] = {}
            elif isinstance(organ_, model.Phloem):
                self.phloem = organ_ # the phloem
            elif isinstance(organ_, model.Roots):
                self.roots = organ_ # the roots
            elif isinstance(organ_, model.Grains):
                self.grains = organ_ # the grains
            self.initial_conditions_mapping[organ_] = {}
            for compartment_name, compartment_initial_value in organ_.initial_conditions.iteritems():
                self.initial_conditions_mapping[organ_][compartment_name] = i
                self.initial_conditions.append(compartment_initial_value)
                i += 1

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
              If `photosynthesis_computation_interval`=0 (the default), then photosynthesis
              is computed for each time step demanded by the solver.

            - `odeint_mxstep` (:class:`int`) - Maximum number of (internally defined) steps allowed for each integration point in time grid.
              `odeint_mxstep` is passed to :func:`scipy.integrate.odeint` as `mxstep`. If `odeint_mxstep` = 0, then `mxstep` is determined by the solver.
              The default value ( `5000` ) normally permits to solve the current model. User should increased this value if a more complex model is defined
              and if this model make the integration failed.

            - `show_progressbar` (:class:`bool`) - True: show the progress bar ; False: DO NOT show the progress bar.

        :Returns:
            Dataframe containing the CN exchanges between organs
            for each desired time.

        :Returns Type:
            :class:`pandas.DataFrame`

        .. warning:: due to the solver of :func:`scipy.integrate.odeint`, `meteo` must provide data for t = `stop_time` + 1.
                     For the same reason, the attribute `PAR` of each organ of `organs` must also provide data for t = `stop_time` + 1.
                     This is automatically checked by the current function.

        .. seealso:: Barillot et al. 2014.

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
        Check the consistency of meteo and photosynthesis data.
        Raise an exception if meteo or photosynthesis data is not valid.
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
        for organ_ in self.PAR_linear_interpolation.iterkeys():
            lowest_t = organ_.PAR.first_valid_index()
            if start_time < lowest_t:
                raise CNWheatInputError('the lowest t ({}) in the PAR of {} is greater than start_time ({}).'.format(lowest_t, organ_.name, start_time))
            highest_t = organ_.PAR.last_valid_index()
            if highest_t < solver_upper_boundary:
                raise CNWheatInputError("""the highest t ({}) in the PAR of {} is lower than stop_time + 1 = {}.
                                scipy.integrate.odeint requires the highest t to be equal or
                                greater than stop_time + 1""".format(highest_t, organ_.name, solver_upper_boundary))


    def _clear_organs(self):
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
            for organ_, compartments in self.initial_conditions_mapping.iteritems():
                for compartment_name, compartment_index in compartments.iteritems():
                    if np.in1d(y_isnan_indices, compartment_index).any():
                        nan_compartments.append(organ_.name + '.' + compartment_name)
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

        total_transpiration = 0.0

        for organ_, PAR_linear_interpolation in self.PAR_linear_interpolation.iteritems():
            organ_photosynthesis_mapping = self.photosynthesis_mapping[organ_]
            # calculate the photosynthesis of organ_ only if it has not been already calculated at t_inf
            if t_inf not in organ_photosynthesis_mapping:
                organ_photosynthesis_mapping[t_inf] = {}
                air_temperature_t_inf = self.meteo_interpolations['air_temperature'](t_inf)
                humidity_t_inf = self.meteo_interpolations['humidity'](t_inf)
                ambient_CO2_t_inf = self.meteo_interpolations['ambient_CO2'](t_inf)
                Wind_top_canopy_inf = self.meteo_interpolations['Wind'](t_inf)
                An, Tr = photosynthesis_model.PhotosynthesisModel.calculate_An(t_inf, organ_, PAR_linear_interpolation, 
                                                                         air_temperature_t_inf, ambient_CO2_t_inf, humidity_t_inf, Wind_top_canopy_inf) # TODO: add dependancy to nitrogen
                #print 't=', t, 'PAR=', PAR_linear_interpolation(t_inf), 'Tr_{} = '.format(organ_.name), Tr
                organ_photosynthesis_mapping[t_inf]['photosynthesis'] = organ_.calculate_photosynthesis(t_inf, An)
                organ_photosynthesis_mapping[t_inf]['transpiration'] = organ_.calculate_transpiration(t_inf, Tr)
            organ_transpiration = organ_photosynthesis_mapping[t_inf]['transpiration']
            total_transpiration += organ_transpiration
        #print 'total_transpiration', total_transpiration

        conc_nitrates_soil = self.roots.calculate_conc_nitrates_soil(t)
        roots_uptake_nitrate, potential_roots_uptake_nitrates = self.roots.calculate_uptake_nitrates(conc_nitrates_soil, nitrates_roots, total_transpiration)
        roots_export_amino_acids = self.roots.calculate_export_amino_acids(amino_acids_roots, total_transpiration)
        #print 'roots_export_amino_acids',roots_export_amino_acids

        for organ_ in self.PAR_linear_interpolation.iterkeys():
            starch = y[self.initial_conditions_mapping[organ_]['starch']]
            sucrose = y[self.initial_conditions_mapping[organ_]['sucrose']]
            triosesP = y[self.initial_conditions_mapping[organ_]['triosesP']]
            fructan = y[self.initial_conditions_mapping[organ_]['fructan']]
            nitrates = y[self.initial_conditions_mapping[organ_]['nitrates']]
            amino_acids = y[self.initial_conditions_mapping[organ_]['amino_acids']]
            proteins = y[self.initial_conditions_mapping[organ_]['proteins']]

            # intermediate variables
            photosynthesis_ = self.photosynthesis_mapping[organ_][t_inf]['photosynthesis']
            organ_transpiration = self.photosynthesis_mapping[organ_][t_inf]['transpiration']

            # flows
            s_starch = organ_.calculate_s_starch(triosesP)
            d_starch = organ_.calculate_d_starch(starch)
            s_sucrose = organ_.calculate_s_sucrose(triosesP)
            organ_.loading_sucrose = organ_.calculate_loading_sucrose(sucrose, sucrose_phloem)
            regul_s_fructan = organ_.calculate_regul_s_fructan(organ_.loading_sucrose)
            d_fructan = organ_.calculate_d_fructan(sucrose, fructan)
            s_fructan = organ_.calculate_s_fructan(sucrose, regul_s_fructan)
            nitrates_import = organ_.calculate_nitrates_import(roots_uptake_nitrate, organ_transpiration, total_transpiration)
            amino_acids_import = organ_.calculate_amino_acids_import(roots_export_amino_acids, organ_transpiration, total_transpiration)
            s_amino_acids = organ_.calculate_s_amino_acids(nitrates, triosesP)
            #print 'Organ SAA=', organ_.name , s_amino_acids, nitrates
            s_proteins = organ_.calculate_s_proteins(amino_acids)
            d_proteins = organ_.calculate_d_proteins(proteins)
            organ_.loading_amino_acids = organ_.calculate_loading_amino_acids(amino_acids, amino_acids_phloem)

            # compartments derivatives
            starch_derivative = organ_.calculate_starch_derivative(s_starch, d_starch)
            sucrose_derivative = organ_.calculate_sucrose_derivative(s_sucrose, d_starch, organ_.loading_sucrose, s_fructan, d_fructan)
            triosesP_derivative = organ_.calculate_triosesP_derivative(photosynthesis_, s_sucrose, s_starch, s_amino_acids)
            fructan_derivative = organ_.calculate_fructan_derivative(s_fructan, d_fructan)
            nitrates_derivative = organ_.calculate_nitrates_derivative(nitrates_import, s_amino_acids)
            #print 'Organ:', organ_.name, 'nitrates_derivative:',nitrates_derivative
            amino_acids_derivative = organ_.calculate_amino_acids_derivative(amino_acids_import, s_amino_acids, s_proteins, d_proteins, organ_.loading_amino_acids)
            #print 'Organ:', organ_.name, 'amino_acids_derivative:',amino_acids_derivative
            proteins_derivative = organ_.calculate_proteins_derivative(s_proteins, d_proteins)
            y_derivatives[self.initial_conditions_mapping[organ_]['starch']] = starch_derivative
            y_derivatives[self.initial_conditions_mapping[organ_]['sucrose']] = sucrose_derivative
            y_derivatives[self.initial_conditions_mapping[organ_]['triosesP']] = triosesP_derivative
            y_derivatives[self.initial_conditions_mapping[organ_]['fructan']] = fructan_derivative
            y_derivatives[self.initial_conditions_mapping[organ_]['nitrates']] = nitrates_derivative
            y_derivatives[self.initial_conditions_mapping[organ_]['amino_acids']] = amino_acids_derivative
            y_derivatives[self.initial_conditions_mapping[organ_]['proteins']] = proteins_derivative


        # Grains
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

        # Roots
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

        # Phloem
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
            * and intermediate and post-processed variables usefull for debug and validation
        """
        solver_output = solver_output.T

        result_items = [('t', t)]

        # phloem
        sucrose_phloem = solver_output[self.initial_conditions_mapping[self.phloem]['sucrose']]
        amino_acids_phloem = solver_output[self.initial_conditions_mapping[self.phloem]['amino_acids']]
        variables = [('Conc_Sucrose_{}'.format(self.phloem.name).rstrip('_'), self.phloem.calculate_conc_sucrose(sucrose_phloem)),
                     ('Conc_C_Sucrose_{}'.format(self.phloem.name).rstrip('_'), self.phloem.calculate_conc_c_sucrose(sucrose_phloem)),
                     ('Conc_Amino_Acids_{}'.format(self.phloem.name).rstrip('_'), self.phloem.calculate_conc_amino_acids(amino_acids_phloem))]

        compartments = [('Sucrose_{}'.format(self.phloem.name).rstrip('_'), sucrose_phloem),
                        ('Amino_Acids_{}'.format(self.phloem.name).rstrip('_'), amino_acids_phloem)]
        result_items.extend(variables + compartments)

        # Photosynthesis and transpiration
        total_transpiration = 0.0

        organs_An = []
        organs_photosynthesis = []
        organs_transpiration = []
        organs_transpiration_DF = pd.DataFrame()
        for organ_, PAR_linear_interpolation in self.PAR_linear_interpolation.iteritems():
            An_Tr_list = []
            for t_ in t:
                An_Tr_list.append(photosynthesis_model.PhotosynthesisModel.calculate_An(
                                      t_,
                                      organ_,
                                      PAR_linear_interpolation,
                                      self.meteo_interpolations['air_temperature'](t_),
                                      self.meteo_interpolations['ambient_CO2'](t_),
                                      self.meteo_interpolations['humidity'](t_),
                                      self.meteo_interpolations['Wind'](t_)))
            An_Tr = np.array(An_Tr_list)

            # Photosynthesis
            photosynthesis_ = map(organ_.calculate_photosynthesis, t, An_Tr[:,0])
            organs_An.append(('Photosynthesis_Surfacic_Rate_{}'.format(organ_.name).rstrip('_'), An_Tr[:,0]))
            organs_photosynthesis.append(('Photosynthesis_{}'.format(organ_.name).rstrip('_'), photosynthesis_))
            # Transpiration
            organ_transpiration = map(organ_.calculate_transpiration, t, An_Tr[:,1])
            organs_transpiration_DF[organ_.name] = organ_transpiration
            organs_transpiration.append(('Transpiration_{}'.format(organ_.name).rstrip('_'), organ_transpiration))

        total_transpiration = organs_transpiration_DF.sum(axis=1)
        total_transpiration_out = [('Total_transpiration', total_transpiration)]
        result_items.extend(organs_An + organs_photosynthesis + organs_transpiration + total_transpiration_out)

        # Roots
        sucrose_roots = solver_output[self.initial_conditions_mapping[self.roots]['sucrose']]
        unloading_sucrose = map(self.roots.calculate_unloading_sucrose, sucrose_phloem)
        nitrates_roots = solver_output[self.initial_conditions_mapping[self.roots]['nitrates']]
        conc_nitrates_soil = map(self.roots.calculate_conc_nitrates_soil,t)
        uptake_nitrates, potential_roots_uptake_nitrates = self.roots.calculate_uptake_nitrates(conc_nitrates_soil, nitrates_roots, total_transpiration)
        amino_acids_roots = solver_output[self.initial_conditions_mapping[self.roots]['amino_acids']]
        roots_export_amino_acids =  map(self.roots.calculate_export_amino_acids, amino_acids_roots, total_transpiration)

        variables = [('Dry_Mass_{}'.format(self.roots.name).rstrip('_'), self.roots.calculate_dry_mass(sucrose_roots)),
                     ('Conc_Nitrates_{}'.format(self.roots.name).rstrip('_'), self.roots.calculate_conc_nitrates(nitrates_roots)),
                     ('Conc_Amino_Acids_{}'.format(self.roots.name).rstrip('_'), self.roots.calculate_conc_amino_acids(amino_acids_roots)),
                     ('Conc_Nitrates_Soil_{}'.format(self.roots.name).rstrip('_'), conc_nitrates_soil)]

        flows = [('Unloading_Sucrose_{}'.format(self.roots.name).rstrip('_'), unloading_sucrose),
                 ('Unloading_Amino_Acids_{}'.format(self.roots.name).rstrip('_'), map(self.roots.calculate_unloading_amino_acids, amino_acids_phloem)),
                 ('Uptake_Nitrates_{}'.format(self.roots.name).rstrip('_'), uptake_nitrates),
                 ('Potential_Uptake_Nitrates_{}'.format(self.roots.name).rstrip('_'), potential_roots_uptake_nitrates),
                 ('Export_Amino_Acids_{}'.format(self.roots.name).rstrip('_'),roots_export_amino_acids),
                 ('S_Amino_Acids_{}'.format(self.roots.name).rstrip('_'), map(self.roots.calculate_s_amino_acids, nitrates_roots, sucrose_roots))]

        compartments = [('Sucrose_{}'.format(self.roots.name).rstrip('_'), sucrose_roots),
                        ('Nitrates_{}'.format(self.roots.name).rstrip('_'), nitrates_roots),
                        ('Amino_Acids_{}'.format(self.roots.name).rstrip('_'), amino_acids_roots)]

        result_items.extend(variables + flows + compartments)


        # Photosynthetic organs
        for organ_ in self.PAR_linear_interpolation.iterkeys():
            triosesP = solver_output[self.initial_conditions_mapping[organ_]['triosesP']]
            starch = solver_output[self.initial_conditions_mapping[organ_]['starch']]
            sucrose = solver_output[self.initial_conditions_mapping[organ_]['sucrose']]
            fructan = solver_output[self.initial_conditions_mapping[organ_]['fructan']]
            loading_sucrose = map(organ_.calculate_loading_sucrose, sucrose, sucrose_phloem)
            regul_s_fructan = map(organ_.calculate_regul_s_fructan, loading_sucrose)

            nitrates = solver_output[self.initial_conditions_mapping[organ_]['nitrates']]
            amino_acids = solver_output[self.initial_conditions_mapping[organ_]['amino_acids']]
            proteins = solver_output[self.initial_conditions_mapping[organ_]['proteins']]

            variables = [('Conc_TriosesP_{}'.format(organ_.name).rstrip('_'), organ_.calculate_conc_triosesP(triosesP)),
                         ('Conc_Starch_{}'.format(organ_.name).rstrip('_'), organ_.calculate_conc_starch(starch)),
                         ('Conc_Sucrose_{}'.format(organ_.name).rstrip('_'), organ_.calculate_conc_sucrose(sucrose)),
                         ('Conc_Fructan_{}'.format(organ_.name).rstrip('_'), organ_.calculate_conc_fructan(fructan)),
                         ('Regul_S_Fructan_{}'.format(organ_.name).rstrip('_'), regul_s_fructan),
                         ('Conc_Nitrates_{}'.format(organ_.name).rstrip('_'), organ_.calculate_conc_nitrates(nitrates)),
                         ('Conc_Amino_Acids_{}'.format(organ_.name).rstrip('_'), organ_.calculate_conc_amino_acids(amino_acids)),
                         ('Conc_Proteins_{}'.format(organ_.name).rstrip('_'), organ_.calculate_conc_proteins(proteins))]


            flows = [('S_Starch_{}'.format(organ_.name).rstrip('_'), map(organ_.calculate_s_starch, triosesP)),
                     ('D_Starch_{}'.format(organ_.name).rstrip('_'), map(organ_.calculate_d_starch, starch)),
                     ('S_Sucrose_{}'.format(organ_.name).rstrip('_'), map(organ_.calculate_s_sucrose, triosesP)),
                     ('Loading_Sucrose_{}'.format(organ_.name).rstrip('_'), loading_sucrose),
                     ('S_Fructan_{}'.format(organ_.name).rstrip('_'), map(organ_.calculate_s_fructan, sucrose, regul_s_fructan)),
                     ('D_Fructan_{}'.format(organ_.name).rstrip('_'), map(organ_.calculate_d_fructan, sucrose, fructan)),
                     ('Nitrates_import_{}'.format(organ_.name).rstrip('_'), map(organ_.calculate_nitrates_import, uptake_nitrates, organs_transpiration_DF[organ_.name], total_transpiration)),
                     ('Amino_Acids_import_{}'.format(organ_.name).rstrip('_'), map(organ_.calculate_amino_acids_import, roots_export_amino_acids, organs_transpiration_DF[organ_.name], total_transpiration)),
                     ('S_Amino_Acids_{}'.format(organ_.name).rstrip('_'), map(organ_.calculate_s_amino_acids, nitrates, triosesP)),
                     ('S_Proteins_{}'.format(organ_.name).rstrip('_'), map(organ_.calculate_s_proteins, amino_acids)),
                     ('D_Proteins_{}'.format(organ_.name).rstrip('_'), map(organ_.calculate_d_proteins, proteins)),
                     ('Loading_Amino_Acids_{}'.format(organ_.name).rstrip('_'), map(organ_.calculate_loading_amino_acids, amino_acids, amino_acids_phloem))]

            compartments = [('TriosesP_{}'.format(organ_.name).rstrip('_'), triosesP),
                            ('Starch_{}'.format(organ_.name).rstrip('_'), starch),
                            ('Sucrose_{}'.format(organ_.name).rstrip('_'), sucrose),
                            ('Fructan_{}'.format(organ_.name).rstrip('_'), fructan),
                            ('Nitrates_{}'.format(organ_.name).rstrip('_'), nitrates),
                            ('Amino_Acids_{}'.format(organ_.name).rstrip('_'), amino_acids),
                            ('Proteins_{}'.format(organ_.name).rstrip('_'), proteins)]

            result_items.extend(variables + flows + compartments)

        # Grains
        structure_grains = solver_output[self.initial_conditions_mapping[self.grains]['structure']]
        starch_grains = solver_output[self.initial_conditions_mapping[self.grains]['starch']]
        proteins_grains = solver_output[self.initial_conditions_mapping[self.grains]['proteins']]

        RGR_structure_grains = map(self.grains.calculate_RGR_structure, sucrose_phloem)

        variables = [('Dry_Mass_{}'.format(self.grains.name).rstrip('_'), self.grains.calculate_dry_mass(structure_grains, starch_grains)),
                     ('RGR_Structure_{}'.format(self.grains.name).rstrip('_'), RGR_structure_grains),
                     ('Proteins_Mass_{}'.format(self.grains.name).rstrip('_'), self.grains.calculate_protein_mass(proteins_grains, structure_grains))]

        s_grain_structure = map(self.grains.calculate_s_grain_structure, t, structure_grains, RGR_structure_grains)
        s_grain_starch = map(self.grains.calculate_s_grain_starch, t, sucrose_phloem)

        flows = [('S_grain_structure{}'.format(self.grains.name).rstrip('_'), s_grain_structure),
                 ('S_grain_starch_{}'.format(self.grains.name).rstrip('_'), s_grain_starch),
                 ('S_Proteins_{}'.format(self.grains.name).rstrip('_'), map(self.grains.calculate_s_proteins, s_grain_structure, s_grain_starch, amino_acids_phloem, sucrose_phloem, structure_grains))]

        compartments = [('Structure_{}'.format(self.grains.name).rstrip('_'), structure_grains),
                        ('Starch_{}'.format(self.grains.name).rstrip('_'), starch_grains),
                        ('Proteins_{}'.format(self.grains.name).rstrip('_'), proteins_grains)]

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

