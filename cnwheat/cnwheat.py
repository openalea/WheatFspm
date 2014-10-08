# -*- coding: latin-1 -*-
"""
    cnwheat.cnwheat
    ~~~~~~~~~~~~~~~

    Function to run the model.

    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.
"""

from __future__ import division # use "//" to do integer division

import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.interpolate import interp1d

import organ
import photosynthesis

def _calculate_all_derivatives(y, t, phloem, organs_without_phloem, meteo_interpolations, photosynthesis_computation_interval):
    """Compute the derivative of `y` at `t`.
    
    :func:`_calculate_all_derivatives` is passed as **func** argument to 
    :func:`odeint(func, y0, t, args=(),...) <scipy.integrate.odeint>`.
    :func:`_calculate_all_derivatives` is called automatically by 
    :func:`scipy.integrate.odeint`.
    
    First call to :func:`_calculate_all_derivatives` uses `y` = **y0** and 
    `t` = **t** [0], where **y0** and **t** are arguments passed to :func:`odeint(func, y0, t, args=(),...) <scipy.integrate.odeint>`.
    
    Following calls to :func:`_calculate_all_derivatives` use `t` in [min( **t** ), max( **t** ) + 1] where 
    **t** is an argument passed to :func:`odeint(func, y0, t, args=(),...) <scipy.integrate.odeint>`. `y` is 
    computed automatically by the solver.
    
    :Parameters:

        - `y` (:class:`list`) - The current y values. The y values of `phloem` 
          must appear first. `y` is automatically set by 
          :func:`scipy.integrate.odeint`. User does not have control over `y`. 
          At first call to :func:`_calculate_all_derivatives` by :func:`scipy.integrate.odeint`, `y` = **y0** 
          where **y0** is one of the arguments passed to :func:`odeint(func, y0, t, args=(),...) <scipy.integrate.odeint>`.
          For each following call to :func:`_calculate_all_derivatives`, `y` is 
          computed automatically by the solver.

        - `t` (:class:`float`) - The current t at which we want to compute the derivatives. 
          `t` is automatically set by :func:`scipy.integrate.odeint`. 
          User does not have control over `t`.
          At first call to :func:`_calculate_all_derivatives` :func:`scipy.integrate.odeint`, 
          `t` = **t** [0], where **t** is one of the arguments passed to :func:`odeint(func, y0, t, args=(),...) <scipy.integrate.odeint>`.
          For each following call to :func:`_calculate_all_derivatives`, `t` belongs 
          to the interval [min( **t** ), max( **t** ) + 1], where **t** is an 
          argument passed to :func:`odeint(func, y0, t, args=(),...) <scipy.integrate.odeint>`.

        - `phloem` (:class:`cnwheat.organ.Phloem`) - The phloem.

        - `organs_without_phloem` (:class:`list`) - List of :class:`cnwheat.organ.Organ`. 
          Does not include the phloem.
        
        - `meteo_interpolations` (:class:`dict`) - A dictionary which contains, 
          for each data of meteo data, a function which permits to find the value 
          of new points using interpolation. Keys are the name of meteo data (:class:`str`). 
          Values are :class:`scipy.interpolate.interpolate.interp1d`.
        
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
    y_iter = iter(y)
    y_derivatives = []

    sucrose_phloem = y_iter.next()

    for organ_ in organs_without_phloem:

        if isinstance(organ_, organ.PhotosyntheticOrgan):
            storage = y_iter.next()
            sucrose = y_iter.next()
            triosesP = y_iter.next()
            fructan = y_iter.next()
            # calculate the time step at which we want to compute the photosynthesis
            if photosynthesis_computation_interval == 0: # i.e. we want to compute the photosynthesis for each time step demanded by the solver
                t_inf = t
            else:
                t_inf  = t // photosynthesis_computation_interval * photosynthesis_computation_interval
            # calculate the photosynthesis of organ_ only if it has not been already calculated at t_inf
            if t_inf not in organ_.photosynthesis_mapping:
                PAR_t_inf = organ_.PAR_linear_interpolation(t_inf)
                Tac_t_inf = meteo_interpolations['Tac'](t_inf)
                hs_t_inf = meteo_interpolations['hs'](t_inf)
                Ca_t_inf = meteo_interpolations['Ca'](t_inf)
                An = photosynthesis.PhotosynthesisModel.calculate_An(t_inf, PAR_t_inf, Tac_t_inf, Ca_t_inf, hs_t_inf)
                organ_.photosynthesis_mapping[t_inf] = organ_.calculate_photosynthesis(t_inf, An)
            photosynthesis_ = organ_.photosynthesis_mapping[t_inf]
            # flows
            d_storage = organ_.calculate_d_storage(storage)
            s_storage = organ_.calculate_s_storage(triosesP)
            s_sucrose = organ_.calculate_s_sucrose(triosesP)
            organ_.loading_sucrose = organ_.calculate_loading_sucrose(sucrose, sucrose_phloem)
            regul_s_fructan = organ_.calculate_regul_s_fructan(organ_.loading_sucrose)
            d_fructan = organ_.calculate_d_fructan(sucrose, fructan)
            s_fructan = organ_.calculate_s_fructan(sucrose, regul_s_fructan)
            # compartments derivatives
            storage_derivative = organ_.calculate_storage_derivative(s_storage, d_storage)
            sucrose_derivative = organ_.calculate_sucrose_derivative(s_sucrose, d_storage, organ_.loading_sucrose, s_fructan, d_fructan)
            triosesP_derivative = organ_.calculate_triosesP_derivative(photosynthesis_, s_sucrose, s_storage)
            fructan_derivative = organ_.calculate_fructan_derivative(s_fructan, d_fructan)
            y_derivatives.extend([storage_derivative, sucrose_derivative, triosesP_derivative, fructan_derivative])

        elif isinstance(organ_, organ.Grains):
            storage = y_iter.next()
            structure = y_iter.next()
            # needed variables
            RGR_structure = organ_.calculate_RGR_structure(sucrose_phloem)
            # flows
            organ_.unloading_sucrose_structure = organ_.calculate_unloading_sucrose_structure(t, structure, RGR_structure)
            organ_.unloading_sucrose_storage = organ_.calculate_unloading_sucrose_storage(t, sucrose_phloem)
            organ_.structure = structure
            # compartments derivatives
            storage_derivative = organ_.calculate_storage_derivative(organ_.unloading_sucrose_storage, organ_.structure)
            structure_derivative = organ_.calculate_structure_derivative(organ_.unloading_sucrose_structure)
            y_derivatives.extend([storage_derivative, structure_derivative])

        elif isinstance(organ_, organ.Roots):
            sucrose = y_iter.next()
            # flows
            organ_.unloading_sucrose = organ_.calculate_unloading_sucrose(sucrose_phloem)
            # compartments derivatives
            sucrose_derivative = organ_.calculate_sucrose_derivative(organ_.unloading_sucrose)
            y_derivatives.extend([sucrose_derivative])

    sucrose_phloem_derivative = phloem.calculate_sucrose_derivative(organs_without_phloem)
    y_derivatives.insert(0, sucrose_phloem_derivative)

    return y_derivatives


def run(start_time, stop_time, number_of_output_steps, organs, meteo, photosynthesis_computation_interval=0, odeint_mxstep=5000):
    """
    Compute CN exchanges in wheat architecture defined by lamina, sheaths, internodes, peduncles, a chaff, a phloem, roots and grains

    The computation is done between `start_time` and `stop_time`, for `number_of_output_steps` steps.

    :Parameters:

        - `start_time` (:class:`int`) - The starting of the time grid.

        - `stop_time` (:class:`int`) - The end of the time grid.

        - `number_of_output_steps` (:class:`int`) - Number of time points for which to compute the CN exchanges in the system.

        - `organs` (:class:`list`) - List of :class:`cnwheat.organ.Organ`.

        - `meteo` (:class:`pandas.DataFrame`) - a :class:`pandas.DataFrame` which index 
          is time in hours and which columns are "Tac", "hs" and "Ca". Due to the 
          solver of :func:`scipy.integrate.odeint`, `meteo` must provide data for t = `stop_time` + 1.

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

    # check the consistency of meteo
    lowest_t = meteo.first_valid_index()
    if start_time < lowest_t:
        raise Exception('Error: the lowest t ({}) in meteo data is greater than start_time ({}).'.format(lowest_t, start_time))

    solver_upper_boundary = stop_time + 1
    highest_t = meteo.last_valid_index()
    if highest_t < solver_upper_boundary:
        raise Exception("""Error: the highest t ({}) in meteo data is lower than stop_time + 1 = {}.
                        scipy.integrate.odeint requires the highest t to be equal or
                        greater than stop_time + 1""".format(highest_t, solver_upper_boundary))

    # check the consistency of the PAR
    for organ_ in organs:
        if isinstance(organ_, organ.PhotosyntheticOrgan):
            lowest_t = organ_.PAR.first_valid_index()
            if start_time < lowest_t:
                raise Exception('Error: the lowest t ({}) in the PAR of {} is greater than start_time ({}).'.format(lowest_t, organ_.name, start_time))
            highest_t = organ_.PAR.last_valid_index()
            if highest_t < solver_upper_boundary:
                raise Exception("""Error: the highest t ({}) in the PAR of {} is lower than stop_time + 1 = {}.
                                scipy.integrate.odeint requires the highest t to be equal or
                                greater than stop_time + 1""".format(highest_t, organ_.name, solver_upper_boundary))

    # interpolate meteo data
    meteo_interpolations = {}
    for column in meteo.columns:
        meteo_interpolations[column] = interp1d(meteo.index, meteo[column])

    # interpolate the PAR and construct the list of initial conditions
    initial_conditions = []
    organs_without_phloem = []

    for organ_ in organs:
        if isinstance(organ_, organ.PhotosyntheticOrgan):
            organ_.PAR_linear_interpolation = interp1d(organ_.PAR.index, organ_.PAR)
        if isinstance(organ_, organ.Phloem):
            phloem = organ_
            initial_conditions = organ_.get_initial_conditions() + initial_conditions
        else:
            organs_without_phloem.append(organ_)
            initial_conditions.extend(organ_.get_initial_conditions())

    t = np.linspace(start_time, stop_time, number_of_output_steps)

    soln, infodict = odeint(_calculate_all_derivatives, initial_conditions, t, (phloem, organs_without_phloem, meteo_interpolations, photosynthesis_computation_interval), full_output=True, mxstep=odeint_mxstep)

    if not set(infodict['mused']).issubset([1,2]): # I'm not sure if this test is robust or not... Beware especially when scipy is updated.
        raise Exception("Error: Integration failed. See the logs of lsoda or try to increase the value of 'mxstep'.")

    soln_iter = iter(soln.T)

    # construct the table of the results
    result_items = [('t', t)]

    sucrose_phloem = soln_iter.next()
    variables = [(('Conc_Sucrose_%s' % phloem.name).rstrip('_'), phloem.calculate_conc_sucrose(sucrose_phloem)),
                 (('Conc_C_Sucrose_%s' % phloem.name).rstrip('_'), phloem.calculate_conc_c_sucrose(sucrose_phloem))]
    compartments = [(('Sucrose_%s' % phloem.name).rstrip('_'), sucrose_phloem)]
    result_items.extend(variables + compartments)

    for organ_ in organs_without_phloem:
        if isinstance(organ_, organ.PhotosyntheticOrgan):
            storage = soln_iter.next()
            sucrose = soln_iter.next()
            triosesP = soln_iter.next()
            fructan = soln_iter.next()
            loading_sucrose = map(organ_.calculate_loading_sucrose, sucrose, sucrose_phloem)
            regul_s_fructan = map(organ_.calculate_regul_s_fructan, loading_sucrose)

            An = np.array(map(photosynthesis.PhotosynthesisModel.calculate_An,
                              t,
                              organ_.PAR_linear_interpolation(t),
                              meteo_interpolations['Tac'](t),
                              meteo_interpolations['Ca'](t),
                              meteo_interpolations['hs'](t)))

            variables = [(('Photosynthesis_%s' % organ_.name).rstrip('_'), map(organ_.calculate_photosynthesis, t, An)),
                         (('Conc_TriosesP_%s' % organ_.name).rstrip('_'), organ_.calculate_conc_triosesP(triosesP)),
                         (('Conc_Storage_%s' % organ_.name).rstrip('_'), organ_.calculate_conc_storage(storage)),
                         (('Conc_Sucrose_%s' % organ_.name).rstrip('_'), organ_.calculate_conc_sucrose(sucrose)),
                         (('Conc_Fructan_%s' % organ_.name).rstrip('_'), organ_.calculate_conc_fructan(fructan)),
                         (('Regul_S_Fructan_%s' % organ_.name).rstrip('_'), regul_s_fructan)]
            flows = [(('D_Storage_%s' % organ_.name).rstrip('_'), map(organ_.calculate_d_storage, storage)),
                     (('S_Storage_%s' % organ_.name).rstrip('_'), map(organ_.calculate_s_storage, triosesP)),
                     (('S_Sucrose_%s' % organ_.name).rstrip('_'), map(organ_.calculate_s_sucrose, triosesP)),
                     (('Loading_Sucrose_%s' % organ_.name).rstrip('_'), loading_sucrose),
                     (('D_Fructan_%s' % organ_.name).rstrip('_'), map(organ_.calculate_d_fructan, sucrose, fructan)),
                     (('S_Fructan_%s' % organ_.name).rstrip('_'), map(organ_.calculate_s_fructan, sucrose, regul_s_fructan))]
            compartments = [(('Storage_%s' % organ_.name).rstrip('_'), storage),
                            (('Sucrose_%s' % organ_.name).rstrip('_'), sucrose),
                            (('TriosesP_%s' % organ_.name).rstrip('_'), triosesP),
                            (('Fructan_%s' % organ_.name).rstrip('_'), fructan)]
        elif isinstance(organ_, organ.Grains):
            storage = soln_iter.next()
            structure = soln_iter.next()
            RGR_structure = map(organ_.calculate_RGR_structure, sucrose_phloem)
            variables = [(('Dry_Mass_%s' % organ_.name).rstrip('_'), organ_.calculate_dry_mass(structure, storage)),
                         (('RGR_Structure_%s' % organ_.name).rstrip('_'), RGR_structure)]
            flows = [(('Unloading_Sucrose_Storage_%s' % organ_.name).rstrip('_'), map(organ_.calculate_unloading_sucrose_storage, t, sucrose_phloem)),
                     (('Unloading_Sucrose_Structure_%s' % organ_.name).rstrip('_'), map(organ_.calculate_unloading_sucrose_structure, t, structure, RGR_structure))]
            compartments = [(('Storage_%s' % organ_.name).rstrip('_'), storage),
                            (('Structure_%s' % organ_.name).rstrip('_'), structure)]
        elif isinstance(organ_, organ.Roots):
            sucrose = soln_iter.next()
            variables = [(('Dry_Mass_%s' % organ_.name).rstrip('_'), organ_.calculate_dry_mass(sucrose))]
            flows = [(('Unloading_Sucrose_%s' % organ_.name).rstrip('_'), map(organ_.calculate_unloading_sucrose, sucrose_phloem))]
            compartments = [(('Sucrose_%s' % organ_.name).rstrip('_'), sucrose)]
        result_items.extend(variables + flows + compartments)

    return pd.DataFrame.from_items(result_items)

