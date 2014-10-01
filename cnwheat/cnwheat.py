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

def _calculate_all_derivatives(y, t, phloem, organs, meteo_interpolations, photosynthesis_computation_interval):
    """Compute the derivative of `y` at `t`.
    """
    y_iter = iter(y)
    y_derivatives = []

    SUCROSE_phloem = y_iter.next()

    def _calculate_Photosynthesis(t, PAR_linear_interpolation, photosynthesis_mapping):
        if photosynthesis_computation_interval == 0:
            t_inf = t
        else:
            t_inf  = t // photosynthesis_computation_interval * photosynthesis_computation_interval
        if t_inf not in photosynthesis_mapping:
            PAR_t_inf = PAR_linear_interpolation(t_inf)
            Tac_t_inf = meteo_interpolations['Tac'](t_inf)
            hs_t_inf = meteo_interpolations['hs'](t_inf)
            Ca_t_inf = meteo_interpolations['Ca'](t_inf)
            An = photosynthesis.PhotosynthesisModel.calculate_An(t_inf, PAR_t_inf, Tac_t_inf, Ca_t_inf, hs_t_inf)
            photosynthesis_mapping[t_inf] = organ_.calculate_Photosynthesis(t, An)
        return photosynthesis_mapping[t_inf]

    for organ_ in organs:

        if isinstance(organ_, (organ.Lamina, organ.Chaff)):
            STORAGE = y_iter.next()
            SUCROSE = y_iter.next()
            TRIOSESP = y_iter.next()
            # needed variables
            Photosynthesis = _calculate_Photosynthesis(t, organ_.PAR_linear_interpolation, organ_.photosynthesis_mapping)
            # flows
            D_storage = organ_.calculate_D_storage(STORAGE)
            S_storage = organ_.calculate_S_storage(TRIOSESP)
            S_sucrose = organ_.calculate_S_sucrose(TRIOSESP)
            organ_.Loading_sucrose = organ_.calculate_Loading_sucrose(SUCROSE, SUCROSE_phloem)
            # compartments derivatives
            STORAGE_derivative = organ_.calculate_STORAGE_derivative(S_storage, D_storage)
            SUCROSE_derivative = organ_.calculate_SUCROSE_derivative(S_sucrose, D_storage, organ_.Loading_sucrose)
            TRIOSESP_derivative = organ_.calculate_TRIOSESP_derivative(Photosynthesis, S_sucrose, S_storage)
            y_derivatives.extend([STORAGE_derivative, SUCROSE_derivative, TRIOSESP_derivative])

        elif isinstance(organ_, (organ.Internode, organ.Peduncle, organ.Sheath)):
            FRUCTAN = y_iter.next()
            STORAGE = y_iter.next()
            SUCROSE = y_iter.next()
            TRIOSESP = y_iter.next()
            # needed variables
            Photosynthesis = _calculate_Photosynthesis(t, organ_.PAR_linear_interpolation, organ_.photosynthesis_mapping)
            # flows
            D_storage = organ_.calculate_D_storage(STORAGE)
            S_storage = organ_.calculate_S_storage(TRIOSESP)
            S_sucrose = organ_.calculate_S_sucrose(TRIOSESP)
            organ_.Loading_sucrose = organ_.calculate_Loading_sucrose(SUCROSE, SUCROSE_phloem)
            Regul_Sfructan = organ_.calculate_Regul_Sfructan(organ_.Loading_sucrose)
            D_fructan = organ_.calculate_D_fructan(SUCROSE, FRUCTAN)
            S_fructan = organ_.calculate_S_fructan(SUCROSE, Regul_Sfructan)
            # compartments derivatives
            FRUCTAN_derivative = organ_.calculate_FRUCTAN_derivative(S_fructan, D_fructan)
            STORAGE_derivative = organ_.calculate_STORAGE_derivative(S_storage, D_storage)
            SUCROSE_derivative = organ_.calculate_SUCROSE_derivative(S_sucrose, D_storage, organ_.Loading_sucrose, S_fructan, D_fructan)
            TRIOSESP_derivative = organ_.calculate_TRIOSESP_derivative(Photosynthesis, S_sucrose, S_storage)
            y_derivatives.extend([FRUCTAN_derivative, STORAGE_derivative, SUCROSE_derivative, TRIOSESP_derivative])

        elif isinstance(organ_, organ.Grains):
            STORAGE = y_iter.next()
            STRUCTURE = y_iter.next()
            # needed variables
            RGR_structure = organ_.calculate_RGR_structure(SUCROSE_phloem)
            # flows
            organ_.Unloading_sucrose_structure = organ_.calculate_Unloading_sucrose_structure(t, STRUCTURE, RGR_structure)
            organ_.Unloading_sucrose_storage = organ_.calculate_Unloading_sucrose_storage(t, SUCROSE_phloem)
            organ_.STRUCTURE = STRUCTURE
            # compartments derivatives
            STORAGE_derivative = organ_.calculate_STORAGE_derivative(organ_.Unloading_sucrose_storage, organ_.STRUCTURE)
            STRUCTURE_derivative = organ_.calculate_STRUCTURE_derivative(organ_.Unloading_sucrose_structure)
            y_derivatives.extend([STORAGE_derivative, STRUCTURE_derivative])

        elif isinstance(organ_, organ.Roots):
            Sucrose = y_iter.next()
            # flows
            organ_.Unloading_sucrose = organ_.calculate_Unloading_sucrose(SUCROSE_phloem)
            # compartments derivatives
            Sucrose_derivative = organ_.calculate_Sucrose_derivative(organ_.Unloading_sucrose)
            y_derivatives.extend([Sucrose_derivative])

    SUCROSE_phloem_derivative = phloem.calculate_SUCROSE_derivative(organs)
    y_derivatives.insert(0, SUCROSE_phloem_derivative)

    return y_derivatives


def run(start_time, stop_time, number_of_output_steps, phloem, organs, meteo, photosynthesis_computation_interval=0, odeint_mxstep=5000):
    """
    Compute CN exchanges in wheat architecture defined by lamina, sheaths, internodes, peduncles, a chaff, a phloem, roots and grains

    The computation is done between `start_time` and `stop_time`, for `number_of_output_steps` steps.

    :Parameters:

        - `start_time` (:class:`int`) - The starting of the time grid.

        - `stop_time` (:class:`int`) - The end of the time grid.

        - `number_of_output_steps` (:class:`int`) - Number of time points for which to compute the CN exchanges in the system.

        - `phloem` (:class:`cnwheat.organ.Phloem`) - The phloem of the system.

        - `organs` (:class:`list`) - List of :class:`cnwheat.organ.Organ`. Does not include the phloem.

        - `meteo` (:class:`pandas.DataFrame`) - The meteo data: 'Tac', 'hs', 'Ca'.

        - `photosynthesis_computation_interval` (:class:`int`) - The interval which defined the time steps
          at which photosynthesis is computed. For example, if `photosynthesis_computation_interval`=4,
          then photosynthesis is computed at t=0, t=4, t=8, ...
          This permits to save computation time for large simulation.
          If `photosynthesis_computation_interval`=0 (the default), then photosynthesis is computed for each time step demanded by the integrator.

        - `odeint_mxstep` (:class:`int`) - Maximum number of (internally defined) steps allowed for each integration point in time grid.
          `odeint_mxstep` is passed to :func:`scipy.integrate.odeint` as `mxstep`. If `odeint_mxstep`==0, then `mxstep` is determined by the solver.
          The default value (`5000`) normally permits to solve the current model. User should increased this value if a more complex model is defined
          and if this model make the integration failed.

    :Returns:
        Dataframe containing the CN exchanges between organs
        for each desired time.

    :Returns Type:
        :class:`pandas.DataFrame`

    .. warning:: due to the integrator of :func:`scipy.integrate.odeint`, `meteo` must provide data for t=`stop_time`+1.
                 For the same reason, the attribute `PAR` of each organ of `organs` must also provide data for t=`stop_time`+1.
                 This is automatically checked by the current function.

    .. seealso:: Barillot et al. 2014.

    """

    # check the consistency of meteo
    lowest_t = meteo.first_valid_index()
    if start_time < lowest_t:
        raise Exception('Error: the lowest t ({}) in meteo data is greater than start_time ({}).'.format(lowest_t, start_time))

    integrator_upper_boundary = stop_time + 1
    highest_t = meteo.last_valid_index()
    if highest_t < integrator_upper_boundary:
        raise Exception("""Error: the highest t ({}) in meteo data is lower than stop_time + 1 = {}.
                        scipy.integrate.odeint requires the highest t to be equal or
                        greater than stop_time + 1""".format(highest_t, integrator_upper_boundary))

    # check the consistency of the PAR
    for organ_ in organs:
        if 'PAR' in organ_.__dict__:
            lowest_t = organ_.PAR.first_valid_index()
            if start_time < lowest_t:
                raise Exception('Error: the lowest t ({}) in the PAR of {} is greater than start_time ({}).'.format(lowest_t, organ_.name, start_time))
            highest_t = organ_.PAR.last_valid_index()
            if highest_t < integrator_upper_boundary:
                raise Exception("""Error: the highest t ({}) in the PAR of {} is lower than stop_time + 1 = {}.
                                scipy.integrate.odeint requires the highest t to be equal or
                                greater than stop_time + 1""".format(highest_t, organ_.name, integrator_upper_boundary))

    # interpolate meteo data
    meteo_interpolations = {}
    for column in meteo.columns:
        meteo_interpolations[column] = interp1d(meteo.index, meteo[column])

    # interpolate the PAR and construct the list of initial conditions
    initial_conditions = []

    initial_conditions.extend(phloem.get_initial_conditions())

    for organ_ in organs:
        try:
            organ_.PAR_linear_interpolation = interp1d(organ_.PAR.index, organ_.PAR)
        except AttributeError:
            pass
        initial_conditions.extend(organ_.get_initial_conditions())

    t = np.linspace(start_time, stop_time, number_of_output_steps)

    soln, infodict = odeint(_calculate_all_derivatives, initial_conditions, t, (phloem, organs, meteo_interpolations, photosynthesis_computation_interval), full_output=True, mxstep=odeint_mxstep)

    if not set(infodict['mused']).issubset([1,2]): # I'm not sure if this test is robust or not... Beware especially when scipy is updated.
        raise Exception("Error: Integration failed. See the logs of lsoda or try to increase the value of 'mxstep'.")

    soln_iter = iter(soln.T)

    # construct the table of the results
    result_items = [('t', t)]

    SUCROSE_phloem = soln_iter.next()
    variables = [(('Conc_Sucrose_%s' % phloem.name).rstrip('_'), phloem.calculate_Conc_Sucrose(SUCROSE_phloem)),
                 (('Conc_C_Sucrose_%s' % phloem.name).rstrip('_'), phloem.calculate_Conc_C_Sucrose(SUCROSE_phloem))]
    compartments = [(('SUCROSE_%s' % phloem.name).rstrip('_'), SUCROSE_phloem)]
    result_items.extend(variables + compartments)

    for organ_ in organs:
        if isinstance(organ_, (organ.Lamina, organ.Chaff)):
            STORAGE = soln_iter.next()
            SUCROSE = soln_iter.next()
            TRIOSESP = soln_iter.next()

            An = np.array(map(photosynthesis.PhotosynthesisModel.calculate_An,
                              t,
                              organ_.PAR_linear_interpolation(t),
                              meteo_interpolations['Tac'](t),
                              meteo_interpolations['Ca'](t),
                              meteo_interpolations['hs'](t)))

            Photosynthesis = map(organ_.calculate_Photosynthesis, t, An)

            variables = [(('Photosynthesis_%s' % organ_.name).rstrip('_'), Photosynthesis),
                         (('Conc_TriosesP_%s' % organ_.name).rstrip('_'), organ_.calculate_Conc_TriosesP(TRIOSESP)),
                         (('Conc_Storage_%s' % organ_.name).rstrip('_'), organ_.calculate_Conc_Storage(STORAGE)),
                         (('Conc_Sucrose_%s' % organ_.name).rstrip('_'), organ_.calculate_Conc_Sucrose(SUCROSE))]
            flows = [(('D_storage_%s' % organ_.name).rstrip('_'), map(organ_.calculate_D_storage, STORAGE)),
                     (('S_storage_%s' % organ_.name).rstrip('_'), map(organ_.calculate_S_storage, TRIOSESP)),
                     (('S_sucrose_%s' % organ_.name).rstrip('_'), map(organ_.calculate_S_sucrose, TRIOSESP)),
                     (('Loading_sucrose_%s' % organ_.name).rstrip('_'), map(organ_.calculate_Loading_sucrose, SUCROSE, SUCROSE_phloem))]
            compartments = [(('STORAGE_%s' % organ_.name).rstrip('_'), STORAGE),
                            (('SUCROSE_%s' % organ_.name).rstrip('_'), SUCROSE),
                            (('TRIOSESP_%s' % organ_.name).rstrip('_'), TRIOSESP)]
        elif isinstance(organ_, (organ.Internode, organ.Peduncle, organ.Sheath)):
            FRUCTAN = soln_iter.next()
            STORAGE = soln_iter.next()
            SUCROSE = soln_iter.next()
            TRIOSESP = soln_iter.next()
            Loading_sucrose = map(organ_.calculate_Loading_sucrose, SUCROSE, SUCROSE_phloem)
            Regul_Sfructan = map(organ_.calculate_Regul_Sfructan, Loading_sucrose)

            An = np.array(map(photosynthesis.PhotosynthesisModel.calculate_An,
                              t,
                              organ_.PAR_linear_interpolation(t),
                              meteo_interpolations['Tac'](t),
                              meteo_interpolations['Ca'](t),
                              meteo_interpolations['hs'](t)))

            variables = [(('Photosynthesis_%s' % organ_.name).rstrip('_'), organ_.calculate_Photosynthesis(t, An)),
                         (('Conc_TriosesP_%s' % organ_.name).rstrip('_'), organ_.calculate_Conc_TriosesP(TRIOSESP)),
                         (('Conc_Storage_%s' % organ_.name).rstrip('_'), organ_.calculate_Conc_Storage(STORAGE)),
                         (('Conc_Sucrose_%s' % organ_.name).rstrip('_'), organ_.calculate_Conc_Sucrose(SUCROSE)),
                         (('Conc_Fructan_%s' % organ_.name).rstrip('_'), organ_.calculate_Conc_Fructan(FRUCTAN)),
                         (('Regul_Sfructan_%s' % organ_.name).rstrip('_'), Regul_Sfructan)]
            flows = [(('D_storage_%s' % organ_.name).rstrip('_'), map(organ_.calculate_D_storage, STORAGE)),
                     (('S_storage_%s' % organ_.name).rstrip('_'), map(organ_.calculate_S_storage, TRIOSESP)),
                     (('S_sucrose_%s' % organ_.name).rstrip('_'), map(organ_.calculate_S_sucrose, TRIOSESP)),
                     (('Loading_sucrose_%s' % organ_.name).rstrip('_'), Loading_sucrose),
                     (('D_fructan_%s' % organ_.name).rstrip('_'), map(organ_.calculate_D_fructan, SUCROSE, FRUCTAN)),
                     (('S_fructan_%s' % organ_.name).rstrip('_'), map(organ_.calculate_S_fructan, SUCROSE, Regul_Sfructan))]
            compartments = [(('STORAGE_%s' % organ_.name).rstrip('_'), STORAGE),
                            (('SUCROSE_%s' % organ_.name).rstrip('_'), SUCROSE),
                            (('TRIOSESP_%s' % organ_.name).rstrip('_'), TRIOSESP),
                            (('FRUCTAN_%s' % organ_.name).rstrip('_'), FRUCTAN)]
        elif isinstance(organ_, organ.Grains):
            STORAGE = soln_iter.next()
            STRUCTURE = soln_iter.next()
            RGR_structure = map(organ_.calculate_RGR_structure, SUCROSE_phloem)
            variables = [(('Dry_mass_%s' % organ_.name).rstrip('_'), organ_.calculate_Dry_mass(STRUCTURE, STORAGE)),
                         (('RGR_structure_%s' % organ_.name).rstrip('_'), RGR_structure)]
            flows = [(('Unloading_sucrose_storage_%s' % organ_.name).rstrip('_'), map(organ_.calculate_Unloading_sucrose_storage, t, SUCROSE_phloem)),
                     (('Unloading_sucrose_structure_%s' % organ_.name).rstrip('_'), map(organ_.calculate_Unloading_sucrose_structure, t, STRUCTURE, RGR_structure))]
            compartments = [(('STORAGE_%s' % organ_.name).rstrip('_'), STORAGE),
                            (('STRUCTURE_%s' % organ_.name).rstrip('_'), STRUCTURE)]
        elif isinstance(organ_, organ.Roots):
            Sucrose = soln_iter.next()
            variables = [(('Dry_mass_%s' % organ_.name).rstrip('_'), organ_.calculate_Dry_mass(Sucrose))]
            flows = [(('Unloading_sucrose_%s' % organ_.name).rstrip('_'), map(organ_.calculate_Unloading_sucrose, SUCROSE_phloem))]
            compartments = [(('Sucrose_%s' % organ_.name).rstrip('_'), Sucrose)]
        result_items.extend(variables + flows + compartments)

    return pd.DataFrame.from_items(result_items)

