# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    elongwheat.model
    ~~~~~~~~~~~~~

    The module :mod:`elongwheat.model` defines the equations of the kinetic of leaf and internode elongation according to C and N status.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2015.
"""

"""
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
"""
import pandas as pd

import parameters

# -------------------------------------------------------------------------------------------------------------------
# --- SAM
# -------------------------------------------------------------------------------------------------------------------

def calculate_growing_temperature(Tair, Tsol, SAM_height):
    """ Return temperature to be used for growth zone

    :Parameters:
        - `Tair` (:class:`float`) - air temperature at t (degree Celsius)
        - `Tsol` (:class:`float`) - soil temperature at t (degree Celsius)
        - `SAM_height` (:class:`float`) - Height of SAM, calculated from internode length (m).
    :Returns:
        Return temperature to be used for growth zone at t (degree Celsius)
    :Returns Type:
        :class:`float`
    """

    if SAM_height > parameters.sowing_depth:
        growth_temperature = Tair
    else:
        growth_temperature = Tsol

    return growth_temperature

def calculate_SAM_sumTT(growth_temperature, sum_TT, delta_t):
    """ Increment thermal time of the SAM

    :Parameters:
        - `growth_temperature` (:class:`float`) - temperature of the SAM at t (degree Celsius)
        - `sum_TT` (:class:`float`) - SAM cumul of temperature (degree day-1)
        - `delta_t` (:class:`float`) - Time step of the simulation (s)
    :Returns:
        Updated thermal time of the SAM at t (degree day-1)
    :Returns Type:
        :class:`float`
    """

    sum_TT += (growth_temperature * delta_t) / (24 * 3600)
    return sum_TT

def calculate_SAM_status(status, sum_TT, nb_leaves):
    """ Update SAM status, leaf number

    :Parameters:
        - `status` (:class:`string`) - SAM status ('vegetative', if emitting leaf primordia or 'reproductive')
        - `sum_TT` (:class:`float`) - SAM cumul of temperature (degree day-1)
        - `nb_leaves` (:class:`float`) - Number of leaves already emited by the SAM.
    :Returns:
        Number of leaf to be initiated (should be 0 or 1), updated leaf number on the SAM, status
    :Returns Type:
        :class:`float`, :class:`int`, :class:`string`
    """

    if status == 'vegetative':
        if sum_TT >= parameters.PLASTOCHRONE*(nb_leaves+1): # first 3/4 leaves are initiated in the seed
            init_leaf = min(parameters.max_nb_leaves, int(sum_TT / parameters.PLASTOCHRONE) - nb_leaves)
            if init_leaf > 1:
                raise ValueError('Error : {} leaf primordia have been created in one time step. Consider shorter timestep'.format(init_leaf))
        else:
            init_leaf = 0

        if (nb_leaves + init_leaf) >= parameters.max_nb_leaves:
            status = 'reproductive'
    else: # Reproductive SAM
        init_leaf = 0

    nb_leaves += init_leaf

    return init_leaf, nb_leaves, status

def calculate_SAM_GA(status, sum_TT):
    """ Synthesis of GA by the SAM according to its stage

    :Parameters:
        - `status` (:class:`string`) - SAM status ('vegetative', if emitting leaf primordia or 'reproductive')
        - `sum_TT` (:class:`float`) - SAM cumul of temperature (degree day-1)
    :Returns:
        GA production
    :Returns Type:
        :class:`bool`
    """
    is_producing_GA = status == 'reproductive' and sum_TT > (parameters.max_nb_leaves * parameters.PLASTOCHRONE) + parameters.delta_TT_GA
    return is_producing_GA

# -------------------------------------------------------------------------------------------------------------------
# --- Leaves
# -------------------------------------------------------------------------------------------------------------------

def calculate_ligule_height(sheath_internode_length, all_element_inputs, SAM_id, all_ligule_height_df):
    """ Calculate ligule heights below each phytomer of an axis from lengths of internodes and sheaths.
        The result is formated into a dataframe which is concatenated with a shared dataframe storing all ligule heights

    :Parameters:
        - `sheath_internode_length` (:class:`dict`) - Dictionary with sheath and internode lengths for a given axis, (m)
        - `all_element_inputs` (:class:`dict`) - Dictionary of all element inputs
        - `SAM_id` (:class:`tuple`) - Id of the SAM (plant_id, axis_id)
        - `all_ligule_height_df` (:cpandas.DataFrame) - Shared dataframe used to store all ligule heights
    :Returns:
        A dataframe with ligule height of a given axis
    :Returns Type:
        :class:`pandas.DataFrame`
    """

    for phytomer_id, lengths in sheath_internode_length.iteritems():
        lamina_id = SAM_id + (phytomer_id, 'blade', 'LeafElement1')
        if all_element_inputs.has_key(lamina_id) and not all_element_inputs[lamina_id]['is_growing']:
            ligule_height = sum(lengths['sheath'] + lengths['cumulated_internode'])
            ligule_height_df = pd.DataFrame([(SAM_id, phytomer_id, ligule_height)], columns=list(all_ligule_height_df))
            all_ligule_height_df = pd.concat((all_ligule_height_df, ligule_height_df))

    return all_ligule_height_df

def calculate_leaf_pseudostem_length(ligule_heights, bottom_hiddenzone_height, phytomer_id):
    """ Distance between the bottom of the hidden zone and the highest previous ligule

    :Parameters:
        - `ligule_heights` (:class:`dict`) - height of all ligules of a given axis (m)
        - `bottom_hiddenzone_height` (:class:`float`) - height of the bottom of the hidden zone given by the cumulated lengths of below internodes (m)
        - `phytomer_id` (:class:`int`) - Id of the SAM (plant_id, axis_id)
    :Returns:
        Distance for the leaf to emerge out of the pseudostem (m)
    :Returns Type:
        :class:`float`
    """
    top_ligule_height = max(ligule_heights[ligule_heights['phytomer'] < phytomer_id]['ligule height']) # highest previous ligule
    leaf_pseudostem_length = top_ligule_height - bottom_hiddenzone_height

    return leaf_pseudostem_length

def calculate_deltaL_preE(sucrose, leaf_L, amino_acids, mstruct, delta_t, leaf_rank):
    """ Delta of leaf length over delta_t as a function of sucrose and amino acids, from initiation to the emergence of the previous leaf.

    :Parameters:
        - `sucrose` (:class:`float`) - Amount of sucrose (µmol C)
        - `leaf_L` (:class:`float`) - Total leaf length (m)
        - `amino_acids` (:class:`float`) - Amount of amino acids (µmol N)
        - `mstruct` (:class:`float`) - Structural mass (g)
    :Returns:
        delta delta_leaf_L (m)
    :Returns Type:
        :class:`float`
    """

##    RER_max = parameters.RERmax
    RER_max = parameters.RERmax_dict[leaf_rank]
    if sucrose > 0:
        delta_leaf_L = leaf_L * RER_max * delta_t# * ((sucrose / mstruct) / (parameters.Kc + (sucrose / mstruct))) * (((amino_acids/mstruct) **3) / (parameters.Kn**3 + (amino_acids / mstruct)**3))
    else:
        delta_leaf_L = 0

    return delta_leaf_L

def calculate_leaf_pseudo_age(leaf_pseudo_age, delta_t):
    """ Pseudo age of the leaf since beginning of automate elongation (s)
    :Parameters:
        - `leaf_pseudo_age` (:class:`float`) - Previous pseudo age of the leaf since beginning of automate elongation (s)
        - `delta_t` (:class:`float`) - Time step (s)
    :Returns:
        Updated leaf_pseudo_age (s)
    :Returns Type:
        :class:`float`
    """

    return leaf_pseudo_age + delta_t

def calculate_L_postE(leaf_pseudo_age, leaf_Lmax):
    """ Leaf length from the emergence of the previous leaf to the end of elongation (automate function depending on leaf pseudo age and final length).

    :Parameters:
        - `leaf_pseudo_age` (:class:`float`) - Pseudo age of the leaf since beginning of automate elongation (s)
        - `leaf_Lmax` (:class:`float`) - Final leaf length (m)
    :Returns:
        leaf_L (m)
    :Returns Type:
        :class:`float`
    """

    if leaf_pseudo_age <= parameters.tb:
        leaf_L = parameters.FITTED_L0 * leaf_Lmax
    elif leaf_pseudo_age < parameters.te :
        leaf_L = min(leaf_Lmax, leaf_Lmax * (abs((1+ (max(0, (parameters.te - leaf_pseudo_age)) / (parameters.te - parameters.tm))) * (min(1, float(leaf_pseudo_age - parameters.tb) / float(parameters.te - parameters.tb))**((parameters.te - parameters.tb)/(parameters.te - parameters.tm)))) + parameters.OFFSET_LEAF))
    else:
        leaf_L = leaf_Lmax

    return leaf_L

def calculate_leaf_emergence(leaf_L, leaf_pseudostem_length):
    """Calculate if a given leaf has emerged from its pseudostem

    :Parameters:
        - `leaf_L` (:class:`float`) - Total leaf length (m)
        - `leaf_pseudostem_length` (:class:`float`) - Length of the pseudostem (m)
    :Returns:
        Specifies if the leaf has emerged (True) or not (False)
    :Returns Type:
        :class:`bool`
    """
    return leaf_L > leaf_pseudostem_length

def calculate_lamina_L(leaf_L, leaf_pseudostem_length, hiddenzone_id):
    """ Emerged lamina length given by the difference between leaf length and leaf pseudostem length.

    :Parameters:
        - `leaf_L` (:class:`float`) - Total leaf length (m)
        - `leaf_pseudostem_length` (:class:`float`) - Length of the pseudostem (m)
        - `hiddenzone_id` (:class:`tuple`) - Id of the hidden zone (plant_id, axis_id, phytomer_id)
    :Returns:
        lamina length (m)
    :Returns Type:
        :class:`float`
    """
    lamina_L = leaf_L - leaf_pseudostem_length
    if lamina_L <=0:
        raise Warning('The leaf is shorther than its pseudostem for {}'.format(hiddenzone_id))
    return max(0, lamina_L)

def calculate_leaf_Lmax(leaf_Lem_prev):
    """ Final leaf length.

    :Parameters:
        - `leaf_Lem_prev` (:class:`float`) - Leaf length at the emergence of the previous leaf (m)
    :Returns:
        Final leaf length (m)
    :Returns Type:
        :class:`float`
    """
    return leaf_Lem_prev * parameters.SCALING_FACTOR_LEAF

def calculate_SL_ratio(phytomer_rank):
    """ Sheath:Lamina final length ratio according to the rank. Parameters from Dornbush (2011).

    :Parameters:
        - `phytomer_rank` (:class:`float`)
    :Returns:
        Sheath:Lamina ratio (dimensionless)
    :Returns Type:
        :class:`float`
    """
    return -0.0021 * phytomer_rank**3 + 0.037 * phytomer_rank**2 - 0.1527 * phytomer_rank + 0.4962

def calculate_lamina_Lmax(leaf_Lmax, sheath_lamina_ratio):
    """ Final lamina length.

    :Parameters:
        - `leaf_Lmax` (:class:`float`) - Final leaf length (m)
        - `sheath_lamina_ratio` (:class:`float`) - Sheath:Lamina ratio (dimensionless)

    :Returns:
        final lamina length (m)
    :Returns Type:
        :class:`float`
    """
    return leaf_Lmax / (1 + sheath_lamina_ratio)

def calculate_sheath_Lmax(leaf_Lmax, lamina_Lmax):
    """ Final sheath length.

    :Parameters:
        - `leaf_Lmax` (:class:`float`) - Final leaf length (m)
        - `lamina_Lmax` (:class:`float`) - Final lamina length (m)

    :Returns:
        final sheath length (m)
    :Returns Type:
        :class:`float`
    """
    return leaf_Lmax - lamina_Lmax

def calculate_leaf_Wmax(lamina_Lmax, fructan, mstruct):
    """ Maximal leaf width.
    0.0575 et 0.12 issu graph Dornbush

    :Parameters:
        - `lamina_Lmax` (:class:`float`) - Maximal lamina length (m)
        - `fructan` (:class:`float`) - Fructan in the hidden zone at the time of the previous leaf emergence (µmol C).
        - `mstruct` (:class:`float`) - Mstruct of the hidden zone at the time of the previous leaf emergence (g).
    :Returns:
        maximal leaf width (m)
    :Returns Type:
        :class:`float`
    """
    return (0.0575 * lamina_Lmax - 0.00012) * (parameters.EC_wmax * 2 * parameters.Ksslw/(parameters.Ksslw + (fructan / mstruct)) + (1-parameters.EC_wmax)) #TODO: a remplacer

def calculate_SSLW(fructan, mstruct):
    """ Structural Specific Lamina Weight.

    :Parameters:
        - `fructan` (:class:`float`) - Fructan in the hidden zone at the time of the previous leaf emergence (µmol C).
        - `mstruct` (:class:`float`) - Mstruct of the hidden zone at the time of the previous leaf emergence (g).
    :Returns:
        Structural Specific Leaf Weight (g m-2)
    :Returns Type:
        :class:`float`
    """
    conc_fructan = fructan / mstruct

    return parameters.min_SSLW + (parameters.max_SSLW - parameters.min_SSLW) * conc_fructan/ (conc_fructan + parameters.Ksslw)

def calculate_SSSW(SSLW):
    """ Structural Specific Sheath Weight.

    :Parameters:
        - `SSLW` (:class:`float`) - Structural Specific Leaf Weight (g m-2).
    :Returns:
        Structural Specific Sheath Weight (g m-2)
    :Returns Type:
        :class:`float`
    """
    return SSLW * parameters.ratio_SSSW_SSLW

def calculate_emerged_sheath_L(leaf_L, leaf_pseudostem_length, lamina_L):
    """ Emerged sheath length. Assumes that leaf_L = leaf_pseudostem_length + emerged_sheath_L + lamina_L

    :Parameters:
        - `leaf_L` (:class:`float`) - Total leaf length (m)
        - `leaf_pseudostem_length` (:class:`float`) - Length of the pseudostem (m)
        - `lamina_L` (:class:`float`) - Lamina length (m)
    :Returns:
        sheath length (m)
    :Returns Type:
        :class:`float`
    """
    return leaf_L - leaf_pseudostem_length - lamina_L


# -------------------------------------------------------------------------------------------------------------------
# --- Internodes
# -------------------------------------------------------------------------------------------------------------------

def calculate_cumulated_internode_length(internode_lengths):
    """ Calculate cumulative internode length

    :Parameters:
        - `internode_lengths` (:class:`list`) - list of internode lengths below a given phytomer index (m)
    :Returns:
        Return the cumulated internode length below a given phytomer index (m)
    :Returns Type:
        :class:`float`
    """

    return sum(internode_lengths)

def calculate_internode_distance_to_emergence(leaf_pseudostem_length, curr_internode_L):
    """ Distance for the internode to be visible given by the pseudostem length.

    :Parameters:
        - `leaf_pseudostem_length` (:class:`float`) - Pseudostem length of the leaf on the same phyomer (m)
        - `curr_internode_L` (:class:`float`) - Internode length of the current phyomer (m).
     :Returns:
        Distance for the internode to be visible (m)
    :Returns Type:
        :class:`float`
    """
    internode_distance_to_emergence = leaf_pseudostem_length - curr_internode_L

    return max(0, internode_distance_to_emergence)

def calculate_internode_Lmax(internode_L_lig):
    """ Final internode length.

    :Parameters:
        - `internode_L_lig` (:class:`float`) - Internode length at the ligulation of the previous leaf (m)
    :Returns:
        Final internode length (m)
    :Returns Type:
        :class:`float`
    """

    internode_Lmax = internode_L_lig * parameters.SCALING_FACTOR_INT

    return internode_Lmax

def calculate_SSIW(SSSW):
    """ Structural Specific Internode Weight.

    :Parameters:
        - `SSSW` (:class:`float`) - Structural Specific Sheath Weight (g m-2).
    :Returns:
        Structural Specific Internode Weight (g m-2)
    :Returns Type:
        :class:`float`
    """
    return SSSW

def calculate_init_internode_elongation(sum_TT, phytomer_rank):
    """Initialize internode elongation.

    :Parameters:
        - `sum_TT` (:class:`float`) - SAM cumul of temperature (degree day-1)
        - `phytomer_rank` (:class:`int`) - Rank of the phytomer.
    :Returns:
        Specifies if the internode has started the elongation (True) or not (False), and initialize internode_L
    :Returns Type:
        :class:`bool`, :class:`float`
    """

    phytomer_thermal_age = sum_TT - phytomer_rank*parameters.PLASTOCHRONE

    if phytomer_thermal_age > parameters.nb_PLASTO_internode_init * parameters.PLASTOCHRONE:
       is_growing = True
       internode_L = parameters.internode_L_init
    else:
        is_growing = False
        internode_L = 0

    return is_growing, internode_L

def calculate_delta_internode_L_preL(internode_rank, sucrose, internode_L, amino_acids, mstruct, delta_t):
    """ delta of internode length over delta_t as a function of sucrose and amino acids, from initiation to the ligulation of the previous leaf.

    :Parameters:
        - 'internode_rank' - Phytomer number.
        - `sucrose` (:class:`float`) - Amount of sucrose (µmol C)
        - `internode_L` (:class:`float`) - Total internode length (m)
        - `amino_acids` (:class:`float`) - Amount of amino acids (µmol N)
        - `mstruct` (:class:`float`) - Structural mass of the hidden zone(g)
    :Returns:
        delta delta_internode_L (m)
    :Returns Type:
        :class:`float`
    """
    #RER_max = parameters.RERmax
    RER_max = parameters.RERmax_dict[internode_rank]

    if sucrose > 0:
        delta_internode_L = internode_L  * RER_max * delta_t# * ((sucrose / mstruct) / (parameters.Kc + (sucrose / mstruct))) * (((amino_acids/mstruct) **3) / (parameters.Kn**3 + (amino_acids / mstruct)**3))
    else:
        delta_internode_L = 0

    return delta_internode_L

def calculate_internode_pseudo_age(internode_pseudo_age, delta_t):
    """ Pseudo age of the internode since beginning of automate elongation (s)
    :Parameters:
        - `internode_pseudo_age` (:class:`float`) - Pseudo age of the leaf since beginning of automate elongation (s)
        - `delta_t` (:class:`float`) - Time step (s)
    :Returns:
        internode_pseudo_age (s)
    :Returns Type:
        :class:`float`
    """
    return internode_pseudo_age + delta_t

def calculate_short_internode_Lmax(internode_L_lig, internode_pseudo_age):
    """ Final internode length.

    :Parameters:
        - `internode_L_lig` (:class:`float`) - Internode length at the ligulation of the previous leaf (m)
        - `internode_pseudo_age` (:class:`float`) - Internode pseudo age since previous leaf ligulation (s)
    :Returns:
        Final internode length (m)
    :Returns Type:
        :class:`float`
    """

    L0 = 1 / calculate_internode_L_postL(internode_pseudo_age, 1) #: Initial relative length of the short internode according to its pseudo age
    internode_Lmax = internode_L_lig * L0

    return internode_Lmax

def calculate_internode_L_postL(internode_pseudo_age, internode_Lmax):
    """ internode length, from the ligulation of the previous leaf to the end of elongation (predefined kinetic depending on internode pseudo age).

    :Parameters:
        - `internode_pseudo_age` (:class:`float`) - Pseudo age of the internode since beginning of internode automate growth (s)
        - `internode_Lmax` (:class:`float`) - Final internode length (m)
    :Returns:
        internode_L (m)
    :Returns Type:
        :class:`float`
    """

    if internode_pseudo_age <= parameters.tb_IN:
        internode_L = (1 / parameters.SCALING_FACTOR_INT) * internode_Lmax
    elif internode_pseudo_age < parameters.te_IN :
        internode_L = min(internode_Lmax, internode_Lmax * (abs((1 + (max(0, (parameters.te_IN - internode_pseudo_age)) / (parameters.te_IN - parameters.tm_IN))) * (min(1, float(internode_pseudo_age - parameters.tb_IN) / float(parameters.te_IN - parameters.tb_IN))**((parameters.te_IN - parameters.tb_IN) / (parameters.te_IN - parameters.tm_IN)))) + parameters.OFFSET_INT))
    else:
        internode_L = internode_Lmax
    return internode_L

def calculate_internode_visibility(internode_L, internode_distance_to_emergence):
    """Calculate if a given internode is visible

    :Parameters:
        - `internode_L` (:class:`float`) - Total internode length (m)
        - `internode_distance_to_emergence` (:class:`float`) - Length of the pseudostem + internode length (m)
    :Returns:
        Specifies if the internode is visible (True) or not (False)
    :Returns Type:
        :class:`bool`
    """
    return internode_L > internode_distance_to_emergence

def calculate_emerged_internode_L(internode_L, internode_distance_to_emergence):
    """ Emerged internode length.

    :Parameters:
        - `internode_L` (:class:`float`) - Total internode length (m)
        - `internode_distance_to_emergence` (:class:`float`) - Length of the pseudostem + internode length (m)
    :Returns:
        Emerged internode length (m)
    :Returns Type:
        :class:`float`
    """
    return internode_L - internode_distance_to_emergence

def calculate_end_internode_elongation(internode_L, internode_Lmax, internode_pseudo_age):
    """Calculate if a given internode has finished elongating

    :Parameters:
        - `internode_L` (:class:`float`) - Total internode length (m)
        - `internode_Lmax` (:class:`float`) - Maximum internode length (m)
        - `internode_pseudo_age` (:class:`float`) - Internode pseudo age (s)
    :Returns:
        Specifies if the internode has completed elongation (True) or not (False)
    :Returns Type:
        :class:`bool`
    """

    return  (internode_L >= internode_Lmax) or (internode_pseudo_age >= parameters.te_IN)