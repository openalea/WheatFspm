# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    elongwheat.model
    ~~~~~~~~~~~~~

    The module :mod:`elongwheat.model` defines the equations of the kinetic of leaf elongation according to CN status.

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

import parameters

# -------------------------------------------------------------------------------------------------------------------
# --- General
# -------------------------------------------------------------------------------------------------------------------

def calculate_cumsum_internode(all_internode_length, all_SAM_ids):
    """ Calculate cumulative sum of the internode lengths for each axis.

    :Parameters:
        - `all_internode_length` (:class:`dict`) - Dictionay with phytomers' keys and internode length as value for all the axis (m)
        - `all_SAM_ids` (:class:`dict`) - List of axis ids
    :Returns:
        Return Dictionay with phytomers' keys and cumulative internode length as value for all the axis (m)
    :Returns Type:
        :class:`float`
    """
    all_internode_cumsum = dict()
    for SAM_id in all_SAM_ids:
        internode_length = {k:v for k,v in all_internode_length.iteritems() if k[:2] == SAM_id }.copy()
        phytomer_nums = range(1,max([k[2] for k in internode_length])+1 )
        for phytomer_num in phytomer_nums:
            all_internode_cumsum.update( { (SAM_id[0],SAM_id[1],phytomer_num) : sum( [internode_length[d] for d in internode_length if d[2] <= phytomer_num ] ) } )

    return all_internode_cumsum

# -------------------------------------------------------------------------------------------------------------------
# --- SAM
# -------------------------------------------------------------------------------------------------------------------

def calculate_SAM_height(internode_cumsum):
    """ Calculate SAM height from lengths of internodes.

    :Parameters:
        - `internode_length` (:class:`dict`) - Dictionay with phytomers' keys and internode cumulative length as value (m)
    :Returns:
        Return updated SAM height (m)
    :Returns Type:
        :class:`float`
    """
    return max(internode_cumsum.values())

def calculate_growing_temperature(Ta, Ts, plant_height):
    """ Return temperature to be used for growth zone

    :Parameters:
        - `Ta` (:class:`float`) - air temperature at t (degree Celsius)
        - `Ts` (:class:`float`) - soil temperature at t (degree Celsius)
        - `plant_height` (:class:`float`) - Sum of all internode lengths (m).
    :Returns:
        Return temperature to be used for growth zone at t (degree Celsius)
    :Returns Type:
        :class:`float`
    """
    if plant_height > parameters.sowing_depth:
        growth_temperature = Ta
    else:
        growth_temperature = Ts
    return growth_temperature

def calculate_SAM_sumTT(T, sum_TT, nb_leaves, status, delta_t):
    """ Return sum of thermal time of the SAM

    :Parameters:
        - `T` (:class:`float`) - growing temperature at t (degree Celsius)
        - `sum_TT` (:class:`float`) - SAM cumul of temperature (degree Celsius)
        - `nb_leaves` (:class:`float`) - Number of leaves already emited by the SAM.
        - `status` (:class:`interger`) - SAM status ('retired' or 'leaf' if emitting leaf primordia).
        - `delta_t` (:class:`float`) - Time step of the simulation (s).
    :Returns:
        Return temperature to be used for growthing zone at t (degree Celsius)
    :Returns Type:
        :class:`float`
    """
    init_leaf = 0
    sum_TT += (T * delta_t) / (24 * 3600)

    if status == 'vegetative': # TODO: peut-etre a deplacer dans le convertisseur
        if sum_TT >= parameters.PLASTO_leaf*(nb_leaves+1): #first 3/4 leaves are initiated in the seed
            init_leaf = min(parameters.max_nb_leaves, int(sum_TT / parameters.PLASTO_leaf) - nb_leaves)
            if init_leaf > 1:
                raise ValueError('Error : {} leaf primordia have been created in one time step.'.format(init_leaf))
    return sum_TT, init_leaf

def calculate_SAM_status(status, nb_leaves, init_leaf):
    if (nb_leaves + init_leaf) >= parameters.max_nb_leaves:
        status = 'reproductive'
    nb_leaves += init_leaf
    return status, nb_leaves

def calculate_SAM_GA(status, sum_TT):
    is_producing_GA = False
    if status == 'reproductive' and sum_TT > (parameters.max_nb_leaves * parameters.PLASTO_leaf) + parameters.delta_TT_GA :
       is_producing_GA = True
    return is_producing_GA


# -------------------------------------------------------------------------------------------------------------------
# --- Leaves
# -------------------------------------------------------------------------------------------------------------------

def calculate_ligule_height(ligulated_sheath_length, all_internode_cumsum):
    """ Calculate ligules height from lengths of internodes and sheaths.

    :Parameters:
        - `ligulated_sheath_length` (:class:`dict`) - Dictionay with phytomers' keys and sheath length as value (m)
        - `internode_cumsum` (:class:`dict`) - Dictionay with phytomers' keys and internode cumulative length as value (m)
    :Returns:
        Return height of the all the ligules (m)
    :Returns Type:
        :class:`float`
    """
    ligule_height = { k : v+all_internode_cumsum[k] for k,v in ligulated_sheath_length.iteritems()}

    return ligule_height

def calculate_leaf_dist_to_emerge(all_ligule_height, bottom_hiddenzone_height, SAM_id, hiddenzone_id):
    """ length of the pseudostem given by the higest ligule

    :Parameters:
        - `all_ligule_height` (:class:`float`) - Elevation of the all the ligules (m)
        - `bottom_hiddenzone_height` (:class:`float`) - Sum of lengths of the internodes from metamer 1 to current metamer N (current metamer included) (m).
        - `SAM_id` (:class:`tuple`) - Axis ID.
        - `hiddenzone_id` (:class:`tuple`) - Hiddenzone ID.
    :Returns:
        Distance for the leaf to emerge out of the pseudostem (m)
    :Returns Type:
        :class:`float`
    """
    ligule_height_below = [v for k,v in all_ligule_height.iteritems() if k[:2] == SAM_id and k[2] < hiddenzone_id[2]]
    top_ligule_height = max(ligule_height_below) if len(ligule_height_below) > 0 else 0
    leaf_dist_to_emerge = top_ligule_height - bottom_hiddenzone_height
    return max(leaf_dist_to_emerge,0)

def calculate_deltaL_preE(leaf_rank, sucrose, leaf_L, amino_acids, mstruct, delta_t):
    """ delta of leaf length over delta_t as a function of sucrose and amino acids, from initiation to the emergence of the previous leaf.

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

    RER_max = parameters.RERmax #parameters.RERmax_dict[leaf_rank]
    if sucrose > 0:
        delta_leaf_L = leaf_L * RER_max * delta_t * ((sucrose / mstruct) / (parameters.Kc + (sucrose / mstruct))) * (((amino_acids/mstruct) **3) / (parameters.Kn**3 + (amino_acids / mstruct)**3))
    else:
        delta_leaf_L = 0

    return delta_leaf_L

def calculate_leaf_psdo_age(leaf_psdo_age, delta_t):
    """ Pseudo age of the leaf since beginning of automate elongation (s)
    :Parameters:
        - `leaf_psdo_age` (:class:`float`) - Pseudo age of the leaf since beginning of automate elongation (s)
        - `delta_t` (:class:`float`) - Time step (s)
    :Returns:
        leaf_psdo_age (s)
    :Returns Type:
        :class:`float`
    """
    return leaf_psdo_age + delta_t

def calculate_L_postE(leaf_psdo_age, leaf_Lmax):
    """ leaf length, from the emergence of the previous leaf to the end of elong (predefined elong kinetic depending on leaf state).

    :Parameters:
        - `leaf_psdo_age` (:class:`float`) - Pseudo age of the leaf since beginning of automate elongation (s)
        - `leaf_Lmax` (:class:`float`) - Final leaf length (m)
    :Returns:
        leaf_L (m)
    :Returns Type:
        :class:`float`
    """

    E = abs( (1+(max(0,parameters.te)/(parameters.te - parameters.tm )))*(min(1,(-parameters.tb)/(parameters.te - parameters.tb))**((parameters.te-parameters.tb)/(parameters.te-parameters.tm)))) # E = beta function at t=0

    if leaf_psdo_age <= parameters.tb:
        leaf_L = (1/parameters.Y0) * leaf_Lmax
    elif leaf_psdo_age < parameters.te :
        leaf_L = min( leaf_Lmax , leaf_Lmax *( abs((1+(max(0,(parameters.te - leaf_psdo_age))/(parameters.te-parameters.tm)))*(min(1,(leaf_psdo_age-parameters.tb)/(parameters.te-parameters.tb))**((parameters.te-parameters.tb)/(parameters.te-parameters.tm)))) + (1/parameters.Y0) - E ) )
    else:
        leaf_L = leaf_Lmax
    return leaf_L

def calculate_leaf_Lmax(leaf_Lem_prev):
    """ Final leaf length.

    :Parameters:
        - `leaf_Lem_prev` (:class:`float`) - Leaf length at the emergence of the previous leaf (m)
    :Returns:
        Final leaf length (m)
    :Returns Type:
        :class:`float`
    """
    return leaf_Lem_prev * parameters.Y0

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

def calculate_leaf_emergence(leaf_L, leaf_dist_to_emerge):
    """Calculate if a given leaf has emerged from the hidden zone

    :Parameters:
        - `leaf_L` (:class:`float`) - Total leaf length (m)
        - `leaf_dist_to_emerge` (:class:`float`) - Length of the pseudostem (m)
    :Returns:
        Specifies if the leaf has emerged (True) or not (False)
    :Returns Type:
        :class:`bool`
    """
    return leaf_L > leaf_dist_to_emerge

def calculate_lamina_L(leaf_L, leaf_dist_to_emerge, hiddenzone_id):
    """ Emerged lamina length given by the difference between leaf length and hidden zone length.

    :Parameters:
        - `leaf_L` (:class:`float`) - Total leaf length (m)
        - `leaf_dist_to_emerge` (:class:`float`) - Length of the pseudostem (m)
    :Returns:
        lamina length (m)
    :Returns Type:
        :class:`float`
    """
    lamina_L = leaf_L - leaf_dist_to_emerge
    if lamina_L <=0:
        raise Warning('The leaf is shorther than its pseudostem for {}'.format(hiddenzone_id) )
    return max(0, lamina_L)

def calculate_sheath_L(leaf_L, leaf_dist_to_emerge, lamina_L):
    """ Emerged sheath length. Assumes that leaf_L = leaf_dist_to_emerge + sheath_L + lamina_L

    :Parameters:
        - `leaf_L` (:class:`float`) - Total leaf length (m)
        - `leaf_dist_to_emerge` (:class:`float`) - Length of the pseudostem (m)
        - `lamina_L` (:class:`float`) - Lamina length (m)
    :Returns:
        sheath length (m)
    :Returns Type:
        :class:`float`
    """
    return leaf_L - leaf_dist_to_emerge - lamina_L


# -------------------------------------------------------------------------------------------------------------------
# --- Internodes
# -------------------------------------------------------------------------------------------------------------------

def calculate_internode_dist_to_emerge(leaf_dist_to_emerge, curr_internode_L):
    """ distance for the internnode to be visible given by the previous sheaths.

    :Parameters:
        - `leaf_dist_to_emerge` (:class:`float`) - Pseudostem length for the leaf on the same metamer (m)
        - `curr_internode_L` (:class:`float`) - Internode length of the current metamer (m).
     :Returns:
        Distance for the internode to be visible (m)
    :Returns Type:
        :class:`float`
    """
    internode_dist_to_emerge = leaf_dist_to_emerge - curr_internode_L
    return max(internode_dist_to_emerge,0)

def calculate_internode_Lmax(internode_L_lig_curr):
    """ Final internode length.

    :Parameters:
        - `internode_L_lig_curr` (:class:`float`) - Internode length at the ligulation of the leaf on the same phytomer (m)
    :Returns:
        Final internode length (m)
    :Returns Type:
        :class:`float`
    """

    internode_Lmax = internode_L_lig_curr * parameters.Y0_IN

    return internode_Lmax

def calculate_short_internode_Lmax(internode_L_GA_prod, internode_psdo_age):
    """ Final internode length.

    :Parameters:
        - `internode_L_GA_prod` (:class:`float`) - Internode length at the ligulation of the leaf on the same phytomer (m)
        - `internode_psdo_age` (:class:`float`) - Internode pseudo age since previous leaf ligulation (s)
    :Returns:
        Final internode length (m)
    :Returns Type:
        :class:`float`
    """

    Y0_short_IN = 1 / calculate_internode_L_postL(internode_psdo_age, 1)
    internode_Lmax = internode_L_GA_prod * Y0_short_IN

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

def calculate_internode_visibility(internode_L, internode_dist_to_emerge):
    """Calculate if a given internode is visible

    :Parameters:
        - `internode_L` (:class:`float`) - Total internode length (m)
        - `internode_dist_to_emerge` (:class:`float`) - Length of the pseudostem + internode length (m)
    :Returns:
        Specifies if the internode is visible (True) or not (False)
    :Returns Type:
        :class:`bool`
    """
    return internode_L > internode_dist_to_emerge

def calculate_delta_internode_L_preL(internode_rank, sucrose, internode_L, amino_acids, mstruct, delta_t):
    """ delta of internode length over delta_t as a function of sucrose and amino acids, from initiation to the ligulation of the leaf on the same phytomer.

    :Parameters:
        - 'internode_rank' - Phytomer number.
        - `sucrose` (:class:`float`) - Amount of sucrose (µmol C)
        - `internode_L` (:class:`float`) - Total internode length (m)
        - `amino_acids` (:class:`float`) - Amount of amino acids (µmol N)
        - `mstruct` (:class:`float`) - Structural mass (g)
    :Returns:
        delta delta_internode_L (m)
    :Returns Type:
        :class:`float`
    """
    RER_max = parameters.RERmax
    ##RER_max = parameters.RERmax_dict[internode_rank]
    if sucrose > 0:
        delta_internode_L = internode_L  * RER_max * delta_t * ((sucrose / mstruct) / (parameters.Kc + (sucrose / mstruct))) * (((amino_acids/mstruct) **3) / (parameters.Kn**3 + (amino_acids / mstruct)**3))
    else:
        delta_internode_L = 0
    return delta_internode_L

def calculate_internode_psdo_age(internode_psdo_age, delta_t):
    """ Pseudo age of the internode since beginning of automate elongation (s)
    :Parameters:
        - `internode_psdo_age` (:class:`float`) - Pseudo age of the leaf since beginning of automate elongation (s)
        - `delta_t` (:class:`float`) - Time step (s)
    :Returns:
        internode_psdo_age (s)
    :Returns Type:
        :class:`float`
    """
    return internode_psdo_age + delta_t

def calculate_internode_L_postL(internode_psdo_age, internode_Lmax):
    """ internode length, from the ligulation of the leaf on the same phytomer to the end of elong (predefined elong kinetic depending on internode state).

    :Parameters:
        - `internode_psdo_age` (:class:`float`) - Pseudo age of the internode since beginning of internode automate growth (s)
        - `internode_Lmax` (:class:`float`) - Final internode length (m)
    :Returns:
        delta internode_L (m)
    :Returns Type:
        :class:`float`
    """

    E = (1+(max(0,parameters.te_IN)/(parameters.te_IN - parameters.tm_IN )))*(min(1,(-parameters.tb_IN)/(parameters.te_IN - parameters.tb_IN))**((parameters.te_IN-parameters.tb_IN)/(parameters.te_IN-parameters.tm_IN))) # E = beta function at t=0

    if internode_psdo_age <= parameters.tb_IN:
        internode_L = (1/parameters.Y0_IN) * internode_Lmax
    elif internode_psdo_age < parameters.te_IN :
        internode_L = min( internode_Lmax , internode_Lmax *( abs((1+(max(0,(parameters.te_IN - internode_psdo_age))/(parameters.te_IN-parameters.tm_IN)))*(min(1,(internode_psdo_age-parameters.tb_IN)/(parameters.te_IN-parameters.tb_IN))**((parameters.te_IN-parameters.tb_IN)/(parameters.te_IN-parameters.tm_IN)))) + (1/parameters.Y0_IN) - E ) )
    else:
        internode_L = internode_Lmax
    return internode_L

def calculate_init_internode_elongation(sum_TT, metamer_rank):
    """Initialize internode elongation.

    :Parameters:
        - `sum_TT` (:class:`float`) - soil temperature at t (degree Celsius)
        - `metamer_rank` (:class:`float`) - Rank of the phytomer.
    :Returns:
        Specifies if the internode has started the elongation (True) or not (False), and initialize internode_L
    :Returns Type:
        :class:`bool`
    """
    is_growing = False
    internode_L = 0
    if (sum_TT - metamer_rank*parameters.PLASTO_leaf) > parameters.nb_PLASTO_internode_init*parameters.PLASTO_leaf:
       is_growing = True
       internode_L = parameters.internode_L_init
    return is_growing, internode_L

def calculate_end_internode_elongation(internode_L, internode_Lmax):
    """Calculate if a given internode has finished elongating

    :Parameters:
        - `internode_L` (:class:`float`) - Total internode length (m)
        - `internode_Lmax` (:class:`float`) - Maximum internode length (m)
    :Returns:
        Specifies if the internode has completed elongation (True) or not (False)
    :Returns Type:
        :class:`bool`
    """
    return  (internode_L >= internode_Lmax)

def calculate_end_internode_elongation_age(internode_psdo_age):
    """Calculate if a given internode has finished elongating

    :Parameters:
        - `internode_psdo_age` (:class:`float`) - Internode psdo age (s)
    :Returns:
        Specifies if the internode has completed elongation (True) or not (False)
    :Returns Type:
        :class:`bool`
    """
    return (internode_psdo_age >= parameters.tm_IN)