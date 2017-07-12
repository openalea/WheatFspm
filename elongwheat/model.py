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
# --- SAM
# -------------------------------------------------------------------------------------------------------------------

def calculate_growing_temp(Ta, Ts, plant_height):
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
        growth_temp = Ta
    else:
        growth_temp = Ts
    return growth_temp

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
        if sum_TT >= parameters.PLASTO_leaf*(nb_leaves+1):
            init_leaf = min(parameters.max_nb_leaves, int(sum_TT / parameters.PLASTO_leaf) - nb_leaves)
            if init_leaf > 1:
                raise ValueError('Error : {} leaf primordia have been created in one time step.'.format(init_leaf))
    return sum_TT, init_leaf

def calculate_SAM_status(status, nb_leaves, init_leaf):
    if (nb_leaves + init_leaf) >= parameters.max_nb_leaves:
        status = 'reproductive'
    nb_leaves += init_leaf
    return status, nb_leaves

def calculate_SAM_GA(status, nb_leaves, sum_TT):
    is_producing_GA = False
    if status == 'reproductive' and sum_TT > (nb_leaves * parameters.PLASTO_leaf) + parameters.delta_TT_GA :
       is_producing_GA = True
    return is_producing_GA


# -------------------------------------------------------------------------------------------------------------------
# --- Leaves
# -------------------------------------------------------------------------------------------------------------------

def calculate_leaf_dist_to_emerge(previous_leaf_dist_to_emerge, internode_L, previous_sheath_visible_L, previous_sheath_final_hidden_L):
    """ length of the pseudostem given by the previous sheaths and internode length

    :Parameters:
        - `previous_leaf_dist_to_emerge` (:class:`float`) - Distance for the previous leaf to emerge out of the pseudostem (m). Could be None is no previous hiddenzone found.
        - `internode_L` (:class:`float`) - Length of the internode (hidden part) (m).
        - `previous_sheath_visible_L` (:class:`float`) - Visible length of the previous sheath (m).
        - `previous_sheath_final_hidden_L` (:class:`float`) - Final hidden length of the previous sheath (m).
    :Returns:
        Distance for the leaf to emerge out of the pseudostem (m)
    :Returns Type:
        :class:`float`
    """
    if previous_leaf_dist_to_emerge:
        leaf_dist_to_emerge = previous_leaf_dist_to_emerge + previous_sheath_visible_L - internode_L
    else:
        leaf_dist_to_emerge = previous_sheath_final_hidden_L + previous_sheath_visible_L - internode_L # here 'previous_sheath_visible_L' is also the final visible length of the previous sheath because the previous leaf has fully elongated.
    return leaf_dist_to_emerge


def calculate_deltaL_preE(sucrose, leaf_L, amino_acids, mstruct, delta_t):
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
    if sucrose > 0:
        delta_leaf_L = leaf_L * ((sucrose / mstruct) / (parameters.Kc + (sucrose / mstruct))) * (((amino_acids/mstruct) **3) / (parameters.Kn**3 + (amino_acids / mstruct)**3)) * parameters.RERmax * delta_t
    else:
        delta_leaf_L = 0
    return delta_leaf_L

def calculate_deltaL_postE(leaf_L, leaf_Lmax, sucrose, delta_t):
    """ delta of leaf length, from the emergence of the previous leaf to the end of elong (predefined elong kinetic depending on leaf state).

    :Parameters:
        - `leaf_L` (:class:`float`) - Total leaf length (m)
        - `leaf_Lmax` (:class:`float`) - Final leaf length (m)
        - `sucrose` (:class:`float`) - Amount of sucrose (µmol C)
    :Returns:
        delta delta_leaf_L (m)
    :Returns Type:
        :class:`float`
    """
    if sucrose > 0:
        delta_leaf_L = parameters.K * leaf_L * max(((leaf_L / leaf_Lmax) - 1)**2, parameters.EPSILON**2)**(parameters.N) * delta_t
    else:
        delta_leaf_L = 0
    return delta_leaf_L

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

def calculate_lamina_L(leaf_L, leaf_dist_to_emerge):
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
        raise Warning('The leaf is shorther than its pseudostem')
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

def calculate_internode_dist_to_emerge(previous_leaf_dist_to_emerge,  previous_sheath_visible_L, previous_sheath_final_hidden_L):
    """ distance for the internnode to be visible given by the previous sheaths.

    :Parameters:
        - `previous_leaf_dist_to_emerge` (:class:`float`) - Distance for the previous leaf to emerge out of the pseudostem (m). Could be None is no previous hiddenzone found.
        - `previous_sheath_visible_L` (:class:`float`) - Visible length of the previous sheath (m).
        - `previous_sheath_final_hidden_L` (:class:`float`) - Final hidden length of the previous sheath (m).
    :Returns:
        Distance for the internode to be visible (m)
    :Returns Type:
        :class:`float`
    """
    if previous_leaf_dist_to_emerge:
        internode_dist_to_emerge = previous_leaf_dist_to_emerge + previous_sheath_visible_L
    else:
        internode_dist_to_emerge = previous_sheath_final_hidden_L + previous_sheath_visible_L # here 'previous_sheath_visible_L' is also the final visible length of the previous sheath because the previous leaf has fully elongated.
    return internode_dist_to_emerge

def calculate_internode_Lmax(internode_L_lig_curr):
    """ Final internode length.

    :Parameters:
        - `internode_L_lig_curr` (:class:`float`) - Internode length at the ligulation of the leaf on the same phytomer (m)
    :Returns:
        Final internode length (m)
    :Returns Type:
        :class:`float`
    """
    return internode_L_lig_curr * parameters.Y0 # So far same parameter than for leaves

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

def calculate_delta_internode_L_preL(sucrose, internode_L, amino_acids, mstruct, delta_t):
    """ delta of internode length over delta_t as a function of sucrose and amino acids, from initiation to the ligulation of the leaf on the same phytomer.

    :Parameters:
        - `sucrose` (:class:`float`) - Amount of sucrose (µmol C)
        - `internode_L` (:class:`float`) - Total internode length (m)
        - `amino_acids` (:class:`float`) - Amount of amino acids (µmol N)
        - `mstruct` (:class:`float`) - Structural mass (g)
    :Returns:
        delta delta_internode_L (m)
    :Returns Type:
        :class:`float`
    """
    if sucrose > 0:
        delta_internode_L = internode_L * ((sucrose / mstruct) / (parameters.Kc + (sucrose / mstruct))) * (((amino_acids/mstruct) **3) / (parameters.Kn**3 + (amino_acids / mstruct)**3)) * parameters.RERmax * delta_t
    else:
        delta_internode_L = 0
    return delta_internode_L

def calculate_delta_internode_L_postL(internode_L, internode_Lmax, sucrose, delta_t):
    """ delta of internode length, from the ligulation of the leaf on the same phytomer to the end of elong (predefined elong kinetic depending on internode state).

    :Parameters:
        - `internode_L` (:class:`float`) - Total internode length (m)
        - `internode_Lmax` (:class:`float`) - Final internode length (m)
        - `sucrose` (:class:`float`) - Amount of sucrose (µmol C)
    :Returns:
        delta delta_internode_L (m)
    :Returns Type:
        :class:`float`
    """
    if sucrose > 0:
        delta_internode_L = parameters.K * internode_L * max(((internode_L / internode_Lmax) - 1)**2, parameters.EPSILON**2)**(parameters.N) * delta_t
    else:
        delta_internode_L = 0
    return delta_internode_L

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
    return internode_L >= internode_Lmax

