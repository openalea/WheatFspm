# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    growthwheat.model
    ~~~~~~~~~~~~~

    The module :mod:`growthwheat.model` defines the equations of the kinetic of leaf growth (mass flows) according to leaf elongation. Also includes root growth.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2015.
"""

"""
    Information about this versioned file:
        $LastChangedBy: mngauthier $
        $LastChangedDate: 2018-07-20 15:10:01 +0200 (ven., 20 juil. 2018) $
        $LastChangedRevision: 40 $
        $URL: https://subversion.renater.fr/growth-wheat/trunk/trunk/growthwheat/model.py $
        $Id: model.py 40 2018-07-20 13:10:01Z mngauthier $
"""

import parameters

def calculate_delta_leaf_enclosed_mstruct(leaf_L, delta_leaf_L):
    """ Relation between length and mstruct for the leaf segment located in the hidden zone.
    Parameters alpha_mass_growth and beta_mass_growth estimated from Williams (1975) and expressed in g of dry mass. #TODO : Check the ref (Williams 1960?)
    Parameter RATIO_MSTRUCT_DM is then used to convert in g of structural dry mass.

    :Parameters:
        - `leaf_L` (:class:`float`) - Total leaf length (m)
        - `delta_leaf_L` (:class:`float`) - delta of leaf length (m)
    :Returns:
        delta_leaf_enclosed_mstruct (g)
    :Returns Type:
        :class:`float`
    """
    return parameters.ALPHA * parameters.BETA * leaf_L**(parameters.BETA-1) * delta_leaf_L * parameters.RATIO_MSTRUCT_DM

def calculate_delta_leaf_enclosed_mstruct_postE(leaf_pseudostem_L, delta_leaf_pseudostem_L, mstruct, LSSW, leaf_pseudo_age, delta_t):
    """ Relation between length and mstruct for the leaf segment located in the hidden zone after leaf emergence.
    Increase mstruct to match sheath mstruct calculation when it is mature.
    #TODO : Hiddenzone mstruct calculation will not work for sheath sorten than previous one.

    :Parameters:
        - `leaf_L` (:class:`float`) - Total leaf length (m)
        - `delta_leaf_L` (:class:`float`) - delta of leaf length (m)
    :Returns:
        delta_leaf_enclosed_mstruct (g)
    :Returns Type:
        :class:`float`
    """
    increment_mstruct = (leaf_pseudostem_L * LSSW - mstruct)/(parameters.te - leaf_pseudo_age)

    return increment_mstruct * delta_t

def calculate_delta_internode_enclosed_mstruct(internode_L, delta_internode_L):
    """ Relation between length and mstruct for the internode segment located in the hidden zone.
    Same relationship than for enclosed leaf corrected by RATIO_ENCLOSED_LEAF_INTERNODE.
    Parameters alpha_mass_growth and beta_mass_growth estimated from Williams (1975) and expressed in g of dry mass.
    Parameter RATIO_MSTRUCT_DM is then used to convert in g of structural dry mass.

    :Parameters:
        - `internode_L` (:class:`float`) - Enclosed internode length (m)
        - `delta_internode_L` (:class:`float`) - delta of enclosed internode length (m)
    :Returns:
        delta_enclosed_internode_mstruct (g)
    :Returns Type:
        :class:`float`
    """

    return parameters.RATIO_ENCLOSED_LEAF_INTERNODE * parameters.ALPHA * parameters.BETA * internode_L**(parameters.BETA-1) * delta_internode_L * parameters.RATIO_MSTRUCT_DM

def calculate_delta_emerged_tissue_mstruct(SSW, previous_mstruct, area):
    """ delta mstruct of emerged tissue (lamina, sheath and internode). Calculated from tissue area.

    :Parameters:
        - `SSW` (:class:`float`) - Structural Specific Weight (g m-2)
        - `previous_mstruct` (:class:`float`) - mstruct at the previous time step i.e. not yet updated (g)
        - `area` (:class:`float`) - Area at the current time step, as updated by the geometrical model (m2)
    :Returns:
        delta mstruct (g)
    :Returns Type:
        :class:`float`
    """
    updated_mstruct = SSW * area
    delta_mstruct = updated_mstruct - previous_mstruct
    return delta_mstruct

def calculate_delta_Nstruct(delta_mstruct):
    """ delta Nstruct of hidden zone and emerged tissue (lamina and sheath).

    :Parameters:
        - `delta_mstruct` (:class:`float`) - delta of mstruct (g)
    :Returns:
        delta Nstruct (g)
    :Returns Type:
        :class:`float`
    """
    return delta_mstruct * parameters.RATIO_AMINO_ACIDS_MSTRUCT

def calculate_export(delta_mstruct, metabolite, hiddenzone_mstruct):
    """Export of metabolite from the hidden zone towards the emerged part of the leaf integrated over delta_t.

    :Parameters:
        - `delta_mstruct` (:class:`float`) - Delta of structural dry mass of the emerged part of the leaf (g)
        - `metabolite` (:class:`float`) - Metabolite amount in the hidden zone (�mol C or N)
        - `hiddenzone_mstruct` (:class:`float`) - Structural mass of the hidden zone (g)

    :Returns:
        metabolite export (�mol N)
    :Returns Type:
        :class:`float`
    """
    return delta_mstruct * max(0, (metabolite / hiddenzone_mstruct))

def calculate_s_Nstruct_amino_acids(delta_hiddenzone_Nstruct, delta_lamina_Nstruct, delta_sheath_Nstruct):
    """Consumption of amino acids for the calculated mstruct growth (µmol N consumed by mstruct growth)

    :Parameters:
        - `delta_hiddenzone_Nstruct` (:class:`float`) - Nstruct growth of the hidden zone (g)
        - `delta_lamina_Nstruct` (:class:`float`) - Nstruct growth of the lamina (g)
        - `delta_sheath_Nstruct` (:class:`float`) - Nstruct growth of the sheath (g)
    :Returns:
        Amino acid consumption (µmol N)
    :Returns Type:
        :class:`float`
    """
    return (delta_hiddenzone_Nstruct + delta_lamina_Nstruct + delta_sheath_Nstruct) / parameters.N_MOLAR_MASS * 1E6

def calculate_s_mstruct_sucrose(delta_hiddenzone_mstruct, delta_lamina_mstruct, delta_sheath_mstruct, s_Nstruct_amino_acids_N):
    """Consumption of sucrose for the calculated mstruct growth (µmol C consumed by mstruct growth)

    :Parameters:
        - `delta_hiddenzone_mstruct` (:class:`float`) - mstruct growth of the hidden zone (g)
        - `delta_lamina_mstruct` (:class:`float`) - mstruct growth of the lamina (g)
        - `delta_sheath_mstruct` (:class:`float`) - mstruct growth of the sheath (g)
        - `s_Nstruct_amino_acids_N` (:class:`float`) - Total amino acid consumption (µmol N) due to Nstruct (µmol N)
    :Returns:
        Sucrose consumption (µmol C)
    :Returns Type:
        :class:`float`
    """
    s_Nstruct_amino_acids = s_Nstruct_amino_acids_N / parameters.AMINO_ACIDS_N_RATIO   #: µmol of AA
    s_mstruct_amino_acids_C = s_Nstruct_amino_acids * parameters.AMINO_ACIDS_C_RATIO   #: µmol of C coming from AA
    s_mstruct_C = (delta_hiddenzone_mstruct + delta_lamina_mstruct + delta_sheath_mstruct) * parameters.RATIO_SUCROSE_MSTRUCT / parameters.C_MOLAR_MASS * 1E6 #: Total C used for mstruct growth (µmol C)
    s_mstruct_sucrose_C = s_mstruct_C - s_mstruct_amino_acids_C                        #: �ol of coming from sucrose

    return s_mstruct_sucrose_C


## Roots
def calculate_roots_mstruct_growth(sucrose, amino_acids, mstruct, delta_t):
    """Root structural dry mass growth integrated over delta_t

    : Parameters:
        - `sucrose` (:class:`float`) - Amount of sucrose in roots (µmol C)
        - `amino_acids` (:class:`float`) - Amount of amino acids in roots (µmol N)
        - `mstruct` (:class:`float`) - Root structural mass (g)

    : Returns:
        mstruct_C_growth (µmol C), mstruct_growth (g), Nstruct_growth (g), Nstruct_N_growth (µmol N)

    :Returns Type:
        :class:`float`
    """
    conc_sucrose = max(0, sucrose/mstruct)

    mstruct_C_growth = (conc_sucrose * parameters.VMAX_ROOTS_GROWTH) / (conc_sucrose + parameters.K_ROOTS_GROWTH) * delta_t * mstruct     #: root growth in C (µmol of C)
    mstruct_growth = (mstruct_C_growth*1E-6 * parameters.C_MOLAR_MASS) / parameters.RATIO_C_MSTRUCT_ROOTS                                 #: root growth (g of structural dry mass)
    Nstruct_growth = mstruct_growth * parameters.RATIO_N_MSTRUCT_ROOTS_                                                                   #: root growth in N (g of structural dry mass)
    Nstruct_N_growth = min(amino_acids, (Nstruct_growth / parameters.N_MOLAR_MASS)*1E6)                                                   #: root growth in nitrogen (µmol N)

    return mstruct_C_growth, mstruct_growth, Nstruct_growth, Nstruct_N_growth
