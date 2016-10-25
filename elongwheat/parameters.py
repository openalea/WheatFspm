# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    elongwheat.parameters
    ~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`elongwheat.parameters` defines the constant parameters.

    :copyright: Copyright 2015 INRA-ECOSYS, see AUTHORS.
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

EC_wmax = 0.3 #: variation de + ou - 15% de maximal leaf width (SU)
Y0 = 137 #: Facteur agrandissement feuille en mode automate (SU)
K = 5.465e-06 #: Parameter of the elong function after previous leaf emergence (s-1)
N = 0.4317 #: Parameter of the elong function after previous leaf emergence
Ksslw = 4160 #10000 #: Affinité SSLW aux fructanes (µmol C g-1)
Kc = 145.6 #350 #: affinité du RER au C (µmol/g)
Kn = 16.64 #40 #: affinité du RER à N (µmol/g)
min_SSLW = 22 #: g m-2
max_SSLW = 50 #: g m-2
ratio_SSSW_SSLW = 5 # ratio gaine/limbe des matieres seches structurales spécifiques (calculé depuis les données de J. Bertheloot, 2004)
RERmax = 4e-06 #: s-1
EPSILON = 0.01 #: A threshold, expressed in relative leaf length that remains to be produced, under which the rate of leaf elongation will be assumed as constant


class HiddenZoneInit:
    """
    Initial values for hidden zones
    """

    leaf_is_growing = True
    hiddenzone_L = 0
    delta_hiddenzone_L = 0
    leaf_L = 4E-08
    delta_leaf_L = 0
    leaf_Lmax = 0
    leaf_Lem_prev = 0
    lamina_Lmax = 0
    sheath_Lmax = 0
    leaf_Wmax = 0
    SSLW = 0
    SSSW = 0
    leaf_is_emerged = False
    sucrose = 1E-3
    amino_acids = 1E-3
    fructan = 0
    mstruct = 2.65E-08

class OrganInit:
    """
    Initial values for organs
    """

    visible_length = 0
    is_growing = True
    final_hidden_length = 0
    length = 0