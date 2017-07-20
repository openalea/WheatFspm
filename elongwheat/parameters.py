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
PLASTO_leaf = 50 # Leaf pastochron (°C d-1)
max_nb_leaves = 10 # Max number of leaves per axis
sowing_depth = 0.05 # Sowing depth (m) used to define plant emergence
internode_L_init = 4E-08
nb_PLASTO_internode_init = 5 # From Malvoisin 1984b, associated with priodia of 1.10-4 m
delta_TT_GA = 0 # Thermal time between floral transition of SAM and Gibberelin production

class HiddenZoneInit(object):
    """
    Initial values for hidden zones
    """
    def __init__(self):
        self.leaf_is_growing = True
        self.internode_is_growing = False
        self.leaf_dist_to_emerge = 4E-5          #: m TODO: ok comme valeur?
        self.delta_leaf_dist_to_emerge = 0       #: m
        self.internode_dist_to_emerge = 0        #: m
        self.delta_internode_dist_to_emerge = 0  #: TODO: needed??
        self.leaf_L = 4E-5                       #: m
        self.delta_leaf_L = 0                    #: TODO: needed??
        self.internode_L = 0                     #: m
        self.delta_internode_L = 0 #: needed??
        self.leaf_Lmax = None                    #: m, no calculation before emergence Ln-1
        self.lamina_Lmax = None                  #: m, no calculation before emergence Ln-1
        self.sheath_Lmax = None                  #: m, no calculation before emergence Ln-1
        self.leaf_Wmax = None                    #: m, no calculation before emergence Ln-1
        self.SSLW = None                         #: g m-2, no calculation before emergence Ln-1
        self.SSSW = None                         #: g m-2, no calculation before emergence Ln-1
        self.leaf_is_emerged = False
        self.internode_Lmax = None               #: m, no calculation before ligulation Ln
        self.SSIW = None                         #: g m-2, no calculation before ligulation Ln
        self.internode_is_visible = False
        self.sucrose = 1E-3                      #: µmol C
        self.amino_acids = 1E-3                  #: µmol N
        self.fructan = 0                         #: µmol C
        self.leaf_enclosed_mstruct = 2.65E-08    #: g TODO: viens d'ou?
        self.internode_mstruct = 0               #: g
        self.mstruct = self.leaf_enclosed_mstruct + self.internode_mstruct #: g
        self.leaf_enclosed_Nstruct = self.leaf_enclosed_mstruct * 0.0322 #: g, parameter value in growth wheat #: g
        self.internode_Nstruct = self.internode_mstruct * 0.0322 #: g, parameter value in growth wheat
        self.proteins = 0                        #: µmol N

class OrganInit(object):
    """
    Initial values for emerged and growing organs
    """
    def __init__(self):
        self.length = 1E-3            #: m
        self.visible_length = 1E-3    #: m
        self.final_hidden_length = 0  #: m
        self.is_growing = True