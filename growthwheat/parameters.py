# -*- coding: latin-1 -*-

from __future__ import division  # use "//" to do integer division

"""
    growthwheat.parameters
    ~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`growthwheat.parameters` defines the constant parameters.

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
C_MOLAR_MASS = 12                     #: Carbon molar mass (g mol-1)
N_MOLAR_MASS = 14                     #: Nitrogen molar mass (g mol-1)

# Shoot
ALPHA = 0.106  # 1.537e-02            #: Parameter of the relation between leaf mass and leaf length (g m^(-BETA))
BETA = 1.28                           #: Parameter of the relation between leaf mass and leaf length (dimensionless)
RATIO_SUCROSE_MSTRUCT = 0.384         #: Mass of C (under carbohydrate form, g) in 1 g of mstruct (Penning de Vries, Witlage and Kremer, 1978)
RATIO_AMINO_ACIDS_MSTRUCT = 0.005     #: Mass of N (under amino acid/protein form, g) in 1 g of mstruct (Penning de Vries, Witlage and Kremer, 1978)
AMINO_ACIDS_C_RATIO = 3.67            #: Mean number of mol of C in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly) # Sources ?? Calculs AS = 4.15
AMINO_ACIDS_N_RATIO = 1.17            #: Mean number of mol of N in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly) # Sources ?? Calculs AS = 1.25
RATIO_MSTRUCT_DM = 0.8                #: Ratio mstruct/dry matter (dimensionless)
RATIO_ENCLOSED_LEAF_INTERNODE = 5     #: We use ratio sheath:lamina of the specific structural dry masses (from data of J. Bertheloot, 2004)

# Leaf Automate elongation
te = 300 * 3600 * 24 / 12  #: end of leaf elongation in automate growth (s at 12°c); fitted from adapted data from Fournier 2005
FITTED_L0 = 0.01557936             #: Fitted value of leaf length at t=0 after rescaling the beta function with L0 (m); Fournier 2005 sur courbe corrigee

# Internode Automate elongation
FITTED_L0_IN = 1/59.0  # 5.2 #: Scaling factor of the internode in automate growth (dimensionless), fitted from Malvoisin 1984 II
te_IN = 210 * 3600 * 24 / 12  #: end of internode elongation in automate growth (s at 12°c) ; fitted from Malvoisin 1984 II

# Roots
VMAX_ROOTS_GROWTH = 0.015             #: Maximal rate of root structural dry matter growth (µmol C s-1 g-1 MS) post flo at 12°C
K_ROOTS_GROWTH = 1250                 #: Affinity coefficient of root structural dry matter growth (µmol C g-1 MS) post flo
RATIO_C_MSTRUCT_ROOTS = 0.384         #: Mean contribution of carbon to root structural dry mass (g C g-1 Mstruct)
RATIO_N_MSTRUCT_ROOTS_ = 0.02         #: Mean contribution of nitrogen to root structural dry mass (g N g-1 Mstruct)


class HiddenZoneInit:
    """
    Initial values for hidden zones
    """
    def __init__(self):

        # MG : I think none of these parameters is usefull since new hiddenzones are only set by elongwheat which already defines all these values
        self.leaf_is_growing = True
        self.internode_is_growing = False
        self.internode_is_mature = False
        self.leaf_dist_to_emerge = 4E-08
        self.delta_leaf_dist_to_emerge = 0
        self.internode_dist_to_emerge = 0
        self.delta_internode_dist_to_emerge = 0
        self.leaf_L = 4E-08
        self.delta_leaf_L = 0
        self.internode_L = 0
        self.delta_internode_L = 0
        self.leaf_Lmax = None  # no calculation before emergence Ln-1
        self.lamina_Lmax = None  # no calculation before emergence Ln-1
        self.sheath_Lmax = None  # no calculation before emergence Ln-1
        self.leaf_Wmax = None  # no calculation before emergence Ln-1
        self.SSLW = None  # no calculation before emergence Ln-1
        self.SSSW = None  # no calculation before emergence Ln-1
        self.leaf_is_emerged = False
        self.leaf_enclosed_mstruct = 2.65E-08
        self.internode_enclosed_mstruct = 0
        self.mstruct = self.leaf_enclosed_mstruct + self.internode_enclosed_mstruct
        self.leaf_enclosed_Nstruct = self.leaf_enclosed_mstruct * 0.0322  # parameter value in growth wheat
        self.internode_enclosed_Nstruct = self.internode_enclosed_mstruct * 0.0322  # parameter value in growth wheat
        self.sucrose = 1E-3
        self.amino_acids = 1E-3
        self.fructan = 0
        self.proteins = 0
        self.cytokinins = 0  #: AU


class OrganInit:
    """
    Initial values for organs
    """
    def __init__(self):
        self.is_growing = True
        self.green_area = 0
        self.sucrose = 0              #: µmol C
        self.amino_acids = 0          #: µmol N
        self.fructan = 0              #: µmol C
        self.proteins = 0             #: µmol N
        self.mstruct = 0              #: g
        self.Nstruct = 0              #: g
        self.cytokinins = 0  #: g
        self.conc_cytokinins = 100.0  #: AU / g mstruct # Not used
