# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

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
        $LastChangedBy: mngauthier $
        $LastChangedDate: 2018-07-20 15:10:01 +0200 (ven., 20 juil. 2018) $
        $LastChangedRevision: 40 $
        $URL: https://subversion.renater.fr/growth-wheat/trunk/trunk/growthwheat/parameters.py $
        $Id: parameters.py 40 2018-07-20 13:10:01Z mngauthier $
"""
C_MOLAR_MASS = 12                     #: Carbon molar mass (g mol-1)
N_MOLAR_MASS = 14                     #: Nitrogen molar mass (g mol-1)

## Shoot
ALPHA = 0.106 # 1.537e-02                     #: Parameter of the relation between leaf mass and leaf length (g m^(-BETA))
BETA = 1.28                           #: Parameter of the relation between leaf mass and leaf length (dimensionless)
RATIO_SUCROSE_MSTRUCT = 0.384        #: Mass of C (under carbohydrate form, g) in 1 g of mstruct (Penning de Vries, Witlage and Kremer, 1978)
RATIO_AMINO_ACIDS_MSTRUCT = 0.005#0.0322   #: Mass of N (under amino acid/protein form, g) in 1 g of mstruct (Penning de Vries, Witlage and Kremer, 1978)
AMINO_ACIDS_C_RATIO = 3.67            #: Mean number of mol of C in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)
AMINO_ACIDS_N_RATIO = 1.17            #: Mean number of mol of N in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)
RATIO_MSTRUCT_DM = 0.8                #: Ratio mstruct/dry matter (dimensionless)
RATIO_ENCLOSED_LEAF_INTERNODE = 5     #: We use ratio sheath:lamina of the specific structural dry masses (from data of J. Bertheloot, 2004)

te = 271 * 3600                       #: Parametre end elongation feuille en mode automate (s); Fournier 2005 sur courbe corrigee - See elongwheat parameters
te_IN = 210 * 3600                    #: Parametre end elongation EN en mode automate (s); Malvoisin 1984 II - See elongwheat parameters


## Roots
VMAX_ROOTS_GROWTH = 0.015             #: Maximal rate of root structural dry matter growth (µmol C s-1 g-1 MS)
K_ROOTS_GROWTH = 1250                 #: Affinity coefficient of root structural dry matter growth (µmol C g-1 MS)
RATIO_C_MSTRUCT_ROOTS = 0.384         #: Mean contribution of carbon to root structural dry mass (g C g-1 Mstruct)
RATIO_N_MSTRUCT_ROOTS_ = 0.02         #: Mean contribution of nitrogen to root structural dry mass (g N g-1 Mstruct)




class HiddenZoneInit:
    """
    Initial values for hidden zones
    """
    def __init__(self):

        ## MG : I think none of these parameters is usefull since new hiddenzones are only set by elongwheat which already defines all these values
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
        self.leaf_Lmax = None # no calculation before emergence Ln-1
        self.lamina_Lmax = None # no calculation before emergence Ln-1
        self.sheath_Lmax = None # no calculation before emergence Ln-1
        self.leaf_Wmax = None # no calculation before emergence Ln-1
        self.SSLW = None # no calculation before emergence Ln-1
        self.SSSW = None # no calculation before emergence Ln-1
        self.leaf_is_emerged = False
        self.leaf_enclosed_mstruct = 2.65E-08
        self.internode_enclosed_mstruct = 0
        self.mstruct = self.leaf_enclosed_mstruct + self.internode_enclosed_mstruct
        self.leaf_enclosed_Nstruct = self.leaf_enclosed_mstruct * 0.0322 # parameter value in growth wheat
        self.internode_enclosed_Nstruct = self.internode_enclosed_mstruct * 0.0322 # parameter value in growth wheat
        self.sucrose = 1E-3
        self.amino_acids = 1E-3
        self.fructan = 0
        self.proteins = 0
        self.cytokinins = 0 #: AU


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
        self.conc_cytokinins = 15     #: AU / g mstruct