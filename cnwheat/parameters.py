# -*- coding: latin-1 -*-

"""
    cnwheat.parameters
    ~~~~~~~~~~~~~~~~~~

    The module :mod:`cnwheat.parameters` defines the constant parameters in a population of plants.

    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2014.
"""

"""
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
"""

class PopulationParameters:
    """
    Internal parameters of populations.
    """
    pass


class PlantParameters:
    """
    Internal parameters of plants.
    """
    pass


class AxisParameters:
    """
    Internal parameters of axes.
    """
    pass


class PhytomerParameters:
    """
    Internal parameters of phytomers.
    """
    pass


class OrganParameters:
    """
    Internal parameters common to all the organs.
    """
    MSTRUCT_AXIS = 2.08                     #: Structural mass  of a plant (g) (Bertheloot, 2011)
    ALPHA_AXIS = 1                          #: Proportion of the structural mass containing the substrates
    DELTA_T = 3600                          #: Timestep of the model (s)

    C_MOLAR_MASS = 12                       #: Molar mass of carbon (g mol-1)
    NB_C_SUCROSE = 12                       #: Number of C in 1 mol of sucrose
    SUCROSE_MOLAR_MASS_C_RATIO = 0.4        #: Contribution of C in sucrose mass

    AMINO_ACIDS_C_RATIO = 3.67              #: Mean number of mol of C in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)
    AMINO_ACIDS_N_RATIO = 1.17              #: Mean number of mol of N in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)
    AMINO_ACIDS_MOLAR_MASS_N_RATIO = 0.145  #: Mean contribution of N in amino acids mass
    N_MOLAR_MASS = 14                       #: Molar mass of nitrogen (g mol-1)
    RATIO_EXPORT_NITRATES_ROOTS = 0.1       #: Proportion of uptaked nitrates actually exported from roots to shoot (1-RATIO_EXPORT_NITRATES_ROOTS = part of nitrates staying in roots)


class PhloemParameters(OrganParameters):
    """
    Internal parameters of phloems.
    """
    ALPHA = 1 #: Proportion of structural mass containing substrate


class GrainsParameters(OrganParameters):
    """
    Internal parameters of grains.
    """
    ALPHA = 1                                   #: Proportion of structural mass containing substrate
    AMINO_ACIDS_MOLAR_MASS_N_RATIO = 0.136      #: Mean contribution of N in amino acids mass contained in gluten (Glu, Gln and Pro)

    # Structure parameters
    VMAX_RGR = 1.9e-06                          #: Maximal value of the Relative Growth Rate of grain structure (dimensionless)
    K_RGR = 300                                 #: Affinity coefficient of the Relative Growth Rate of grain structure (dimensionless)

    # Starch parameters
    VMAX_STARCH = 0.75                          #: Maximal rate of grain filling of starch (µmol C s-1 g-1 MS)
    K_STARCH = 100                              #: Affinity coefficient of grain filling of starch (µmol C g-1 MS)

    Y_GRAINS = 0.75                             #: Proportion of C loaded from phloem actually used for grain structure and starch (1 - Y_GRAINS is a kind of growth respiration)
    FILLING_INIT = 360                          #: Time (h) at which phloem loading switch from grain structure to accumulation of starch


class RootsParameters(OrganParameters):
    """
    Internal parameters of roots.
    """
    ALPHA = 1                                 #: Proportion of structural mass containing substrate

    VMAX_SUCROSE_UNLOADING = 0.02             #: Maximal rate of sucrose unloading from phloem to roots (µmol C sucrose s-1 g-1 MS)
    K_SUCROSE_UNLOADING = 100                 #: Affinity coefficient of sucrose unloading from phloem to roots (µmol C sucrose g-1 MS)
    VMAX_AMINO_ACIDS_UNLOADING = 0.001        #: Maximal rate of amino acids unloading from phloem to roots (µmol N amino acids s-1 g-1 MS)
    K_AMINO_ACIDS_UNLOADING = 100             #: Affinity coefficient of amino acids unloading from phloem to roots (µmol N amino acids g-1 MS)

    # Nitrate uptake
    K_TR_UPTAKE_NITRATES = 1E-2     #: Affinity coefficient of nitrate uptake by roots (mm H20)
    A_VMAX_HATS = 0.1333            #: Parameter for estimating the maximal rate of nitrates uptake at saturating soil N concentration;HATS (dimsensionless)
    LAMBDA_VMAX_HATS = 0.0025       #: Parameter for estimating the maximal rate of nitrates uptake at saturating soil N concentration;HATS (s-1)
    A_K_HATS = 211812               #: Parameter for estimating the affinity coefficient of nitrates uptake at saturating soil N concentration;HATS (dimsensionless)
    LAMBDA_K_HATS = 0.0018          #: Parameter for estimating the affinity coefficient of nitrates uptake at saturating soil N concentration;HATS (g m-3)
    A_LATS = 4.614E-09              #: Parameter for estimating the rate of nitrates uptake at low soil N concentration; LATS (dimensionless)
    LAMBDA_LATS = 1.6517E-03        #: Parameter for estimating the rate of nitrates uptake at low soil N concentration; LATS (m3 µmol-1 s-1)

    # Amino acids
    VMAX_AMINO_ACIDS = 0.001        #: Maximal rate of amino acid synthesis (µmol N s-1 g-1 MS)
    K_AMINO_ACIDS_NITRATES = 3      #: Affinity coefficient of amino acid synthesis from nitrates (µmol N g-1 MS)
    K_AMINO_ACIDS_SUCROSE = 3       #: Affinity coefficient of amino acid synthesis from triosesP (µmol C g-1 MS)
    K_TR_EXPORT_AMINO_ACIDS = 1E-8  #: Affinity coefficient of amino acids export from roots to shoot (mm H20)


class PhotosyntheticOrganParameters(OrganParameters):
    """
    Internal parameters of photosynthetic organs.
    """
    # Sucrose
    VMAX_SUCROSE = 1                #: Maximal rate of sucrose synthesis (µmol C s-1 g-1 MS)
    K_SUCROSE = 0.66                #: Affinity coefficient of sucrose synthesis (µmol C g-1 MS)

    # Starch
    VMAX_STARCH = 2                 #: Maximal rate of starch synthesis (µmol C s-1 g-1 MS)
    K_STARCH = 20                   #: Affinity coefficient of starch synthesis (µmol C g-1 MS)
    DELTA_DSTARCH = 0.0001          #: Relative rate of starch degradation (s-1)

    # Fructans
    VMAX_SFRUCTAN = 0.2             #: Maximal rate of fructan synthesis (µmol C s-1 g-1 MS)
    K_SFRUCTAN = 5000               #: Affinity coefficient of fructan synthesis (µmol C g-1 MS)
    N_SFRUCTAN = 3                  #: Number of "substrates" for fructan synthesis (dimensionless)
    VMAX_REGUL_SFRUCTAN = 1         #: Maximal value of the regulation function of fructan synthesis (dimensionless)
    K_REGUL_SFRUCTAN = 8            #: Affinity coefficient of the regulation function of fructan synthesis (dimensionless)
    N_REGUL_SFRUCTAN = 15           #: Parameter of the regulation function of fructan synthesis (dimensionless)
    VMAX_DFRUCTAN = 0.035           #: Maximal rate of fructan degradation (µmol C s-1 g-1 MS)
    K_DFRUCTAN = 100                #: Affinity coefficient of fructan degradation (µmol C g-1 MS)

    # Loading sucrose
    SIGMA_SUCROSE = 1.85e-07        #: Conductivity of an organ-phloem pathway (g mol-1 m-2 s-1) ; used to compute the sucrose loaded to the phloem
    SIGMA_AMINO_ACIDS = 1e-06       #: Conductivity of an organ-phloem pathway (g mol-1 m-2 s-1) ; used to compute the amino acids loaded to the phloem
    BETA = 1                        #: Kind of volumetric mass density at power -2/3 ((g m-3)**(-2/3))

    # Amino acids
    VMAX_AMINO_ACIDS = 1            #: Maximal rate of amino acid synthesis (µmol N s-1 g-1 MS)
    K_AMINO_ACIDS_NITRATES = 3      #: Affinity coefficient of amino acid synthesis from nitrates (µmol N g-1 MS)
    K_AMINO_ACIDS_TRIOSESP = 0.2    #: Affinity coefficient of amino acid synthesis from triosesP (µmol C g-1 MS)

    # Proteins
    VMAX_SPROTEINS = 0.002          #: Maximal rate of protein synthesis (µmol N s-1 g-1 MS)
    K_SPROTEINS = 100               #: Affinity coefficient of protein synthesis (µmol N g-1 MS)
    DELTA_DPROTEINS = 1.85E-6       #: Relative rate of protein degradation (s-1)


class ChaffParameters(PhotosyntheticOrganParameters):
    """
    Internal parameters of chaffs.
    """
    pass


class LaminaParameters(PhotosyntheticOrganParameters):
    """
    Internal parameters of laminae.
    """
    pass


class InternodeParameters(PhotosyntheticOrganParameters):
    """
    Internal parameters of internodes.
    """
    pass


class PeduncleParameters(PhotosyntheticOrganParameters):
    """
    Internal parameters of peduncles.
    """
    pass


class SheathParameters(PhotosyntheticOrganParameters):
    """
    Internal parameters of sheaths.
    """
    pass


class PhotosyntheticOrganElementParameters:
    """
    Internal parameters of photosynthetic organs elements
    """
    pass


class ChaffElementParameters(PhotosyntheticOrganElementParameters):
    """
    Internal parameters of chaffs elements.
    """
    ALPHA = 1 #: Proportion of structural mass containing substrate


class LaminaElementParameters(PhotosyntheticOrganElementParameters):
    """
    Internal parameters of laminae elements.
    """
    ALPHA = 1 #: Proportion of structural mass containing substrate

    #: Temporary estimation of lamina senescence ({lamina_index: (time of senescence beginning (h), offset of the linear regression)})
    INFLEXION_POINTS = {1: (585, 8.57E-3),
                        2: (465, 7.61E-3),
                        3: (345, 6.64E-3)}


class InternodeElementParameters(PhotosyntheticOrganElementParameters):
    """
    Internal parameters of internodes elements.
    """
    ALPHA = 1 #: Proportion of structural mass containing substrate


class PeduncleElementParameters(PhotosyntheticOrganElementParameters):
    """
    Internal parameters of peduncles elements.
    """
    ALPHA = 1 #: Proportion of structural mass containing substrate


class SheathElementParameters(PhotosyntheticOrganElementParameters):
    """
    Internal parameters of sheaths elements.
    """
    ALPHA = 1 #: Proportion of structural mass containing substrate





