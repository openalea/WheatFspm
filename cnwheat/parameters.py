# -*- coding: latin-1 -*-

"""
    cnwheat.parameters
    ~~~~~~~~~~~~~~~~~~

    The module :mod:`cnwheat.parameters` defines the constant parameters in a population of plants.

    :copyright: Copyright 2014-2017 INRA-ECOSYS, see AUTHORS.
    :license: CeCILL-C, see LICENSE for details.
    
    **Acknowledgments**: The research leading these results has received funding through the 
    Investment for the Future programme managed by the Research National Agency 
    (BreedWheat project ANR-10-BTBR-03).
    
    .. seealso:: Barillot et al. 2016.
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
    ALPHA = 1                          #: Proportion of the structural mass containing the substrates


class PhytomerParameters:
    """
    Internal parameters of phytomers.
    """
    pass


class OrganParameters:
    """
    Internal parameters common to all the organs.
    """

    # class model.Organ
    C_MOLAR_MASS = 12                       #: Molar mass of carbon (g mol-1)
    NB_C_TRIOSEP = 3                        #: Number of C in 1 mol of trioseP
    NB_C_HEXOSES = 6                        #: Number of C in 1 mol of hexoses (glucose, fructose)
    NB_C_SUCROSE = 12                       #: Number of C in 1 mol of sucrose
    HEXOSE_MOLAR_MASS_C_RATIO = 0.4         #: Contribution of C in hexose mass
    RATIO_C_MSTRUCT = 0.384                 #: Mean contribution of carbon to structural dry mass (g C g-1 Mstruct)

    AMINO_ACIDS_C_RATIO = 3.67              #: Mean number of mol of C in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)
    AMINO_ACIDS_N_RATIO = 1.17              #: Mean number of mol of N in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)
    AMINO_ACIDS_MOLAR_MASS_N_RATIO = 0.145  #: Mean contribution of N in amino acids mass
    N_MOLAR_MASS = 14                       #: Molar mass of nitrogen (g mol-1)


class HiddenZoneParameters(OrganParameters):
    SIGMA = 5E-2                          #: Coefficient de diffusion surfacique. Utilisé dans la loi de Fick (g m-2 s-1)
    Vmax_Regul_Sfructans = 1
    K_Regul_Sfructans = 0.5
    n_Regul_Sfructans = 15
    Vmax_Sfructans = 0.2 # µmol/g/s
    delta_Dproteins = 1.85e-06

class HiddenZoneInitCompartments(object):
    """
    Initial values for hidden zones
    """
    def __init__(self):
        self.sucrose = 1E-3      #: µmol C
        self.fructan = 0         #: µmol C
        self.amino_acids = 1E-3  #: µmol N
        self.proteins = 0        #: µmol N
        self.mstruct = 6.39E-08  #: g
        self.Nstruct = 2.06E-09  #: g

class PhloemParameters(OrganParameters):
    """
    Internal parameters of phloems.
    """
    ALPHA = 1 #: Proportion of structural mass containing substrate

class PhloemInitCompartments(object):
    """
    Initial values for phloem
    """
    def __init__(self):
        self.sucrose = 500       #: µmol C
        self.amino_acids = 100   #: µmol N

class GrainsParameters(OrganParameters):
    """
    Internal parameters of grains.
    """
    ALPHA = 1                                   #: Proportion of structural mass containing substrate
    AMINO_ACIDS_MOLAR_MASS_N_RATIO = 0.136      #: Mean contribution of N in amino acids mass contained in gluten (Glu, Gln and Pro)

    # Structure parameters
    VMAX_RGR = 1.5e-06                          #: Maximal value of the Relative Growth Rate of grain structure (s-1)
    K_RGR = 300                                 #: Affinity coefficient of the Relative Growth Rate of grain structure (µmol C)

    # Starch parameters
    VMAX_STARCH = 0.35                          #: Maximal rate of grain filling of starch (µmol C s-1 g-1 MS)
    K_STARCH = 400                              #: Affinity coefficient of grain filling of starch (µmol C g-1 MS)

    FILLING_INIT = 360 * 3600                   #: Time (s) at which phloem loading switch from grain structure to accumulation of starch
    FILLING_END = 900 * 3600                    #: Time (s) at which grains filling stops. (Bertheloot et al., 2011)

class GrainsInitCompartments(object):
    """
    Initial values for grains
    """
    def __init__(self):
        self.age_from_flowering = 0 #: second
        self.starch = 0             #: µmol C
        self.structure = 1          #: µmol C
        self.proteins = 0           #: µmol N

class RootsParameters(OrganParameters):
    """
    Internal parameters of roots.
    """
    ALPHA = 1                       #: Proportion of structural mass containing substrate

    VMAX_SUCROSE_UNLOADING = 0.03   #: Maximal rate of sucrose unloading from phloem to roots (µmol C sucrose s-1 g-1 MS)
    K_SUCROSE_UNLOADING = 1000      #: Affinity coefficient of sucrose unloading from phloem to roots (µmol C sucrose g-1 MS)

    # Regulation function by transpiration of nitrate uptake
    K_TRANSPIRATION = 1             #: Affinity coefficient for the regulation function by culm transpiration (mmol H20 m-2 s-1)

    # Regulation function by C in roots of nitrate uptake
    K_C = 4000                      #: Affinity coefficient for the regulation function by root C (µmol C sucrose g-1 MS)

    # Nitrate uptake
    NET_INFLUX_UPTAKE_RATIO = 0.6   #: ratio (net uptake : nitrate influx)
    A_VMAX_HATS = 0.1333            #: Parameter for estimating the maximal rate of nitrates uptake at saturating soil N concentration;HATS (dimensionless)
    LAMBDA_VMAX_HATS = 0.0025       #: Parameter for estimating the maximal rate of nitrates uptake at saturating soil N concentration;HATS (s-1)
    A_K_HATS = 211812               #: Parameter for estimating the affinity coefficient of nitrates uptake at saturating soil N concentration;HATS (dimsensionless)
    LAMBDA_K_HATS = 0.0018          #: Parameter for estimating the affinity coefficient of nitrates uptake at saturating soil N concentration;HATS (g m-3)
    A_LATS = 4.614E-09              #: Parameter for estimating the rate of nitrates uptake at low soil N concentration; LATS (dimensionless)
    LAMBDA_LATS = 1.6517E-03        #: Parameter for estimating the rate of nitrates uptake at low soil N concentration; LATS (m3 µmol-1 s-1)


    # Nitrate export
    K_NITRATE_EXPORT = 1E-6         #: Relative rate of nitrate export from roots (s-1)

    # Amino acids
    VMAX_AMINO_ACIDS = 0.001        #: Maximal rate of amino acid synthesis (µmol N s-1 g-1 MS)
    K_AMINO_ACIDS_NITRATES = 3      #: Affinity coefficient of amino acid synthesis from nitrates (µmol N g-1 MS)
    K_AMINO_ACIDS_SUCROSE = 350     #: Affinity coefficient of amino acid synthesis from triosesP (µmol C g-1 MS)
    K_AMINO_ACIDS_EXPORT = 3E-5     #: Relative rate of amino acids export from roots (s-1)

    # Exudation
    C_EXUDATION = 0.20              #: Proportion of C exudated from C sucrose unloaded to roots (Keith et al., 1986)

    # Cytokinins
    VMAX_S_CYTOKININS = 4.5E-04     #: Maximal rate of cytokinins synthesis (UA g-1 mstruct s-1)
    K_NITRATES_CYTOKININS = 200     #: Affinity coefficient of cytokinins synthesis for nitrates (µmol N nitrates g-1 mstruct)
    K_SUCROSE_CYTOKININS = 1250     #: Affinity coefficient of cytokinins synthesis for sucrose (µmol C sucrose g-1 mstruct)
    N_SUC_CYTOKININS = 10           #: A parameter for cytokinins synthesis (dimensionless)
    N_NIT_CYTOKININS = 0.7          #: A parameter for cytokinins synthesis (dimensionless)
    K_CYTOKININS_EXPORT = 2E-4      #: Relative rate of cytokinins export from roots (s-1)

class RootsInitCompartments(object):
    """
    Initial values for roots
    """
    def __init__(self):
        self.sucrose = 0       #: µmol C
        self.nitrates = 0      #: µmol C
        self.amino_acids = 0   #: µmol N
        self.cytokinins = 0    #: AU
        self.mstruct = 0.15    #: g
        self.Nstruct = 0.0045  #: g

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
    VMAX_SFRUCTAN_POT = 0.015       #: Potential maximal rate of fructan synthesis (µmol C s-1 g-1 MS)
    K_SFRUCTAN = 5000               #: Affinity coefficient of fructan synthesis (µmol C g-1 MS)
    K_REGUL_SFRUCTAN = 0.001        #: Affinity coefficient of the regulation function of fructan synthesis (µmol g-1 MS)
    N_REGUL_SFRUCTAN = 3            #: Parameter of the regulation function of fructan synthesis (dimensionless)
    VMAX_DFRUCTAN = 0.035           #: Maximal rate of fructan degradation (µmol C s-1 g-1 MS)
    K_DFRUCTAN = 100                #: Affinity coefficient of fructan degradation (µmol C g-1 MS)

    # Loading sucrose and amino acids
    SIGMA_SUCROSE = 1e-08           #: Conductivity of an organ-phloem pathway (g2 µmol-1 m-2 s-1) ; used to compute the sucrose loaded to the phloem
    SIGMA_AMINO_ACIDS = 1e-07       #: Conductivity of an organ-phloem pathway (g2 µmol-1 m-2 s-1) ; used to compute the amino acids loaded to the phloem
    BETA = 1                        #: Kind of volumetric mass density at power -2/3 ((g m-3)**(-2/3))

    # Amino acids
    VMAX_AMINO_ACIDS = 1            #: Maximal rate of amino acid synthesis (µmol N s-1 g-1 MS)
    K_AMINO_ACIDS_NITRATES = 3      #: Affinity coefficient of amino acid synthesis from nitrates (µmol N g-1 MS)
    K_AMINO_ACIDS_TRIOSESP = 0.2    #: Affinity coefficient of amino acid synthesis from triosesP (µmol C g-1 MS)

    # Proteins
    VMAX_SPROTEINS = 0.0015         #: Maximal rate of protein synthesis (µmol N s-1 g-1 MS)
    K_SPROTEINS = 100               #: Affinity coefficient of protein synthesis (µmol N g-1 MS)
    VMAX_DPROTEINS = 2.5E-6         #: Maximal rate of protein degradation (µmol g-1 mstruct s-1)
    K_DPROTEINS = 50                #: Affinity coefficient with cytokinins for protein degradation (UA g-1 mstruct)
    N_DPROTEINS = 2.1               #: A coefficient for the regulation of protein degradation by cytokines (dimensionless)

    # cytokinins
    DELTA_D_CYTOKININS = 3.E-6      #: Relative rate of cytokinins degradation (s-1)


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

class PhotosyntheticOrganElementInitCompartments(object):
    """
    Initial values for photosynthetic organ elements
    """
    def __init__(self):
        self.green_area = 1E-4  #: m2
        self.mstruct = 5E-3     #: g
        self.Nstruct = 1E-3     #: g

        self.triosesP = 0       #: µmol C
        self.starch = 0         #: µmol C
        self.sucrose = 0        #: µmol C
        self.fructan = 0        #: µmol C
        self.nitrates = 0       #: µmol N
        self.amino_acids = 0    #: µmol N
        self.proteins = 0       #: µmol N
        self.cytokinins = 15    #: AU

        self.is_growing = False #: Flag indicating if the element is growing or not (:class:`bool`)
        self.Tr = 0             #: Transpiration rate (mmol m-2 s-1)
        self.Ag = 0             #: Gross assimilation (µmol m-2 s-1)
        self.Ts = 15            #: Organ temperature (degree Celsius)

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


class SoilParameters:
    """
    Internal parameters of soil.
    """
    MINERALISATION_RATE = 2.05E-6           # Mineralisation rate (µmol N nitrates m-3 s-1)