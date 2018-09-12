# -*- coding: latin-1 -*-
import pandas as pd

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

SECOND_TO_HOUR_RATE_CONVERSION = 3600


def from_dataframe(object_, dataframe_):
    """Set attributes of *object_* from data in *dataframe_*.

    :Parameters:
        - `object_` (:class:`object`) - The object to set.
        - `dataframe_` (:class:`pandas.DataFrame`) - The dataframe used to set the attribute(s)
          of *object_*.
          *dataframe_* must have only 2 rows:

              * one row is for the header and contains the name of each attribute,
              * and one row contains the value of each attribute.
    """
    object_.__dict__.update(dataframe_.to_dict(orient='index')[dataframe_.first_valid_index()])


def to_dataframe(object_):
    """Create and return a dataframe from attributes of *object_*.

    :Parameters:
        - `object_` (:class:`object`) - The object used to create the dataframe.

    :Returns:
        A dataframe which contains the attributes of *object_*, with only 2 rows:

          * one row is for the header and contains the name of each attribute,
          * and one row contains the value of each attribute.

    :Returns Type:
        :class:`pandas.DataFrame`
    """
    return pd.DataFrame(object_.__dict__, index=[0]).sort_index(axis=1)


class PopulationParameters(object):
    """
    Internal parameters of populations.
    """
    def __init__(self):
        pass


#: The instance of class :class:`cnwheat.parameters.PopulationParameters` for current process
POPULATION_PARAMETERS = PopulationParameters()


class PlantParameters(object):
    """
    Internal parameters of plants.
    """
    def __init__(self):
        pass


#: The instance of class :class:`cnwheat.parameters.PlantParameters` for current process
PLANT_PARAMETERS = PlantParameters()


class AxisParameters(object):
    """
    Internal parameters of axes.
    """
    def __init__(self):
        self.ALPHA = 1                          #: Proportion of the structural mass containing the substrates


#: The instance of class :class:`cnwheat.parameters.AxisParameters` for current process
AXIS_PARAMETERS = AxisParameters()


class PhytomerParameters(object):
    """
    Internal parameters of phytomers.
    """
    def __init__(self):
        pass


#: The instance of class :class:`cnwheat.parameters.PhytomerParameters` for current process
PHYTOMER_PARAMETERS = PhytomerParameters()


class HiddenZoneParameters(object):
    """
    Internal parameters of hidden growing zones.
    """
    def __init__(self):
        self.SIGMA = 5E-2                #: Coefficient of surface diffusion. Used in Fick's law (g m-2 s-1).
        self.Vmax_Regul_Sfructans = 1    #: Maximal rate of fructan synthesis in the regulation function of fructan synthesis (:math:`\mu` mol C s-1 g-1 MS)
        self.K_Regul_Sfructans = 0.5     #: Affinity coefficient of the regulation function of fructan synthesis (:math:`\mu` mol g-1 MS)
        self.n_Regul_Sfructans = 15      #: Parameter of the regulation function of fructan synthesis (dimensionless)
        self.Vmax_Sfructans = 0.2        #: Maximal rate of fructan synthesis (:math:`\mu` mol C s-1 g-1 MS)
        self.delta_Dproteins = 0.25e-6  #1.85e-06  #: Relative rate of proteins degradation (s-1)


#: The instance of class :class:`cnwheat.parameters.HiddenZoneParameters` for current process
HIDDEN_ZONE_PARAMETERS = HiddenZoneParameters()


class HiddenZoneInitCompartments(object):
    """
    Initial values for compartments of hidden zones.
    """
    def __init__(self):
        self.sucrose = 1E-3      #: initial value of sucrose (:math:`\mu` mol C)
        self.fructan = 0         #: initial value of fructan (:math:`\mu` mol C)
        self.amino_acids = 1E-3  #: initial value of amino_acids (:math:`\mu` mol N)
        self.proteins = 0        #: initial value of proteins (:math:`\mu` mol N)
        self.mstruct = 6.39E-08  #: initial value of mstruct (g)
        self.Nstruct = 2.06E-09  #: initial value of Nstruct (g)


#: The instance of class :class:`cnwheat.parameters.HiddenZoneInitCompartments` for current process
HIDDEN_ZONE_INIT_COMPARTMENTS = HiddenZoneInitCompartments()


class PhloemParameters(object):
    """
    Internal parameters of phloems.
    """
    def __init__(self):
        self.ALPHA = 1  #: Proportion of structural mass containing substrate


#: The instance of class :class:`cnwheat.parameters.PhloemParameters` for current process
PHLOEM_PARAMETERS = PhloemParameters()


class PhloemInitCompartments(object):
    """
    Initial values for compartments of phloem.
    """
    def __init__(self):
        self.sucrose = 500       #: initial value of sucrose (:math:`\mu` mol C)
        self.amino_acids = 100   #: initial value of amino_acids (:math:`\mu` mol N)


#: The instance of class :class:`cnwheat.parameters.PhloemInitCompartments` for current process
PHLOEM_INIT_COMPARTMENTS = PhloemInitCompartments()


class GrainsParameters(object):
    """
    Internal parameters of grains.
    """
    def __init__(self):
        self.ALPHA = 1                                   #: Proportion of structural mass containing substrate

        # Structure parameters
        self.VMAX_RGR = 1.5e-06                          #: Maximal value of the Relative Growth Rate of grain structure (s-1)
        self.K_RGR = 300                                 #: Affinity coefficient of the Relative Growth Rate of grain structure (:math:`\mu` mol C)

        # Starch parameters
        self.VMAX_STARCH = 0.35                          #: Maximal rate of grain filling of starch (:math:`\mu` mol C s-1 g-1 MS)
        self.K_STARCH = 400                              #: Affinity coefficient of grain filling of starch (:math:`\mu` mol C g-1 MS)

        self.FILLING_INIT = 360 * 3600                   #: Time (s) at which phloem loading switch from grain structure to accumulation of starch
        self.FILLING_END = 900 * 3600                    #: Time (s) at which grains filling stops. (Bertheloot et al., 2011)


#: The instance of class :class:`cnwheat.parameters.GrainsParameters` for current process
GRAINS_PARAMETERS = GrainsParameters()


class GrainsInitCompartments(object):
    """
    Initial values for compartments of grains.
    """
    def __init__(self):
        self.age_from_flowering = 0  #: initial value of age_from_flowering (s)
        self.starch = 0              #: initial value of starch (:math:`\mu` mol C)
        self.structure = 1           #: initial value of structure (:math:`\mu` mol C)
        self.proteins = 0            #: initial value of proteins (:math:`\mu` mol N)


#: The instance of class :class:`cnwheat.parameters.GrainsInitCompartments` for current process
GRAINS_INIT_COMPARTMENTS = GrainsInitCompartments()


class RootsParameters(object):
    """
    Internal parameters of roots.
    """
    def __init__(self):
        self.ALPHA = 1                       #: Proportion of structural mass containing substrate

        self.VMAX_SUCROSE_UNLOADING = 0.03   #: Maximal rate of sucrose unloading from phloem to roots (:math:`\mu` mol C sucrose s-1 g-1 MS)
        self.K_SUCROSE_UNLOADING = 1000      #: Affinity coefficient of sucrose unloading from phloem to roots (:math:`\mu` mol C sucrose g-1 MS)

        # Regulation function by transpiration of nitrate uptake
        self.K_TRANSPIRATION = 1             #: Affinity coefficient for the regulation function by culm transpiration (mmol H20 m-2 s-1)

        # Regulation function by C in roots of nitrate uptake
        self.K_C = 4000                      #: Affinity coefficient for the regulation function by root C (:math:`\mu` mol C sucrose g-1 MS)

        # Nitrate uptake
        self.NET_INFLUX_UPTAKE_RATIO = 0.6   #: ratio (net uptake : nitrate influx)
        self.A_VMAX_HATS = 0.1333            #: Parameter for estimating the maximal rate of nitrates uptake at saturating soil N concentration;HATS (:math:`\mu` mol g-1 s-1)
        self.LAMBDA_VMAX_HATS = 0.0025       #: Parameter for estimating the maximal rate of nitrates uptake at saturating soil N concentration;HATS (g :math:`\mu` mol-1)
        self.A_K_HATS = 211812               #: Parameter for estimating the affinity coefficient of nitrates uptake at saturating soil N concentration;HATS (:math:`\mu` mol m-3)
        self.LAMBDA_K_HATS = 0.0018          #: Parameter for estimating the affinity coefficient of nitrates uptake at saturating soil N concentration;HATS (g :math:`\mu` mol-1)
        self.A_LATS = 4.614E-09              #: Parameter for estimating the rate of nitrates uptake at low soil N concentration; LATS (m3 g-1 s-1)
        self.LAMBDA_LATS = 1.6517E-03        #: Parameter for estimating the rate of nitrates uptake at low soil N concentration; LATS (g :math:`\mu` mol-1)

        # Nitrate export
        self.K_NITRATE_EXPORT = 1E-6         #: Relative rate of nitrate export from roots (s-1)

        # Amino acids
        self.VMAX_AMINO_ACIDS = 0.001        #: Maximal rate of amino acid synthesis (:math:`\mu` mol N s-1 g-1 MS)
        self.K_AMINO_ACIDS_NITRATES = 3      #: Affinity coefficient of amino acid synthesis from nitrates (:math:`\mu` mol N g-1 MS)
        self.K_AMINO_ACIDS_SUCROSE = 350     #: Affinity coefficient of amino acid synthesis from triosesP (:math:`\mu` mol C g-1 MS)
        self.K_AMINO_ACIDS_EXPORT = 3E-5     #: Relative rate of amino acids export from roots (s-1)

        # Exudation
        self.C_EXUDATION = 0.20              #: Proportion of C exudated from C sucrose unloaded to roots (Keith et al., 1986)

        # Cytokinins
        self.VMAX_S_CYTOKININS = 4.5E-04     #: Maximal rate of cytokinins synthesis (UA g-1 mstruct s-1)
        self.K_NITRATES_CYTOKININS = 200     #: Affinity coefficient of cytokinins synthesis for nitrates (:math:`\mu` mol N nitrates g-1 mstruct)
        self.K_SUCROSE_CYTOKININS = 1250     #: Affinity coefficient of cytokinins synthesis for sucrose (:math:`\mu` mol C sucrose g-1 mstruct)
        self.N_SUC_CYTOKININS = 10           #: A parameter for cytokinins synthesis (dimensionless)
        self.N_NIT_CYTOKININS = 0.7          #: A parameter for cytokinins synthesis (dimensionless)
        self.K_CYTOKININS_EXPORT = 2E-4      #: Relative rate of cytokinins export from roots (s-1)


#: The instance of class :class:`cnwheat.parameters.RootsParameters` for current process
ROOTS_PARAMETERS = RootsParameters()


class RootsInitCompartments(object):
    """
    Initial values for compartments of roots.
    """
    def __init__(self):
        self.sucrose = 0       #: initial value of sucrose (:math:`\mu` mol C)
        self.nitrates = 0      #: initial value of nitrates (:math:`\mu` mol C)
        self.amino_acids = 0   #: initial value of amino_acids (:math:`\mu` mol N)
        self.cytokinins = 0    #: initial value of cytokinins (AU)
        self.mstruct = 0.15    #: initial value of mstruct (g)
        self.Nstruct = 0.0045  #: initial value of Nstruct (g)


#: The instance of class :class:`cnwheat.parameters.RootsInitCompartments` for current process
ROOTS_INIT_COMPARTMENTS = RootsInitCompartments()


class PhotosyntheticOrganParameters(object):
    """
    Internal parameters of photosynthetic organs.
    """
    def __init__(self):
        # Sucrose
        self.VMAX_SUCROSE = 1                #: Maximal rate of sucrose synthesis (:math:`\mu` mol C s-1 g-1 MS)
        self.K_SUCROSE = 0.66                #: Affinity coefficient of sucrose synthesis (:math:`\mu` mol C g-1 MS)

        # Starch
        self.VMAX_STARCH = 2                 #: Maximal rate of starch synthesis (:math:`\mu` mol C s-1 g-1 MS)
        self.K_STARCH = 20                   #: Affinity coefficient of starch synthesis (:math:`\mu` mol C g-1 MS)
        self.DELTA_DSTARCH = 0.0001          #: Relative rate of starch degradation (s-1)

        # Fructans
        self.VMAX_SFRUCTAN_POT = 0.015       #: Potential maximal rate of fructan synthesis (:math:`\mu` mol C s-1 g-1 MS)
        self.K_SFRUCTAN = 5000               #: Affinity coefficient of fructan synthesis (:math:`\mu` mol C g-1 MS)
        self.K_REGUL_SFRUCTAN = 0.001        #: Affinity coefficient of the regulation function of fructan synthesis (:math:`\mu` mol g-1 MS)
        self.N_REGUL_SFRUCTAN = 3            #: Parameter of the regulation function of fructan synthesis (dimensionless)
        self.VMAX_DFRUCTAN = 0.035           #: Maximal rate of fructan degradation (:math:`\mu` mol C s-1 g-1 MS)
        self.K_DFRUCTAN = 100                #: Affinity coefficient of fructan degradation (:math:`\mu` mol C g-1 MS)

        # Loading sucrose and amino acids
        self.SIGMA_SUCROSE = 1e-08           #: Conductivity of an organ-phloem pathway (g2 :math:`\mu` mol-1 m-2 s-1) ; used to compute the sucrose loaded to the phloem
        self.SIGMA_AMINO_ACIDS = 1e-07       #: Conductivity of an organ-phloem pathway (g2 :math:`\mu` mol-1 m-2 s-1) ; used to compute the amino acids loaded to the phloem
        self.BETA = 1                        #: Kind of volumetric mass density at power -2/3 ((g m-3)**(-2/3))

        # Amino acids
        self.VMAX_AMINO_ACIDS = 1            #: Maximal rate of amino acid synthesis (:math:`\mu` mol N s-1 g-1 MS)
        self.K_AMINO_ACIDS_NITRATES = 3      #: Affinity coefficient of amino acid synthesis from nitrates (:math:`\mu` mol N g-1 MS)
        self.K_AMINO_ACIDS_TRIOSESP = 0.2    #: Affinity coefficient of amino acid synthesis from triosesP (:math:`\mu` mol C g-1 MS)

        # Proteins
        self.VMAX_SPROTEINS = 0.0015         #: Maximal rate of protein synthesis (:math:`\mu` mol N s-1 g-1 MS)
        self.K_SPROTEINS = 100               #: Affinity coefficient of protein synthesis (:math:`\mu` mol N g-1 MS)
        self.VMAX_DPROTEINS = 2.5E-6         #: Maximal rate of protein degradation (:math:`\mu` mol g-1 mstruct s-1)
        self.K_DPROTEINS = 50                #: Affinity coefficient with cytokinins for protein degradation (UA g-1 mstruct)
        self.N_DPROTEINS = 2.1               #: A coefficient for the regulation of protein degradation by cytokines (dimensionless)

        # cytokinins
        self.DELTA_D_CYTOKININS = 3.E-6      #: Relative rate of cytokinins degradation (s-1)


#: The instance of class :class:`cnwheat.parameters.PhotosyntheticOrganParameters` for current process
PHOTOSYNTHETIC_ORGAN_PARAMETERS = PhotosyntheticOrganParameters()


class ChaffParameters(PhotosyntheticOrganParameters):
    """
    Internal parameters of chaffs.
    """
    def __init__(self):
        pass


#: The instance of class :class:`cnwheat.parameters.ChaffParameters` for current process
CHAFF_PARAMETERS = ChaffParameters()


class LaminaParameters(PhotosyntheticOrganParameters):
    """
    Internal parameters of laminae.
    """
    def __init__(self):
        pass


#: The instance of class :class:`cnwheat.parameters.LaminaParameters` for current process
LAMINA_PARAMETERS = LaminaParameters()


class InternodeParameters(PhotosyntheticOrganParameters):
    """
    Internal parameters of internodes.
    """
    def __init__(self):
        pass


#: The instance of class :class:`cnwheat.parameters.InternodeParameters` for current process
INTERNODE_PARAMETERS = InternodeParameters()


class PeduncleParameters(PhotosyntheticOrganParameters):
    """
    Internal parameters of peduncles.
    """
    def __init__(self):
        pass


#: The instance of class :class:`cnwheat.parameters.PeduncleParameters` for current process
PEDUNCLE_PARAMETERS = PeduncleParameters()


class SheathParameters(PhotosyntheticOrganParameters):
    """
    Internal parameters of sheaths.
    """
    def __init__(self):
        pass


#: The instance of class :class:`cnwheat.parameters.SheathParameters` for current process
SHEATH_PARAMETERS = SheathParameters()


class PhotosyntheticOrganElementParameters(object):
    """
    Internal parameters of photosynthetic organs elements.
    """
    def __init__(self):
        pass


#: The instance of class :class:`cnwheat.parameters.PhotosyntheticOrganElementParameters` for current process
PHOTOSYNTHETIC_ORGAN_ELEMENT_PARAMETERS = PhotosyntheticOrganElementParameters()


class PhotosyntheticOrganElementInitCompartments(object):
    """
    Initial values for compartments of photosynthetic organ elements.
    """
    def __init__(self):
        self.green_area = 1E-4   #: initial value of green_area (m2)
        self.mstruct = 5E-3      #: initial value of mstruct (g)
        self.Nstruct = 1E-3      #: initial value of Nstruct (g)

        self.triosesP = 0        #: initial value of triosesP (:math:`\mu` mol C)
        self.starch = 0          #: initial value of starch (:math:`\mu` mol C)
        self.sucrose = 0         #: initial value of sucrose (:math:`\mu` mol C)
        self.fructan = 0         #: initial value of fructan (:math:`\mu` mol C)
        self.nitrates = 0        #: initial value of nitrates (:math:`\mu` mol N)
        self.amino_acids = 0     #: initial value of amino_acids (:math:`\mu` mol N)
        self.proteins = 0        #: initial value of proteins (:math:`\mu` mol N)
        self.cytokinins = 0     #: initial value of cytokinins (AU)

        self.is_growing = False  #: initial value of is_growing (Flag indicating if the element is growing or not (:class:`bool`)
        self.Tr = 0              #: initial value of Tr (Transpiration rate (mmol m-2 s-1)
        self.Ag = 0              #: initial value of Ag (Gross assimilation (:math:`\mu` mol m-2 s-1)
        self.Ts = 15             #: initial value of Ts (Organ temperature (degree Celsius)


#: The instance of class :class:`cnwheat.parameters.PhotosyntheticOrganElementInitCompartments` for current process
PHOTOSYNTHETIC_ORGAN_ELEMENT_INIT_COMPARTMENTS = PhotosyntheticOrganElementInitCompartments()


class ChaffElementParameters(PhotosyntheticOrganElementParameters):
    """
    Internal parameters of chaffs elements.
    """
    def __init__(self):
        self.ALPHA = 1  #: Proportion of structural mass containing substrate


#: The instance of class :class:`cnwheat.parameters.ChaffElementParameters` for current process
CHAFF_ELEMENT_PARAMETERS = ChaffElementParameters()


class LaminaElementParameters(PhotosyntheticOrganElementParameters):
    """
    Internal parameters of laminae elements.
    """
    def __init__(self):
        self.ALPHA = 1  #: Proportion of structural mass containing substrate


#: The instance of class :class:`cnwheat.parameters.LaminaElementParameters` for current process
LAMINA_ELEMENT_PARAMETERS = LaminaElementParameters()


class InternodeElementParameters(PhotosyntheticOrganElementParameters):
    """
    Internal parameters of internodes elements.
    """
    def __init__(self):
        self.ALPHA = 1  #: Proportion of structural mass containing substrate


#: The instance of class :class:`cnwheat.parameters.InternodeElementParameters` for current process
INTERNODE_ELEMENT_PARAMETERS = InternodeElementParameters()


class PeduncleElementParameters(PhotosyntheticOrganElementParameters):
    """
    Internal parameters of peduncles elements.
    """
    def __init__(self):
        self.ALPHA = 1  #: Proportion of structural mass containing substrate


#: The instance of class :class:`cnwheat.parameters.PeduncleElementParameters` for current process
PEDUNCLE_ELEMENT_PARAMETERS = PeduncleElementParameters()


class SheathElementParameters(PhotosyntheticOrganElementParameters):
    """
    Internal parameters of sheaths elements.
    """
    def __init__(self):
        self.ALPHA = 1  #: Proportion of structural mass containing substrate


#: The instance of class :class:`cnwheat.parameters.SheathElementParameters` for current process
SHEATH_ELEMENT_PARAMETERS = SheathElementParameters()


class SoilParameters(object):
    """
    Internal parameters of soil.
    """
    def __init__(self):
        self.MINERALISATION_RATE = 2.05E-6  # Mineralisation rate (:math:`\mu` mol N nitrates m-3 s-1)


#: The instance of class :class:`cnwheat.parameters.SoilParameters` for current process
SOIL_PARAMETERS = SoilParameters()
