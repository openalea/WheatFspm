# -*- coding: latin-1 -*-

from __future__ import division  # use "//" to do integer division
import numpy as np

import parameters


"""
    cnwheat.model
    ~~~~~~~~~~~~~

    The module :mod:`cnwheat.model` defines the equations of the CN exchanges in a population of plants.

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


class EcophysiologicalConstants:
    """
    Ecophysiological constants.
    """
    C_MOLAR_MASS = 12                       #: Molar mass of carbon (g mol-1)
    NB_C_TRIOSEP = 3                        #: Number of C in 1 mol of trioseP
    NB_C_HEXOSES = 6                        #: Number of C in 1 mol of hexoses (glucose, fructose)
    NB_C_SUCROSE = 12                       #: Number of C in 1 mol of sucrose
    HEXOSE_MOLAR_MASS_C_RATIO = 0.4         #: Contribution of C in hexose mass
    RATIO_C_mstruct = 0.384                 #: Mean contribution of carbon to structural dry mass (g C g-1 mstruct)

    AMINO_ACIDS_C_RATIO = 3.67              #: Mean number of mol of C in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)
    AMINO_ACIDS_N_RATIO = 1.17              #: Mean number of mol of N in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)
    AMINO_ACIDS_MOLAR_MASS_N_RATIO = 0.145  #: Mean contribution of N in amino acids mass
    N_MOLAR_MASS = 14                       #: Molar mass of nitrogen (g mol-1)


class Population(object):
    """
    The class :class:`Population` defines the CN exchanges at population scale.

    A :class:`population <Population>` must have at least one :class:`plant <Plant>`.
    """

    PARAMETERS = parameters.POPULATION_PARAMETERS  #: the internal parameters of the population

    def __init__(self, plants=None):
        if plants is None:
            plants = []
        self.plants = plants  #: the list of plants

    def calculate_aggregated_variables(self):
        """Calculate the integrative variables of the population recursively.
        """
        for plant in self.plants:
            plant.calculate_aggregated_variables()


class Plant(object):
    """
    The class :class:`Plant` defines the CN exchanges at plant scale.

    A :class:`plant <Plant>` must have at least one :class:`axis <Axis>`.
    """

    PARAMETERS = parameters.PLANT_PARAMETERS  #: the internal parameters of the plants

    def __init__(self, index=None, axes=None):
        self.index = index  #: the index of the plant
        if axes is None:
            axes = []
        self.axes = axes  #: the list of axes

    def calculate_aggregated_variables(self):
        """Calculate the integrative variables of the plant recursively.
        """
        for axis in self.axes:
            axis.calculate_aggregated_variables()


class Axis(object):
    """
    The class :class:`Axis` defines the CN exchanges at axis scale.

    An :class:`axis <Axis>` must have:
        * one :class:`set of roots <Roots>`,
        * one :class:`phloem <Phloem>`,
        * zero or one :class:`set of grains <Grains>`,
        * at least one :class:`phytomer<Phytomer>`.
    """

    PARAMETERS = parameters.AXIS_PARAMETERS  #: the internal parameters of the axes

    def __init__(self, label=None, roots=None, phloem=None, grains=None, phytomers=None):
        self.label = label  #: the label of the axis
        self.roots = roots  #: the roots
        self.phloem = phloem  #: the phloem
        self.grains = grains  #: the grains
        if phytomers is None:
            phytomers = []
        self.phytomers = phytomers  #: the list of phytomers
        self.mstruct = None  #: structural mass of the axis (g)
        # integrative variables
        self.Total_Transpiration = None  #: the total transpiration (mmol s-1)

    def calculate_aggregated_variables(self):
        """Calculate the integrative variables of the axis recursively.
        """
        self.mstruct = 0
        if self.roots is not None:
            self.roots.calculate_aggregated_variables()
            self.mstruct += self.roots.mstruct
        if self.phloem is not None:
            self.phloem.calculate_aggregated_variables()
        if self.grains is not None:
            self.grains.calculate_aggregated_variables()
            self.mstruct += self.grains.structural_dry_mass
        for phytomer in self.phytomers:
            phytomer.calculate_aggregated_variables()
            self.mstruct += phytomer.mstruct


class Phytomer(object):
    """
    The class :class:`Phytomer` defines the CN exchanges at phytomer scale.

    A :class:`phytomer <Phytomer>` must have at least:
        * 1 photosynthetic organ: :class:`chaff <Chaff>`, :class:`peduncle <Peduncle>`,
                                  :class:`lamina <Lamina>`, :class:`internode <Internode>`,
                                  or :class:`sheath <Sheath>`.
        * or 1 :class:`hiddenzone <HiddenZone>`.
    """

    PARAMETERS = parameters.PHYTOMER_PARAMETERS  #: the internal parameters of the phytomers

    def __init__(self, index=None, chaff=None, peduncle=None, lamina=None, internode=None, sheath=None, hiddenzone=None):
        self.index = index  #: the index of the phytomer
        self.chaff = chaff  #: the chaff
        self.peduncle = peduncle  #: the peduncle
        self.lamina = lamina  #: the lamina
        self.internode = internode  #: the internode
        self.sheath = sheath  #: the sheath
        self.hiddenzone = hiddenzone  #: the hidden zone
        self.mstruct = None  #: the structural mass of the phytomer (g)

    def calculate_aggregated_variables(self):
        """Calculate the integrative variables of the phytomer recursively.
        """
        self.mstruct = 0
        for organ_ in (self.chaff, self.peduncle, self.lamina, self.internode, self.sheath, self.hiddenzone):
            if organ_ is not None:
                organ_.calculate_aggregated_variables()
                self.mstruct += organ_.mstruct


class Organ(object):
    """
    The class :class:`Organ` defines the CN exchanges at organ scale.

    :class:`Organ` is the base class of all organs. DO NOT INSTANTIATE IT.
    """

    def __init__(self, label):
        self.label = label  #: the label of the organ

    def initialize(self):
        """Initialize the derived attributes of the organ.
        """
        pass

    def calculate_aggregated_variables(self):
        """Calculate the integrative variables of the organ recursively.
        """
        pass


class HiddenZone(Organ):
    """
    The class :class:`HiddenZone` defines the CN exchanges in an hidden zone.
    """

    PARAMETERS = parameters.HIDDEN_ZONE_PARAMETERS                #: the internal parameters of the hidden zone
    INIT_COMPARTMENTS = parameters.HIDDEN_ZONE_INIT_COMPARTMENTS  #: the initial values of compartments and state parameters

    def __init__(self, label='hiddenzone', mstruct=INIT_COMPARTMENTS.mstruct, Nstruct=INIT_COMPARTMENTS.Nstruct,
                 sucrose=INIT_COMPARTMENTS.sucrose, fructan=INIT_COMPARTMENTS.fructan, amino_acids=INIT_COMPARTMENTS.amino_acids, proteins=INIT_COMPARTMENTS.proteins):

        super(HiddenZone, self).__init__(label)

        # state parameters
        self.mstruct = mstruct                      #: g
        self.Nstruct = Nstruct                      #: g

        # state variables
        self.sucrose = sucrose                      #: :math:`\mu mol` C
        self.fructan = fructan                      #: :math:`\mu mol` C
        self.amino_acids = amino_acids              #: :math:`\mu mol` N
        self.proteins = proteins                    #: :math:`\mu mol` N

        # fluxes from phloem
        self.Unloading_Sucrose = None               #: current Unloading of sucrose from phloem to hiddenzone integrated over delta t (:math:`\mu mol` C)
        self.Unloading_Amino_Acids = None           #: current Unloading of amino acids from phloem to hiddenzone integrated over delta t (:math:`\mu mol` N)

        # other fluxes
        self.S_Proteins = None                      #: protein synthesis (:math:`\mu mol` N g-1 mstruct)
        self.S_Fructan = None                       #: fructan synthesis (:math:`\mu mol` C g-1 mstruct)
        self.D_Fructan = None                       #: fructan degradation (:math:`\mu mol` C g-1 mstruct)
        self.D_Proteins = None                      #: protein degradation (:math:`\mu mol` N g-1 mstruct)

    # FLUXES

    def calculate_Unloading_Sucrose(self, sucrose, sucrose_phloem, mstruct_axis):
        """Rate of sucrose Unloading from phloem to the hidden zone (:math:`\mu mol` C sucrose unloaded g-1 mstruct h-1).
        Transport-resistance equation

        :Parameters:
            - `sucrose` (:class:`float`) - Sucrose amount in the hidden zone (:math:`\mu mol` C)
            - `sucrose_phloem` (:class:`float`) - Sucrose amount in phloem (:math:`\mu mol` C)
            - `mstruct_axis` (:class:`float`) -The structural dry mass of the axis (g)
        :Returns:
            Rate of Sucrose Unloading (:math:`\mu mol` C h-1)
        :Returns Type:
            :class:`float`
        """
        conductance = parameters.HIDDEN_ZONE_PARAMETERS.SIGMA * parameters.PHOTOSYNTHETIC_ORGAN_PARAMETERS.BETA * self.mstruct**(2/3)  # TODO: choix valeurs param / flux phloem-hgz
        return ((sucrose_phloem / mstruct_axis) - (sucrose / self.mstruct)) * conductance * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    def calculate_Unloading_Amino_Acids(self, amino_acids, amino_acids_phloem, mstruct_axis):
        """Rate of amino acids Unloading from phloem to the hidden zone (:math:`\mu mol` N amino acids unloaded g-1 mstruct h-1).
        Transport-resistance equation

        :Parameters:
            - `amino_acids` (:class:`float`) - Amino_acids amount in the hidden zone (:math:`\mu mol` N)
            - `amino_acids_phloem` (:class:`float`) - Amino_acids amount in phloem (:math:`\mu mol` N)
            - `mstruct_axis` (:class:`float`) -The structural dry mass of the axis (g)
        :Returns:
            Rate of Amino_acids Unloading (:math:`\mu mol` N h-1)
        :Returns Type:
            :class:`float`
        """
        conductance = parameters.HIDDEN_ZONE_PARAMETERS.SIGMA * parameters.PHOTOSYNTHETIC_ORGAN_PARAMETERS.BETA * self.mstruct**(2/3)
        return ((amino_acids_phloem / mstruct_axis) - (amino_acids / self.mstruct)) * conductance * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    def calculate_S_proteins(self, amino_acids):
        """Rate of protein synthesis (:math:`\mu mol` N proteins h-1 g-1 MS).
        Michaelis-Menten function of amino acids.

        :Parameters:
            - `amino_acids` (:class:`float`) - Amino acid amount in the hidden zone (:math:`\mu mol` N).
        :Returns:
            Rate of Protein synthesis (:math:`\mu mol` N g-1 mstruct h-1)
        :Returns Type:
            :class:`float`
        """
        return 2*(parameters.PHOTOSYNTHETIC_ORGAN_PARAMETERS.VMAX_SPROTEINS * max(0, (amino_acids / self.mstruct))) / (parameters.PHOTOSYNTHETIC_ORGAN_PARAMETERS.K_SPROTEINS +
                                                                                                                     max(0, (amino_acids / self.mstruct))) * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    def calculate_D_Proteins(self, proteins):
        """Rate of protein degradation (:math:`\mu mol` N proteins h-1 g-1 MS).
        First order kinetic

        :Parameters:
            - `proteins` (:class:`float`) - Protein amount in the hidden zone (:math:`\mu mol` N).
        :Returns:
            Rate of Protein degradation (:math:`\mu mol` N g-1 mstruct h-1)
        :Returns Type:
            :class:`float`
        """
        return max(0, (parameters.HIDDEN_ZONE_PARAMETERS.delta_Dproteins * (proteins / self.mstruct))) * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    def calculate_Regul_S_Fructan(self, Loading_Sucrose):
        """Regulating function for fructan maximal rate of synthesis.
        Negative regulation by the loading of sucrose from the phloem ("swith-off" sigmoïdal kinetic).

        :Parameters:
            - `Loading_Sucrose` (:class:`float`) - Sucrose loading (:math:`\mu mol` C)
        :Returns:
            Maximal rate of fructan synthesis (:math:`\mu mol` C g-1 mstruct)
        :Returns Type:
        """
        if Loading_Sucrose <= 0:
            Vmax_Sfructans = PhotosyntheticOrgan.PARAMETERS.VMAX_SFRUCTAN_POT
        else:  # Regulation by sucrose loading
            rate_Loading_Sucrose_massic = Loading_Sucrose/self.mstruct/parameters.SECOND_TO_HOUR_RATE_CONVERSION
            Vmax_Sfructans = ((PhotosyntheticOrgan.PARAMETERS.VMAX_SFRUCTAN_POT * PhotosyntheticOrgan.PARAMETERS.K_REGUL_SFRUCTAN**PhotosyntheticOrgan.PARAMETERS.N_REGUL_SFRUCTAN) /
                              (max(0, rate_Loading_Sucrose_massic**PhotosyntheticOrgan.PARAMETERS.N_REGUL_SFRUCTAN) +
                               PhotosyntheticOrgan.PARAMETERS.K_REGUL_SFRUCTAN**PhotosyntheticOrgan.PARAMETERS.N_REGUL_SFRUCTAN))
        return Vmax_Sfructans

    def calculate_S_Fructan(self, sucrose, Regul_S_Fructan):
        """Rate of fructan synthesis (:math:`\mu mol` C fructan g-1 mstruct h-1).
        Sigmoïdal function of sucrose.

        :Parameters:
            - `sucrose` (:class:`float`) - Amount of sucrose (:math:`\mu mol` C)
            - `Regul_S_Fructan` (:class:`float`) - Maximal rate of fructan synthesis regulated by sucrose loading (:math:`\mu mol` C g-1 mstruct)
        :Returns:
            Rate of Fructan synthesis (:math:`\mu mol` C g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return ((max(0, sucrose)/self.mstruct) * Regul_S_Fructan)/((max(0, sucrose)/self.mstruct) + PhotosyntheticOrgan.PARAMETERS.K_SFRUCTAN) * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    def calculate_D_Fructan(self, sucrose, fructan):
        """Rate of fructan degradation (:math:`\mu mol` C fructan g-1 mstruct h-1).
        Inhibition function by the end product i.e. sucrose (Bancal et al., 2012).

        :Parameters:
            - `sucrose` (:class:`float`) - Amount of sucrose (:math:`\mu mol` C)
            - `fructan` (:class:`float`) - Amount of fructan (:math:`\mu mol` C)
        :Returns:
            Rate of Fructan degradation (:math:`\mu mol` C g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        d_potential = ((PhotosyntheticOrgan.PARAMETERS.K_DFRUCTAN * PhotosyntheticOrgan.PARAMETERS.VMAX_DFRUCTAN) /
                       ((max(0, sucrose) / self.mstruct) + PhotosyntheticOrgan.PARAMETERS.K_DFRUCTAN)) * parameters.SECOND_TO_HOUR_RATE_CONVERSION
        d_actual = min(d_potential, max(0, fructan))
        return d_actual

    # COMPARTMENTS

    def calculate_sucrose_derivative(self, Unloading_Sucrose, S_Fructan, D_Fructan, hiddenzone_Loading_Sucrose_contribution):
        """delta sucrose of hidden zone.

        :Parameters:
            - `Unloading_Sucrose` (:class:`float`) - Sucrose unloaded (:math:`\mu mol` C)
            - `S_Fructan` (:class:`float`) - Fructan synthesis (:math:`\mu mol` C g-1 mstruct)
            - `D_Fructan` (:class:`float`) - Fructan degradation (:math:`\mu mol` C g-1 mstruct)
            - `hiddenzone_Loading_Sucrose_contribution` (:class:`float`) - Sucrose imported from the emerged tissues (:math:`\mu mol` C)
        :Returns:
            delta sucrose (:math:`\mu mol` C sucrose)
        :Returns Type:
            :class:`float`
        """
        return Unloading_Sucrose + (D_Fructan - S_Fructan) * self.mstruct + hiddenzone_Loading_Sucrose_contribution

    def calculate_amino_acids_derivative(self, Unloading_Amino_Acids, S_Proteins, D_Proteins, hiddenzone_Loading_Amino_Acids_contribution):
        """delta amino acids of hidden zone.

        :Parameters:
            - `Unloading_Amino_Acids` (:class:`float`) - Amino acids unloaded (:math:`\mu mol` N)
            - `S_Proteins` (:class:`float`) - Protein synthesis (:math:`\mu mol` N g-1 mstruct)
            - `D_Proteins` (:class:`float`) - Protein degradation (:math:`\mu mol` N g-1 mstruct)
            - `hiddenzone_Loading_Amino_Acids_contribution` (:class:`float`) - Amino acids imported from the emerged tissues (:math:`\mu mol` N)
        :Returns:
            delta amino acids (:math:`\mu mol` N amino acids)
        :Returns Type:
            :class:`float`
        """
        return Unloading_Amino_Acids + (D_Proteins - S_Proteins) * self.mstruct + hiddenzone_Loading_Amino_Acids_contribution

    def calculate_fructan_derivative(self, S_Fructan, D_Fructan):
        """delta fructans of hidden zone.

        :Parameters:
            - `S_Fructan` (:class:`float`) - Fructans synthesis (:math:`\mu mol` C g-1 mstruct)
            - `D_Fructan` (:class:`float`) - Fructans degradation (:math:`\mu mol` C g-1 mstruct)
        :Returns:
            delta fructans (:math:`\mu mol` C fructans)
        :Returns Type:
            :class:`float`
        """
        return (S_Fructan - D_Fructan) * self.mstruct

    def calculate_proteins_derivative(self, S_Proteins, D_Proteins):
        """delta proteins of hidden zone.

        :Parameters:
            - `S_Proteins` (:class:`float`) - Protein synthesis (:math:`\mu mol` N g-1 mstruct)
            - `D_Proteins` (:class:`float`) - Protein degradation (:math:`\mu mol` N g-1 mstruct)
        :Returns:
            delta proteins (:math:`\mu mol` N proteins)
        :Returns Type:
            :class:`float`
        """
        return (S_Proteins - D_Proteins) * self.mstruct


class Phloem(Organ):
    """
    The class :class:`Phloem` defines the CN exchanges in a phloem.
    """

    PARAMETERS = parameters.PHLOEM_PARAMETERS                #: the internal parameters of the phloem
    INIT_COMPARTMENTS = parameters.PHLOEM_INIT_COMPARTMENTS  #: the initial values of compartments and state parameters

    def __init__(self, label='phloem', sucrose=INIT_COMPARTMENTS.sucrose, amino_acids=INIT_COMPARTMENTS.amino_acids):

        super(Phloem, self).__init__(label)

        # state variables
        self.sucrose = sucrose          #: :math:`\mu mol` C sucrose
        self.amino_acids = amino_acids  #: :math:`\mu mol` N amino acids

    # COMPARTMENTS

    def calculate_sucrose_derivative(self, contributors):
        """delta sucrose

        :Parameters:
            - `contributors` (:class:`object`) - Organs exchanging C with the phloem
        :Returns:
            delta sucrose (:math:`\mu mol` C sucrose)
        :Returns Type:
            :class:`float`
        """
        sucrose_derivative = 0
        for contributor in contributors:
            if isinstance(contributor, PhotosyntheticOrganElement):
                sucrose_derivative += contributor.Loading_Sucrose
            elif isinstance(contributor, Grains):
                sucrose_derivative -= contributor.S_grain_structure + (contributor.S_grain_starch * contributor.structural_dry_mass)
            elif isinstance(contributor, Roots):
                sucrose_derivative -= contributor.Unloading_Sucrose * contributor.mstruct * contributor.__class__.PARAMETERS.ALPHA
            elif isinstance(contributor, HiddenZone):
                sucrose_derivative -= contributor.Unloading_Sucrose

        return sucrose_derivative

    def calculate_amino_acids_derivative(self, contributors):
        """delta amino acids

        :Parameters:
            - `contributors` (:class:`object`) - Organs exchanging N with the phloem
        :Returns:
            delta amino acids (:math:`\mu mol` N amino acids)
        :Returns Type:
            :class:`float`
        """
        amino_acids_derivative = 0
        for contributor in contributors:
            if isinstance(contributor, PhotosyntheticOrganElement):
                amino_acids_derivative += contributor.Loading_Amino_Acids
            elif isinstance(contributor, Grains):
                amino_acids_derivative -= contributor.S_Proteins
            elif isinstance(contributor, Roots):
                amino_acids_derivative -= contributor.Unloading_Amino_Acids * contributor.mstruct * contributor.__class__.PARAMETERS.ALPHA
            elif isinstance(contributor, HiddenZone):
                amino_acids_derivative -= contributor.Unloading_Amino_Acids

        return amino_acids_derivative


class Grains(Organ):
    """
    The class :class:`Grains` defines the CN exchanges in a set of grains.
    """

    AMINO_ACIDS_MOLAR_MASS_N_RATIO = 0.136                   #: Mean contribution of N in amino acids mass contained in gluten (Glu, Gln and Pro)

    PARAMETERS = parameters.GRAINS_PARAMETERS                #: the internal parameters of the grains
    INIT_COMPARTMENTS = parameters.GRAINS_INIT_COMPARTMENTS  #: the initial values of compartments and state parameters

    def __init__(self, label='grains', age_from_flowering=INIT_COMPARTMENTS.age_from_flowering, starch=INIT_COMPARTMENTS.starch, structure=INIT_COMPARTMENTS.structure,
                 proteins=INIT_COMPARTMENTS.proteins):

        super(Grains, self).__init__(label)

        # state variables
        self.age_from_flowering = age_from_flowering    #: seconds
        self.starch = starch                            #: :math:`\mu mol` of C starch
        self.structure = structure                      #: :math:`\mu mol` of C sucrose
        self.proteins = proteins                        #: :math:`\mu mol` of N proteins

        # derived attributes
        self.structural_dry_mass = None                 #: g of MS

        # fluxes from phloem
        self.S_grain_structure = None                   #: current synthesis of grain structure integrated over a delta t (:math:`\mu mol` C)
        self.S_grain_starch = None                      #: current synthesis of grain starch integrated over a delta t (:math:`\mu mol` C g-1 mstruct)
        self.S_Proteins = None                          #: current synthesis of grain proteins integrated over a delta t (:math:`\mu mol` N)

        # intermediate variables
        self.RGR_Structure = None                       #: RGR of grain structure (dimensionless?)
        self.R_grain_growth_struct = None               #: grain struct  respiration (:math:`\mu mol` C respired)
        self.R_grain_growth_starch = None               #: grain starch growth respiration (:math:`\mu mol` C respired)

    def initialize(self):
        """Initialize the derived attributes of the organ.
        """
        self.structural_dry_mass = self.calculate_structural_dry_mass(self.structure)

    # VARIABLES

    def calculate_structural_dry_mass(self, structure):
        """Grain structural dry mass.

        :Parameters:
            - `structure` (:class:`float`) - Grain structural C mass (:math:`\mu mol` C)
        :Returns:
            Grain structural dry mass (g)
        :Returns Type:
            :class:`float`
        """
        return (structure*1E-6*EcophysiologicalConstants.C_MOLAR_MASS) / EcophysiologicalConstants.RATIO_C_mstruct

    def calculate_RGR_Structure(self, sucrose_phloem, mstruct_axis):
        """Relative Growth Rate of grain structure, regulated by sucrose concentration in phloem.

        :Parameters:
            - `sucrose_phloem` (:class:`float`) - Sucrose amount in phloem (:math:`\mu mol` C)
            - `mstruct_axis` (:class:`float`) -The structural dry mass of the axis (g)
        :Returns:
            RGR of grain structure (dimensionless?)
        :Returns Type:
            :class:`float`
        """
        return ((max(0, sucrose_phloem) / (mstruct_axis * Axis.PARAMETERS.ALPHA)) * Grains.PARAMETERS.VMAX_RGR) / ((max(0, sucrose_phloem)/(mstruct_axis * Axis.PARAMETERS.ALPHA))
                                                                                                                   + Grains.PARAMETERS.K_RGR)

    # FLUXES

    def calculate_S_grain_structure(self, prec_structure, RGR_Structure):
        """Rate of grain structure synthesis (:math:`\mu mol` C structure h-1).
        Exponential function, RGR regulated by sucrose concentration in the phloem.

        :Parameters:
            - `prec_structure` (:class:`float`) - Grain structure at t-1 (:math:`\mu mol` C)
            - `RGR_Structure` (:class:`float`) - Relative Growth Rate of grain structure (dimensionless?)
        :Returns:
            Rate of Synthesis of grain structure (:math:`\mu mol` C h-1)
        :Returns Type:
            :class:`float`
        """
        if self.age_from_flowering <= Grains.PARAMETERS.FILLING_INIT:  #: Grain enlargment
            S_grain_structure = prec_structure * RGR_Structure * parameters.SECOND_TO_HOUR_RATE_CONVERSION
        else:                                                          #: Grain filling
            S_grain_structure = 0
        return S_grain_structure

    def calculate_S_grain_starch(self, sucrose_phloem, mstruct_axis):
        """Rate of starch synthesis in grains (i.e. grain filling) (:math:`\mu mol` C starch g-1 mstruct h-1).
        Michaelis-Menten function of sucrose concentration in the phloem.

        :Parameters:
            - `sucrose_phloem` (:class:`float`) - Sucrose concentration in phloem (:math:`\mu mol` C g-1 mstruct)
            - `mstruct_axis` (:class:`float`) -The structural dry mass of the axis (g)
        :Returns:
            Rate of Synthesis of grain starch (:math:`\mu mol` C g-1 mstruct h-1)
        :Returns Type:
            :class:`float`
        """
        if self.age_from_flowering <= Grains.PARAMETERS.FILLING_INIT:   #: Grain enlargment
            S_grain_starch = 0
        elif self.age_from_flowering > Grains.PARAMETERS.FILLING_END:   #: Grain maturity
            S_grain_starch = 0
        else:                                                           #: Grain filling
            S_grain_starch = (((max(0, sucrose_phloem)/(mstruct_axis * Axis.PARAMETERS.ALPHA)) * Grains.PARAMETERS.VMAX_STARCH) /
                              ((max(0, sucrose_phloem)/(mstruct_axis * Axis.PARAMETERS.ALPHA)) + Grains.PARAMETERS.K_STARCH)) * parameters.SECOND_TO_HOUR_RATE_CONVERSION
        return S_grain_starch

    def calculate_S_proteins(self, S_grain_structure, S_grain_starch, amino_acids_phloem, sucrose_phloem, structural_dry_mass):
        """Protein synthesis in grains.
        N is assumed to be co-transported along with the unloaded sucrose from phloem (using the ratio amino acids:sucrose of phloem).

        :Parameters:
            - `S_grain_structure` (:class:`float`) - Synthesis of grain structure (:math:`\mu mol` C)
            - `S_grain_starch` (:class:`float`) - Synthesis of grain starch (:math:`\mu mol` C g-1 mstruct)
            - `amino_acids_phloem` (:class:`float`) - Amino acids concentration in phloem (:math:`\mu mol` N g-1 mstruct)
            - `sucrose_phloem` (:class:`float`) - Sucrose concentration in phloem (:math:`\mu mol` C g-1 mstruct)
            - `structural_dry_mass` (:class:`float`) - Grain structural dry mass (g)
        :Returns:
            Synthesis of grain proteins (:math:`\mu mol` N)
        :Returns Type:
            :class:`float`
        """
        if sucrose_phloem > 0:
            S_Proteins = (S_grain_structure + S_grain_starch*structural_dry_mass) * (amino_acids_phloem / sucrose_phloem)
        else:
            S_Proteins = 0
        return S_Proteins

    # COMPARTMENTS

    def calculate_structure_derivative(self, S_grain_structure, R_growth):
        """delta grain structure.

        :Parameters:
            - `S_grain_structure` (:class:`float`) - Synthesis of grain structure (:math:`\mu mol` C)
            - `R_growth` (:class:`float`) - Grain growth respiration (:math:`\mu mol` C respired)
        :Returns:
            delta grain structure (:math:`\mu mol` C structure)
        :Returns Type:
            :class:`float`
        """
        return S_grain_structure - R_growth

    def calculate_starch_derivative(self, S_grain_starch, structural_dry_mass, R_growth):
        """delta grain starch.

        :Parameters:
            - `S_grain_starch` (:class:`float`) - Synthesis of grain starch (:math:`\mu mol` C g-1 mstruct)
            - `structural_dry_mass` (:class:`float`) - Grain structural dry mass (g)
            - `R_growth` (:class:`float`) - Grain growth respiration (:math:`\mu mol` C respired)
        :Returns:
            delta grain starch (:math:`\mu mol` C starch)
        :Returns Type:
            :class:`float`
        """
        return (S_grain_starch * structural_dry_mass) - R_growth

    def calculate_proteins_derivative(self, S_Proteins):
        """delta grain proteins.

        :Parameters:
            - `S_Proteins` (:class:`float`) - Synthesis of grain proteins (:math:`\mu mol` N)
        :Returns:
            delta grain proteins (:math:`\mu mol` N proteins)
        :Returns Type:
            :class:`float`
        """
        return S_Proteins


class Roots(Organ):
    """
    The class :class:`Roots` defines the CN exchanges in a set of roots.
    """

    PARAMETERS = parameters.ROOTS_PARAMETERS                #: the internal parameters of the roots
    INIT_COMPARTMENTS = parameters.ROOTS_INIT_COMPARTMENTS  #: the initial values of compartments and state parameters

    def __init__(self, label='roots', mstruct=INIT_COMPARTMENTS.mstruct, Nstruct=INIT_COMPARTMENTS.Nstruct,
                                      sucrose=INIT_COMPARTMENTS.sucrose, nitrates=INIT_COMPARTMENTS.nitrates, amino_acids=INIT_COMPARTMENTS.amino_acids,
                                      cytokinins=INIT_COMPARTMENTS.cytokinins):

        super(Roots, self).__init__(label)

        # state parameters
        self.mstruct = mstruct                 #: Structural mass (g)
        self.Nstruct = Nstruct                 #: Structural N mass (g)

        # state variables
        self.sucrose = sucrose                 #: :math:`\mu mol` C sucrose
        self.nitrates = nitrates               #: :math:`\mu mol` N nitrates
        self.amino_acids = amino_acids         #: :math:`\mu mol` N amino acids
        self.cytokinins = cytokinins           #: AU cytokinins

        # fluxes from phloem
        self.Unloading_Sucrose = None          #: current Unloading of sucrose from phloem to roots
        self.Unloading_Amino_Acids = None      #: current Unloading of amino acids from phloem to roots

        # other fluxes
        self.Export_Nitrates = None            #: Total export of nitrates from roots to shoot organs integrated over a delta t (:math:`\mu mol` N)
        self.Export_Amino_Acids = None         #: Total export of amino acids from roots to shoot organs integrated over a delta t (:math:`\mu mol` N)
        self.S_Amino_Acids = None              #: Rate of amino acid synthesis in roots integrated over a delta t (:math:`\mu mol` N g-1 mstruct)
        self.Uptake_Nitrates = None            #: Rate of nitrate uptake by roots integrated over a delta t (:math:`\mu mol` N nitrates)
        self.S_cytokinins = None               #: Rate of cytokinin synthesis integrated over a delta t (AU g-1 mstruct)
        self.Export_cytokinins = None          #: Total export of cytokinin from roots to shoot organs integrated over a delta t (AU)

        # Integrated variables
        self.Total_Organic_Nitrogen = None     #: current amount of organic N (:math:`\mu mol` N)

        # intermediate variables
        self.R_Nnit_upt = None                 #: Nitrate uptake respiration (:math:`\mu mol` C respired)
        self.R_Nnit_red = None                 #: Nitrate reduction-linked respiration (:math:`\mu mol` C respired)
        self.R_residual = None                 #: Residual maintenance respiration (cost from protein turn-over, cell ion gradients, futile cycles...) (:math:`\mu mol` C respired)
        self.C_exudation = None                #: C sucrose lost by root exudation integrated over a delta t (:math:`\mu mol` C g-1 mstruct)
        self.N_exudation = None                #: N amino acids lost by root exudation integrated over a delta t (:math:`\mu mol` N g-1 mstruct)
        self.regul_transpiration = None        #: Dimensionless regulating factor of metabolite exports from roots by shoot transpiration
        self.HATS_LATS = None                  #: Nitrate influx (:math:`\mu mol` N)
        self.sum_respi = None                  #: Sum of respirations for roots i.e. related to N uptake, amino acids synthesis and residual (:math:`\mu mol` C)

    def calculate_aggregated_variables(self):
        self.Total_Organic_Nitrogen = self.calculate_Total_Organic_Nitrogen(self.amino_acids, self.Nstruct)

    # VARIABLES

    def calculate_Total_Organic_Nitrogen(self, amino_acids, Nstruct):
        """Total amount of organic N (amino acids + Nstruct).
        Used to calculate residual respiration.

        :Parameters:
            - `amino_acids` (:class:`float`) - Amount of amino acids (:math:`\mu mol` N)
            - `Nstruct` (:class:`float`) -Structural N mass (g)
        :Returns:
            Total amount of organic N (:math:`\mu mol` N)
        :Returns Type:
            :class:`float`
        """
        return amino_acids + (Nstruct / EcophysiologicalConstants.N_MOLAR_MASS)*1E6

    def calculate_regul_transpiration(self, total_surfacic_transpiration):
        """A function to regulate metabolite exports from roots by shoot transpiration

        :Parameters:
            - `total_surfacic_transpiration` (:class:`float`) - Surfacic rate of total transpiration (mmol m-2 s-1)
        :Returns:
            Dimensionless regulating factor
        :Returns Type:
            :class:`float`
        """
        return total_surfacic_transpiration / (total_surfacic_transpiration + Roots.PARAMETERS.K_TRANSPIRATION)

    # FLUXES

    def calculate_Unloading_Sucrose(self, sucrose_phloem, mstruct_axis):
        """Rate of sucrose Unloading from phloem to roots (:math:`\mu mol` C sucrose unloaded g-1 mstruct h-1).
        Michaelis-Menten function of the sucrose concentration in phloem.

        :Parameters:
            - `sucrose_phloem` (:class:`float`) - Sucrose concentration in phloem (:math:`\mu mol` C g-1 mstruct)
            - `mstruct_axis` (:class:`float`) -The structural dry mass of the axis (g)
        :Returns:
            Rate of Sucrose Unloading (:math:`\mu mol` C g-1 mstruct h-1)
        :Returns Type:
            :class:`float`
        """
        return (((max(0, sucrose_phloem)/(mstruct_axis * Axis.PARAMETERS.ALPHA)) * Roots.PARAMETERS.VMAX_SUCROSE_UNLOADING) /
                ((max(0, sucrose_phloem)/(mstruct_axis * Axis.PARAMETERS.ALPHA)) + Roots.PARAMETERS.K_SUCROSE_UNLOADING)) * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    def calculate_Unloading_Amino_Acids(self, Unloading_Sucrose, sucrose_phloem, amino_acids_phloem):
        """Unloading of amino_acids from phloem to roots.
        Amino acids are assumed to be co-transported along with the unloaded sucrose from phloem (using the ratio amino acids:sucrose of phloem).

        :Parameters:
            - `Unloading_Sucrose` (:class:`float`) - Sucrose Unloading (:math:`\mu mol` C g-1 mstruct)
            - `sucrose_phloem` (:class:`float`) - Sucrose concentration in phloem (:math:`\mu mol` C g-1 mstruct)
            - `amino_acids_phloem` (:class:`float`) - Amino acids concentration in phloem (:math:`\mu mol` N g-1 mstruct)
        :Returns:
            Amino acids Unloading (:math:`\mu mol` N g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        if amino_acids_phloem <= 0:
            Unloading_Amino_Acids = 0
        else:
            Unloading_Amino_Acids = Unloading_Sucrose * (amino_acids_phloem / sucrose_phloem)
        return Unloading_Amino_Acids

    def calculate_Uptake_Nitrates(self, Conc_Nitrates_Soil, nitrates_roots, sucrose_roots):
        """Rate of nitrate uptake by roots
            - Nitrate uptake is calculated as the sum of the 2 transport systems: HATS and LATS
            - HATS and LATS parameters are calculated as a function of root nitrate concentration (negative regulation)
            - Nitrate uptake is finally regulated by the total culm transpiration and sucrose concentration (positive regulation)

        :Parameters:
            - `Conc_Nitrates_Soil` (:class:`float`) - Soil nitrate concentration Unloading (:math:`\mu mol` N m-3 soil)
            - `nitrates_roots` (:class:`float`) - Amount of nitrates in roots (:math:`\mu mol` N)
            - `sucrose_roots` (:class:`float`) - Amount of sucrose in roots (:math:`\mu mol` C)
        :Returns:
            Nitrate uptake (:math:`\mu mol` N nitrates) and nitrate influxes HATS and LATS (:math:`\mu mol` N h-1)
        :Returns Type:
            :class:`tuple` of 2 :class:`float`
        """
        conc_nitrates_roots = nitrates_roots / self.mstruct

        #: High Affinity Transport System (HATS)
        VMAX_HATS_MAX = Roots.PARAMETERS.A_VMAX_HATS * np.exp(-Roots.PARAMETERS.LAMBDA_VMAX_HATS * conc_nitrates_roots)  #: Maximal rate of nitrates influx at saturating soil N concentration;HATS (:math:`\mu mol` N nitrates g-1 mstruct s-1)
        K_HATS = Roots.PARAMETERS.A_K_HATS * np.exp(-Roots.PARAMETERS.LAMBDA_K_HATS * conc_nitrates_roots)               #: Affinity coefficient of nitrates influx at saturating soil N concentration;HATS (:math:`\mu mol` m-3)
        HATS = (VMAX_HATS_MAX * Conc_Nitrates_Soil)/ (K_HATS + Conc_Nitrates_Soil)                                       #: Rate of nitrate influx by HATS (:math:`\mu mol` N nitrates uptaked s-1 g-1 mstruct)

        #: Low Affinity Transport System (LATS)
        K_LATS = Roots.PARAMETERS.A_LATS * np.exp(-Roots.PARAMETERS.LAMBDA_LATS * conc_nitrates_roots)                   #: Rate constant for nitrates influx at low soil N concentration; LATS (m3 g-1 mstruct s-1)
        LATS = (K_LATS * Conc_Nitrates_Soil)                                                                             #: Rate of nitrate influx by LATS (:math:`\mu mol` N nitrates uptaked s-1 g-1 mstruct)

        #: Nitrate influx (:math:`\mu mol` N)
        HATS_LATS = (HATS + LATS) * self.mstruct * parameters.SECOND_TO_HOUR_RATE_CONVERSION

        # Regulations
        regul_C = (sucrose_roots/self.mstruct) / ((sucrose_roots/self.mstruct) + Roots.PARAMETERS.K_C)                   #: Nitrate uptake regulation by root C
        net_nitrate_uptake = HATS_LATS * Roots.PARAMETERS.NET_INFLUX_UPTAKE_RATIO * regul_C                              #: Net nitrate uptake (:math:`\mu mol` N nitrates uptaked by roots)
        return net_nitrate_uptake, HATS_LATS

    def calculate_S_amino_acids(self, nitrates, sucrose):
        """Rate of amino acid synthesis in roots (:math:`\mu mol` N amino acids g-1 mstruct h-1).
        Bi-substrate Michaelis-Menten function of nitrates and sucrose.

        :Parameters:
            - `nitrates` (:class:`float`) - Amount of nitrates in roots (:math:`\mu mol` N)
            - `sucrose` (:class:`float`) - Amount of sucrose in roots (:math:`\mu mol` C)
        :Returns:
            Amino acids synthesis (:math:`\mu mol` N g-1 mstruct h-1)
        :Returns Type:
            :class:`float`
        """
        return Roots.PARAMETERS.VMAX_AMINO_ACIDS / ((1 + Roots.PARAMETERS.K_AMINO_ACIDS_NITRATES/(nitrates/(self.mstruct*Roots.PARAMETERS.ALPHA))) *
                                                    (1 + Roots.PARAMETERS.K_AMINO_ACIDS_SUCROSE/(sucrose/(self.mstruct*Roots.PARAMETERS.ALPHA)))) * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    def calculate_Export_Nitrates(self, nitrates, regul_transpiration):
        """Total export of nitrates from roots to shoot organs
        Export is calculated as a function on nitrate concentration and culm transpiration.

        :Parameters:
            - `nitrates` (:class:`float`) - Amount of nitrates in roots (:math:`\mu mol` N)
            - `regul_transpiration` (:class:`float`) - Regulating factor by transpiration (mmol H2O m-2 s-1)
        :Returns:
            Rate of Export of nitrates (:math:`\mu mol` N h-1)
        :Returns Type:
            :class:`float`
        """
        if nitrates <= 0 or regul_transpiration <= 0:
            Export_Nitrates = 0
        else:
            f_nitrates = (nitrates / (self.mstruct * Roots.PARAMETERS.ALPHA)) * Roots.PARAMETERS.K_NITRATE_EXPORT           #: :math:`\mu mol` g-1 s-1
            Export_Nitrates = f_nitrates * self.mstruct * regul_transpiration * parameters.SECOND_TO_HOUR_RATE_CONVERSION   #: Nitrate export regulation by transpiration (:math:`\mu mol` N)
        return Export_Nitrates

    def calculate_Export_Amino_Acids(self, amino_acids, regul_transpiration):
        """Total export of amino acids from roots to shoot organs
        Amino acids export is calculated as a function of nitrate export using the ratio amino acids:nitrates in roots.

        :Parameters:
            - `amino_acids` (:class:`float`) - Amount of amino acids in roots (:math:`\mu mol` N)
            - `regul_transpiration` (:class:`float`) - Regulating factor by transpiration (mmol H2O m-2 s-1)

        :Returns:
            Rate of Export of amino acids (:math:`\mu mol` N h-1)
        :Returns Type:
            :class:`float`
        """
        if amino_acids <= 0:
            Export_Amino_Acids = 0
        else:
            f_amino_acids = (amino_acids/(self.mstruct * Roots.PARAMETERS.ALPHA)) * Roots.PARAMETERS.K_AMINO_ACIDS_EXPORT
            Export_Amino_Acids = f_amino_acids * self.mstruct * regul_transpiration * parameters.SECOND_TO_HOUR_RATE_CONVERSION  #: Amino acids export regulation by plant transpiration (:math:`\mu mol` N)

        return Export_Amino_Acids

    def calculate_exudation(self, Unloading_Sucrose, sucrose_roots, amino_acids_roots, amino_acids_phloem):
        """C sucrose and N amino acids lost by root exudation (:math:`\mu mol` C or N g-1 mstruct).
            - C exudation is calculated as a fraction of C Unloading from phloem
            - N exudation is calculated from C exudation using the ratio amino acids:sucrose of the phloem

        :Parameters:
            - `Unloading_Sucrose` (:class:`float`) - Sucrose Unloading (:math:`\mu mol` C g-1 mstruct h-1)
            - `sucrose_roots` (:class:`float`) - Amount of sucrose in roots (:math:`\mu mol` C)
            - `amino_acids_roots` (:class:`float`) - Amount of amino acids in roots (:math:`\mu mol` N)
            - `amino_acids_phloem` (:class:`float`) - Amount of amino acids in phloem (:math:`\mu mol` N)
        :Returns:
            Rates of C exudated (:math:`\mu mol` C g-1 mstruct h-1) and N_exudation (:math:`\mu mol` N g-1 mstruct h-1)
        :Returns Type:
            :class:`tuple` of 2 :class:`float`
        """
        C_exudation = min(sucrose_roots, Unloading_Sucrose * Roots.PARAMETERS.C_EXUDATION)  #: C exudated (:math:`\mu mol` g-1 mstruct)
        if amino_acids_phloem <= 0 or amino_acids_roots <= 0 or sucrose_roots <= 0:
            N_exudation = 0
        else:
            N_exudation = (amino_acids_roots/sucrose_roots) * C_exudation
        return C_exudation, N_exudation

    def calculate_S_cytokinins(self, sucrose_roots, nitrates_roots):
        """ Rate of cytokinin synthesis (AU cytokinins g-1 mstruct h-1).
        Cytokinin synthesis regulated by both root sucrose and nitrates. As a signal molecule, cytokinins are assumed have a neglected effect on sucrose.
        Thus, no cost in C is applied to the sucrose pool.

        :Parameters:
            - `sucrose_roots` (:class:`float`) - Amount of sucrose in roots (:math:`\mu mol` C)
            - `nitrates_roots` (:class:`float`) - Amount of nitrates in roots (:math:`\mu mol` N)
        :Returns:
            Rate of Cytokinin synthesis (AU g-1 mstruct h-1)
        :Returns Type:
            :class:`float`
        """
        conc_sucrose = max(0, (sucrose_roots/self.mstruct))
        conc_Nitrates = max(0, (nitrates_roots/self.mstruct))
        f_sucrose = conc_sucrose**Roots.PARAMETERS.N_SUC_CYTOKININS/(conc_sucrose**Roots.PARAMETERS.N_SUC_CYTOKININS + Roots.PARAMETERS.K_SUCROSE_CYTOKININS**Roots.PARAMETERS.N_SUC_CYTOKININS)
        f_nitrates = conc_Nitrates**Roots.PARAMETERS.N_NIT_CYTOKININS/(conc_Nitrates**Roots.PARAMETERS.N_NIT_CYTOKININS + Roots.PARAMETERS.K_NITRATES_CYTOKININS**Roots.PARAMETERS.N_NIT_CYTOKININS)
        S_cytokinins = Roots.PARAMETERS.VMAX_S_CYTOKININS * f_sucrose * f_nitrates * parameters.SECOND_TO_HOUR_RATE_CONVERSION
        return S_cytokinins

    def calculate_Export_cytokinins(self, cytokinins, regul_transpiration):
        """Total export of cytokinin from roots to shoot organs
        Cytokinin export is calculated as a function of cytokinin concentration and culm transpiration.

        :Parameters:
            - `cytokinins` (:class:`float`) - Amount of cytokinins in roots (AU)
            - `regul_transpiration` (:class:`float`) - Regulating factor by transpiration (mmol H2O m-2 s-1)
        :Returns:
            Rate of Cytokinin export (AU h-1)
        :Returns Type:
            :class:`float`
        """
        if cytokinins <= 0:
            Export_cytokinins = 0
        else:
            f_cytokinins = (cytokinins / (self.mstruct*Roots.PARAMETERS.ALPHA)) * Roots.PARAMETERS.K_CYTOKININS_EXPORT
            Export_cytokinins = f_cytokinins * self.mstruct * regul_transpiration * parameters.SECOND_TO_HOUR_RATE_CONVERSION  #: Cytokinin export regulation by plant transpiration (AU)

        return Export_cytokinins

    # COMPARTMENTS

    def calculate_sucrose_derivative(self, Unloading_Sucrose, S_Amino_Acids, C_exudation, sum_respi):
        """delta root sucrose.

        :Parameters:
            - `Unloading_Sucrose` (:class:`float`) - Sucrose Unloading (:math:`\mu mol` C g-1 mstruct)
            - `S_Amino_Acids` (:class:`float`) - Amino acids synthesis (:math:`\mu mol` N g-1 mstruct)
            - `C_exudation` (:class:`float`) - C exudation (:math:`\mu mol` C g-1 mstruct)
            - `sum_respi` (:class:`float`) - Sum of respirations for roots i.e. related to N uptake, amino acids synthesis and residual (:math:`\mu mol` C)
        :Returns:
            delta root sucrose (:math:`\mu mol` C sucrose)
        :Returns Type:
            :class:`float`
        """
        sucrose_consumption_AA = (S_Amino_Acids / EcophysiologicalConstants.AMINO_ACIDS_N_RATIO) * EcophysiologicalConstants.AMINO_ACIDS_C_RATIO  #: Contribution of sucrose to the synthesis of amino_acids
        return (Unloading_Sucrose - sucrose_consumption_AA - C_exudation) * self.mstruct - sum_respi

    def calculate_nitrates_derivative(self, Uptake_Nitrates, Export_Nitrates, S_Amino_Acids):
        """delta root nitrates.

        :Parameters:
            - `Uptake_Nitrates` (:class:`float`) - Nitrate uptake (:math:`\mu mol` N nitrates)
            - `Export_Nitrates` (:class:`float`) - Export of nitrates (:math:`\mu mol` N)
            - `S_Amino_Acids` (:class:`float`) - Amino acids synthesis (:math:`\mu mol` N g-1 mstruct)
        :Returns:
            delta root nitrates (:math:`\mu mol` N nitrates)
        :Returns Type:
            :class:`float`
        """
        import_nitrates_roots = Uptake_Nitrates
        nitrate_reduction_AA = S_Amino_Acids                                                #: Contribution of nitrates to the synthesis of amino_acids
        return import_nitrates_roots - Export_Nitrates - nitrate_reduction_AA*self.mstruct

    def calculate_amino_acids_derivative(self, Unloading_Amino_Acids, S_Amino_Acids, Export_Amino_Acids, N_exudation):
        """delta root amino acids.

        :Parameters:
            - `Unloading_Amino_Acids` (:class:`float`) - Amino acids Unloading (:math:`\mu mol` N g-1 mstruct)
            - `S_Amino_Acids` (:class:`float`) - Amino acids synthesis (:math:`\mu mol` N g-1 mstruct)
            - `Export_Amino_Acids` (:class:`float`) - Export of amino acids (:math:`\mu mol` N)
            - `N_exudation` (:class:`float`) - N exudated (:math:`\mu mol` g-1 mstruct)
        :Returns:
            delta root amino acids (:math:`\mu mol` N amino acids)
        :Returns Type:
            :class:`float`
        """
        return (Unloading_Amino_Acids + S_Amino_Acids - N_exudation)*self.mstruct - Export_Amino_Acids

    def calculate_cytokinins_derivative(self, S_cytokinins, Export_cytokinins):
        """delta root cytokinins.

        :Parameters:
            - `S_cytokinins` (:class:`float`) - Cytokinin synthesis (AU g-1 mstruct)
            - `Export_cytokinins` (:class:`float`) - Cytokinin export (AU)
        :Returns:
            delta root cytokinins (AU cytokinins)
        :Returns Type:
            :class:`float`
        """
        return S_cytokinins * self.mstruct - Export_cytokinins


class PhotosyntheticOrgan(Organ):
    """
    The class :class:`PhotosyntheticOrgan` defines the CN exchanges in a photosynthetic organ.

    A :class:`photosynthetic organ <PhotosyntheticOrgan>` must have at least 1
    :class:`photosynthetic organ element <PhotosyntheticOrganElement>`:
    :class:`chaff element <ChaffElement>`, :class:`lamina element <LaminaElement>`,
    :class:`internode element <InternodeElement>`, :class:`peduncle element <PeduncleElement>`,
    or :class:`sheath element <SheathElement>`.

    :class:`PhotosyntheticOrgan` is the base class of all photosynthetic organs. DO NOT INSTANTIATE IT.
    """

    PARAMETERS = parameters.PHOTOSYNTHETIC_ORGAN_PARAMETERS  #: the internal parameters of the photosynthetic organs

    def __init__(self, label, exposed_element, enclosed_element):

        super(PhotosyntheticOrgan, self).__init__(label)

        self.exposed_element = exposed_element    #: the exposed element
        self.enclosed_element = enclosed_element  #: the enclosed element
        self.mstruct = None                       #: the structural dry mass

    def calculate_aggregated_variables(self):
        self.mstruct = 0
        for element in (self.exposed_element, self.enclosed_element):
            if element is not None:
                element.calculate_aggregated_variables()
                self.mstruct += element.mstruct


class Chaff(PhotosyntheticOrgan):
    """
    The class :class:`Chaff` defines the CN exchanges in a chaff.
    """

    PARAMETERS = parameters.CHAFF_PARAMETERS  #: the internal parameters of the chaffs

    def __init__(self, label=None, exposed_element=None, enclosed_element=None):
        super(Chaff, self).__init__(label, exposed_element, enclosed_element)


class Lamina(PhotosyntheticOrgan):
    """
    The class :class:`Lamina` defines the CN exchanges in a lamina.
    """

    PARAMETERS = parameters.LAMINA_PARAMETERS  #: the internal parameters of the laminae

    def __init__(self, label=None, exposed_element=None, enclosed_element=None):
        super(Lamina, self).__init__(label, exposed_element, enclosed_element)


class Internode(PhotosyntheticOrgan):
    """
    The class :class:`Internode` defines the CN exchanges in an internode.
    """

    PARAMETERS = parameters.INTERNODE_PARAMETERS  #: the internal parameters of the internodes

    def __init__(self, label=None, exposed_element=None, enclosed_element=None):
        super(Internode, self).__init__(label, exposed_element, enclosed_element)


class Peduncle(PhotosyntheticOrgan):
    """
    The class :class:`Peduncle` defines the CN exchanges in a peduncle.
    """

    PARAMETERS = parameters.PEDUNCLE_PARAMETERS  #: the internal parameters of the peduncles

    def __init__(self, label=None, exposed_element=None, enclosed_element=None):
        super(Peduncle, self).__init__(label, exposed_element, enclosed_element)


class Sheath(PhotosyntheticOrgan):
    """
    The class :class:`Sheath` defines the CN exchanges in a sheath.
    """

    PARAMETERS = parameters.SHEATH_PARAMETERS  #: the internal parameters of the sheaths

    def __init__(self, label=None, exposed_element=None, enclosed_element=None):
        super(Sheath, self).__init__(label, exposed_element, enclosed_element)


class PhotosyntheticOrganElement(object):
    """
    The class :class:`PhotosyntheticOrganElement` defines the CN exchanges in a photosynthetic organ element.

    An element must belong to an organ of the same type (e.g. a class:`LaminaElement` must belong to a class:`Lamina`).

    :class:`PhotosyntheticOrganElement` is the base class of all photosynthetic organs elements. DO NOT INSTANTIATE IT.
    """

    PARAMETERS = parameters.PHOTOSYNTHETIC_ORGAN_ELEMENT_PARAMETERS                #: the internal parameters of the photosynthetic organs elements
    INIT_COMPARTMENTS = parameters.PHOTOSYNTHETIC_ORGAN_ELEMENT_INIT_COMPARTMENTS  #: the initial values of compartments and state parameters

    def __init__(self, label=None, green_area=INIT_COMPARTMENTS.green_area, mstruct=INIT_COMPARTMENTS.mstruct, Nstruct=INIT_COMPARTMENTS.Nstruct,
                       triosesP=INIT_COMPARTMENTS.triosesP, starch=INIT_COMPARTMENTS.starch, sucrose=INIT_COMPARTMENTS.sucrose, fructan=INIT_COMPARTMENTS.fructan,
                       nitrates=INIT_COMPARTMENTS.nitrates, amino_acids=INIT_COMPARTMENTS.amino_acids, proteins=INIT_COMPARTMENTS.proteins, cytokinins=INIT_COMPARTMENTS.cytokinins,
                       Tr=INIT_COMPARTMENTS.Tr, Ag=INIT_COMPARTMENTS.Ag, Ts=INIT_COMPARTMENTS.Ts, is_growing=INIT_COMPARTMENTS.is_growing):

        self.label = label                   #: the label of the element

        # state parameters
        self.mstruct = mstruct               #: Structural dry mass (g)
        self.Nstruct = Nstruct               #: Structural N mass (g)
        self.is_growing = is_growing         #: Flag indicating if the element is growing or not (:class:`bool`)
        self.green_area = green_area         #: green area (m-2)
        self.Tr = Tr                         #: Transpiration rate (mmol m-2 s-1)
        self.Ag = Ag                         #: Gross assimilation (:math:`\mu mol` m-2 s-1)
        self.Ts = Ts                         #: Organ temperature (:math:`^{\circ}C`)

        # state variables
        self.triosesP = triosesP             #: :math:`\mu mol` C
        self.starch = starch                 #: :math:`\mu mol` C
        self.sucrose = sucrose               #: :math:`\mu mol` C
        self.fructan = fructan               #: :math:`\mu mol` C
        self.nitrates = nitrates             #: :math:`\mu mol` N
        self.amino_acids = amino_acids       #: :math:`\mu mol` N
        self.proteins = proteins             #: :math:`\mu mol` N
        self.cytokinins = cytokinins         #: AU

        # fluxes to phloem
        self.Loading_Sucrose = None          #: Rate of sucrose loading to phloem (:math:`\mu mol` C)
        self.Loading_Amino_Acids = None      #: Rate of amino acids loading to phloem (:math:`\mu mol` N)

        # other fluxes
        self.S_Proteins = None              #: Rate of protein synthesis (:math:`\mu mol` N g-1 mstruct)
        self.S_Amino_Acids = None           #: Rate of amino acids synthesis (:math:`\mu mol` N g-1 mstruct)
        self.Regul_S_Fructan = None         #: Maximal rate of fructan synthesis (:math:`\mu mol` C g-1 mstruct)
        self.S_Starch = None                #: Rate of starch synthesis (:math:`\mu mol` C g-1 mstruct)
        self.D_Starch = None                #: Rate of starch degradation (:math:`\mu mol` C g-1 mstruct)
        self.S_Sucrose = None               #: Rate of sucrose synthesis (:math:`\mu mol` C g-1 mstruct)
        self.S_Fructan = None               #: Rate of fructan synthesis (:math:`\mu mol` C g-1 mstruct)
        self.D_Fructan = None               #: Rate of fructan degradation ((:math:`\mu mol` C g-1 mstruct)
        self.Nitrates_import = None         #: Total nitrates imported from roots (:math:`\mu mol` N nitrates)
        self.Amino_Acids_import = None      #: Total amino acids imported from roots (:math:`\mu mol` N amino acids)
        self.k_proteins = None              #: First order kinetic regulated by cytokinins concentration
        self.D_Proteins = None              #: Rate of protein degradation (:math:`\mu mol` N g-1 mstruct)
        self.cytokinins_import = None       #: Import of cytokinins (AU)
        self.D_cytokinins = None            #: Rate of cytokinins degradation (AU g-1 mstruct)

        # Integrated variables
        self.Total_Organic_Nitrogen = None  #: current total nitrogen amount (:math:`\mu mol` N)

        # intermediate variables
        self.R_Nnit_red = None              #: Nitrate reduction-linked respiration (:math:`\mu mol` C respired)
        self.R_residual = None              #: Residual maintenance respiration (cost from protein turn-over, cell ion gradients, futile cycles...) (:math:`\mu mol` C respired)
        self.Transpiration = None           #: Surfacic transpiration rate of an element (mmol H2O s-1)
        self.R_phloem_loading = None        #: Phloem loading respiration (:math:`\mu mol` C respired)
        self.Photosynthesis = None          #: Total Photosynthesis of an element integrated over a delta t (:math:`\mu mol` C)
        self.sum_respi = None               #: Sum of respirations for the element i.e. related to C loading to phloem, amino acids synthesis and residual (:math:`\mu mol` C)

    def calculate_aggregated_variables(self):
        """Calculate the integrative variables of the element.
        """
        self.Total_Organic_Nitrogen = self.calculate_Total_Organic_Nitrogen(self.amino_acids, self.proteins, self.Nstruct)

    # VARIABLES

    def calculate_total_Photosynthesis(self, Ag, green_area):
        """Total Photosynthesis of an element (:math:`\mu mol` C m-2 h-1 * m2).

        :Parameters:
            - `Ag` (:class:`float`) - Gross Photosynthesis rate (:math:`\mu mol` C m-2 s-1)
            - `green_area` (:class:`float`) - Green area (m2)
        :Returns:
            Rate of Total Photosynthesis (:math:`\mu mol` C h-1)
        :Returns Type:
            :class:`float`
        """
        return Ag * green_area * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    def calculate_Total_Transpiration(self, Tr, green_area):
        """Surfacic transpiration rate of an element

        :Parameters:
            - `Tr` (:class:`float`) - Transpiration rate (mmol H2O m-2 s-1)
            - `green_area` (:class:`float`) - Green area (m2)
        :Returns:
            Total transpiration (mmol H2O s-1)
        :Returns Type:
            :class:`float`
        """
        return Tr * green_area

    def calculate_Regul_S_Fructan(self, Loading_Sucrose):
        """Regulating function for fructan maximal rate of synthesis.
        Negative regulation by the loading of sucrose from the phloem ("swith-off" sigmoïdal kinetic).

        :Parameters:
            - `Loading_Sucrose` (:class:`float`) - Sucrose loading (:math:`\mu mol` C)
        :Returns:
            Maximal rate of fructan synthesis (:math:`\mu mol` C g-1 mstruct h-1)
        :Returns Type:
            :class:`float`
        """
        if Loading_Sucrose <= 0:
            Vmax_Sfructans = PhotosyntheticOrgan.PARAMETERS.VMAX_SFRUCTAN_POT
        else:  # Regulation by sucrose loading
            rate_Loading_Sucrose_massic = Loading_Sucrose/self.mstruct/parameters.SECOND_TO_HOUR_RATE_CONVERSION
            Vmax_Sfructans = ((PhotosyntheticOrgan.PARAMETERS.VMAX_SFRUCTAN_POT * PhotosyntheticOrgan.PARAMETERS.K_REGUL_SFRUCTAN**PhotosyntheticOrgan.PARAMETERS.N_REGUL_SFRUCTAN) /
                              (max(0, rate_Loading_Sucrose_massic**PhotosyntheticOrgan.PARAMETERS.N_REGUL_SFRUCTAN) +
                               PhotosyntheticOrgan.PARAMETERS.K_REGUL_SFRUCTAN**PhotosyntheticOrgan.PARAMETERS.N_REGUL_SFRUCTAN))
        return Vmax_Sfructans

    def calculate_Total_Organic_Nitrogen(self, amino_acids, proteins, Nstruct):
        """Total amount of organic N (amino acids + proteins + Nstruct).
        Used to calculate residual respiration.

        :Parameters:
            - `amino_acids` (:class:`float`) - Amount of amino acids (:math:`\mu mol` N)
            - `proteins` (:class:`float`) - Amount of proteins (:math:`\mu mol` N)
            - `Nstruct` (:class:`float`) - Structural N mass (g)
        :Returns:
            Total amount of organic N (:math:`\mu mol` N)
        :Returns Type:
            :class:`float`
        """
        return amino_acids + proteins + (Nstruct / EcophysiologicalConstants.N_MOLAR_MASS)*1E6

    # FLUXES

    def calculate_S_Starch(self, triosesP):
        """Rate of starch synthesis (:math:`\mu mol` C starch g-1 mstruct h-1).
        Michaelis-Menten function of triose phosphates.

        :Parameters:
            - `triosesP` (:class:`float`) - Amount of triose phosphates (:math:`\mu mol` C)
        :Returns:
            Rate of Starch synthesis (:math:`\mu mol` C g-1 mstruct h-1)
        :Returns Type:
            :class:`float`
        """
        return (((max(0, triosesP)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) * PhotosyntheticOrgan.PARAMETERS.VMAX_STARCH) /
                ((max(0, triosesP)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) + PhotosyntheticOrgan.PARAMETERS.K_STARCH)) * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    def calculate_D_Starch(self, starch):
        """Rate of starch degradation (:math:`\mu mol` C starch g-1 mstruct h-1).
        First order kinetic.

        :Parameters:
            - `starch` (:class:`float`) - Amount of starch (:math:`\mu mol` C)
        :Returns:
            Starch degradation (:math:`\mu mol` C g-1 mstruct h-1)
        :Returns Type:
            :class:`float`
        """
        return max(0, PhotosyntheticOrgan.PARAMETERS.DELTA_DSTARCH * (starch/(self.mstruct*self.__class__.PARAMETERS.ALPHA))) * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    def calculate_S_Sucrose(self, triosesP):
        """Rate of sucrose synthesis (:math:`\mu mol` C sucrose g-1 mstruct h-1).
        Michaelis-Menten function of triose phosphates.

        :Parameters:
            - `triosesP` (:class:`float`) - Amount of triose phosphates (:math:`\mu mol` C)
        :Returns:
            Rate of Sucrose synthesis (:math:`\mu mol` C g-1 mstruct h-1)
        :Returns Type:
            :class:`float`
        """
        return (((max(0, triosesP)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) * PhotosyntheticOrgan.PARAMETERS.VMAX_SUCROSE) /
                ((max(0, triosesP)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) + PhotosyntheticOrgan.PARAMETERS.K_SUCROSE)) * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    def calculate_Loading_Sucrose(self, sucrose, sucrose_phloem, mstruct_axis):
        """Rate of sucrose loading to phloem (:math:`\mu mol` C sucrose h-1).
        Transport-resistance model.

        :Parameters:
            - `sucrose` (:class:`float`) - Amount of sucrose in the element (:math:`\mu mol` C)
            - `sucrose_phloem` (:class:`float`) - Amount of sucrose in the phloem (:math:`\mu mol` C)
            - `mstruct_axis` (:class:`float`) - Structural dry mass of the axis (g)
        :Returns:
            Rate of Sucrose loading (:math:`\mu mol` C h-1)
        :Returns Type:
            :class:`float`
        """
        conc_sucrose_element = sucrose / (self.mstruct*self.__class__.PARAMETERS.ALPHA)
        conc_sucrose_phloem = sucrose_phloem / (mstruct_axis * parameters.AXIS_PARAMETERS.ALPHA)
        #: Driving compartment (:math:`\mu mol` C g-1 mstruct)
        driving_sucrose_compartment = max(conc_sucrose_element, conc_sucrose_phloem)
        #: Gradient of sucrose between the element and the phloem (:math:`\mu mol` C g-1 mstruct)
        diff_sucrose = conc_sucrose_element - conc_sucrose_phloem
        #: Conductance depending on mstruct (g2 :math:`\mu mol`-1 s-1)
        conductance = PhotosyntheticOrgan.PARAMETERS.SIGMA_SUCROSE * PhotosyntheticOrgan.PARAMETERS.BETA * self.mstruct**(2/3)

        return driving_sucrose_compartment * diff_sucrose * conductance * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    def calculate_export_sucrose(self, sucrose, sucrose_hiddenzone, mstruct_hiddenzone):
        """Rate of sucrose exportation to hidden zone (:math:`\mu mol` C sucrose h-1).
        Transport-resistance model.

        :Parameters:
            - `sucrose` (:class:`float`) - Amount of sucrose in the element (:math:`\mu mol` C)
            - `sucrose_hiddenzone` (:class:`float`) - Sucrose amount in the hidden zone (:math:`\mu mol` C)
            - `mstruct_hiddenzone` (:class:`float`) - mstruct of the hidden zone (g)
        :Returns:
            Rate of Sucrose export (:math:`\mu mol` C h-1)
        :Returns Type:
            :class:`float`
        """
        conc_sucrose_element = sucrose / (self.mstruct*self.__class__.PARAMETERS.ALPHA)
        conc_sucrose_hiddenzone = sucrose_hiddenzone / mstruct_hiddenzone
        #: Gradient of sucrose between the element and the hidden zone (:math:`\mu mol` C g-1 mstruct)
        diff_sucrose = conc_sucrose_element - conc_sucrose_hiddenzone
        #: Conductance depending on mstruct
        conductance = HiddenZone.PARAMETERS.SIGMA * PhotosyntheticOrgan.PARAMETERS.BETA * mstruct_hiddenzone**(2/3)

        return diff_sucrose * conductance * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    def calculate_S_Fructan(self, sucrose, Regul_S_Fructan):
        """Rate of fructan synthesis (:math:`\mu mol` C fructan g-1 mstruct h-1).
        Sigmoïdal function of sucrose.

        :Parameters:
            - `sucrose` (:class:`float`) - Amount of sucrose (:math:`\mu mol` C)
            - `Regul_S_Fructan` (:class:`float`) - Maximal rate of fructan synthesis regulated by sucrose loading (:math:`\mu mol` C g-1 mstruct)
        :Returns:
            Rate of Fructan synthesis (:math:`\mu mol` C g-1 mstruct h-1)
        :Returns Type:
            :class:`float`
        """
        return ((max(0, sucrose)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) * Regul_S_Fructan) /\
               ((max(0, sucrose) / (self.mstruct*self.__class__.PARAMETERS.ALPHA)) + PhotosyntheticOrgan.PARAMETERS.K_SFRUCTAN) * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    def calculate_D_Fructan(self, sucrose, fructan):
        """Rate of fructan degradation (:math:`\mu mol` C fructan g-1 mstruct h-1).
        Inhibition function by the end product i.e. sucrose (Bancal et al., 2012).

        :Parameters:
            - `sucrose` (:class:`float`) - Amount of sucrose (:math:`\mu mol` C)
            - `fructan` (:class:`float`) - Amount of fructan (:math:`\mu mol` C)
        :Returns:
            Rate of Fructan degradation (:math:`\mu mol` C g-1 mstruct h-1)
        :Returns Type:
            :class:`float`
        """
        d_potential = ((PhotosyntheticOrgan.PARAMETERS.K_DFRUCTAN * PhotosyntheticOrgan.PARAMETERS.VMAX_DFRUCTAN) /
                       ((max(0, sucrose)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) + PhotosyntheticOrgan.PARAMETERS.K_DFRUCTAN)) * parameters.SECOND_TO_HOUR_RATE_CONVERSION
        d_actual = min(d_potential, max(0, fructan))
        return d_actual

    def calculate_Nitrates_import(self, Export_Nitrates, element_transpiration, Total_Transpiration):
        """Total nitrates imported from roots (:math:`\mu mol` N nitrates).
        Nitrates coming from roots (fraction of uptake + direct export) are distributed according to the contribution of the element to culm transpiration.

        :Parameters:
            - `Export_Nitrates` (:class:`float`) - Exported nitrates by roots (:math:`\mu mol` N)
            - `element_transpiration` (:class:`float`) - Element transpiration (mmol H2O s-1)
            - `Total_Transpiration` (:class:`float`) - Culm transpiration (mmol H2O s-1)
        :Returns:
            Total nitrates import (:math:`\mu mol` N nitrates)
        :Returns Type:
            :class:`float`
        """
        if Total_Transpiration > 0:
            Nitrates_import = Export_Nitrates * (element_transpiration/Total_Transpiration)      #: Proportion of exported nitrates from roots to element
        else:  # Avoids further float division by zero error
            Nitrates_import = 0
        return Nitrates_import

    def calculate_Amino_Acids_import(self, roots_exported_amino_acids, element_transpiration, Total_Transpiration):
        """Total amino acids imported from roots  (:math:`\mu mol` N amino acids).
        Amino acids exported by roots are distributed according to the contribution of the element to culm transpiration.

        :Parameters:
            - `roots_exported_amino_acids` (:class:`float`) - Exported amino acids by roots (:math:`\mu mol` N)
            - `element_transpiration` (:class:`float`) - Element transpiration (mmol H2O s-1)
            - `Total_Transpiration` (:class:`float`) - Culm transpiration (mmol H2O s-1)
        :Returns:
            Total amino acids import (:math:`\mu mol` N amino acids)
        :Returns Type:
            :class:`float`
        """
        if Total_Transpiration > 0:
            Amino_Acids_import = roots_exported_amino_acids * (element_transpiration/Total_Transpiration)  #: Proportion of exported amino acids from roots to organ
        else:
            Amino_Acids_import = 0
        return Amino_Acids_import

    def calculate_S_amino_acids(self, nitrates, triosesP):
        """Rate of amino acids synthesis (:math:`\mu mol` N amino acids h-1 g-1 MS).
        Bi-substrate Michaelis-Menten function of nitrates and triose phosphates.

        :Parameters:
            - `nitrates` (:class:`float`) - Amount of nitrates (:math:`\mu mol` N)
            - `triosesP` (:class:`float`) - Amount of triosesP (:math:`\mu mol` C)
        :Returns:
            Rate of Amino acids synthesis (:math:`\mu mol` N h-1 g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        if nitrates <= 0 or triosesP <= 0:
            calculate_S_amino_acids = 0
        else:
            calculate_S_amino_acids = PhotosyntheticOrgan.PARAMETERS.VMAX_AMINO_ACIDS /\
                                      ((1 + PhotosyntheticOrgan.PARAMETERS.K_AMINO_ACIDS_NITRATES / (nitrates/(self.mstruct*self.__class__.PARAMETERS.ALPHA))) *
                                       (1 + PhotosyntheticOrgan.PARAMETERS.K_AMINO_ACIDS_TRIOSESP/(triosesP/(self.mstruct*self.__class__.PARAMETERS.ALPHA)))) * parameters.SECOND_TO_HOUR_RATE_CONVERSION
        return calculate_S_amino_acids

    def calculate_S_proteins(self, amino_acids):
        """Rate of protein synthesis (:math:`\mu mol` N proteins h-1 g-1 MS).
        Michaelis-Menten function of amino acids.

        :Parameters:
            - `amino_acids` (:class:`float`) - Amount of amino acids (:math:`\mu mol` N).
        :Returns:
            Protein synthesis (:math:`\mu mol` N h-1 g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        calculate_S_proteins = (((max(0, amino_acids) / (self.mstruct*self.__class__.PARAMETERS.ALPHA)) * PhotosyntheticOrgan.PARAMETERS.VMAX_SPROTEINS) /
                                ((max(0, amino_acids) / (self.mstruct*self.__class__.PARAMETERS.ALPHA)) + PhotosyntheticOrgan.PARAMETERS.K_SPROTEINS)) * parameters.SECOND_TO_HOUR_RATE_CONVERSION
        return calculate_S_proteins

    def calculate_D_Proteins(self, proteins, cytokinins):
        """Rate of protein degradation (:math:`\mu mol` N proteins s-1 g-1 MS h-1).
        First order kinetic regulated by cytokinins concentration.

        :Parameters:
            - `proteins` (:class:`float`) - Amount of proteins (:math:`\mu mol` N)
            - `cytokinins` (:class:`float`) - Amount of cytokinins (AU)
        :Returns:
            First order kinetic regulated by cytokinins concentration,
            and rate of protein degradation (:math:`\mu mol` N g-1 mstruct)
        :Returns Type:
            :class:`tuple` of 2 :class:`float`
        """
        conc_cytokinins = max(0, cytokinins / self.mstruct)
        k_proteins = (PhotosyntheticOrgan.PARAMETERS.VMAX_DPROTEINS * PhotosyntheticOrgan.PARAMETERS.K_DPROTEINS**PhotosyntheticOrgan.PARAMETERS.N_DPROTEINS) /\
                     (conc_cytokinins**PhotosyntheticOrgan.PARAMETERS.N_DPROTEINS + PhotosyntheticOrgan.PARAMETERS.K_DPROTEINS**PhotosyntheticOrgan.PARAMETERS.N_DPROTEINS)

        return k_proteins, max(0, k_proteins * (proteins / (self.mstruct*self.__class__.PARAMETERS.ALPHA))) * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    def calculate_Loading_Amino_Acids(self, amino_acids, amino_acids_phloem, mstruct_axis):
        """Rate of amino acids loading to phloem (:math:`\mu mol` N amino acids h-1).
        Transport-resistance model.

        :Parameters:
            - `amino_acids` (:class:`float`) - Amount of amino acids in the element (:math:`\mu mol` N)
            - `amino_acids_phloem` (:class:`float`) - Amount of amino acids in the phloem (:math:`\mu mol` N)
            - `mstruct_axis` (:class:`float`) - Structural dry mass of the axis (g)
        :Returns:
            Amino acids loading (:math:`\mu mol` N h-1)
        :Returns Type:
            :class:`float`
        """
        Conc_Amino_Acids_element = amino_acids / (self.mstruct*self.__class__.PARAMETERS.ALPHA)
        Conc_Amino_Acids_phloem = amino_acids_phloem / (mstruct_axis * parameters.AXIS_PARAMETERS.ALPHA)
        #: Driving compartment (:math:`\mu mol` N g-1 mstruct)
        driving_amino_acids_compartment = max(Conc_Amino_Acids_element, Conc_Amino_Acids_phloem)
        #: Gradient of amino acids between the element and the phloem (:math:`\mu mol` N g-1 mstruct)
        diff_amino_acids = Conc_Amino_Acids_element - Conc_Amino_Acids_phloem
        #: Conductance depending on mstruct (g2 :math:`\mu mol`-1 s-1)
        conductance = PhotosyntheticOrgan.PARAMETERS.SIGMA_AMINO_ACIDS * PhotosyntheticOrgan.PARAMETERS.BETA * self.mstruct**(2/3)

        return driving_amino_acids_compartment * diff_amino_acids * conductance * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    def calculate_Export_Amino_Acids(self, amino_acids, amino_acids_hiddenzone, mstruct_hiddenzone):
        """Rate of amino acids exportation to hidden zone (:math:`\mu mol` N amino acids h-1).
        Transport-resistance model.

        :Parameters:
            - `amino_acids` (:class:`float`) - Amount of amino acids in the element (:math:`\mu mol` N)
            - `amino_acids_hiddenzone` (:class:`float`) - Amino acids amount in the hidden zone (:math:`\mu mol` N)
            - `mstruct_hiddenzone` (:class:`float`) - mstruct of the hidden zone (g)
        :Returns:
            Rate of Amino acids export (:math:`\mu mol` N h-1)
        :Returns Type:
            :class:`float`
        """
        Conc_Amino_Acids_element = amino_acids / (self.mstruct*self.__class__.PARAMETERS.ALPHA)
        Conc_Amino_Acids_hiddenzone = amino_acids_hiddenzone / mstruct_hiddenzone
        #: Gradient of amino acids between the element and the hidden zone (:math:`\mu mol` N g-1 mstruct)
        diff_amino_acids = Conc_Amino_Acids_element - Conc_Amino_Acids_hiddenzone
        #: Conductance depending on mstruct
        conductance = HiddenZone.PARAMETERS.SIGMA * PhotosyntheticOrgan.PARAMETERS.BETA * mstruct_hiddenzone**(2/3)

        return diff_amino_acids * conductance * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    def calculate_cytokinins_import(self, roots_exporteD_cytokinins, element_transpiration, Total_Transpiration):
        """Import of cytokinins (AU).
        Cytokinin exported by roots are distributed according to the contribution of the element to culm transpiration.

        :Parameters:
            - `roots_exporteD_cytokinins` (:class:`float`) - Exported cytokinins from roots (AU)
            - `element_transpiration` (:class:`float`) - Element transpiration (mmol H2O s-1)
            - `Total_Transpiration` (:class:`float`) - Culm transpiration (mmol H2O s-1)
        :Returns:
            Cytokinin import (AU)
        :Returns Type:
            :class:`float`
        """
        if Total_Transpiration > 0:
            cytokinins_import = roots_exporteD_cytokinins * (element_transpiration / Total_Transpiration)
        else:
            cytokinins_import = 0
        return cytokinins_import

    def calculate_D_cytokinins(self, cytokinins):
        """Rate of cytokinins degradation (AU g-1 mstruct h-1).
        First order kinetic.

        :Parameters:
            - `cytokinins` (:class:`float`) - Amount of cytokinins (AU)
        :Returns:
            Rate of Cytokinin degradation (AU h-1 g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return max(0, PhotosyntheticOrgan.PARAMETERS.DELTA_D_CYTOKININS * (cytokinins/(self.mstruct*self.__class__.PARAMETERS.ALPHA))) * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    # COMPARTMENTS

    def calculate_triosesP_derivative(self, Photosynthesis, S_Sucrose, S_Starch, S_Amino_Acids):
        """ delta triose phosphates of element.

        :Parameters:
            - `Photosynthesis` (:class:`float`) - Total gross Photosynthesis (:math:`\mu mol` C)
            - `S_Sucrose` (:class:`float`) - Sucrose synthesis (:math:`\mu mol` C g-1 mstruct)
            - `S_Starch` (:class:`float`) - Starch synthesis (:math:`\mu mol` C g-1 mstruct)
            - `S_Amino_Acids` (:class:`float`) - Amino acids synthesis (:math:`\mu mol` N g-1 mstruct)
        :Returns:
            delta triose phosphates (:math:`\mu mol` C triose phosphates)
        :Returns Type:
            :class:`float`
        """
        #: Contribution of triosesP to the synthesis of amino_acids
        triosesP_consumption_AA = (S_Amino_Acids / EcophysiologicalConstants.AMINO_ACIDS_N_RATIO) * EcophysiologicalConstants.AMINO_ACIDS_C_RATIO
        return Photosynthesis - (S_Sucrose + S_Starch + triosesP_consumption_AA) * (self.mstruct*self.__class__.PARAMETERS.ALPHA)

    def calculate_starch_derivative(self, S_Starch, D_Starch):
        """delta starch of element.

        :Parameters:
            - `S_Starch` (:class:`float`) - Starch synthesis (:math:`\mu mol` C g-1 mstruct)
            - `D_Starch` (:class:`float`) - Starch degradation (:math:`\mu mol` C g-1 mstruct)
        :Returns:
            delta starch (:math:`\mu mol` C starch)
        :Returns Type:
            :class:`float`
        """
        return (S_Starch - D_Starch) * (self.mstruct*self.__class__.PARAMETERS.ALPHA)

    def calculate_sucrose_derivative(self, S_Sucrose, D_Starch, Loading_Sucrose, S_Fructan, D_Fructan, sum_respi):
        """delta sucrose of element.

        :Parameters:
            - `S_Sucrose` (:class:`float`) - Sucrose synthesis (:math:`\mu mol` C g-1 mstruct)
            - `D_Starch` (:class:`float`) - Starch degradation (:math:`\mu mol` C g-1 mstruct)
            - `Loading_Sucrose` (:class:`float`) - Sucrose loading (:math:`\mu mol` C)
            - `S_Fructan` (:class:`float`) - Fructan synthesis (:math:`\mu mol` C g-1 mstruct)
            - `D_Fructan` (:class:`float`) - Fructan degradation (:math:`\mu mol` C g-1 mstruct)
            - `sum_respi` (:class:`float`) - Sum of respirations for the element i.e. related to C loading to phloem, amino acids synthesis and residual (:math:`\mu mol` C)
        :Returns:
            delta sucrose (:math:`\mu mol` C sucrose)
        :Returns Type:
            :class:`float`
        """
        return (S_Sucrose + D_Starch + D_Fructan - S_Fructan) * self.mstruct - sum_respi - Loading_Sucrose

    def calculate_fructan_derivative(self, S_Fructan, D_Fructan):
        """delta fructan of element.

        :Parameters:
            - `S_Fructan` (:class:`float`) - Fructan synthesis (:math:`\mu mol` C g-1 mstruct)
            - `D_Fructan` (:class:`float`) - Fructan degradation (:math:`\mu mol` C g-1 mstruct)
        :Returns:
            delta fructan (:math:`\mu mol` C fructan)
        :Returns Type:
            :class:`float`
        """
        return (S_Fructan - D_Fructan) * (self.mstruct*self.__class__.PARAMETERS.ALPHA)

    def calculate_nitrates_derivative(self, Nitrates_import, S_Amino_Acids):
        """delta nitrates of element.

        :Parameters:
            - `Nitrates_import` (:class:`float`) - Nitrate import from roots (:math:`\mu mol` N)
            - `S_Amino_Acids` (:class:`float`) - Amino acids synthesis (:math:`\mu mol` N g-1 mstruct)
        :Returns:
            delta nitrates (:math:`\mu mol` N nitrates)
        :Returns Type:
            :class:`float`
        """
        nitrate_reduction_AA = S_Amino_Acids  #: Contribution of nitrates to the synthesis of amino_acids
        return Nitrates_import - (nitrate_reduction_AA * self.mstruct*self.__class__.PARAMETERS.ALPHA)

    def calculate_amino_acids_derivative(self, Amino_Acids_import, S_Amino_Acids, S_Proteins, D_Proteins, Loading_Amino_Acids):
        """delta amino acids of element.

        :Parameters:
            - `Amino_Acids_import` (:class:`float`) - Amino acids import from roots (:math:`\mu mol` N)
            - `S_Amino_Acids` (:class:`float`) - Amino acids synthesis (:math:`\mu mol` N g-1 mstruct)
            - `S_Proteins` (:class:`float`) - Protein synthesis (:math:`\mu mol` N g-1 mstruct)
            - `D_Proteins` (:class:`float`) - Protein degradation (:math:`\mu mol` N g-1 mstruct)
            - `Loading_Amino_Acids` (:class:`float`) - Amino acids loading (:math:`\mu mol` N)
        :Returns:
            delta amino acids (:math:`\mu mol` N amino acids)
        :Returns Type:
            :class:`float`
        """
        return Amino_Acids_import - Loading_Amino_Acids + (S_Amino_Acids + D_Proteins - S_Proteins) * (self.mstruct * self.__class__.PARAMETERS.ALPHA)

    def calculate_proteins_derivative(self, S_Proteins, D_Proteins):
        """delta proteins of element.

        :Parameters:
            - `S_Proteins` (:class:`float`) - Protein synthesis (:math:`\mu mol` N g-1 mstruct)
            - `D_Proteins` (:class:`float`) - Protein degradation (:math:`\mu mol` N g-1 mstruct)
        :Returns:
            delta proteins (:math:`\mu mol` N proteins)
        :Returns Type:
            :class:`float`
        """
        return (S_Proteins - D_Proteins) * (self.mstruct*self.__class__.PARAMETERS.ALPHA)

    def calculate_cytokinins_derivative(self, import_cytokinins, D_cytokinins):
        """delta cytokinins of element.

        :Parameters:
            - `import_cytokinins` (:class:`float`) - Cytokinin import (AU)
            - `D_cytokinins` (:class:`float`) - Cytokinin degradation (AU g-1 mstruct)
        :Returns:
            delta cytokinins (AU cytokinins)
        :Returns Type:
            :class:`float`
        """
        return import_cytokinins - D_cytokinins * (self.mstruct*self.__class__.PARAMETERS.ALPHA)


class ChaffElement(PhotosyntheticOrganElement):
    """
    The class :class:`ChaffElement` defines the CN exchanges in a chaff element.
    """

    PARAMETERS = parameters.CHAFF_ELEMENT_PARAMETERS  #: the internal parameters of the chaffs elements


class LaminaElement(PhotosyntheticOrganElement):
    """
    The class :class:`LaminaElement` defines the CN exchanges in a lamina element.
    """

    PARAMETERS = parameters.LAMINA_ELEMENT_PARAMETERS  #: the internal parameters of the laminae elements


class InternodeElement(PhotosyntheticOrganElement):
    """
    The class :class:`InternodeElement` defines the CN exchanges in an internode element.
    """

    PARAMETERS = parameters.INTERNODE_ELEMENT_PARAMETERS  #: the internal parameters of the internodes elements


class PeduncleElement(PhotosyntheticOrganElement):
    """
    The class :class:`PeduncleElement` defines the CN exchanges in a peduncle element.
    """

    PARAMETERS = parameters.PEDUNCLE_ELEMENT_PARAMETERS  #: the internal parameters of the peduncles elements


class SheathElement(PhotosyntheticOrganElement):
    """
    The class :class:`SheathElement` defines the CN exchanges in a sheath element.
    """

    PARAMETERS = parameters.SHEATH_ELEMENT_PARAMETERS  #: the internal parameters of the sheaths elements


class Soil(object):
    """
    The class :class:`Soil` defines the amount of nitrogen in the volume of soil explored by roots.
    """

    PARAMETERS = parameters.SOIL_PARAMETERS  #: the internal parameters of the soil

    def __init__(self, volume=None, nitrates=None, Tsoil=None):

        # state parameters
        self.volume = volume                   #: volume of soil explored by roots (m3)
        self.Tsoil = Tsoil                     #: soil temperature (:math:`^{\circ}C`)

        # state variables
        self.nitrates = nitrates               #: :math:`\mu mol` N nitrates

        # intermediate variables
        self.Conc_Nitrates_Soil = None         #: soil nitrate concentration Unloading (:math:`\mu mol` N m-3 soil)
        self.mineralisation = None             #: mineralisation on organic N into nitrates in soil (:math:`\mu mol`)

    # VARIABLES

    def calculate_Conc_Nitrates(self, nitrates):
        """Nitrate concentration in soil.

        :Parameters:
            - `nitrates` (:class:`float`) - Amount of nitrates (:math:`\mu mol` N)
        :Returns:
            Nitrate concentration (:math:`\mu mol` nitrates m-3)
        :Returns Type:
            :class:`float`
        """
        return max(0, (nitrates / self.volume))

    # FLUX
    def calculate_mineralisation(self):
        """Mineralisation on organic N into nitrates in soil.

        :Parameters:

        :Returns:
            Rate of Nitrate mineralisation (:math:`\mu mol` h-1)
        :Returns Type:
            :class:`float`
        """
        return parameters.SOIL_PARAMETERS.MINERALISATION_RATE * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    # COMPARTMENTS

    def calculate_nitrates_derivative(self, mineralisation, soil_contributors, culm_density):
        """delta soil nitrates.

        :Parameters:
            - `mineralisation` (:class:`float`) - N mineralisation in soil (:math:`\mu mol` m-2 N nitrates)
            - `soil_contributors` (:class:`tuple`) - A tuple with (Nitrate uptake per axis (:math:`\mu mol` N nitrates), the plant id)
            - `culm_density` (:class:`dict`) - A dictionary of culm density (culm_density = {plant_id: culm_density, ...})
        :Returns:
            delta nitrates (:math:`\mu mol` N nitrates)
        :Returns Type:
            :class:`float`
        """
        Uptake_Nitrates = 0
        for root_uptake, plant_id in soil_contributors:
            Uptake_Nitrates += root_uptake * culm_density[plant_id]  # TODO: temporary, will be removed in next version
        return mineralisation - Uptake_Nitrates
