# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    cnwheat.model
    ~~~~~~~~~~~~~

    The module :mod:`cnwheat.model` defines the equations of the CN exchanges in a population of plants.

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

import warnings
import logging

import numpy as np

import parameters

class ModelWarning(UserWarning): pass
class ModelInputWarning(ModelWarning): pass

class ModelError(Exception): pass
class ModelInputError(ModelError): pass

warnings.simplefilter('always', ModelInputWarning)


def enum(**enums):
    return type('Enum', (), enums)


class Population(object):
    """
    The class :class:`Population` defines the CN exchanges at the population scale.

    A :class:`population <Population>` must have one or several :class:`plants <Plant>`.
    """

    PARAMETERS = parameters.PopulationParameters #: the internal parameters of the population

    def __init__(self, t=0, plants=None):
        if plants is None:
            plants = []
        self.plants = plants #: the list of plants
        self.t = t #: Time (h)

    def calculate_integrative_variables(self):
        """Calculate the integrative variables of the population recursively.
        """
        for plant in self.plants:
            plant.calculate_integrative_variables(self.t)


class Plant(object):
    """
    The class :class:`Plant` defines the CN exchanges at the plants scale.

    A :class:`plant <Plant>` must have one or several :class:`axes <Axis>`.
    """

    PARAMETERS = parameters.PlantParameters #: the internal parameters of the plants

    def __init__(self, axes=None, index=1):
        if axes is None:
            axes = []
        self.axes = axes #: the list of axes
        self.index = index #: the index of the plant, from 1 to n.

    def calculate_integrative_variables(self, t):
        """Calculate the integrative variables of the plant recursively.
        """
        for axis in self.axes:
            axis.calculate_integrative_variables(t)


class Axis(object):
    """
    The class :class:`Axis` defines the CN exchanges at the axes scale.

    An :class:`axis <Axis>` must have:
        * one :class:`set of roots <Roots>`,
        * one :class:`phloem <Phloem>`,
        * zero or one :class:`set of grains <Grains>`,
        * one or several :class:`phytomer<Phytomer>.
    """

    PARAMETERS = parameters.AxisParameters #: the internal parameters of the axes

    Types = enum(MAIN_STEM='MS', TILLER='T') #: the authorized types of the axes

    TYPES_STRINGS = Types.__dict__.values() #: the string values of the authorized types of the axes

    def __init__(self, roots=None, phloem=None, grains=None, phytomers=None, axis_type=Types.MAIN_STEM, index=0):
        self.roots = roots #: the roots
        self.phloem = phloem #: the phloem
        self.grains = grains #: the grains
        if phytomers is None:
            phytomers = []
        self.phytomers = phytomers #: the list of phytomers

        self.type = axis_type #: the type of the axis ; can be either Axis.Types.MAIN_STEM or Axis.Types.TILLER
        self.index = index #: the index of the axis ; 0: main stem, 1..n: tiller.

        if axis_type not in Axis.TYPES_STRINGS:
            logger = logging.getLogger(__name__)
            message = 'axis_type not in {}.'.format(Axis.TYPES_STRINGS)
            logger.exception(message)
            raise ModelInputError(message)

        if axis_type == Axis.Types.MAIN_STEM and index != 0:
            logger = logging.getLogger(__name__)
            message = 'non-zero index for {}'.format(Axis.Types.MAIN_STEM)
            logger.exception(message)
            warnings.warn(message, ModelInputWarning)

        self.id = axis_type #: the id of the axis ; the id is built from axis_type and index
        if axis_type != Axis.Types.MAIN_STEM:
            self.id += str(index)

    def calculate_integrative_variables(self, t):
        """Calculate the integrative variables of the axis recursively.
        """
        if self.roots is not None:
            self.roots.calculate_integrative_variables()
        if self.phloem is not None:
            self.phloem.calculate_integrative_variables(t)
        if self.grains is not None:
            self.grains.calculate_integrative_variables(t)
        for phytomer in self.phytomers:
            phytomer.calculate_integrative_variables(t, phytomer.index)


class Phytomer(object):
    """
    The class :class:`Phytomer` defines the CN exchanges at the phytomers scale.

    A :class:`phytomer <Phytomer>` must have either:
        * one :class:`chaff <Chaff>`,
        * OR one :class:`peduncle <Peduncle>`,
        * OR one :class:`lamina <Lamina>`, one :class:`internode <Internode>` and one :class:`sheath <Sheath>`.
    """

    PARAMETERS = parameters.PhytomerParameters #: the internal parameters of the phytomers

    def __init__(self, chaff=None, peduncle=None, lamina=None, internode=None, sheath=None, index=1):
        self.chaff = chaff #: the chaff
        self.peduncle = peduncle #: the peduncle
        self.lamina = lamina #: the lamina
        self.internode = internode #: the internode
        self.sheath = sheath #: the sheath
        self.index = index #: the index of the phytomer, from 1 to n.

    def calculate_integrative_variables(self, t, phytomer_index):
        """Calculate the integrative variables of the phytomer recursively.
        """
        if self.chaff is not None:
            self.chaff.calculate_integrative_variables(t, phytomer_index)
        if self.peduncle is not None:
            self.peduncle.calculate_integrative_variables(t, phytomer_index)
        if self.lamina is not None:
            self.lamina.calculate_integrative_variables(t, phytomer_index)
        if self.internode is not None:
            self.internode.calculate_integrative_variables(t, phytomer_index)
        if self.sheath is not None:
            self.sheath.calculate_integrative_variables(t, phytomer_index)


class Organ(object):
    """
    The class :class:`Organ` defines the CN exchanges at the organs scale.

    :class:`Organ` is the base class of all organs. DO NOT INSTANTIATE IT.
    """

    PARAMETERS = parameters.OrganParameters #: the internal parameters of the organs

    def calculate_integrative_variables(self, t, phytomer_index=None):
        """Calculate the integrative variables of the organ recursively.
        """
        pass


class Phloem(Organ):
    """
    The class :class:`Phloem` defines the CN exchanges in a phloem.
    """

    PARAMETERS = parameters.PhloemParameters #: the internal parameters of the phloems

    def __init__(self, sucrose, amino_acids):

        # variables
        self.sucrose = sucrose          #: µmol C sucrose
        self.amino_acids = amino_acids  #: µmol N amino acids


    # VARIABLES

    def calculate_conc_sucrose(self, sucrose):
        """sucrose concentration (µmol sucrose g-1 MS)
        """
        return (sucrose/Organ.PARAMETERS.NB_C_SUCROSE)/Organ.PARAMETERS.MSTRUCT_AXIS

    def calculate_conc_amino_acids(self, amino_acids):
        """Amino_acid concentration (µmol amino_acids g-1 MS)
        """
        return (amino_acids/Organ.PARAMETERS.AMINO_ACIDS_N_RATIO) / Organ.PARAMETERS.MSTRUCT_AXIS

    # COMPARTMENTS

    def calculate_sucrose_derivative(self, contributors):
        """delta sucrose of phloem integrated over delta_t (µmol C sucrose)
        """
        sucrose_derivative = 0
        for contributor in contributors:
            if isinstance(contributor, PhotosyntheticOrganElement):
                sucrose_derivative += contributor.loading_sucrose * contributor.mstruct * contributor.__class__.PARAMETERS.ALPHA
            elif isinstance(contributor, Grains):
                sucrose_derivative -= contributor.s_grain_structure + (contributor.s_grain_starch * ((contributor.structure*1E-6) * Organ.PARAMETERS.C_MOLAR_MASS)) #: Conversion of structure from umol of C to g of C
            elif isinstance(contributor, Roots):
                sucrose_derivative -= contributor.unloading_sucrose * contributor.mstruct * contributor.__class__.PARAMETERS.ALPHA
        return sucrose_derivative

    def calculate_amino_acids_derivative(self, contributors):
        """delta amino acids of phloem integrated over delat_t (µmol N amino acids)
        """
        amino_acids_derivative = 0
        for contributor in contributors:
            if isinstance(contributor, PhotosyntheticOrganElement):
                amino_acids_derivative += contributor.loading_amino_acids * contributor.mstruct * contributor.__class__.PARAMETERS.ALPHA
            elif isinstance(contributor, Grains):
                amino_acids_derivative -= contributor.s_proteins
            elif isinstance(contributor, Roots):
                amino_acids_derivative -= contributor.unloading_amino_acids * contributor.mstruct * contributor.__class__.PARAMETERS.ALPHA
        return amino_acids_derivative


class Grains(Organ):
    """
    The class :class:`Grains` defines the CN exchanges in a set of grains.
    """

    PARAMETERS = parameters.GrainsParameters #: the internal parameters of the grains

    def __init__(self, starch, structure, proteins):

        # variables
        self.starch = starch                     #: µmol of C starch
        self.structure = structure               #: µmol of C sucrose
        self.proteins = proteins                 #: µmol of N proteins

        # fluxes from phloem
        self.s_grain_structure = None            #: current rate of grain structure synthesis
        self.s_grain_starch = None               #: current rate of grain starch C synthesis
        self.s_proteins = None                   #: current rate of grain protein synthesis


    # VARIABLES

    def calculate_dry_mass(self, structure, starch, proteins):
        """Grain total dry mass (g)
        """
        #: Carbohydrates mass, grain carbohydrates supposed to be mainly starch i.e. glucose polymers (C6 H12 O6)
        C_mass = ((structure + starch)*1E-6*Organ.PARAMETERS.C_MOLAR_MASS) / Organ.PARAMETERS.HEXOSE_MOLAR_MASS_C_RATIO

        #: N mass, grain proteins were supposed to be gluten mainly composed of Glu, Gln and Pro
        N_mass = (proteins*1E-6*Organ.PARAMETERS.N_MOLAR_MASS) / Grains.PARAMETERS.AMINO_ACIDS_MOLAR_MASS_N_RATIO

        return  C_mass + N_mass

    def calculate_structural_dry_mass(self, structure):
        """Grain structural mass (g)
        """
        return (structure*1E-6*Organ.PARAMETERS.C_MOLAR_MASS) / Organ.PARAMETERS.HEXOSE_MOLAR_MASS_C_RATIO

    def calculate_protein_mass(self, proteins):
        """Grain total protein mass (g)
        """
        mass_N_proteins = proteins*1E-6 * Organ.PARAMETERS.N_MOLAR_MASS                        #: Mass of nitrogen in proteins (g)
        #mass_proteins = mass_N_proteins / Organ.PARAMETERS.AMINO_ACIDS_MOLAR_MASS_N_RATIO     #: Total mass of proteins (g)
        return mass_N_proteins

    def calculate_RGR_structure(self, sucrose_phloem):
        """Relative Growth Rate of grain structure, regulated by phloem concentrations
        """
        return ((max(0, sucrose_phloem)/(Organ.PARAMETERS.MSTRUCT_AXIS*Organ.PARAMETERS.ALPHA_AXIS)) * Grains.PARAMETERS.VMAX_RGR) / ((max(0, sucrose_phloem)/(Organ.PARAMETERS.MSTRUCT_AXIS*Organ.PARAMETERS.ALPHA_AXIS)) + Grains.PARAMETERS.K_RGR)

    def calculate_unloading_sucrose(self, s_grain_structure, s_grain_starch, structural_dry_mass):
        """Unloading of sucrose from phloem to grains integrated * delta_t (µmol sucrose)
        """
        return (s_grain_structure + (s_grain_starch * structural_dry_mass)) / Organ.PARAMETERS.NB_C_SUCROSE


    # FLUXES

    def calculate_s_grain_structure(self, t, prec_structure, RGR_structure):
        """Synthesis of grain structure integrated over delta_t (µmol C structure s-1 * DELTA_T). RGR_structure is regulated by phloem concentrations
        """
        if t<=Grains.PARAMETERS.FILLING_INIT: #: Grain enlargment
            s_grain_structure = prec_structure * RGR_structure * Organ.PARAMETERS.DELTA_T
        else:                                 #: Grain filling
            s_grain_structure = 0
        return s_grain_structure

    def calculate_s_grain_starch(self, t, sucrose_phloem):
        """Synthesis of grain C starch integrated over delta_t (µmol C starch s-1 g-1 MS * DELTA_T). Rate regulated by phloem concentrations and unloading
        """
        if t<=Grains.PARAMETERS.FILLING_INIT: #: Grain enlargment
            s_grain_starch = 0
        else:                                 #: Grain filling
            s_grain_starch = (((max(0, sucrose_phloem)/(Organ.PARAMETERS.MSTRUCT_AXIS*Organ.PARAMETERS.ALPHA_AXIS)) * Grains.PARAMETERS.VMAX_STARCH) / ((max(0, sucrose_phloem)/(Organ.PARAMETERS.MSTRUCT_AXIS*Organ.PARAMETERS.ALPHA_AXIS)) + Grains.PARAMETERS.K_STARCH)) * Organ.PARAMETERS.DELTA_T
        return s_grain_starch

    def calculate_s_proteins(self, s_grain_structure, s_grain_starch, amino_acids_phloem, sucrose_phloem, structural_dry_mass):
        """Synthesis of grain proteins over delta_t (µmol N proteins). Rate regulated by phloem concentrations and unloading. Co-transported with sucrose relatively to the ratio amino acids:sucrose in phloem
        """
        if sucrose_phloem >0:
            s_proteins = (s_grain_structure + s_grain_starch*structural_dry_mass) * (amino_acids_phloem / sucrose_phloem)
        else:
            s_proteins = 0
        return s_proteins

    # COMPARTMENTS

    def calculate_structure_derivative(self, s_grain_structure, R_growth):
        """delta grain structure integrated over delat_t (µmol C structure)
        """
        return s_grain_structure - R_growth

    def calculate_starch_derivative(self, s_grain_starch, structural_dry_mass, R_growth):
        """delta grain starch integrated over delat_t (µmol C starch)
        """
        return (s_grain_starch * structural_dry_mass) - R_growth

    def calculate_proteins_derivative(self, s_proteins):
        """delta grain proteins integrated over delat_t (µmol N proteins)
        """
        return s_proteins


class Roots(Organ):
    """
    The class :class:`Roots` defines the CN exchanges in a set of roots.
    """

    PARAMETERS = parameters.RootsParameters #: the internal parameters of the roots

    def __init__(self, mstruct, Nstruct, sucrose, nitrates, amino_acids):

        # variables
        self.mstruct = mstruct                 #: Structural mass (g)
        self.Nstruct = Nstruct                 #: Structural nitrogen (g)
        self.sucrose = sucrose                 #: µmol C sucrose
        self.nitrates = nitrates               #: µmol N nitrates
        self.amino_acids = amino_acids         #: µmol N amino acids

        # fluxes from phloem
        self.unloading_sucrose = None          #: current unloading of sucrose from phloem to roots
        self.unloading_amino_acids = None      #: current unloading of amino acids from phloem to roots

        # Integrated variables
        self.total_nitrogen = None            #: current total nitrogen amount (µmol N)

    def calculate_integrative_variables(self):
        self.total_nitrogen = self.calculate_total_nitrogen(self.nitrates, self.amino_acids, self.Nstruct)

    # VARIABLES

    def calculate_conc_sucrose(self, sucrose):
        """Sucrose concentration (µmol sucrose g-1 MS).
        This is a concentration output (expressed in amount of substance g-1 MS).
        """
        return (sucrose/self.mstruct)/Organ.PARAMETERS.NB_C_SUCROSE

    def calculate_conc_nitrates_soil(self, t):
        """Nitrate concetration in soil (µmol nitrates m-3)
        """
        return -500*t + 5E+05 # TODO: Temporary

    def calculate_conc_nitrates(self, nitrates):
        """Nitrate concentration (µmol nitrates g-1 MS)
        """
        return (nitrates/self.mstruct)

    def calculate_conc_amino_acids(self, amino_acids):
        """Amino_acid concentration (µmol amino_acids g-1 MS)
        """
        return (amino_acids/Organ.PARAMETERS.AMINO_ACIDS_N_RATIO)/self.mstruct

    def calculate_total_nitrogen(self, nitrates, amino_acids, Nstruct):
        return nitrates + amino_acids + (Nstruct / Roots.PARAMETERS.N_MOLAR_MASS)*1E6

    # FLUXES

    def calculate_unloading_sucrose(self, sucrose_phloem):
        """Unloading of sucrose from phloem to roots (µmol C sucrose unloaded s-1 g-1 MS * DELTA_T)
        """
        return (((max(0, sucrose_phloem)/(Organ.PARAMETERS.MSTRUCT_AXIS*Organ.PARAMETERS.ALPHA_AXIS)) * Roots.PARAMETERS.VMAX_SUCROSE_UNLOADING) / ((max(0, sucrose_phloem)/(Organ.PARAMETERS.MSTRUCT_AXIS*Organ.PARAMETERS.ALPHA_AXIS)) + Roots.PARAMETERS.K_SUCROSE_UNLOADING)) * Organ.PARAMETERS.DELTA_T

    def calculate_unloading_amino_acids(self, amino_acids_phloem):
        """Unloading of amino_acids from phloem to roots over delta_t (µmol N amino_acids unloaded s-1 g-1 MS)
        """
        return (((max(0, amino_acids_phloem)/(Organ.PARAMETERS.MSTRUCT_AXIS*Organ.PARAMETERS.ALPHA_AXIS)) * Roots.PARAMETERS.VMAX_AMINO_ACIDS_UNLOADING) / ((max(0, amino_acids_phloem)/(Organ.PARAMETERS.MSTRUCT_AXIS*Organ.PARAMETERS.ALPHA_AXIS)) + Roots.PARAMETERS.K_AMINO_ACIDS_UNLOADING)) * Organ.PARAMETERS.DELTA_T # TODO: Temporary

    def calculate_uptake_nitrates(self, conc_nitrates_soil, nitrates_roots, total_transpiration):
        """Uptake of nitrates by roots (µmol N nitrates imported s-1 * DELTA_T)
        """
        VMAX_HATS_MAX = Roots.PARAMETERS.A_VMAX_HATS * np.exp(-Roots.PARAMETERS.LAMBDA_VMAX_HATS*(nitrates_roots/self.mstruct))        #: Maximal rate of nitrates uptake at saturating soil N concentration;HATS (µmol N nitrates g-1 s-1)
        K_HATS = Roots.PARAMETERS.A_K_HATS * np.exp(-Roots.PARAMETERS.LAMBDA_K_HATS*(nitrates_roots/self.mstruct))                     #: Affinity coefficient of nitrates uptake at saturating soil N concentration;HATS (µmol m-3)
        HATS = (VMAX_HATS_MAX * conc_nitrates_soil)/ (K_HATS + conc_nitrates_soil)                                                     #: High Affinity Transport System (µmol N nitrates uptaked s-1 g-1 MS roots)
        K_LATS = Roots.PARAMETERS.A_LATS * np.exp(-Roots.PARAMETERS.LAMBDA_LATS*(nitrates_roots/self.mstruct))                         #: Rate of nitrates uptake at low soil N concentration; LATS (m3 g-1 s-1)
        LATS = (K_LATS * conc_nitrates_soil)                                                                                           #: Low Affinity Transport System (µmol N nitrates uptaked s-1 g-1 MS roots)

        potential_uptake = (HATS + LATS) * self.mstruct * Organ.PARAMETERS.DELTA_T                                                     #: Potential nitrate uptake (µmol N nitrates uptaked by roots integrated over delta_t)
        actual_uptake = potential_uptake * (total_transpiration/(total_transpiration + Roots.PARAMETERS.K_TR_UPTAKE_NITRATES))         #: Nitrate uptake regulated by plant transpiration (µmol N nitrates uptaked by roots)
        return actual_uptake, potential_uptake

    def calculate_s_amino_acids(self, nitrates, sucrose):
        """Rate of amino acid synthesis in roots(µmol N amino acids s-1 g-1 MS * DELTA_T)
        """
        return Roots.PARAMETERS.VMAX_AMINO_ACIDS / ((1 + Roots.PARAMETERS.K_AMINO_ACIDS_NITRATES/(nitrates/(self.mstruct*Roots.PARAMETERS.ALPHA))) * (1 + Roots.PARAMETERS.K_AMINO_ACIDS_SUCROSE/(sucrose/(self.mstruct*Roots.PARAMETERS.ALPHA)))) * Organ.PARAMETERS.DELTA_T

    def calculate_export_amino_acids(self, amino_acids, total_transpiration):
        """Total export of amino acids from roots to shoot organs (abstraction of the xylem compartment) (µmol N amino acids exported during DELTA_T (already accounted in Transpiration))
        """
        return (amino_acids/(self.mstruct * Roots.PARAMETERS.ALPHA)) * (total_transpiration/(total_transpiration + Roots.PARAMETERS.K_TR_EXPORT_AMINO_ACIDS))

    # COMPARTMENTS

    def calculate_sucrose_derivative(self, unloading_sucrose, s_amino_acids, R_Nnit_upt, R_Nnit_red, R_residual):
        """delta root sucrose integrated over delat_t (µmol C sucrose)
        """
        sucrose_consumption_AA = (s_amino_acids / Organ.PARAMETERS.AMINO_ACIDS_N_RATIO) * Organ.PARAMETERS.AMINO_ACIDS_C_RATIO      #: Contribution of sucrose to the synthesis of amino_acids

        return (unloading_sucrose - sucrose_consumption_AA) * self.mstruct - R_Nnit_upt - R_Nnit_red - R_residual

    def calculate_nitrates_derivative(self, uptake_nitrates, s_amino_acids):
        """delta root nitrates integrated over delat_t (µmol N nitrates)
        """
        import_nitrates_roots = uptake_nitrates * (1-Organ.PARAMETERS.RATIO_EXPORT_NITRATES_ROOTS)                                  #: Proportion of uptaked nitrates staying in roots
        nitrate_reduction_AA = s_amino_acids                                                                                        #: Contribution of nitrates to the synthesis of amino_acids
        return import_nitrates_roots - (nitrate_reduction_AA*self.mstruct)

    def calculate_amino_acids_derivative(self, unloading_amino_acids, s_amino_acids, export_amino_acids):
        """delta root amino acids integrated over delat_t (µmol N amino acids)
        """
        return (unloading_amino_acids + s_amino_acids)*self.mstruct  - export_amino_acids


class PhotosyntheticOrgan(Organ):
    """
    The class :class:`PhotosyntheticOrgan` defines the CN exchanges in a photosynthetic organ.

    :class:`PhotosyntheticOrgan` is the base class of all photosynthetic organs. DO NOT INSTANTIATE IT.
    """

    PARAMETERS = parameters.PhotosyntheticOrganParameters #: the internal parameters of the photosynthetic organs

    def __init__(self, exposed_element=None, enclosed_element=None):
        # variables
        self.exposed_element = exposed_element #: the exposed element
        self.enclosed_element = enclosed_element #: the enclosed element

    def calculate_integrative_variables(self, t, phytomer_index):
        for element in (self.exposed_element, self.enclosed_element):
            if element is not None:
                element.calculate_integrative_variables(t, phytomer_index)


class Chaff(PhotosyntheticOrgan):
    """
    The class :class:`Chaff` defines the CN exchanges in a chaff.
    """

    PARAMETERS = parameters.ChaffParameters #: the internal parameters of the chaffs

    def __init__(self, exposed_element=None, enclosed_element=None):
        super(Chaff, self).__init__(exposed_element, enclosed_element)


class Lamina(PhotosyntheticOrgan):
    """
    The class :class:`Lamina` defines the CN exchanges in a lamina.
    """

    PARAMETERS = parameters.LaminaParameters #: the internal parameters of the laminae

    def __init__(self, exposed_element=None, enclosed_element=None):
        super(Lamina, self).__init__(exposed_element, enclosed_element)


class Internode(PhotosyntheticOrgan):
    """
    The class :class:`Internode` defines the CN exchanges in an internode.
    """

    PARAMETERS = parameters.InternodeParameters #: the internal parameters of the internodes

    def __init__(self, exposed_element=None, enclosed_element=None):
        super(Internode, self).__init__(exposed_element, enclosed_element)


class Peduncle(PhotosyntheticOrgan):
    """
    The class :class:`Peduncle` defines the CN exchanges in a peduncle.
    """

    PARAMETERS = parameters.PeduncleParameters #: the internal parameters of the peduncles

    def __init__(self, exposed_element=None, enclosed_element=None):
        super(Peduncle, self).__init__(exposed_element, enclosed_element)


class Sheath(PhotosyntheticOrgan):
    """
    The class :class:`Sheath` defines the CN exchanges in a sheath.
    """

    PARAMETERS = parameters.SheathParameters #: the internal parameters of the sheaths

    def __init__(self, exposed_element=None, enclosed_element=None):
        super(Sheath, self).__init__(exposed_element, enclosed_element)


class PhotosyntheticOrganElement(object):
    """
    The class :class:`PhotosyntheticOrganElement` defines the CN exchanges in a photosynthetic organ element.

    :class:`PhotosyntheticOrganElement` is the base class of all photosynthetic organs elements. DO NOT INSTANTIATE IT.
    """

    PARAMETERS = parameters.PhotosyntheticOrganElementParameters #: the internal parameters of the photosynthetic organs elements

    def __init__(self, area, mstruct, Nstruct, width, height, triosesP, starch,
                 sucrose, fructan, nitrates, amino_acids, proteins,
                 An=None, Tr=None):

        self.area = area                     #: area (m-2)
        self.mstruct = mstruct               #: Structural mass (g)
        self.Nstruct = Nstruct               #: Structural nitrogen (g)
        self.width = width                   #: Width (or diameter for stem organ elements) (m)
        self.height = height                 #: Height of the element from soil (m)
        self.An = An                         #: Net assimilation (µmol m-2 s-1)
        self.Tr = Tr                         #: Transpiration (mm s-1)

        self.triosesP = triosesP
        self.starch = starch
        self.sucrose = sucrose
        self.fructan = fructan
        self.nitrates = nitrates
        self.amino_acids = amino_acids
        self.proteins = proteins

        # fluxes to phloem
        self.loading_sucrose = None           #: current rate of sucrose loading to phloem
        self.loading_amino_acids = None       #: current rate of amino acids loading to phloem

        # Integrated variables
        self.surfacic_nitrogen = None         #: current surfacic nitrogen (g m-2)
        self.total_nitrogen = None            #: current total nitrogen amount (µmol N)

    def calculate_integrative_variables(self, t, phytomer_index):
        """Calculate the integrative variables of the element.
        """
        self.surfacic_nitrogen = self.calculate_surfacic_nitrogen(t, self.nitrates, self.amino_acids, self.proteins, phytomer_index)
        self.total_nitrogen = self.calculate_total_nitrogen(self.nitrates, self.amino_acids, self.proteins, self.Nstruct)

    # VARIABLES
    def calculate_photosynthesis(self, t, An, phytomer_index):
        """Total photosynthesis of an element integrated over DELTA_T (µmol CO2 on element area integrated over delat_t)
        """
        return An * self._calculate_green_area(t, phytomer_index) * PhotosyntheticOrgan.PARAMETERS.DELTA_T

    def calculate_total_Rd(self, t, Rd, phytomer_index):
        """Total respiration of an element integrated over DELTA_T (µmol CO2 on element area integrated over delat_t)
        """
        return Rd * self._calculate_green_area(t, phytomer_index) * PhotosyntheticOrgan.PARAMETERS.DELTA_T

    def calculate_transpiration(self, t, Tr, phytomer_index):
        """Total transpiration of an element integrated over DELTA_T (mm of H2O on element area integrated over delat_t)
        """
        return Tr * self._calculate_green_area(t, phytomer_index) * PhotosyntheticOrgan.PARAMETERS.DELTA_T

    def _calculate_green_area(self, t, phytomer_index):
        """Compute green area of the element.
        """
        return self.area

    def calculate_conc_triosesP(self, triosesP):
        """Trioses Phosphate concentration (µmol triosesP g-1 MS).
        This is a concentration output (expressed in amount of substance g-1 MS).
        """
        return (triosesP/self.mstruct)/Organ.PARAMETERS.NB_C_TRIOSEP

    def calculate_conc_sucrose(self, sucrose):
        """Sucrose concentration (µmol sucrose g-1 MS).
        This is a concentration output (expressed in amount of substance g-1 MS).
        """
        return (sucrose/self.mstruct)/Organ.PARAMETERS.NB_C_SUCROSE

    def calculate_conc_starch(self, starch):
        """Starch concentration (µmol starch g-1 MS (eq glucose)).
        This is a concentration output (expressed in amount of substance g-1 MS).
        """
        return (starch/self.mstruct)/Organ.PARAMETERS.NB_C_HEXOSES

    def calculate_conc_fructan(self, fructan):
        """fructan concentration (µmol fructan g-1 MS (eq glucose))
        """
        return (fructan/self.mstruct)/Organ.PARAMETERS.NB_C_HEXOSES

    def calculate_regul_s_fructan(self, loading_sucrose):
        """Inhibition of fructan synthesis by the loading of sucrose to phloem
        """
        return ((PhotosyntheticOrgan.PARAMETERS.VMAX_REGUL_SFRUCTAN * PhotosyntheticOrgan.PARAMETERS.K_REGUL_SFRUCTAN**(PhotosyntheticOrgan.PARAMETERS.N_REGUL_SFRUCTAN)) / (max(0, loading_sucrose**(PhotosyntheticOrgan.PARAMETERS.N_REGUL_SFRUCTAN)) + PhotosyntheticOrgan.PARAMETERS.K_REGUL_SFRUCTAN**(PhotosyntheticOrgan.PARAMETERS.N_REGUL_SFRUCTAN)))

    def calculate_conc_nitrates(self, nitrates):
        """Nitrate concentration (µmol nitrates g-1 MS)
        """
        return (nitrates/self.mstruct)

    def calculate_conc_amino_acids(self, amino_acids):
        """Amino_acid concentration (µmol amino acids g-1 MS)
        """
        return (amino_acids/PhotosyntheticOrgan.PARAMETERS.AMINO_ACIDS_N_RATIO) / self.mstruct

    def calculate_conc_proteins(self, proteins):
        """Protein concentration (g proteins g-1 MS)
        """
        mass_N_proteins = proteins*1E-6 * PhotosyntheticOrgan.PARAMETERS.N_MOLAR_MASS                        #: Mass of nitrogen in proteins (g)
        mass_proteins = mass_N_proteins / PhotosyntheticOrgan.PARAMETERS.AMINO_ACIDS_MOLAR_MASS_N_RATIO      #: Total mass of proteins (g)
        return (mass_proteins / self.mstruct)

    def calculate_surfacic_nitrogen(self, t, nitrates, amino_acids, proteins, phytomer_index):
        """Surfacic content of nitrogen (g m-2)
        """
        mass_N_tot = (nitrates + amino_acids + proteins)*1E-6 * PhotosyntheticOrgan.PARAMETERS.N_MOLAR_MASS + self.Nstruct
        green_area = self._calculate_green_area(t, phytomer_index)
        return (mass_N_tot / green_area)

    def calculate_total_nitrogen(self, nitrates, amino_acids, proteins, Nstruct):
        return nitrates + amino_acids + proteins + (Nstruct / PhotosyntheticOrgan.PARAMETERS.N_MOLAR_MASS)*1E6

    # FLUXES

    def calculate_s_starch(self, triosesP):
        """Rate of starch synthesis from triosesP (µmol C starch s-1 g-1 MS * DELTA_T).
        """
        return (((max(0, triosesP)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) * PhotosyntheticOrgan.PARAMETERS.VMAX_STARCH) / ((max(0, triosesP)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) + PhotosyntheticOrgan.PARAMETERS.K_STARCH)) * PhotosyntheticOrgan.PARAMETERS.DELTA_T

    def calculate_d_starch(self, starch):
        """Rate of starch degradation from triosesP (µmol C starch s-1 g-1 MS * DELTA_T).
        """
        return max(0, PhotosyntheticOrgan.PARAMETERS.DELTA_DSTARCH * (starch/(self.mstruct*self.__class__.PARAMETERS.ALPHA))) * PhotosyntheticOrgan.PARAMETERS.DELTA_T

    def calculate_s_sucrose(self, triosesP):
        """Rate of sucrose synthesis from triosesP (µmol C sucrose s-1 g-1 MS * DELTA_T).
        """
        return (((max(0,triosesP)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) * PhotosyntheticOrgan.PARAMETERS.VMAX_SUCROSE) / ((max(0, triosesP)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) + PhotosyntheticOrgan.PARAMETERS.K_SUCROSE)) * PhotosyntheticOrgan.PARAMETERS.DELTA_T

    def calculate_loading_sucrose(self, sucrose, sucrose_phloem):
        """Rate of sucrose loading to phloem (µmol C sucrose s-1 g-1 MS * DELTA_T).
        """
        driving_sucrose_compartment = max(sucrose / (self.mstruct*self.__class__.PARAMETERS.ALPHA), sucrose_phloem/(PhotosyntheticOrgan.PARAMETERS.MSTRUCT_AXIS*parameters.OrganParameters.ALPHA_AXIS))
        diff_sucrose = sucrose/(self.mstruct*self.__class__.PARAMETERS.ALPHA) - sucrose_phloem/(PhotosyntheticOrgan.PARAMETERS.MSTRUCT_AXIS*parameters.OrganParameters.ALPHA_AXIS)
        conductance = PhotosyntheticOrgan.PARAMETERS.SIGMA_SUCROSE * PhotosyntheticOrgan.PARAMETERS.BETA * self.mstruct**(2/3)
        return driving_sucrose_compartment * diff_sucrose * conductance * PhotosyntheticOrgan.PARAMETERS.DELTA_T

    def calculate_s_fructan(self, sucrose, regul_s_fructan):
        """Rate of fructan synthesis (µmol C fructan s-1 g-1 MS * DELTA_T)
        """
        return (((max(0, sucrose)/(self.mstruct*self.__class__.PARAMETERS.ALPHA))**(PhotosyntheticOrgan.PARAMETERS.N_SFRUCTAN) * PhotosyntheticOrgan.PARAMETERS.VMAX_SFRUCTAN) / ((max(0, sucrose)/(self.mstruct*self.__class__.PARAMETERS.ALPHA))**(PhotosyntheticOrgan.PARAMETERS.N_SFRUCTAN) + PhotosyntheticOrgan.PARAMETERS.K_SFRUCTAN**(PhotosyntheticOrgan.PARAMETERS.N_SFRUCTAN))) * regul_s_fructan * PhotosyntheticOrgan.PARAMETERS.DELTA_T

    def calculate_d_fructan(self, sucrose, fructan):
        """Rate of fructan degradation (µmol C fructan s-1 g-1 MS)
        """
        return min((PhotosyntheticOrgan.PARAMETERS.K_DFRUCTAN * PhotosyntheticOrgan.PARAMETERS.VMAX_DFRUCTAN) / ((max(0, sucrose)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) + PhotosyntheticOrgan.PARAMETERS.K_DFRUCTAN) , max(0, fructan)) * PhotosyntheticOrgan.PARAMETERS.DELTA_T

    def calculate_nitrates_import(self, roots_uptake_nitrate, organ_transpiration, total_transpiration):
        """Total nitrates imported from roots (through xylem) distributed relatively to element transpiration (µmol N nitrates integrated over delta_t [already accounted in transpiration])
        """
        if organ_transpiration>0:
            nitrates_import = roots_uptake_nitrate * (organ_transpiration/total_transpiration)* PhotosyntheticOrgan.PARAMETERS.RATIO_EXPORT_NITRATES_ROOTS # Proportion of uptaked nitrates exported from roots to shoot
        else: # Avoids further float division by zero error
            nitrates_import = 0
        return nitrates_import

    def calculate_amino_acids_import(self, roots_exported_amino_acids, organ_transpiration, total_transpiration):
        """Total Amino acids imported from roots (through xylem) distributed relatively to element transpiration (µmol N Amino acids integrated over delta_t [already accounted in transpiration])
        """
        if organ_transpiration>0:
            amino_acids_import = roots_exported_amino_acids * (organ_transpiration/total_transpiration)
        else:
            amino_acids_import = 0
        return amino_acids_import

    def calculate_s_amino_acids(self, nitrates, triosesP):
        """Rate of amino acid synthesis (µmol N amino acids s-1 g-1 MS * DELTA_T)
        """
        if nitrates <=0 or triosesP <=0:
            calculate_s_amino_acids = 0
        else:
            calculate_s_amino_acids = PhotosyntheticOrgan.PARAMETERS.VMAX_AMINO_ACIDS / ((1 + PhotosyntheticOrgan.PARAMETERS.K_AMINO_ACIDS_NITRATES/(nitrates/(self.mstruct*self.__class__.PARAMETERS.ALPHA))) * (1 + PhotosyntheticOrgan.PARAMETERS.K_AMINO_ACIDS_TRIOSESP/(triosesP/(self.mstruct*self.__class__.PARAMETERS.ALPHA)))) * PhotosyntheticOrgan.PARAMETERS.DELTA_T
        return calculate_s_amino_acids

    def calculate_s_proteins(self, amino_acids):
        """Rate of protein synthesis (µmol N proteins s-1 g-1 MS * DELTA_T)
        """
        calculate_s_proteins = (((max(0,amino_acids)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) * PhotosyntheticOrgan.PARAMETERS.VMAX_SPROTEINS) / ((max(0, amino_acids)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) + PhotosyntheticOrgan.PARAMETERS.K_SPROTEINS)) * PhotosyntheticOrgan.PARAMETERS.DELTA_T
        return calculate_s_proteins

    def calculate_d_proteins(self, proteins):
        """Rate of protein degradation (µmol N proteins s-1 g-1 MS * DELTA_T)
        """
        return max(0, PhotosyntheticOrgan.PARAMETERS.DELTA_DPROTEINS * (proteins/(self.mstruct*self.__class__.PARAMETERS.ALPHA))) * PhotosyntheticOrgan.PARAMETERS.DELTA_T

    def calculate_loading_amino_acids(self, amino_acids, amino_acids_phloem):
        """Rate of amino acids loading to phloem (µmol N amino acids s-1 g-1 MS * DELTA_T)
        """
        driving_amino_acids_compartment = max(amino_acids / (self.mstruct*self.__class__.PARAMETERS.ALPHA), amino_acids_phloem/(PhotosyntheticOrgan.PARAMETERS.MSTRUCT_AXIS*parameters.OrganParameters.ALPHA_AXIS))
        diff_amino_acids = amino_acids/(self.mstruct*self.__class__.PARAMETERS.ALPHA) - amino_acids_phloem/(PhotosyntheticOrgan.PARAMETERS.MSTRUCT_AXIS*parameters.OrganParameters.ALPHA_AXIS)
        conductance = PhotosyntheticOrgan.PARAMETERS.SIGMA_AMINO_ACIDS * PhotosyntheticOrgan.PARAMETERS.BETA * self.mstruct**(2/3)
        return driving_amino_acids_compartment * diff_amino_acids * conductance * PhotosyntheticOrgan.PARAMETERS.DELTA_T

    # COMPARTMENTS

    def calculate_triosesP_derivative(self, photosynthesis, s_sucrose, s_starch, s_amino_acids):
        """ delta triosesP of element integrated over delat_t (µmol C triosesP).
        """
        triosesP_consumption_AA = (s_amino_acids / PhotosyntheticOrgan.PARAMETERS.AMINO_ACIDS_N_RATIO) * PhotosyntheticOrgan.PARAMETERS.AMINO_ACIDS_C_RATIO #: Contribution of triosesP to the synthesis of amino_acids
        return max(0, photosynthesis) - (s_sucrose + s_starch + triosesP_consumption_AA) * (self.mstruct*self.__class__.PARAMETERS.ALPHA)

    def calculate_starch_derivative(self, s_starch, d_starch):
        """delta starch of element integrated over delat_t (µmol C starch).
        """
        return (s_starch - d_starch) * (self.mstruct*self.__class__.PARAMETERS.ALPHA)

    def calculate_sucrose_derivative(self, s_sucrose, d_starch, loading_sucrose, s_fructan, d_fructan, R_phloem_loading, R_Nnit_red, R_residual):
        """delta sucrose of element integrated over delat_t (µmol C sucrose)
        """
        return (s_sucrose + d_starch + d_fructan - s_fructan - loading_sucrose) * (self.mstruct*self.__class__.PARAMETERS.ALPHA) - R_phloem_loading - R_Nnit_red - R_residual

    def calculate_fructan_derivative(self, s_fructan, d_fructan):
        """delta fructan integrated over delat_t (µmol C fructan)
        """
        return (s_fructan - d_fructan)* (self.mstruct*self.__class__.PARAMETERS.ALPHA)

    def calculate_nitrates_derivative(self, nitrates_import, s_amino_acids):
        """delta nitrates integrated over delat_t (µmol N nitrates)
        """
        nitrate_reduction_AA = s_amino_acids  #: Contribution of nitrates to the synthesis of amino_acids
        return nitrates_import - (nitrate_reduction_AA*self.mstruct*self.__class__.PARAMETERS.ALPHA)

    def calculate_amino_acids_derivative(self, amino_acids_import, s_amino_acids, s_proteins, d_proteins, loading_amino_acids):
        """delta amino acids integrated over delat_t (µmol N amino acids)
        """
        return amino_acids_import + (s_amino_acids + d_proteins - s_proteins - loading_amino_acids) * (self.mstruct*self.__class__.PARAMETERS.ALPHA)

    def calculate_proteins_derivative(self, s_proteins, d_proteins):
        """delta proteins integrated over delat_t (µmol N proteins)
        """
        return (s_proteins - d_proteins) * (self.mstruct*self.__class__.PARAMETERS.ALPHA)


class ChaffElement(PhotosyntheticOrganElement):
    """
    The class :class:`ChaffElement` defines the CN exchanges in a chaff element.
    """

    PARAMETERS = parameters.ChaffElementParameters #: the internal parameters of the chaffs elements


class LaminaElement(PhotosyntheticOrganElement):
    """
    The class :class:`LaminaElement` defines the CN exchanges in a lamina element.
    """

    PARAMETERS = parameters.LaminaElementParameters #: the internal parameters of the laminae elements

    # VARIABLES

    def _calculate_green_area(self, t, phytomer_index):
        """Compute green area of the lamina element.
        """
        value_inflexion = LaminaElement.PARAMETERS.INFLEXION_POINTS[phytomer_index]
        green_area = max(0, min(self.area, -8.04E-6*t + value_inflexion))

        return green_area


class InternodeElement(PhotosyntheticOrganElement):
    """
    The class :class:`InternodeElement` defines the CN exchanges in an internode element.
    """

    PARAMETERS = parameters.InternodeElementParameters #: the internal parameters of the internodes elements


class PeduncleElement(PhotosyntheticOrganElement):
    """
    The class :class:`PeduncleElement` defines the CN exchanges in a peduncle element.
    """

    PARAMETERS = parameters.PeduncleElementParameters #: the internal parameters of the peduncles elements


class SheathElement(PhotosyntheticOrganElement):
    """
    The class :class:`SheathElement` defines the CN exchanges in a sheath element.
    """

    PARAMETERS = parameters.SheathElementParameters #: the internal parameters of the sheaths elements

