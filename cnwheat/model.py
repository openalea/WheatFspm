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

import numpy as np

import parameters


class Population(object):
    """
    The class :class:`Population` defines the CN exchanges at the population scale.

    A :class:`population <Population>` must have one or several :class:`plants <Plant>`.
    """

    PARAMETERS = parameters.PopulationParameters #: the internal parameters of the population

    def __init__(self, plants=None):
        if plants is None:
            plants = []
        self.plants = plants #: the list of plants

    def calculate_integrative_variables(self):
        """Calculate the integrative variables of the population recursively.
        """
        for plant in self.plants:
            plant.calculate_integrative_variables()


class Plant(object):
    """
    The class :class:`Plant` defines the CN exchanges at the plants scale.

    A :class:`plant <Plant>` must have one or several :class:`axes <Axis>`.
    """

    PARAMETERS = parameters.PlantParameters #: the internal parameters of the plants

    def __init__(self, axes=None):
        if axes is None:
            axes = []
        self.axes = axes #: the list of axes
        
    def calculate_integrative_variables(self):
        """Calculate the integrative variables of the plant recursively.
        """
        for axis in self.axes:
            axis.calculate_integrative_variables()


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

    def __init__(self, roots=None, phloem=None, grains=None, soil=None, phytomers=None):
        self.roots = roots #: the roots
        self.phloem = phloem #: the phloem
        self.grains = grains #: the grains
        self.soil = soil #: the soil
        if phytomers is None:
            phytomers = []
        self.phytomers = phytomers #: the list of phytomers

    def calculate_integrative_variables(self):
        """Calculate the integrative variables of the axis recursively.
        """
        if self.roots is not None:
            self.roots.calculate_integrative_variables()
        if self.phloem is not None:
            self.phloem.calculate_integrative_variables()
        if self.grains is not None:
            self.grains.calculate_integrative_variables()
        if self.soil is not None:
            self.soil.calculate_integrative_variables()
        for phytomer in self.phytomers:
            phytomer.calculate_integrative_variables()
            
    @classmethod
    def get_axis_id(cls, axis_index):
        if axis_index == 0:
            axis_id = 'MS'
        else:
            axis_id = 'T' + str(axis_index)
        return axis_id


class Phytomer(object):
    """
    The class :class:`Phytomer` defines the CN exchanges at the phytomers scale.

    A :class:`phytomer <Phytomer>` must have either:
        * one :class:`chaff <Chaff>`,
        * OR one :class:`peduncle <Peduncle>`,
        * OR one :class:`lamina <Lamina>`, one :class:`internode <Internode>` and one :class:`sheath <Sheath>`.
    """

    PARAMETERS = parameters.PhytomerParameters #: the internal parameters of the phytomers

    def __init__(self, chaff=None, peduncle=None, lamina=None, internode=None, sheath=None):
        self.chaff = chaff #: the chaff
        self.peduncle = peduncle #: the peduncle
        self.lamina = lamina #: the lamina
        self.internode = internode #: the internode
        self.sheath = sheath #: the sheath

    def calculate_integrative_variables(self):
        """Calculate the integrative variables of the phytomer recursively.
        """
        if self.chaff is not None:
            self.chaff.calculate_integrative_variables()
        if self.peduncle is not None:
            self.peduncle.calculate_integrative_variables()
        if self.lamina is not None:
            self.lamina.calculate_integrative_variables()
        if self.internode is not None:
            self.internode.calculate_integrative_variables()
        if self.sheath is not None:
            self.sheath.calculate_integrative_variables()


class Organ(object):
    """
    The class :class:`Organ` defines the CN exchanges at the organs scale.

    :class:`Organ` is the base class of all organs. DO NOT INSTANTIATE IT.
    """

    PARAMETERS = parameters.OrganParameters #: the internal parameters of the organs
    
    def initialize(self):
        """Initialize the derived attributes of the organ.
        """
        pass

    def calculate_integrative_variables(self):
        """Calculate the integrative variables of the organ recursively.
        """
        pass
    

class Phloem(Organ):
    """
    The class :class:`Phloem` defines the CN exchanges in a phloem.
    """

    PARAMETERS = parameters.PhloemParameters #: the internal parameters of the phloems

    def __init__(self, sucrose=None, amino_acids=None):

        # state variables
        self.sucrose = sucrose          #: µmol C sucrose
        self.amino_acids = amino_acids  #: µmol N amino acids


    # VARIABLES

    def calculate_conc_sucrose(self, sucrose):
        """Sucrose concentration. Related to the structural dry mass of the culm

        :Parameters:
            - `sucrose` (:class:`float`) - Amount of sucrose in phloem (µmol C)
        :Returns:
            Sucrose concentration (µmol sucrose g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return (sucrose/Organ.PARAMETERS.NB_C_SUCROSE)/Organ.PARAMETERS.MSTRUCT_AXIS

    def calculate_conc_amino_acids(self, amino_acids):
        """Amino_acids concentration. Related to the structural dry mass of the culm.

        :Parameters:
            - `amino_acids` (:class:`float`) - Amount of amino_acids in phloem (µmol N)
        :Returns:
            Amino_acids concentration (µmol amino_acids g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return (amino_acids/Organ.PARAMETERS.AMINO_ACIDS_N_RATIO) / Organ.PARAMETERS.MSTRUCT_AXIS

    # COMPARTMENTS

    def calculate_sucrose_derivative(self, contributors):
        """delta sucrose

        :Parameters:
            - `contributors` (:class:`object`) - Organs exchanging C with the phloem
        :Returns:
            delta sucrose (µmol C sucrose)
        :Returns Type:
            :class:`float`
        """
        sucrose_derivative = 0
        for contributor in contributors:
            if isinstance(contributor, PhotosyntheticOrganElement):
                sucrose_derivative += contributor.loading_sucrose
            elif isinstance(contributor, Grains):
                sucrose_derivative -= contributor.s_grain_structure + (contributor.s_grain_starch * contributor.structural_dry_mass)
            elif isinstance(contributor, Roots):
                sucrose_derivative -= contributor.unloading_sucrose * contributor.mstruct * contributor.__class__.PARAMETERS.ALPHA

        return sucrose_derivative

    def calculate_amino_acids_derivative(self, contributors):
        """delta amino acids

        :Parameters:
            - `contributors` (:class:`object`) - Organs exchanging N with the phloem
        :Returns:
            delta amino acids (µmol N amino acids)
        :Returns Type:
            :class:`float`
        """
        amino_acids_derivative = 0
        for contributor in contributors:
            if isinstance(contributor, PhotosyntheticOrganElement):
                amino_acids_derivative += contributor.loading_amino_acids
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

    def __init__(self, starch=None, structure=None, proteins=None):

        # state variables
        self.starch = starch                     #: µmol of C starch
        self.structure = structure               #: µmol of C sucrose
        self.proteins = proteins                 #: µmol of N proteins
        
        # derived attributes
        self.structural_dry_mass = None          #: g of MS

        # fluxes from phloem
        self.s_grain_structure = None            #: current synthesis of grain structure
        self.s_grain_starch = None               #: current synthesis of grain starch
        self.s_proteins = None                   #: current synthesis of grain proteins

    def initialize(self):
        """Initialize the derived attributes of the organ.
        """
        self.structural_dry_mass = self.calculate_structural_dry_mass(self.structure)

    # VARIABLES

    def calculate_dry_mass(self, structure, starch, proteins):
        """Grain total dry mass.

        :Parameters:
            - `structure` (:class:`float`) - Grain structural C mass (µmol C)
            - `starch` (:class:`float`) - Grain starch content (µmol C)
            - `proteins` (:class:`float`) - Grain protein content (µmol N)
        :Returns:
            Grain total dry mass (g)
        :Returns Type:
            :class:`float`
        """
        #: Carbohydrates mass, grain carbohydrates supposed to be mainly starch i.e. glucose polymers (C6 H12 O6)
        C_mass = ((structure + starch)*1E-6*Organ.PARAMETERS.C_MOLAR_MASS) / Organ.PARAMETERS.HEXOSE_MOLAR_MASS_C_RATIO

        #: N mass, grain proteins were supposed to be gluten mainly composed of Glu, Gln and Pro
        N_mass = (proteins*1E-6*Organ.PARAMETERS.N_MOLAR_MASS) / Grains.PARAMETERS.AMINO_ACIDS_MOLAR_MASS_N_RATIO

        return  C_mass + N_mass

    def calculate_structural_dry_mass(self, structure):
        """Grain structural dry mass.

        :Parameters:
            - `structure` (:class:`float`) - Grain structural C mass (µmol C)
        :Returns:
            Grain structural dry mass (g)
        :Returns Type:
            :class:`float`
        """
        return (structure*1E-6*Organ.PARAMETERS.C_MOLAR_MASS) / Organ.PARAMETERS.RATIO_C_MSTRUCT

    def calculate_protein_mass(self, proteins):
        """Grain total protein mass.

        :Parameters:
            - `proteins` (:class:`float`) - Grain protein content (µmol N)
        :Returns:
            Grain total protein mass (g)
        :Returns Type:
            :class:`float`
        """
        mass_N_proteins = proteins*1E-6 * Organ.PARAMETERS.N_MOLAR_MASS                        #: Mass of nitrogen in proteins (g)
        #mass_proteins = mass_N_proteins / Organ.PARAMETERS.AMINO_ACIDS_MOLAR_MASS_N_RATIO     #: Total mass of proteins (g)
        return mass_N_proteins

    def calculate_RGR_structure(self, sucrose_phloem):
        """Relative Growth Rate of grain structure, regulated by sucrose concentration in phloem.

        :Parameters:
            - `sucrose_phloem` (:class:`float`) - Sucrose concentration in phloem (µmol C g-1 mstruct)
        :Returns:
            RGR of grain structure (dimensionless?)
        :Returns Type:
            :class:`float`
        """
        return ((max(0, sucrose_phloem)/(Organ.PARAMETERS.MSTRUCT_AXIS*Organ.PARAMETERS.ALPHA_AXIS)) * Grains.PARAMETERS.VMAX_RGR) / ((max(0, sucrose_phloem)/(Organ.PARAMETERS.MSTRUCT_AXIS*Organ.PARAMETERS.ALPHA_AXIS)) + Grains.PARAMETERS.K_RGR)


    # FLUXES

    def calculate_s_grain_structure(self, t, prec_structure, RGR_structure):
        """Rate of grain structure synthesis integrated over DELTA_T (µmol C structure s-1 * DELTA_T).
        Exponential function, RGR regulated by sucrose concentration in the phloem.

        :Parameters:
            - `t` (:class:`float`) - Time of the simulation (s)
            - `prec_structure` (:class:`float`) - Grain structure at t-1 (µmol C)
            - `RGR_structure` (:class:`float`) - Relative Growth Rate of grain structure (dimensionless?)
        :Returns:
            Synthesis of grain structure (µmol C)
        :Returns Type:
            :class:`float`
        """
        if t<=Grains.PARAMETERS.FILLING_INIT: #: Grain enlargment
            s_grain_structure = prec_structure * RGR_structure * Organ.PARAMETERS.DELTA_T
        else:                                 #: Grain filling
            s_grain_structure = 0
        return s_grain_structure

    def calculate_s_grain_starch(self, t, sucrose_phloem):
        """Rate of starch synthesis in grains (i.e. grain filling) integrated over DELTA_T (µmol C starch g-1 mstruct s-1 * DELTA_T).
        Michaelis-Menten function of sucrose concentration in the phloem.

        :Parameters:
            - `t` (:class:`float`) - Time of the simulation (s)
            - `sucrose_phloem` (:class:`float`) - Sucrose concentration in phloem (µmol C g-1 mstruct)
        :Returns:
            Synthesis of grain starch (µmol C g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        if t<= Grains.PARAMETERS.FILLING_INIT:   #: Grain enlargment
            s_grain_starch = 0
        elif t> Grains.PARAMETERS.FILLING_END:   #: Grain maturity
            s_grain_starch = 0
        else:                                    #: Grain filling
            s_grain_starch = (((max(0, sucrose_phloem)/(Organ.PARAMETERS.MSTRUCT_AXIS*Organ.PARAMETERS.ALPHA_AXIS)) * Grains.PARAMETERS.VMAX_STARCH) / ((max(0, sucrose_phloem)/(Organ.PARAMETERS.MSTRUCT_AXIS*Organ.PARAMETERS.ALPHA_AXIS)) + Grains.PARAMETERS.K_STARCH)) * Organ.PARAMETERS.DELTA_T
        return s_grain_starch

    def calculate_s_proteins(self, s_grain_structure, s_grain_starch, amino_acids_phloem, sucrose_phloem, structural_dry_mass):
        """Protein synthesis in grains integrated over DELTA_T.
        N is assumed to be co-transported along with the unloaded sucrose from phloem (using the ratio amino acids:sucrose of phloem).

        :Parameters:
            - `s_grain_structure` (:class:`float`) - Synthesis of grain structure (µmol C)
            - `s_grain_starch` (:class:`float`) - Synthesis of grain starch (µmol C g-1 mstruct)
            - `amino_acids_phloem` (:class:`float`) - Amino acids concentration in phloem (µmol N g-1 mstruct)
            - `sucrose_phloem` (:class:`float`) - Sucrose concentration in phloem (µmol C g-1 mstruct)
            - `structural_dry_mass` (:class:`float`) - Grain structural dry mass (g)
        :Returns:
            Synthesis of grain proteins (µmol N)
        :Returns Type:
            :class:`float`
        """
        if sucrose_phloem >0:
            s_proteins = (s_grain_structure + s_grain_starch*structural_dry_mass) * (amino_acids_phloem / sucrose_phloem)
        else:
            s_proteins = 0
        return s_proteins

    # COMPARTMENTS

    def calculate_structure_derivative(self, s_grain_structure, R_growth):
        """delta grain structure.

        :Parameters:
            - `s_grain_structure` (:class:`float`) - Synthesis of grain structure (µmol C)
            - `R_growth` (:class:`float`) - Grain growth respiration (µmol C respired)
        :Returns:
            delta grain structure (µmol C structure)
        :Returns Type:
            :class:`float`
        """
        return s_grain_structure - R_growth

    def calculate_starch_derivative(self, s_grain_starch, structural_dry_mass, R_growth):
        """delta grain starch.

        :Parameters:
            - `s_grain_starch` (:class:`float`) - Synthesis of grain starch (µmol C g-1 mstruct)
            - `structural_dry_mass` (:class:`float`) - Grain structural dry mass (g)
            - `R_growth` (:class:`float`) - Grain growth respiration (µmol C respired)
        :Returns:
            delta grain starch (µmol C starch)
        :Returns Type:
            :class:`float`
        """
        return (s_grain_starch * structural_dry_mass) - R_growth

    def calculate_proteins_derivative(self, s_proteins):
        """delta grain proteins.

        :Parameters:
            - `s_proteins` (:class:`float`) - Synthesis of grain proteins (µmol N)
        :Returns:
            delta grain proteins (µmol N proteins)
        :Returns Type:
            :class:`float`
        """
        return s_proteins


class Roots(Organ):
    """
    The class :class:`Roots` defines the CN exchanges in a set of roots.
    """

    PARAMETERS = parameters.RootsParameters #: the internal parameters of the roots

    def __init__(self, mstruct=None, Nstruct=None, sucrose=None, nitrates=None, amino_acids=None, cytokinines=None, mstruct_C_growth=None,
                 Nstruct_N_growth=None):

        # state parameters
        self.mstruct_C_growth = mstruct_C_growth #: Growth of root structural dry mass (µmol C)
        self.Nstruct_N_growth = Nstruct_N_growth #: Growth of root structural N mass (µmol N)

        # state variables
        self.mstruct = mstruct                 #: Structural mass (g)
        self.Nstruct = Nstruct                 #: Structural N mass (g)
        self.sucrose = sucrose                 #: µmol C sucrose
        self.nitrates = nitrates               #: µmol N nitrates
        self.amino_acids = amino_acids         #: µmol N amino acids
        self.cytokinines = cytokinines         #: UA cytokinines

        # fluxes from phloem
        self.unloading_sucrose = None          #: current unloading of sucrose from phloem to roots
        self.unloading_amino_acids = None      #: current unloading of amino acids from phloem to roots

        # Integrated variables
        self.total_organic_nitrogen = None     #: current amount of organic N (µmol N)

    def calculate_integrative_variables(self):
        self.total_organic_nitrogen = self.calculate_total_organic_nitrogen(self.amino_acids, self.Nstruct)

    # VARIABLES

    def calculate_conc_sucrose(self, sucrose):
        """Sucrose concentration.

        :Parameters:
            - `sucrose` (:class:`float`) - Amount of sucrose (µmol C)
        :Returns:
            Sucrose concentration (µmol sucrose g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return (sucrose/self.mstruct)/Organ.PARAMETERS.NB_C_SUCROSE

    def calculate_conc_nitrates(self, nitrates):
        """Nitrate concentration.

        :Parameters:
            - `nitrates` (:class:`float`) - Amount of nitrates (µmol N)
        :Returns:
            Nitrate concentration (µmol nitrates g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return (nitrates/self.mstruct)

    def calculate_conc_amino_acids(self, amino_acids):
        """Amino_acid concentration.

        :Parameters:
            - `amino_acids` (:class:`float`) - Amount of amino acids (µmol N)
        :Returns:
            Amino_acid concentration (µmol amino_acids g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return (amino_acids/Organ.PARAMETERS.AMINO_ACIDS_N_RATIO)/self.mstruct

    def calculate_conc_cytokinines(self, cytokinines):
        """Cytokinines concentration.

        :Parameters:
            - `cytokinines` (:class:`float`) - Amount of cytokinines (UA)
        :Returns:
            cytokinines concentration (UA cytokinines g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return cytokinines/self.mstruct

    def calculate_total_organic_nitrogen(self, amino_acids, Nstruct):
        """Total amount of organic N (amino acids + Nstruct).
        Used to calculate residual respiration.

        :Parameters:
            - `amino_acids` (:class:`float`) - Amount of amino acids (µmol N)
            - `Nstruct` (:class:`float`) -Structural N mass (g)
        :Returns:
            Total amount of organic N (µmol N)
        :Returns Type:
            :class:`float`
        """
        return amino_acids + (Nstruct / Roots.PARAMETERS.N_MOLAR_MASS)*1E6

    # FLUXES

    def calculate_unloading_sucrose(self, sucrose_phloem):
        """Rate of sucrose unloading from phloem to roots integrated over DELTA_T (µmol C sucrose unloaded g-1 mstruct s-1 * DELTA_T).
        Michaelis-Menten function of the sucrose concentration in phloem.

        :Parameters:
            - `sucrose_phloem` (:class:`float`) - Sucrose concentration in phloem (µmol C g-1 mstruct)
        :Returns:
            Sucrose unloading (µmol C g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return (((max(0, sucrose_phloem)/(Organ.PARAMETERS.MSTRUCT_AXIS*Organ.PARAMETERS.ALPHA_AXIS)) * Roots.PARAMETERS.VMAX_SUCROSE_UNLOADING) / ((max(0, sucrose_phloem)/(Organ.PARAMETERS.MSTRUCT_AXIS*Organ.PARAMETERS.ALPHA_AXIS)) + Roots.PARAMETERS.K_SUCROSE_UNLOADING)) * Organ.PARAMETERS.DELTA_T

    def calculate_unloading_amino_acids(self, unloading_sucrose, sucrose_phloem, amino_acids_phloem):
        """Unloading of amino_acids from phloem to roots integrated over DELTA_T.
        Amino acids are assumed to be co-transported along with the unloaded sucrose from phloem (using the ratio amino acids:sucrose of phloem).

        :Parameters:
            - `unloading_sucrose` (:class:`float`) - Sucrose unloading (µmol C g-1 mstruct)
            - `sucrose_phloem` (:class:`float`) - Sucrose concentration in phloem (µmol C g-1 mstruct)
            - `amino_acids_phloem` (:class:`float`) - Amino acids concentration in phloem (µmol N g-1 mstruct)
        :Returns:
            Amino acids unloading (µmol N g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        if amino_acids_phloem <= 0:
            unloading_amino_acids = 0
        else:
            unloading_amino_acids =  unloading_sucrose * (amino_acids_phloem / sucrose_phloem)
        return unloading_amino_acids

    def calculate_uptake_nitrates(self, conc_nitrates_soil, nitrates_roots, total_surfacic_transpiration, sucrose_roots):
        """Rate of nitrate uptake by roots integrated over DELTA_T.
            - Nitrate uptake is calculated as the sum of the 2 transport systems: HATS and LATS
            - HATS and LATS parameters are calculated as a function of root nitrate concentration (negative regulation)
            - Nitrate uptake is finally regulated by the total culm transpiration and sucrose concentration (positive regulation)

        :Parameters:
            - `conc_nitrates_soil` (:class:`float`) - Soil nitrate concentration unloading (µmol N m-3 soil)
            - `nitrates_roots` (:class:`float`) - Amount of nitrates in roots (µmol N)
            - `total_surfacic_transpiration` (:class:`float`) - Total culm transpiration (mmol H2O m-2 s-1)
            - `sucrose_roots` (:class:`float`) - Amount of sucrose in roots (µmol C)
        :Returns:
            Nitrate uptake (µmol N nitrates)
        :Returns Type:
            :class:`float`
        """
        #: High Affinity Transport System (HATS)
        VMAX_HATS_MAX = Roots.PARAMETERS.A_VMAX_HATS * np.exp(-Roots.PARAMETERS.LAMBDA_VMAX_HATS*(nitrates_roots/self.mstruct))        #: Maximal rate of nitrates influx at saturating soil N concentration;HATS (µmol N nitrates g-1 mstruct s-1)
        K_HATS = Roots.PARAMETERS.A_K_HATS * np.exp(-Roots.PARAMETERS.LAMBDA_K_HATS*(nitrates_roots/self.mstruct))                     #: Affinity coefficient of nitrates influx at saturating soil N concentration;HATS (µmol m-3)
        HATS = (VMAX_HATS_MAX * conc_nitrates_soil)/ (K_HATS + conc_nitrates_soil)                                                     #: Rate of nitrate influx by HATS (µmol N nitrates uptaked s-1 g-1 mstruct)

        #: Low Affinity Transport System (LATS)
        K_LATS = Roots.PARAMETERS.A_LATS * np.exp(-Roots.PARAMETERS.LAMBDA_LATS*(nitrates_roots/self.mstruct))                         #: Rate constant for nitrates influx at low soil N concentration; LATS (m3 g-1 mstruct s-1)
        LATS = (K_LATS * conc_nitrates_soil)                                                                                           #: Rate of nitrate influx by LATS (µmol N nitrates uptaked s-1 g-1 mstruct)
        HATS_LATS = (HATS + LATS) * self.mstruct * Organ.PARAMETERS.DELTA_T                                                            #: Nitrate influx (µmol N)

        # Regulations
        regul_transpiration = (total_surfacic_transpiration/(total_surfacic_transpiration + Roots.PARAMETERS.K_TRANSPIRATION))         #: Nitrate uptake regulation by plant transpiration
        regul_C = (sucrose_roots/self.mstruct) / ((sucrose_roots/self.mstruct) + 2000)                                                  #: Nitrate uptake regulation by root C

        net_nitrate_uptake = HATS_LATS * Roots.PARAMETERS.NET_INFLUX_UPTAKE_RATIO * regul_transpiration * regul_C                      #: Net nitrate uptake (µmol N nitrates uptaked by roots)
        return net_nitrate_uptake, HATS_LATS, regul_transpiration

    def calculate_s_amino_acids(self, nitrates, sucrose):
        """Rate of amino acid synthesis in roots integrated over DELTA_T (µmol N amino acids g-1 mstruct s-1 * DELTA_T).
        Bi-substrate Michaelis-Menten function of nitrates and sucrose.

        :Parameters:
            - `nitrates` (:class:`float`) - Amount of nitrates in roots (µmol N)
            - `sucrose` (:class:`float`) - Amount of sucrose in roots (µmol C)
        :Returns:
            Amino acids synthesis (µmol N g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return Roots.PARAMETERS.VMAX_AMINO_ACIDS / ((1 + Roots.PARAMETERS.K_AMINO_ACIDS_NITRATES/(nitrates/(self.mstruct*Roots.PARAMETERS.ALPHA))) * (1 + Roots.PARAMETERS.K_AMINO_ACIDS_SUCROSE/(sucrose/(self.mstruct*Roots.PARAMETERS.ALPHA)))) * Organ.PARAMETERS.DELTA_T

    def calculate_export_nitrates(self, nitrates, total_surfacic_transpiration):
        """Total export of nitrates from roots to shoot organs integrated over DELTA_T.
        Export is calculated as a function on nitrate concentration and culm transpiration.

        :Parameters:
            - `nitrates` (:class:`float`) - Amount of nitrates in roots (µmol N)
            - `total_surfacic_transpiration` (:class:`float`) - Total culm transpiration (mmol H2O m-2 s-1)
        :Returns:
            Export of nitrates (µmol N)
        :Returns Type:
            :class:`float`
        """
        if nitrates <=0 or total_surfacic_transpiration<=0:
            export_nitrates = 0
        else:
            f_nitrates = (nitrates/(self.mstruct*Roots.PARAMETERS.ALPHA))*Roots.PARAMETERS.K_NITRATE_EXPORT
            regul_transpiration = (total_surfacic_transpiration/(total_surfacic_transpiration + Roots.PARAMETERS.K_TRANSPIRATION)) #: Nitrate export regulation by plant transpiration
            export_nitrates = f_nitrates* self.mstruct * regul_transpiration * Organ.PARAMETERS.DELTA_T                                 #: Actual nitrate export (µmol N)
        return export_nitrates

    def calculate_export_amino_acids(self, amino_acids, total_surfacic_transpiration):
        """Total export of amino acids from roots to shoot organs integrated over DELTA_T.
        Amino acids export is calculated as a function of nitrate export using the ratio amino acids:nitrates in roots.

        :Parameters:
            - `amino_acids` (:class:`float`) - Amount of amino acids in roots (µmol N)
            - `total_surfacic_transpiration` (:class:`float`) - Total culm transpiration (mmol H2O m-2 s-1)

        :Returns:
            Export of amino acids (µmol N)
        :Returns Type:
            :class:`float`
        """
        if amino_acids <=0:
            export_amino_acids = 0
        else:
            f_amino_acids = (amino_acids/(self.mstruct*Roots.PARAMETERS.ALPHA))*Roots.PARAMETERS.K_AMINO_ACIDS_EXPORT
            regul_transpiration = (total_surfacic_transpiration/(total_surfacic_transpiration + Roots.PARAMETERS.K_TRANSPIRATION)) #: Amino acids export regulation by plant transpiration
            export_amino_acids = f_amino_acids* self.mstruct * regul_transpiration * Organ.PARAMETERS.DELTA_T                           #: Actual nitrate export (µmol N)

        return export_amino_acids

    def calculate_exudation(self, unloading_sucrose, sucrose_roots, sucrose_phloem, amino_acids_roots, amino_acids_phloem):
        """C sucrose and N amino acids lost by root exudation integrated over DELTA_T (µmol C or N g-1 mstruct).
            - C exudation is calculated as a fraction of C unloading from phloem
            - N exudation is calculated from C exudation using the ratio amino acids:sucrose of the phloem

        :Parameters:
            - `unloading_sucrose` (:class:`float`) - Sucrose unloading (µmol C g-1 mstruct)
            - `sucrose_roots` (:class:`float`) - Amount of sucrose in roots (µmol C)
            - `sucrose_phloem` (:class:`float`) - Amount of sucrose in phloem (µmol C)
            - `amino_acids_roots` (:class:`float`) - Amount of amino acids in roots (µmol N)
            - `amino_acids_phloem` (:class:`float`) - Amount of amino acids in phloem (µmol N)
        :Returns:
            C exudated, N_exudated (µmol C g-1 mstruct, µmol N g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        C_exudated = min(sucrose_roots, unloading_sucrose * Roots.PARAMETERS.C_EXUDATION)           #: C exudated (µmol g-1 mstruct)
        if amino_acids_phloem <= 0 or amino_acids_roots <= 0:
            N_exudated = 0
        else:
            N_exudated = min(amino_acids_roots, C_exudated * (amino_acids_phloem / sucrose_phloem )) #: N exudated (µmol g-1 mstruct)
        return C_exudated, N_exudated

    def calculate_s_cytokinines(self, sucrose_roots):
        """ Rate of cytokinines synthesis integrated over DELTA_T (UA cytokinines g-1 mstruct s-1 * DELTA_T).
        Michaelis-Menten function of sucrose concentration in roots. As a signal molecule, cytokinines are assumed have a neglected effect on sucrose. Thus, no cost in C is applied to the sucrose pool.

        :Parameters:
            - `sucrose_roots` (:class:`float`) - Amount of sucrose in roots (µmol C)
        :Returns:
            Cytokinines synthesis (UA g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return (Roots.PARAMETERS.VMAX_S_CYTOKININES * max(0, (sucrose_roots/self.mstruct))) / (Roots.PARAMETERS.K_S_CYTOKININES + max(0, (sucrose_roots/self.mstruct))) * Organ.PARAMETERS.DELTA_T

    def calculate_export_cytokinines(self, cytokinines, total_surfacic_transpiration):
        """Total export of cytokinines from roots to shoot organs integrated over DELTA_T.
        Cytokinines export is calculated as a function of cytokinines concentration and culm transpiration.

        :Parameters:
            - `cytokinines` (:class:`float`) - Amount of cytokinines in roots (UA)
            - `total_surfacic_transpiration` (:class:`float`) - Total culm transpiration (mmol H2O m-2 s-1)
        :Returns:
            Cytokinines export (UA)
        :Returns Type:
            :class:`float`
        """
        if cytokinines <=0:
            export_cytokinines = 0
        else:
            f_cytokinines = (cytokinines/(self.mstruct*Roots.PARAMETERS.ALPHA)) * Roots.PARAMETERS.K_CYTOKININES_EXPORT                   #: Potential export of cytokinines (UA g-1 mstruct s-1)
            regul_transpiration = (total_surfacic_transpiration/(total_surfacic_transpiration + Roots.PARAMETERS.K_TRANSPIRATION))        #: Cytokinines export regulation by plant transpiration
            export_cytokinines = f_cytokinines* self.mstruct * regul_transpiration * Organ.PARAMETERS.DELTA_T                             #: Actual cytokinines export (UA)

        return export_cytokinines

    def calculate_d_cytokinines(self, cytokinines):
        """Rate of cytokinines degradation integrated over DELTA_T (UA g-1 mstruct s-1 * DELTA_T).
        Degradation calculated as a first order kinetic.

        :Parameters:
            - `cytokinines` (:class:`float`) - Amount of cytokinines in roots (UA)
        :Returns:
            Cytokinines degradation (UA g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return max(0, Roots.PARAMETERS.DELTA_D_CYTOKININES * (cytokinines/self.mstruct)) * PhotosyntheticOrgan.PARAMETERS.DELTA_T

    # COMPARTMENTS

    def calculate_sucrose_derivative(self, unloading_sucrose, s_amino_acids, mstruct_C_growth, C_exudated, sum_respi):
        """delta root sucrose.

        :Parameters:
            - `unloading_sucrose` (:class:`float`) - Sucrose unloading (µmol C g-1 mstruct)
            - `s_amino_acids` (:class:`float`) - Amino acids synthesis (µmol N g-1 mstruct)
            - `mstruct_C_growth` (:class:`float`) - Growth of root structural dry mass (µmol C)
            - `C_exudated` (:class:`float`) - C exudation (µmol C g-1 mstruct)
            - `sum_respi` (:class:`float`) - Sum of respirations for roots i.e. related to N uptake, amino acids synthesis, mstruct growth and residual (µmol C)
        :Returns:
            delta root sucrose (µmol C sucrose)
        :Returns Type:
            :class:`float`
        """
        sucrose_consumption_AA = (s_amino_acids / Organ.PARAMETERS.AMINO_ACIDS_N_RATIO) * Organ.PARAMETERS.AMINO_ACIDS_C_RATIO      #: Contribution of sucrose to the synthesis of amino_acids
        return (unloading_sucrose - sucrose_consumption_AA - C_exudated) * self.mstruct - sum_respi - mstruct_C_growth

    def calculate_nitrates_derivative(self, uptake_nitrates, export_nitrates, s_amino_acids):
        """delta root nitrates.

        :Parameters:
            - `uptake_nitrates` (:class:`float`) - Nitrate uptake (µmol N nitrates)
            - `export_nitrates` (:class:`float`) - Export of nitrates (µmol N)
            - `s_amino_acids` (:class:`float`) - Amino acids synthesis (µmol N g-1 mstruct)
        :Returns:
            delta root nitrates (µmol N nitrates)
        :Returns Type:
            :class:`float`
        """
        import_nitrates_roots = uptake_nitrates * (1-Organ.PARAMETERS.RATIO_EXPORT_NITRATES_ROOTS)                                  #: Proportion of uptaked nitrates staying in roots
        nitrate_reduction_AA = s_amino_acids                                                                                        #: Contribution of nitrates to the synthesis of amino_acids
        return import_nitrates_roots - export_nitrates - nitrate_reduction_AA*self.mstruct

    def calculate_amino_acids_derivative(self, unloading_amino_acids, s_amino_acids, export_amino_acids, Nstruct_N_growth, N_exudated):
        """delta root amino acids.

        :Parameters:
            - `unloading_amino_acids` (:class:`float`) - Amino acids unloading (µmol N g-1 mstruct)
            - `s_amino_acids` (:class:`float`) - Amino acids synthesis (µmol N g-1 mstruct)
            - `export_amino_acids` (:class:`float`) - Export of amino acids (µmol N)
            - `Nstruct_N_growth` (:class:`float`) - Growth of root structural N mass (µmol N)
            - `N_exudated` (:class:`float`) - N exudated (µmol g-1 mstruct)
        :Returns:
            delta root amino acids (µmol N amino acids)
        :Returns Type:
            :class:`float`
        """
        return (unloading_amino_acids + s_amino_acids- N_exudated)*self.mstruct - export_amino_acids - Nstruct_N_growth

    # VARIABLES

    def calculate_cytokinines_derivative(self, s_cytokinines, d_cytokinines, export_cytokinines):
        """delta root cytokinines.

        :Parameters:
            - `s_cytokinines` (:class:`float`) - Cytokinines synthesis (UA g-1 mstruct)
            - `d_cytokinines` (:class:`float`) - Cytokinines degradation (UA g-1 mstruct)
            - `export_cytokinines` (:class:`float`) - Cytokinines export (UA)
        :Returns:
            delta root cytokinines (UA cytokinines)
        :Returns Type:
            :class:`float`
        """
        return (s_cytokinines - d_cytokinines)*self.mstruct - export_cytokinines


class Soil(Organ):
    """
    The class :class:`Soil` defines the amount of nitrogen in the volume of soil explored by roots.
    """

    PARAMETERS = parameters.SoilParameters #: the internal parameters of the soil

    def __init__(self, volume=None, nitrates=None, Tsoil=None):

        # variables
        self.volume = volume                   #: Volume of soil explored by roots (m3)
        self.nitrates = nitrates               #: µmol N nitrates
        self.Tsoil = Tsoil                     #: soil temperature (°C)

    # VARIABLES

    def calculate_conc_nitrates(self, nitrates):
        """Nitrate concentration in soil.

        :Parameters:
            - `nitrates` (:class:`float`) - Amount of nitrates (µmol N)
        :Returns:
            Nitrate concentration (µmol nitrates m-3)
        :Returns Type:
            :class:`float`
        """
        return (nitrates/self.volume)

    # COMPARTMENTS

    def calculate_nitrates_derivative(self, uptake_nitrates):
        """delta soil nitrates.

        :Parameters:
            - `uptake_nitrates` (:class:`float`) - Nitrate uptake (µmol N nitrates)
        :Returns:
            delta nitrates (µmol N nitrates)
        :Returns Type:
            :class:`float`
        """
        return -uptake_nitrates


class PhotosyntheticOrgan(Organ):
    """
    The class :class:`PhotosyntheticOrgan` defines the CN exchanges in a photosynthetic organ.

    :class:`PhotosyntheticOrgan` is the base class of all photosynthetic organs. DO NOT INSTANTIATE IT.
    """

    PARAMETERS = parameters.PhotosyntheticOrganParameters #: the internal parameters of the photosynthetic organs

    def __init__(self, exposed_element=None, enclosed_element=None):
        self.exposed_element = exposed_element #: the exposed element
        self.enclosed_element = enclosed_element #: the enclosed element

    def calculate_integrative_variables(self):
        for element in (self.exposed_element, self.enclosed_element):
            if element is not None:
                element.calculate_integrative_variables()

    def calculate_total_green_area(self):
        """Calculate the sum of the element green area belonging to the organ.
        """
        total_green_area = 0
        for element in (self.exposed_element, self.enclosed_element):
            if element is not None:
                total_green_area += element.green_area
        return total_green_area


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

    def __init__(self, green_area=None, mstruct=None, Nstruct=None, triosesP=None, starch=None,
                 sucrose=None, fructan=None, nitrates=None, amino_acids=None, proteins=None, cytokinines=None,
                 Tr=None, Ag=None, Rd=None, Ts=None):

        # state parameters
        self.green_area = green_area         #: green area (m-2)
        self.Tr = Tr                         #: Transpiration rate (mmol m-2 s-1)
        self.Ag = Ag                         #: Gross assimilation (µmol m-2 s-1)
        self.Rd = Rd                         #: Mitochondrial respiration rate (µmol m-2 s-1)
        self.Ts = Ts                         #: Organ temperature (°C)

        # state variables
        self.mstruct = mstruct               #: Structural dry mass (g)
        self.Nstruct = Nstruct               #: Structural N mass (g)
        self.triosesP = triosesP             #: µmol C
        self.starch = starch                 #: µmol C
        self.sucrose = sucrose               #: µmol C
        self.fructan = fructan               #: µmol C
        self.nitrates = nitrates             #: µmol N
        self.amino_acids = amino_acids       #: µmol N
        self.proteins = proteins             #: µmol N
        self.cytokinines = cytokinines       #: UA

        # fluxes to phloem
        self.loading_sucrose = None          #: current rate of sucrose loading to phloem
        self.loading_amino_acids = None      #: current rate of amino acids loading to phloem

        # Integrated variables
        self.total_organic_nitrogen = None           #: current total nitrogen amount (µmol N)

    def calculate_integrative_variables(self):
        """Calculate the integrative variables of the element.
        """
        self.total_organic_nitrogen = self.calculate_total_organic_nitrogen(self.amino_acids, self.proteins, self.Nstruct)

    # VARIABLES
    def calculate_photosynthesis(self, Ag, green_area):
        """Total photosynthesis of an element integrated over DELTA_T (µmol C m-2 s-1 * m2 * DELTA_T).

        :Parameters:
            - `Ag` (:class:`float`) - Gross photosynthesis rate (µmol C m-2 s-1)
            - `green_area` (:class:`float`) - Green area (m2)
        :Returns:
            Total photosynthesis (µmol C)
        :Returns Type:
            :class:`float`
        """
        return Ag * green_area * PhotosyntheticOrgan.PARAMETERS.DELTA_T

    def calculate_total_Rd(self, Rd, green_area):
        """Total respiration of an element integrated over DELTA_T (µmol C m-2 s-1 * m2 * DELTA_T).

        :Parameters:
            - `Rd` (:class:`float`) - Respiration rate (µmol C m-2 s-1)
            - `green_area` (:class:`float`) - Green area (m2)
        :Returns:
            Total respiration (µmol C)
        :Returns Type:
            :class:`float`
        """
        return Rd * green_area * PhotosyntheticOrgan.PARAMETERS.DELTA_T

    def calculate_total_transpiration(self, Tr, green_area):
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

    def calculate_conc_triosesP(self, triosesP):
        """Triose Phosphates concentration.

        :Parameters:
            - `triosesP` (:class:`float`) - Amount of triose phosphates (µmol C)
        :Returns:
            Triose phosphates concentration (µmol triosesP g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return (triosesP/self.mstruct)/Organ.PARAMETERS.NB_C_TRIOSEP

    def calculate_conc_sucrose(self, sucrose):
        """Sucrose concentration.

        :Parameters:
            - `sucrose` (:class:`float`) - Amount of sucrose (µmol C)
        :Returns:
            Sucrose concentration (µmol sucrose g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return (sucrose/self.mstruct)/Organ.PARAMETERS.NB_C_SUCROSE

    def calculate_conc_starch(self, starch):
        """Starch concentration.

        :Parameters:
            - `starch` (:class:`float`) - Amount of sucrose (µmol C)
        :Returns:
            Starch concentration (µmol starch g-1 mstruct, eq. glucose).
        :Returns Type:
            :class:`float`
        """
        return (starch/self.mstruct)/Organ.PARAMETERS.NB_C_HEXOSES

    def calculate_conc_fructan(self, fructan):
        """Fructan concentration.

        :Parameters:
            - `fructan` (:class:`float`) - Amount of fructan (µmol C)
        :Returns:
            Fructan concentration (µmol fructan g-1 mstruct, eq. glucose).
        :Returns Type:
            :class:`float`
        """
        return (fructan/self.mstruct)/Organ.PARAMETERS.NB_C_HEXOSES

    def calculate_regul_s_fructan(self, loading_sucrose):
        """Regulating function for fructan synthesis.
        Negative regulation by the loading of sucrose from the phloem ("swith-off" sigmoïdal kinetic).

        :Parameters:
            - `loading_sucrose` (:class:`float`) - Sucrose loading (µmol C)
        :Returns:
            Regulating function for fructan synthesis (dimensionless, from 0 to 1)
        :Returns Type:
        """
        if loading_sucrose <=0:
            #regul_s_fructan = PhotosyntheticOrgan.PARAMETERS.VMAX_REGUL_SFRUCTAN
            Vmax_Sfructans = PhotosyntheticOrgan.PARAMETERS.VMAX_SFRUCTAN
        else:
            rate_loading_sucrose_massic = loading_sucrose/self.mstruct/PhotosyntheticOrgan.PARAMETERS.DELTA_T
            #regul_s_fructan = ((PhotosyntheticOrgan.PARAMETERS.VMAX_REGUL_SFRUCTAN * PhotosyntheticOrgan.PARAMETERS.K_REGUL_SFRUCTAN**(PhotosyntheticOrgan.PARAMETERS.N_REGUL_SFRUCTAN)) / (max(0, rate_loading_sucrose_massic**(PhotosyntheticOrgan.PARAMETERS.N_REGUL_SFRUCTAN)) + PhotosyntheticOrgan.PARAMETERS.K_REGUL_SFRUCTAN**(PhotosyntheticOrgan.PARAMETERS.N_REGUL_SFRUCTAN)))
            Vmax_Sfructans = ((PhotosyntheticOrgan.PARAMETERS.VMAX_SFRUCTAN * PhotosyntheticOrgan.PARAMETERS.K_REGUL_SFRUCTAN**(PhotosyntheticOrgan.PARAMETERS.N_REGUL_SFRUCTAN)) / (max(0, rate_loading_sucrose_massic**(PhotosyntheticOrgan.PARAMETERS.N_REGUL_SFRUCTAN)) + PhotosyntheticOrgan.PARAMETERS.K_REGUL_SFRUCTAN**(PhotosyntheticOrgan.PARAMETERS.N_REGUL_SFRUCTAN)))
        return Vmax_Sfructans#regul_s_fructan

    def calculate_conc_nitrates(self, nitrates):
        """Nitrate concentration.

        :Parameters:
            - `nitrates` (:class:`float`) - Amount of nitrates (µmol N)
        :Returns:
            Nitrate concentration (µmol nitrates g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return (nitrates/self.mstruct)

    def calculate_conc_amino_acids(self, amino_acids):
        """Amino_acid concentration.

        :Parameters:
            - `amino_acids` (:class:`float`) - Amount of amino acids (µmol N)
        :Returns:
            Amino_acid concentration (µmol amino acids g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return (amino_acids/PhotosyntheticOrgan.PARAMETERS.AMINO_ACIDS_N_RATIO) / self.mstruct

    def calculate_conc_proteins(self, proteins):
        """Protein concentration.

        :Parameters:
            - `proteins` (:class:`float`) - Amount of proteins (µmol N)
        :Returns:
            Protein concentration (g proteins g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        mass_N_proteins = proteins*1E-6 * PhotosyntheticOrgan.PARAMETERS.N_MOLAR_MASS                        #: Mass of N in proteins (g)
        mass_proteins = mass_N_proteins / PhotosyntheticOrgan.PARAMETERS.AMINO_ACIDS_MOLAR_MASS_N_RATIO      #: Total mass of proteins (g)
        return (mass_proteins / self.mstruct)

    def calculate_total_organic_nitrogen(self, amino_acids, proteins, Nstruct):
        """Total amount of organic N (amino acids + proteins + Nstruct).
        Used to calculate residual respiration.

        :Parameters:
            - `amino_acids` (:class:`float`) - Amount of amino acids (µmol N)
            - `proteins` (:class:`float`) - Amount of proteins (µmol N)
            - `Nstruct` (:class:`float`) - Structural N mass (g)
        :Returns:
            Total amount of organic N (µmol N)
        :Returns Type:
            :class:`float`
        """
        return amino_acids + proteins + (Nstruct / PhotosyntheticOrgan.PARAMETERS.N_MOLAR_MASS)*1E6

    def calculate_surfacic_nitrogen(self, nitrates, amino_acids, proteins, Nstruct, green_area):
        """Surfacic content of nitrogen.

        :Parameters:
            - `nitrates` (:class:`float`) - Amount of nitrates (µmol N)
            - `amino_acids` (:class:`float`) - Amount of amino acids (µmol N)
            - `proteins` (:class:`float`) - Amount of proteins (µmol N)
            - `Nstruct` (:class:`float`) -  Structural N mass (g)
            - `green_area` (:class:`float`) -  Green area (m2)
        :Returns:
            Surfacic content of nitrogen (g m-2)
        :Returns Type:
            :class:`float`
        """
        mass_N_tot = (nitrates + amino_acids + proteins)*1E-6 * PhotosyntheticOrgan.PARAMETERS.N_MOLAR_MASS + Nstruct
        return (mass_N_tot / green_area)

    def calculate_conc_cytokinines(self, cytokinines):
        """Cytokinines concentration.

        :Parameters:
            - `cytokinines` (:class:`float`) - Amount of cytokinines (UA)
        :Returns:
            Cytokinines concentration (UA g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return cytokinines/self.mstruct

    # FLUXES

    def calculate_s_starch(self, triosesP):
        """Rate of starch synthesis (µmol C starch g-1 mstruct s-1 * DELTA_T).
        Michaelis-Menten function of triose phosphates.

        :Parameters:
            - `triosesP` (:class:`float`) - Amount of triose phosphates (µmol C)
        :Returns:
            Starch synthesis (µmol C g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return (((max(0, triosesP)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) * PhotosyntheticOrgan.PARAMETERS.VMAX_STARCH) / ((max(0, triosesP)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) + PhotosyntheticOrgan.PARAMETERS.K_STARCH)) * PhotosyntheticOrgan.PARAMETERS.DELTA_T

    def calculate_d_starch(self, starch):
        """Rate of starch degradation (µmol C starch g-1 mstruct s-1 * DELTA_T).
        First order kinetic.

        :Parameters:
            - `starch` (:class:`float`) - Amount of starch (µmol C)
        :Returns:
            Starch degradation (µmol C g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return max(0, PhotosyntheticOrgan.PARAMETERS.DELTA_DSTARCH * (starch/(self.mstruct*self.__class__.PARAMETERS.ALPHA))) * PhotosyntheticOrgan.PARAMETERS.DELTA_T

    def calculate_s_sucrose(self, triosesP):
        """Rate of sucrose synthesis (µmol C sucrose g-1 mstruct s-1 * DELTA_T).
        Michaelis-Menten function of triose phosphates.

        :Parameters:
            - `triosesP` (:class:`float`) - Amount of triose phosphates (µmol C)
        :Returns:
            Sucrose synthesis (µmol C g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return (((max(0,triosesP)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) * PhotosyntheticOrgan.PARAMETERS.VMAX_SUCROSE) / ((max(0, triosesP)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) + PhotosyntheticOrgan.PARAMETERS.K_SUCROSE)) * PhotosyntheticOrgan.PARAMETERS.DELTA_T

    def calculate_loading_sucrose(self, sucrose, sucrose_phloem):
        """Rate of sucrose loading to phloem (µmol C sucrose s-1 * DELTA_T).
        Transport-resistance model.

        :Parameters:
            - `sucrose` (:class:`float`) - Amount of sucrose in the element (µmol C)
            - `sucrose_phloem` (:class:`float`) - Amount of sucrose in the phloem (µmol C)
        :Returns:
            Sucrose loading (µmol C)
        :Returns Type:
            :class:`float`
        """
        #: Driving compartment (µmol C g-1 mstruct)
        driving_sucrose_compartment = max(sucrose / (self.mstruct*self.__class__.PARAMETERS.ALPHA), sucrose_phloem/(PhotosyntheticOrgan.PARAMETERS.MSTRUCT_AXIS*parameters.OrganParameters.ALPHA_AXIS))
        #: Gradient of sucrose between the element and the phloem (µmol C g-1 mstruct)
        diff_sucrose = sucrose/(self.mstruct*self.__class__.PARAMETERS.ALPHA) - sucrose_phloem/(PhotosyntheticOrgan.PARAMETERS.MSTRUCT_AXIS*parameters.OrganParameters.ALPHA_AXIS)
        #: Conductance depending on mstruct (g2 µmol-1 s-1)
        conductance = PhotosyntheticOrgan.PARAMETERS.SIGMA_SUCROSE * PhotosyntheticOrgan.PARAMETERS.BETA * self.mstruct**(2/3)

        return driving_sucrose_compartment * diff_sucrose * conductance * PhotosyntheticOrgan.PARAMETERS.DELTA_T

    def calculate_s_fructan(self, sucrose, regul_s_fructan):
        """Rate of fructan synthesis (µmol C fructan g-1 mstruct s-1 * DELTA_T).
        Sigmoïdal function of sucrose.

        :Parameters:
            - `sucrose` (:class:`float`) - Amount of sucrose (µmol C)
            - `regul_s_fructan` (:class:`float`) - Regulating function for fructan synthesis (dimensionless)
        :Returns:
            Fructan synthesis (µmol C g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        #return (((max(0, sucrose)/(self.mstruct*self.__class__.PARAMETERS.ALPHA))**(PhotosyntheticOrgan.PARAMETERS.N_SFRUCTAN) * PhotosyntheticOrgan.PARAMETERS.VMAX_SFRUCTAN) / ((max(0, sucrose)/(self.mstruct*self.__class__.PARAMETERS.ALPHA))**(PhotosyntheticOrgan.PARAMETERS.N_SFRUCTAN) + PhotosyntheticOrgan.PARAMETERS.K_SFRUCTAN**(PhotosyntheticOrgan.PARAMETERS.N_SFRUCTAN))) * regul_s_fructan * PhotosyntheticOrgan.PARAMETERS.DELTA_T
        return ((max(0, sucrose)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) * regul_s_fructan)/((max(0, sucrose)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) + PhotosyntheticOrgan.PARAMETERS.K_SFRUCTAN) * PhotosyntheticOrgan.PARAMETERS.DELTA_T


    def calculate_d_fructan(self, sucrose, fructan):
        """Rate of fructan degradation (µmol C fructan g-1 mstruct s-1 * DELTA_T).
        Inhibition function by the end product i.e. sucrose (Bancal et al., 2012).

        :Parameters:
            - `sucrose` (:class:`float`) - Amount of sucrose (µmol C)
            - `fructan` (:class:`float`) - Amount of fructan (µmol C)
        :Returns:
            Fructan degradation (µmol C g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        d_potential = ((PhotosyntheticOrgan.PARAMETERS.K_DFRUCTAN * PhotosyntheticOrgan.PARAMETERS.VMAX_DFRUCTAN) / ((max(0, sucrose)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) + PhotosyntheticOrgan.PARAMETERS.K_DFRUCTAN)) * PhotosyntheticOrgan.PARAMETERS.DELTA_T
        d_actual = min(d_potential , max(0, fructan))
        return d_actual

    def calculate_nitrates_import(self, roots_uptake_nitrates, roots_export_nitrates, element_transpiration, total_transpiration):
        """Total nitrates imported from roots (µmol N nitrates integrated over DELTA_T).
        Nitrates coming from roots (fraction of uptake + direct export) are distributed according to the contribution of the element to culm transpiration.

        :Parameters:
            - `roots_uptake_nitrates` (:class:`float`) - Nitrate uptake by roots (µmol N)
            - `roots_export_nitrates` (:class:`float`) - Exported nitrates by roots (µmol N)
            - `element_transpiration` (:class:`float`) - Element transpiration (mmol H2O s-1)
            - `total_transpiration` (:class:`float`) - Culm transpiration (mmol H2O s-1)
        :Returns:
            Total nitrates import (µmol N nitrates)
        :Returns Type:
            :class:`float`
        """
        if total_transpiration>0:
            nitrates_import = roots_uptake_nitrates * PhotosyntheticOrgan.PARAMETERS.RATIO_EXPORT_NITRATES_ROOTS * (element_transpiration/total_transpiration) #: Proportion of uptaked nitrates transfered to element
            nitrates_import += roots_export_nitrates * (element_transpiration/total_transpiration)                                                             #: Proportion of exported nitrates from roots to element
        else: # Avoids further float division by zero error
            nitrates_import = 0
        return nitrates_import

    def calculate_amino_acids_import(self, roots_exported_amino_acids, element_transpiration, total_transpiration):
        """Total amino acids imported from roots  (µmol N amino acids integrated over DELTA_T).
        Amino acids exported by roots are distributed according to the contribution of the element to culm transpiration.

        :Parameters:
            - `roots_exported_amino_acids` (:class:`float`) - Exported amino acids by roots (µmol N)
            - `element_transpiration` (:class:`float`) - Element transpiration (mmol H2O s-1)
            - `total_transpiration` (:class:`float`) - Culm transpiration (mmol H2O s-1)
        :Returns:
            Total amino acids import (µmol N amino acids)
        :Returns Type:
            :class:`float`
        """
        if total_transpiration>0:
            amino_acids_import = roots_exported_amino_acids * (element_transpiration/total_transpiration) #: Proportion of exported amino acids from roots to organ
        else:
            amino_acids_import = 0
        return amino_acids_import

    def calculate_s_amino_acids(self, nitrates, triosesP):
        """Rate of amino acids synthesis (µmol N amino acids s-1 g-1 MS * DELTA_T).
        Bi-substrate Michaelis-Menten function of nitrates and triose phosphates.

        :Parameters:
            - `nitrates` (:class:`float`) - Amount of nitrates (µmol N)
            - `triosesP` (:class:`float`) - Amount of triosesP (µmol C)
        :Returns:
            Amino acids synthesis (µmol N g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        if nitrates <=0 or triosesP <=0:
            calculate_s_amino_acids = 0
        else:
            calculate_s_amino_acids = PhotosyntheticOrgan.PARAMETERS.VMAX_AMINO_ACIDS / ((1 + PhotosyntheticOrgan.PARAMETERS.K_AMINO_ACIDS_NITRATES/(nitrates/(self.mstruct*self.__class__.PARAMETERS.ALPHA))) * (1 + PhotosyntheticOrgan.PARAMETERS.K_AMINO_ACIDS_TRIOSESP/(triosesP/(self.mstruct*self.__class__.PARAMETERS.ALPHA)))) * PhotosyntheticOrgan.PARAMETERS.DELTA_T
        return calculate_s_amino_acids

    def calculate_s_proteins(self, amino_acids):
        """Rate of protein synthesis (µmol N proteins s-1 g-1 MS * DELTA_T).
        Michaelis-Menten function of amino acids.

        :Parameters:
            - `amino_acids` (:class:`float`) - Amount of amino acids (µmol N).
        :Returns:
            Protein synthesis (µmol N g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        calculate_s_proteins = (((max(0,amino_acids)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) * PhotosyntheticOrgan.PARAMETERS.VMAX_SPROTEINS) / ((max(0, amino_acids)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) + PhotosyntheticOrgan.PARAMETERS.K_SPROTEINS)) * PhotosyntheticOrgan.PARAMETERS.DELTA_T
        return calculate_s_proteins

    def calculate_d_proteins(self, proteins, cytokinines):
        """Rate of protein degradation (µmol N proteins s-1 g-1 MS * DELTA_T).
        First order kinetic regulated by cytokinines concentration.

        :Parameters:
            - `proteins` (:class:`float`) - Amount of proteins (µmol N)
            - `cytokinines` (:class:`float`) - Amount of cytokinines (UA)
        :Returns:
            Protein degradation (µmol N g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        conc_cytokinines = max(0,cytokinines/self.mstruct)
        k = (PhotosyntheticOrgan.PARAMETERS.VMAX_DPROTEINS * PhotosyntheticOrgan.PARAMETERS.K_DPROTEINS**PhotosyntheticOrgan.PARAMETERS.N_DPROTEINS) / (conc_cytokinines**PhotosyntheticOrgan.PARAMETERS.N_DPROTEINS + PhotosyntheticOrgan.PARAMETERS.K_DPROTEINS**PhotosyntheticOrgan.PARAMETERS.N_DPROTEINS)

        return k, max(0, k * (proteins/(self.mstruct*self.__class__.PARAMETERS.ALPHA))) * PhotosyntheticOrgan.PARAMETERS.DELTA_T

    def calculate_loading_amino_acids(self, amino_acids, amino_acids_phloem):
        """Rate of amino acids loading to phloem (µmol N amino acids s-1 * DELTA_T).
        Transport-resistance model.

        :Parameters:
            - `amino_acids` (:class:`float`) - Amount of amino acids in the element (µmol N)
            - `amino_acids_phloem` (:class:`float`) - Amount of amino acids in the phloem (µmol N)
        :Returns:
            Amino acids loading (µmol N)
        :Returns Type:
            :class:`float`
        """
        #: Driving compartment (µmol N g-1 mstruct)
        driving_amino_acids_compartment = max(amino_acids / (self.mstruct*self.__class__.PARAMETERS.ALPHA), amino_acids_phloem/(PhotosyntheticOrgan.PARAMETERS.MSTRUCT_AXIS*parameters.OrganParameters.ALPHA_AXIS))
        #: Gradient of amino acids between the element and the phloem (µmol N g-1 mstruct)
        diff_amino_acids = amino_acids/(self.mstruct*self.__class__.PARAMETERS.ALPHA) - amino_acids_phloem/(PhotosyntheticOrgan.PARAMETERS.MSTRUCT_AXIS*parameters.OrganParameters.ALPHA_AXIS)
        #: Conductance depending on mstruct (g2 µmol-1 s-1)
        conductance = PhotosyntheticOrgan.PARAMETERS.SIGMA_AMINO_ACIDS * PhotosyntheticOrgan.PARAMETERS.BETA * self.mstruct**(2/3)

        return driving_amino_acids_compartment * diff_amino_acids * conductance * PhotosyntheticOrgan.PARAMETERS.DELTA_T

    # test
    def calculate_cytokinines_import(self, roots_exported_cytokinines, element_transpiration, total_transpiration):
        """Import of cytokinines (UA cytokinines integrated over DELTA_T).
        Cytokinines exported by roots are distributed according to the contribution of the element to culm transpiration.

        :Parameters:
            - `roots_exported_cytokinines` (:class:`float`) - Exported cytokinines from roots (UA)
            - `element_transpiration` (:class:`float`) - Element transpiration (mmol H2O s-1)
            - `total_transpiration` (:class:`float`) - Culm transpiration (mmol H2O s-1)
        :Returns:
            Cytokinines import (UA)
        :Returns Type:
            :class:`float`
        """
        if total_transpiration>0:
            cytokinines_import = roots_exported_cytokinines * (element_transpiration/total_transpiration)
        else:
            cytokinines_import = 0
        return cytokinines_import

    def calculate_d_cytokinines(self, cytokinines):
        """Rate of cytokinines degradation integrated over DELTA_T (UA g-1 mstruct s-1 * DELTA_T).
        First order kinetic.

        :Parameters:
            - `cytokinines` (:class:`float`) - Amount of cytokinines (UA)
        :Returns:
            Cytokinines degradation (UA g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return max(0, PhotosyntheticOrgan.PARAMETERS.DELTA_D_CYTOKININES * (cytokinines/(self.mstruct*self.__class__.PARAMETERS.ALPHA))) * PhotosyntheticOrgan.PARAMETERS.DELTA_T

    # COMPARTMENTS

    def calculate_triosesP_derivative(self, photosynthesis, s_sucrose, s_starch, s_amino_acids):
        """ delta triose phosphates of element.

        :Parameters:
            - `photosynthesis` (:class:`float`) - Total gross photosynthesis (µmol C)
            - `s_sucrose` (:class:`float`) - Sucrose synthesis (µmol C g-1 mstruct)
            - `s_starch` (:class:`float`) - Starch synthesis (µmol C g-1 mstruct)
            - `s_amino_acids` (:class:`float`) - Amino acids synthesis (µmol N g-1 mstruct)
        :Returns:
            delta triose phosphates (µmol C triose phosphates)
        :Returns Type:
            :class:`float`
        """
        triosesP_consumption_AA = (s_amino_acids / PhotosyntheticOrgan.PARAMETERS.AMINO_ACIDS_N_RATIO) * PhotosyntheticOrgan.PARAMETERS.AMINO_ACIDS_C_RATIO #: Contribution of triosesP to the synthesis of amino_acids
        return photosynthesis - (s_sucrose + s_starch + triosesP_consumption_AA) * (self.mstruct*self.__class__.PARAMETERS.ALPHA)

    def calculate_starch_derivative(self, s_starch, d_starch):
        """delta starch of element.

        :Parameters:
            - `s_starch` (:class:`float`) - Starch synthesis (µmol C g-1 mstruct)
            - `d_starch` (:class:`float`) - Starch degradation (µmol C g-1 mstruct)
        :Returns:
            delta starch (µmol C starch)
        :Returns Type:
            :class:`float`
        """
        return (s_starch - d_starch) * (self.mstruct*self.__class__.PARAMETERS.ALPHA)

    def calculate_sucrose_derivative(self, s_sucrose, d_starch, loading_sucrose, s_fructan, d_fructan, sum_respi):
        """delta sucrose of element.

        :Parameters:
            - `s_sucrose` (:class:`float`) - Sucrose synthesis (µmol C g-1 mstruct)
            - `d_starch` (:class:`float`) - Starch degradation (µmol C g-1 mstruct)
            - `loading_sucrose` (:class:`float`) - Sucrose loading (µmol C)
            - `s_fructan` (:class:`float`) - Fructan synthesis (µmol C g-1 mstruct)
            - `d_fructan` (:class:`float`) - Fructan degradation (µmol C g-1 mstruct)
            - `sum_respi` (:class:`float`) - Sum of respirations for the element i.e. related to C loading to phloem, amino acids synthesis and residual (µmol C)
        :Returns:
            delta sucrose (µmol C sucrose)
        :Returns Type:
            :class:`float`
        """
        return (s_sucrose + d_starch + d_fructan - s_fructan) * (self.mstruct) - sum_respi - loading_sucrose

    def calculate_fructan_derivative(self, s_fructan, d_fructan):
        """delta fructan of element.

        :Parameters:
            - `s_fructan` (:class:`float`) - Fructan synthesis (µmol C g-1 mstruct)
            - `d_fructan` (:class:`float`) - Fructan degradation (µmol C g-1 mstruct)
        :Returns:
            delta fructan (µmol C fructan)
        :Returns Type:
            :class:`float`
        """
        return (s_fructan - d_fructan)* (self.mstruct*self.__class__.PARAMETERS.ALPHA)

    def calculate_nitrates_derivative(self, nitrates_import, s_amino_acids):
        """delta nitrates of element.

        :Parameters:
            - `nitrates_import` (:class:`float`) - Nitrate import from roots (µmol N)
            - `s_amino_acids` (:class:`float`) - Amino acids synthesis (µmol N g-1 mstruct)
        :Returns:
            delta nitrates (µmol N nitrates)
        :Returns Type:
            :class:`float`
        """
        nitrate_reduction_AA = s_amino_acids  #: Contribution of nitrates to the synthesis of amino_acids
        return nitrates_import - (nitrate_reduction_AA*self.mstruct*self.__class__.PARAMETERS.ALPHA)

    def calculate_amino_acids_derivative(self, amino_acids_import, s_amino_acids, s_proteins, d_proteins, loading_amino_acids):
        """delta amino acids of element.

        :Parameters:
            - `amino_acids_import` (:class:`float`) - Amino acids import from roots (µmol N)
            - `s_amino_acids` (:class:`float`) - Amino acids synthesis (µmol N g-1 mstruct)
            - `s_proteins` (:class:`float`) - Protein synthesis (µmol N g-1 mstruct)
            - `d_proteins` (:class:`float`) - Protein degradation (µmol N g-1 mstruct)
            - `loading_amino_acids` (:class:`float`) - Amino acids loading (µmol N)
        :Returns:
            delta amino acids (µmol N amino acids)
        :Returns Type:
            :class:`float`
        """
        return amino_acids_import - loading_amino_acids + (s_amino_acids + d_proteins - s_proteins) * (self.mstruct*self.__class__.PARAMETERS.ALPHA)

    def calculate_proteins_derivative(self, s_proteins, d_proteins):
        """delta proteins of element.

        :Parameters:
            - `s_proteins` (:class:`float`) - Protein synthesis (µmol N g-1 mstruct)
            - `d_proteins` (:class:`float`) - Protein degradation (µmol N g-1 mstruct)
        :Returns:
            delta proteins (µmol N proteins)
        :Returns Type:
            :class:`float`
        """
        return (s_proteins - d_proteins) * (self.mstruct*self.__class__.PARAMETERS.ALPHA)

    def calculate_cytokinines_derivative(self, import_cytokinines, d_cytokinines):
        """delta cytokinines of element.

        :Parameters:
            - `import_cytokinines` (:class:`float`) - Cytokinines import (UA)
            - `d_cytokinines` (:class:`float`) - Cytokinines degradation (UA g-1 mstruct)
        :Returns:
            delta cytokinines (UA cytokinines)
        :Returns Type:
            :class:`float`
        """
        return import_cytokinines - d_cytokinines * (self.mstruct*self.__class__.PARAMETERS.ALPHA)

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

