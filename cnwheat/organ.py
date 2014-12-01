# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division
import numpy as np

"""
    cnwheat.organ
    ~~~~~~~~~~~~~

    The classes of the organs.

    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.
"""

"""
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
"""

class Organ(object):
    """
    Base class for any organ. DO NOT INSTANTIATE IT TO USE IT DIRECTLY.
    """

    MSTRUCT_AXIS = 2.08                     #: Structural mass  of a plant (g) (Bertheloot, 2011)
    ALPHA_AXIS = 1                          #: Proportion of the structural mass containing the substrates
    DELTA_T = 3600                          #: Timestep of the model (s)

    AMINO_ACIDS_C_RATIO = 3.67              #: Mean ratio C:N in the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)
    AMINO_ACIDS_N_RATIO = 1.17              #: Mean number of mol of N in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)
    AMINO_ACIDS_MOLAR_MASS_N_RATIO = 0.145  #: Mean contribution of N in amino acids mass
    N_MOLAR_MASS = 14                       #: Molar mass of nitrogen (g mol-1)
    RATIO_EXPORT_NITRATES_ROOTS = 0.75      #: Proportion of uptaked nitrates actually exported from roots to shoot (1-RATIO_EXPORT_NITRATES_ROOTS = part of nitrates staying in roots)


    def __init__(self, name):
        if name is None:
            name = self.__class__.__name__
        self.name = name #: the name of the organ
        self.initial_conditions = {} #: the initial value of each compartment of the organ


class PhotosyntheticOrgan(Organ):
    """
    Base class for any photosynthetic organ. DO NOT INSTANTIATE IT TO USE IT DIRECTLY.
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
    SIGMA = 1.85e-07                #: Conductivity of an organ-phloem pathway (g mol-1 m-2 s-1) ; used to compute the sucrose loaded to the phloem
    BETA = 1                        #: Kind of volumetric mass density at power -2/3 ((g m-3)**(-2/3)) ; used to compute the sucrose loaded to the phloem

    # Amino acids
    VMAX_AMINO_ACIDS = 0.001        #: Maximal rate of amino acid synthesis (µmol N s-1 g-1 MS)
    K_AMINO_ACIDS_NITRATES = 5      #: Affinity coefficient of amino acid synthesis from nitrates (µmol N g-1 MS)
    K_AMINO_ACIDS_TRIOSESP = 0.1    #: Affinity coefficient of amino acid synthesis from triosesP (µmol C g-1 MS)

    # Proteins
    VMAX_SPROTEINS = 1.E-12         #: Maximal rate of protein synthesis (µmol N s-1 g-1 MS)
    K_SPROTEINS = 20                #: Affinity coefficient of protein synthesis (µmol N g-1 MS)
    DELTA_DPROTEINS = 1.85E-6       #: Relative rate of protein degradation (s-1)


    def __init__(self, area, mstruct, width, height, PAR, triosesP_0, starch_0, sucrose_0, fructan_0,
                 nitrates_0, amino_acids_0, proteins_0, name=None):

        super(PhotosyntheticOrgan, self).__init__(name)

        # parameters
        self.area = area                     #: area (m-2)
        self.mstruct = mstruct               #: Structural mass (g)
        self.width = width                   #: Width (or diameter for stem organs) (m)
        self.height = height                 #: Height of organ from soil (m)
        self.PAR = PAR                       #: PAR. Must be a :class:`pandas.Series` which index is time in hours

        self.loading_sucrose = 0             #: current rate of sucrose loading to phloem
        self.loading_amino_acids = 0         #: current rate of amino acids loading to phloem

        # initialize the compartments
        self.initial_conditions = {'triosesP':triosesP_0, 'starch':starch_0, 'sucrose':sucrose_0 , 'fructan':fructan_0,
                                   'nitrates':nitrates_0 , 'amino_acids':amino_acids_0, 'proteins':proteins_0}

    # VARIABLES

    def calculate_photosynthesis(self, t, An):
        """Total photosynthesis of an organ integrated over DELTA_T (µmol CO2 on organ area integrated over delat_t)
        """
        return An * self._calculate_green_area(t) * Organ.DELTA_T

    def calculate_transpiration(self, t, Tr):
        """Total transpiration of an organ integrated over DELTA_T (mm of H2O on organ area integrated over delat_t)
        """
        return Tr * self._calculate_green_area(t) * Organ.DELTA_T

    def _calculate_green_area(self, t):
        """Compute green area of the organ.
        """
        return self.area

    def calculate_conc_triosesP(self, triosesP):
        """Trioses Phosphate concentration (µmol triosesP g-1 MS).
        This is a concentration output (expressed in amount of substance g-1 MS).
        """
        return (triosesP/self.mstruct)/3

    def calculate_conc_sucrose(self, sucrose):
        """Sucrose concentration (µmol sucrose g-1 MS).
        This is a concentration output (expressed in amount of substance g-1 MS).
        """
        return (sucrose/self.mstruct)/12

    def calculate_conc_starch(self, starch):
        """Starch concentration (µmol starch g-1 MS (eq glucose)).
        This is a concentration output (expressed in amount of substance g-1 MS).
        """
        return (starch/self.mstruct)/6

    def calculate_conc_fructan(self, fructan):
        """fructan concentration (µmol fructan g-1 MS (eq glucose))
        """
        return (fructan/self.mstruct)/6

    def calculate_regul_s_fructan(self, loading_sucrose):
        """Inhibition of fructan synthesis by the loading of sucrose to phloem
        """
        return ((PhotosyntheticOrgan.VMAX_REGUL_SFRUCTAN * PhotosyntheticOrgan.K_REGUL_SFRUCTAN**(PhotosyntheticOrgan.N_REGUL_SFRUCTAN)) / (max(0, loading_sucrose**(PhotosyntheticOrgan.N_REGUL_SFRUCTAN)) + PhotosyntheticOrgan.K_REGUL_SFRUCTAN**(PhotosyntheticOrgan.N_REGUL_SFRUCTAN)))

    def calculate_conc_nitrates(self, nitrates):
        """Nitrate concentration (µmol nitrates g-1 MS)
        """
        return (nitrates/self.mstruct)

    def calculate_conc_amino_acids(self, amino_acids):
        """Amino_acid concentration (µmol amino acids g-1 MS)
        """
        return (amino_acids/Organ.AMINO_ACIDS_N_RATIO) / self.mstruct

    def calculate_conc_proteins(self, proteins):
        """Protein concentration (g proteins g-1 MS)
        """
        mass_N_proteins = proteins*1E6 * Organ.N_MOLAR_MASS                         #: Mass of nitrogen in proteins (g)
        mass_proteins = mass_N_proteins / Organ.AMINO_ACIDS_MOLAR_MASS_N_RATIO      #: Total mass of proteins (g)
        return (mass_proteins / self.mstruct)

    # FLOWS

    def calculate_s_starch(self, triosesP):
        """Rate of starch synthesis from triosesP (µmol C starch s-1 g-1 MS * DELTA_T).
        """
        return (((max(0, triosesP)/(self.mstruct*self.__class__.ALPHA)) * PhotosyntheticOrgan.VMAX_STARCH) / ((max(0, triosesP)/(self.mstruct*self.__class__.ALPHA)) + PhotosyntheticOrgan.K_STARCH)) * Organ.DELTA_T

    def calculate_d_starch(self, starch):
        """Rate of starch degradation from triosesP (µmol C starch s-1 g-1 MS * DELTA_T).
        """
        return max(0, PhotosyntheticOrgan.DELTA_DSTARCH * (starch/(self.mstruct*self.__class__.ALPHA))) * Organ.DELTA_T

    def calculate_s_sucrose(self, triosesP):
        """Rate of sucrose synthesis from triosesP (µmol C sucrose s-1 g-1 MS * DELTA_T).
        """
        return (((max(0,triosesP)/(self.mstruct*self.__class__.ALPHA)) * PhotosyntheticOrgan.VMAX_SUCROSE) / ((max(0, triosesP)/(self.mstruct*self.__class__.ALPHA)) + PhotosyntheticOrgan.K_SUCROSE)) * Organ.DELTA_T

    def calculate_loading_sucrose(self, sucrose, sucrose_phloem):
        """Rate of sucrose loading to phloem (µmol C sucrose s-1 g-1 MS * DELTA_T).
        """
        driving_sucrose_compartment = max(sucrose / (self.mstruct*self.__class__.ALPHA), sucrose_phloem/(Organ.MSTRUCT_AXIS*self.__class__.ALPHA_AXIS))
        diff_sucrose = sucrose/(self.mstruct*self.__class__.ALPHA) - sucrose_phloem/(Organ.MSTRUCT_AXIS*self.__class__.ALPHA_AXIS)
        conductance = PhotosyntheticOrgan.SIGMA * PhotosyntheticOrgan.BETA * self.mstruct**(2/3)
        return driving_sucrose_compartment * diff_sucrose * conductance * Organ.DELTA_T

    def calculate_s_fructan(self, sucrose, regul_s_fructan):
        """Rate of fructan synthesis (µmol C fructan s-1 g-1 MS * DELTA_T)
        """
        return (((max(0, sucrose)/(self.mstruct*self.__class__.ALPHA))**(PhotosyntheticOrgan.N_SFRUCTAN) * PhotosyntheticOrgan.VMAX_SFRUCTAN) / ((max(0, sucrose)/(self.mstruct*self.__class__.ALPHA))**(PhotosyntheticOrgan.N_SFRUCTAN) + PhotosyntheticOrgan.K_SFRUCTAN**(PhotosyntheticOrgan.N_SFRUCTAN))) * regul_s_fructan * Organ.DELTA_T

    def calculate_d_fructan(self, sucrose, fructan):
        """Rate of fructan degradation (µmol C fructan s-1 g-1 MS)
        """
        return min((PhotosyntheticOrgan.K_DFRUCTAN * PhotosyntheticOrgan.VMAX_DFRUCTAN) / ((max(0, sucrose)/(self.mstruct*self.__class__.ALPHA)) + PhotosyntheticOrgan.K_DFRUCTAN) , max(0, fructan)) * Organ.DELTA_T

    def calculate_nitrates_import(self, roots_uptake_nitrate, organ_transpiration, total_transpiration):
        """Total nitrates imported from roots (through xylem) distributed relatively to organ transpiration (µmol N nitrates integrated over delta_t [already accounted in transpiration])
        """
        if total_transpiration>0:
            nitrates_import = roots_uptake_nitrate * (organ_transpiration/total_transpiration)* Organ.RATIO_EXPORT_NITRATES_ROOTS # Proportion of uptaked nitrates exported from roots to shoot
        else:
            nitrates_import = 0
        #print self.name, 'Ratio transpiration', nitrates_import
        return nitrates_import

    def calculate_amino_acids_import(self, roots_exported_amino_acids, organ_transpiration, total_transpiration):
        """Total Amino acids imported from roots (through xylem) distributed relatively to organ transpiration (µmol N Amino acids integrated over delta_t [already accounted in transpiration)
        """
        if total_transpiration>0:
            amino_acids_import = roots_exported_amino_acids * (organ_transpiration/total_transpiration)
        else:
            amino_acids_import = 0
        return amino_acids_import

    def calculate_s_amino_acids(self, nitrates, triosesP):
        """Rate of amino acid synthesis (µmol N amino acids s-1 g-1 MS * DELTA_T)
        """
        calculate_s_amino_acids = PhotosyntheticOrgan.VMAX_AMINO_ACIDS / ((1 + PhotosyntheticOrgan.K_AMINO_ACIDS_NITRATES/(nitrates/(self.mstruct*self.__class__.ALPHA))) * (1 + PhotosyntheticOrgan.K_AMINO_ACIDS_TRIOSESP/(triosesP/(self.mstruct*self.__class__.ALPHA)))) * Organ.DELTA_T
        #print calculate_s_amino_acids
        return calculate_s_amino_acids

    def calculate_s_proteins(self, amino_acids):
        """Rate of protein synthesis (µmol N proteins s-1 g-1 MS * DELTA_T)
        """
        calculate_s_proteins = (((max(0,amino_acids)/(self.mstruct*self.__class__.ALPHA)) * PhotosyntheticOrgan.VMAX_SPROTEINS) / ((max(0, amino_acids)/(self.mstruct*self.__class__.ALPHA)) + PhotosyntheticOrgan.K_SPROTEINS)) * Organ.DELTA_T
        #print calculate_s_proteins
        return calculate_s_proteins

    def calculate_d_proteins(self, proteins):
        """Rate of protein degradation (µmol N proteins s-1 g-1 MS * DELTA_T)
        """
        return max(0, PhotosyntheticOrgan.DELTA_DPROTEINS * (proteins/(self.mstruct*self.__class__.ALPHA))) * Organ.DELTA_T

    def calculate_loading_amino_acids(self, amino_acids, amino_acids_phloem):
        """Rate of amino acids loading to phloem (µmol N amino acids s-1 g-1 MS * DELTA_T) # TODO: formalism to be tested
        """
        driving_amino_acids_compartment = max(amino_acids / (self.mstruct*self.__class__.ALPHA), amino_acids_phloem/(Organ.MSTRUCT_AXIS*self.__class__.ALPHA_AXIS))
        diff_amino_acids = amino_acids/(self.mstruct*self.__class__.ALPHA) - amino_acids_phloem/(Organ.MSTRUCT_AXIS*self.__class__.ALPHA_AXIS)
        conductance = PhotosyntheticOrgan.SIGMA * PhotosyntheticOrgan.BETA * self.mstruct**(2/3)
        return driving_amino_acids_compartment * diff_amino_acids * conductance * Organ.DELTA_T

    # COMPARTMENTS

    def calculate_triosesP_derivative(self, photosynthesis, s_sucrose, s_starch, s_amino_acids):
        """ delta triosesP of organ integrated over delat_t (µmol C triosesP).
        """
        triosesP_consumption_AA = (s_amino_acids / Organ.AMINO_ACIDS_N_RATIO) * Organ.AMINO_ACIDS_C_RATIO #: Contribution of triosesP to the synthesis of amino_acids
        return max(0, photosynthesis) - (s_sucrose + s_starch + triosesP_consumption_AA) * (self.mstruct*self.__class__.ALPHA)

    def calculate_starch_derivative(self, s_starch, d_starch):
        """delta starch of organ integrated over delat_t (µmol C starch).
        """
        return (s_starch - d_starch) * (self.mstruct*self.__class__.ALPHA)

    def calculate_sucrose_derivative(self, s_sucrose, d_starch, loading_sucrose, s_fructan, d_fructan):
        """delta sucrose of organ integrated over delat_t (µmol C sucrose)
        """
        return (s_sucrose + d_starch + d_fructan - s_fructan - loading_sucrose) * (self.mstruct*self.__class__.ALPHA)

    def calculate_fructan_derivative(self, s_fructan, d_fructan):
        """delta fructan integrated over delat_t (µmol C fructan)
        """
        return (s_fructan - d_fructan)* (self.mstruct*self.__class__.ALPHA)

    def calculate_nitrates_derivative(self, nitrates_import, s_amino_acids):
        """delta nitrates integrated over delat_t (µmol N nitrates)
        """
        nitrate_reduction_AA = s_amino_acids  #: Contribution of nitrates to the synthesis of amino_acids
        return nitrates_import - (nitrate_reduction_AA*self.mstruct*self.__class__.ALPHA)

    def calculate_amino_acids_derivative(self, amino_acids_import, s_amino_acids, s_proteins, d_proteins, loading_amino_acids):
        """delta amino acids integrated over delat_t (µmol N amino acids)
        """
##        if self. name == 'lamina1':
##            print 'AA deriv+', amino_acids_import, s_amino_acids, d_proteins,'AA deriv-',s_proteins, loading_amino_acids
        return amino_acids_import + (s_amino_acids + d_proteins - s_proteins - loading_amino_acids) * (self.mstruct*self.__class__.ALPHA)

    def calculate_proteins_derivative(self, s_proteins, d_proteins):
        """delta proteins integrated over delat_t (µmol N proteins)
        """
        return (s_proteins - d_proteins) * (self.mstruct*self.__class__.ALPHA)


class Chaff(PhotosyntheticOrgan):
    """
    Class for organ chaff.
    """

    ALPHA = 1 #: Proportion of structural mass containing substrate


class Lamina(PhotosyntheticOrgan):
    """
    Class for organ lamina.
    """

    ALPHA = 1 #: Proportion of structural mass containing substrate

    #: Temporary estimation of lamina senescence ({'lamina_order': (time of senescence beginning (h), offset of the linear regression)})
    LAMINAE_INFLEXION_POINTS = {'lamina1': (600, 78.75),
                                'lamina2': (480, 68.61),
                                'lamina3': (360, 48.76)}

    # VARIABLES

    def _calculate_green_area(self, t):
        """Compute green area of the organ.
        """
        t_inflexion, value_inflexion = Lamina.LAMINAE_INFLEXION_POINTS.get(self.name, (float("inf"), None))
        if t <= t_inflexion: # Non-senescent lamina
            green_area = self.area
        else: # Senescent lamina
            green_area = ((-0.0721*t + value_inflexion)/10000)
        return green_area


class Internode(PhotosyntheticOrgan):
    """
    Class for organ internode.
    """

    ALPHA = 1 #: Proportion of structural mass containing substrate


class Peduncle(PhotosyntheticOrgan):
    """
    Class for organ peduncle.
    """

    ALPHA = 1 #: Proportion of structural mass containing substrate


class Sheath(PhotosyntheticOrgan):
    """
    Class for organ sheath.
    """

    ALPHA = 1 #: Proportion of structural mass containing substrate


class Phloem(Organ):
    """
    Class for organ phloem.
    """

    ALPHA = 1 #: Proportion of structural mass containing substrate

    def __init__(self, sucrose_0, amino_acids_0, name=None):
        super(Phloem, self).__init__(name)

        # initialize the compartment
        self.initial_conditions = {'sucrose': sucrose_0, 'amino_acids': amino_acids_0}

    # VARIABLES

    def calculate_conc_sucrose(self, sucrose):
        """sucrose concentration (µmol sucrose g-1 MS)
        """
        return (sucrose/Organ.MSTRUCT_AXIS)/12

    def calculate_conc_c_sucrose(self, sucrose):
        """sucrose concentration expressed in C (µmol C sucrose g-1 MS)
        """
        return sucrose/(Organ.MSTRUCT_AXIS)

    def calculate_conc_amino_acids(self, amino_acids):
        """Amino_acid concentration (µmol amino_acids g-1 MS)
        """
        return (amino_acids/Organ.AMINO_ACIDS_N_RATIO) / Organ.MSTRUCT_AXIS

    # COMPARTMENTS

    def calculate_sucrose_derivative(self, organs):
        """delta sucrose of phloem integrated over delat_t (µmol C sucrose)
        """
        sucrose_derivative = 0
        for organ_ in organs:
            if isinstance(organ_, PhotosyntheticOrgan):
                sucrose_derivative += organ_.loading_sucrose * organ_.mstruct * organ_.__class__.ALPHA
            elif isinstance(organ_, Grains):
                sucrose_derivative -= (organ_.s_grain_structure + (organ_.s_grain_starch * ((organ_.structure/1E6)*12))) #: Conversion of structure from umol of C to g of C
            elif isinstance(organ_, Roots):
                sucrose_derivative -= organ_.unloading_sucrose * organ_.mstruct * organ_.__class__.ALPHA
        return sucrose_derivative

    def calculate_amino_acids_derivative(self, organs):
        """delta amino acids of phloem integrated over delat_t (µmol N amino acids)
        """
        amino_acids_derivative = 0
        for organ_ in organs:
            if isinstance(organ_, PhotosyntheticOrgan):
                amino_acids_derivative += organ_.loading_amino_acids * organ_.mstruct * organ_.__class__.ALPHA
            elif isinstance(organ_, Grains):
                amino_acids_derivative -= (organ_.s_proteins * ((organ_.structure/1E6)*12)) #: Conversion of structure from umol of C to g of C
            elif isinstance(organ_, Roots):
                amino_acids_derivative -= organ_.unloading_amino_acids * organ_.mstruct * organ_.__class__.ALPHA
        return amino_acids_derivative

class Grains(Organ):
    """
    Class for organ grains.
    """

    ALPHA = 1 #: Proportion of structural mass containing substrate

    # Structure parameters
    VMAX_RGR = 1.9e-06      #: Maximal value of the Relative Growth Rate of grain structure (dimensionless)
    K_RGR = 300             #: Affinity coefficient of the Relative Growth Rate of grain structure (dimensionless)

    # Starch parameters
    VMAX_STARCH = 0.5       #: Maximal rate of grain filling of starch (µmol C s-1 g-1 MS)
    K_STARCH = 100          #: Affinity coefficient of grain filling of starch (µmol C g-1 MS)

    Y_GRAINS = 0.75         #: Proportion of C loaded from phloem actually used for grain structure and starch (1 - Y_GRAINS is a kind of growth respiration)
    FILLING_INIT = 360      #: Time (h) at which phloem loading switch from grain structure to accumulation of starch

    def __init__(self, starch_0, structure_0, proteins_0, name=None):
        super(Grains, self).__init__(name)

        # flow from phloem
        self.s_grain_structure = 0            #: current rate of grain structure synthesis
        self.s_grain_starch_0 = 0             #: current rate of grain starch C synthesis
        self.s_proteins = 0                   #: current rate of grain protein synthesis

        self.structure = 0

        # initialize the compartments
        self.initial_conditions = {'starch': starch_0, 'structure': structure_0, 'proteins': proteins_0}

    # VARIABLES

    def calculate_dry_mass(self, structure, starch):
        """Grain total dry mass (g) # TODO: ajouter la masse des prot?
        """
        return ((structure + starch)/1000000)*12

    def calculate_protein_mass(self, proteins, structure):
        """Grain total protein mass                                                 # TODO trouver stoechiometrie proteines grains
        """
        mass_N_proteins = proteins*1E6 * Organ.N_MOLAR_MASS                         #: Mass of nitrogen in proteins (g)
        mass_proteins = mass_N_proteins / Organ.AMINO_ACIDS_MOLAR_MASS_N_RATIO      #: Total mass of proteins (g)
        return (mass_proteins / structure)

    def calculate_RGR_structure(self, sucrose_phloem):
        """Relative Growth Rate of grain structure, regulated by phloem concentrations
        """
        return ((max(0, sucrose_phloem)/(Organ.MSTRUCT_AXIS*Organ.ALPHA_AXIS)) * Grains.VMAX_RGR) / ((max(0, sucrose_phloem)/(Organ.MSTRUCT_AXIS*Organ.ALPHA_AXIS)) + Grains.K_RGR)

    # FLOWS

    def calculate_s_grain_structure(self, t, prec_structure, RGR_structure):
        """Synthesis of grain structure integrated over delta_t (µmol C structure s-1 * DELTA_T). Rate regulated by phloem concentrations
        """
        if t<=Grains.FILLING_INIT: #: Grain enlargment
            s_grain_structure = prec_structure * RGR_structure * Organ.DELTA_T
        else:                      #: Grain filling
            s_grain_structure = 0
        return s_grain_structure

    def calculate_s_grain_starch(self, t, sucrose_phloem):
        """Synthesis of grain C starch integrated over delta_t (µmol C starch s-1 g-1 MS * DELTA_T). Rate regulated by phloem concentrations and unloading
        """
        if t<=Grains.FILLING_INIT: #: Grain enlargment
            s_grain_starch = 0
        else:                      #: Grain filling
            s_grain_starch = (((max(0, sucrose_phloem)/(Organ.MSTRUCT_AXIS*Organ.ALPHA_AXIS)) * Grains.VMAX_STARCH) / ((max(0, sucrose_phloem)/(Organ.MSTRUCT_AXIS*Organ.ALPHA_AXIS)) + Grains.K_STARCH)) * Organ.DELTA_T
        return s_grain_starch

    def calculate_s_proteins(self, s_grain_structure, s_grain_starch, amino_acids_phloem, sucrose_phloem, structure):
        """Synthesis of grain proteins over delta_t (µmol N proteins). Rate regulated by phloem concentrations and unloading. Co-transported with sucrose relatively to the ratio amino acids:sucrose in phloem
        """
        if sucrose_phloem >0:
            s_proteins = (s_grain_structure + s_grain_starch*((structure/1E6)*12)) * (amino_acids_phloem / sucrose_phloem)
        else:
            s_proteins = 0
        return s_proteins

    # COMPARTMENTS

    def calculate_structure_derivative(self, s_grain_structure):
        """delta grain structure integrated over delat_t (µmol C structure)
        """
        return s_grain_structure * Grains.Y_GRAINS

    def calculate_starch_derivative(self, s_grain_starch, structure):
        """delta grain starch integrated over delat_t (µmol C starch)
        """
        return s_grain_starch * Grains.Y_GRAINS * ((structure/1E6)*12) #: Conversion of grain structure from µmol of C to g of C

    def calculate_proteins_derivative(self, s_proteins):
        """delta grain proteins integrated over delat_t (µmol N proteins)
        """
        return s_proteins


class Roots(Organ):
    """
    Class for organ roots.
    """

    ALPHA = 1                                 #: Proportion of structural mass containing substrate

    VMAX_SUCROSE_UNLOADING = 0.015            #: Maximal rate of sucrose unloading from phloem to roots (µmol C sucrose s-1 g-1 MS)
    K_SUCROSE_UNLOADING = 100                 #: Affinity coefficient of sucrose unloading from phloem to roots (µmol C sucrose g-1 MS)
    VMAX_AMINO_ACIDS_UNLOADING = 0.015        #: Maximal rate of amino acids unloading from phloem to roots (µmol N amino acids s-1 g-1 MS)
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
    VMAX_AMINO_ACIDS = 1            #: Maximal rate of amino acid synthesis (µmol N s-1 g-1 MS)
    K_AMINO_ACIDS_NITRATES = 5      #: Affinity coefficient of amino acid synthesis from nitrates (µmol N g-1 MS)
    K_AMINO_ACIDS_SUCROSE = 0.01    #: Affinity coefficient of amino acid synthesis from triosesP (µmol C g-1 MS)
    K_TR_EXPORT_AMINO_ACIDS = 1E-2  #: Affinity coefficient of amino acids export from roots to shoot (mm H20)


    def __init__(self, mstruct, sucrose_0, nitrates_0, amino_acids_0, name=None):
        super(Roots, self).__init__(name)

        # parameters
        self.mstruct = mstruct  #: Structural mass (g)

        self.unloading_sucrose = 0          #: current unloading of sucrose from phloem to roots
        self.unloading_amino_acids = 0      #: current unloading of amino acids from phloem to roots

        # initialize the compartment
        self.initial_conditions = {'sucrose': sucrose_0, 'nitrates': nitrates_0, 'amino_acids': amino_acids_0}

    # VARIABLES

    def calculate_dry_mass(self, sucrose):
        """Dry mass of roots (g)
        """
        return (sucrose*12)/1000000

    def calculate_conc_nitrates_soil(self, t):
        """Nitrate concetration in soil (µmol nitrates m-3)
        """
        return -52083*t + 5E+07 # TODO: Temporary

    def calculate_conc_nitrates(self, nitrates):
        """Nitrate concentration (µmol nitrates g-1 MS)
        """
        return (nitrates/self.mstruct)

    def calculate_conc_amino_acids(self, amino_acids):
        """Amino_acid concentration (µmol amino_acids g-1 MS)
        """
        return (amino_acids/Organ.AMINO_ACIDS_N_RATIO)/self.mstruct

    # FLOWS

    def calculate_unloading_sucrose(self, sucrose_phloem):
        """Unloading of sucrose from phloem to roots (µmol C sucrose unloaded s-1 g-1 MS * DELTA_T)
        """
        return (((max(0, sucrose_phloem)/(Organ.MSTRUCT_AXIS*Organ.ALPHA_AXIS)) * Roots.VMAX_SUCROSE_UNLOADING) / ((max(0, sucrose_phloem)/(Organ.MSTRUCT_AXIS*Organ.ALPHA_AXIS)) + Roots.K_SUCROSE_UNLOADING)) * Organ.DELTA_T

    def calculate_unloading_amino_acids(self, amino_acids_phloem):
        """Unloading of amino_acids from phloem to roots over delta_t (µmol N amino_acids unloaded s-1 g-1 MS)
        """
        return (((max(0, amino_acids_phloem)/(Organ.MSTRUCT_AXIS*Organ.ALPHA_AXIS)) * Roots.VMAX_AMINO_ACIDS_UNLOADING) / ((max(0, amino_acids_phloem)/(Organ.MSTRUCT_AXIS*Organ.ALPHA_AXIS)) + Roots.K_AMINO_ACIDS_UNLOADING)) * Organ.DELTA_T # TODO: Temporary

    def calculate_uptake_nitrates(self, conc_nitrates_soil, nitrates_roots, total_transpiration):
        """Uptake of nitrates by roots (µmol N nitrates imported s-1 * DELTA_T)
        """
        VMAX_HATS_MAX = Roots.A_VMAX_HATS * np.exp(-Roots.LAMBDA_VMAX_HATS*(nitrates_roots/self.mstruct))        #: Maximal rate of nitrates uptake at saturating soil N concentration;HATS (µmol N nitrates g-1 s-1)
        #print 'VMAX_HATS_MAX', VMAX_HATS_MAX
        K_HATS = Roots.A_K_HATS * np.exp(-Roots.LAMBDA_K_HATS*(nitrates_roots/self.mstruct))                     #: Affinity coefficient of nitrates uptake at saturating soil N concentration;HATS (µmol m-3)
        #print 'K_HATS', K_HATS
        HATS = (VMAX_HATS_MAX * conc_nitrates_soil)/ (K_HATS + conc_nitrates_soil)                               #: High Affinity Transport System (µmol N nitrates uptaked s-1 g-1 MS roots)
        #print 'HATS', HATS
        #print 'conc_nitrates_soil', conc_nitrates_soil
        K_LATS = Roots.A_LATS * np.exp(-Roots.LAMBDA_LATS*(nitrates_roots/self.mstruct))                         #: Rate of nitrates uptake at low soil N concentration; LATS (m3 g-1 s-1)
        LATS = (K_LATS * conc_nitrates_soil)                                                                     #: Low Affinity Transport System (µmol N nitrates uptaked s-1 g-1 MS roots)

        potential_uptake = (HATS + LATS) * self.mstruct * Organ.DELTA_T                                          #: Potential nitrate uptake (µmol N nitrates uptaked by roots integrated over delta_t)
        #print 'potential_uptake', potential_uptake
        #print 'Tr', total_transpiration
        #print 'f(Tr)', total_transpiration/(total_transpiration + Roots.K_TR_UPTAKE_NITRATES)
        actual_uptake = potential_uptake * (total_transpiration/(total_transpiration + Roots.K_TR_UPTAKE_NITRATES)) #: Nitrate uptake regulated by plant transpiration (µmol N nitrates uptaked by roots)
        #print 'actual_uptake',actual_uptake
        return actual_uptake, potential_uptake

    def calculate_s_amino_acids(self, nitrates, sucrose):
        """Rate of amino acid synthesis in roots(µmol N amino acids s-1 g-1 MS * DELTA_T)
        """
        return Roots.VMAX_AMINO_ACIDS / ((1 + Roots.K_AMINO_ACIDS_NITRATES/(nitrates/(self.mstruct*Roots.ALPHA))) * (1 + Roots.K_AMINO_ACIDS_SUCROSE/(sucrose/(self.mstruct*Roots.ALPHA)))) * Organ.DELTA_T

    def calculate_export_amino_acids(self, amino_acids, total_transpiration):
        """Total export of amino acids from roots to shoot organs (abstraction of the xylem compartment) (µmol N amino acids exported during DELTA_T (already accounted in Transpiration))
        """
        return (amino_acids/(self.mstruct * Roots.ALPHA)) * (total_transpiration/(total_transpiration + Roots.K_TR_EXPORT_AMINO_ACIDS))

    # COMPARTMENTS

    def calculate_sucrose_derivative(self, unloading_sucrose, s_amino_acids):
        """delta root sucrose integrated over delat_t (µmol C sucrose)
        """
        sucrose_consumption_AA = (s_amino_acids / Organ.AMINO_ACIDS_N_RATIO) * Organ.AMINO_ACIDS_C_RATIO        #: Contribution of sucrose to the synthesis of amino_acids

        return (unloading_sucrose - sucrose_consumption_AA) * self.mstruct

    def calculate_nitrates_derivative(self, uptake_nitrates, s_amino_acids):
        """delta root nitrates integrated over delat_t (µmol N nitrates)
        """
        import_nitrates_roots = uptake_nitrates * (1-Organ.RATIO_EXPORT_NITRATES_ROOTS)                         #: Proportion of uptaked nitrates staying in roots
        nitrate_reduction_AA = s_amino_acids                                                                    #: Contribution of nitrates to the synthesis of amino_acids
        return import_nitrates_roots - (nitrate_reduction_AA*self.mstruct)

    def calculate_amino_acids_derivative(self, unloading_amino_acids, s_amino_acids, export_amino_acids):
        """delta root amino acids integrated over delat_t (µmol N amino acids)
        """
        #print 'unloading_amino_acids, s_amino_acids, export_amino_acids', unloading_amino_acids, s_amino_acids, export_amino_acids
        return (unloading_amino_acids + s_amino_acids)*self.mstruct  - export_amino_acids # TODO: verifier apres modif