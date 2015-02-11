# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division
import numpy as np

import logging

import parameters

"""
    cnwheat.model
    ~~~~~~~~~~~~~

    The classes of the organs.

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

class Organ(object):
    """
    The Organ class defines the CN exchanges in an organ.

    The Organ class is the base class of all organs. DO NOT INSTANTIATE IT.
    """

    PARAMETERS = parameters.OrganParameters #: the parameters of the organ

    def __init__(self, t, name):
        self.t = t #: the thermal time
        self.name = name #: the name of the organ


class PhotosyntheticOrgan(Organ):
    """
    The PhotosyntheticOrgan class defines the CN exchanges in a photosynthetic organ.

    The PhotosyntheticOrgan class is the base class of all photosynthetic organs. DO NOT INSTANTIATE IT.
    """

    PARAMETERS = parameters.PhotosyntheticOrganParameters #: the parameters of the organ

    def __init__(self, t, area, mstruct, width, height, triosesP, starch,
                 sucrose, fructan, nitrates, amino_acids, proteins,  
                 An=None, Tr=None, name=''):

        super(PhotosyntheticOrgan, self).__init__(t, name)

        # variables
        self.area = area                     #: area (m-2)
        self.mstruct = mstruct               #: Structural mass (g)
        self.width = width                   #: Width (or diameter for stem organs) (m)
        self.height = height                 #: Height of organ from soil (m)
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
        

    # VARIABLES

    def calculate_photosynthesis(self, t, An):
        """Total photosynthesis of an organ integrated over DELTA_T (µmol CO2 on organ area integrated over delat_t)
        """
        return An * self._calculate_green_area(t) * Organ.PARAMETERS.DELTA_T

    def calculate_transpiration(self, t, Tr):
        """Total transpiration of an organ integrated over DELTA_T (mm of H2O on organ area integrated over delat_t)
        """
        return Tr * self._calculate_green_area(t) * Organ.PARAMETERS.DELTA_T

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
        return ((PhotosyntheticOrgan.PARAMETERS.VMAX_REGUL_SFRUCTAN * PhotosyntheticOrgan.PARAMETERS.K_REGUL_SFRUCTAN**(PhotosyntheticOrgan.PARAMETERS.N_REGUL_SFRUCTAN)) / (max(0, loading_sucrose**(PhotosyntheticOrgan.PARAMETERS.N_REGUL_SFRUCTAN)) + PhotosyntheticOrgan.PARAMETERS.K_REGUL_SFRUCTAN**(PhotosyntheticOrgan.PARAMETERS.N_REGUL_SFRUCTAN)))

    def calculate_conc_nitrates(self, nitrates):
        """Nitrate concentration (µmol nitrates g-1 MS)
        """
        return (nitrates/self.mstruct)

    def calculate_conc_amino_acids(self, amino_acids):
        """Amino_acid concentration (µmol amino acids g-1 MS)
        """
        return (amino_acids/Organ.PARAMETERS.AMINO_ACIDS_N_RATIO) / self.mstruct

    def calculate_conc_proteins(self, proteins):
        """Protein concentration (g proteins g-1 MS)
        """
        mass_N_proteins = proteins*1E-6 * Organ.PARAMETERS.N_MOLAR_MASS                        #: Mass of nitrogen in proteins (g)
        mass_proteins = mass_N_proteins / Organ.PARAMETERS.AMINO_ACIDS_MOLAR_MASS_N_RATIO      #: Total mass of proteins (g)
        return (mass_proteins / self.mstruct)

    # FLUXES

    def calculate_s_starch(self, triosesP):
        """Rate of starch synthesis from triosesP (µmol C starch s-1 g-1 MS * DELTA_T).
        """
        return (((max(0, triosesP)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) * PhotosyntheticOrgan.PARAMETERS.VMAX_STARCH) / ((max(0, triosesP)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) + PhotosyntheticOrgan.PARAMETERS.K_STARCH)) * Organ.PARAMETERS.DELTA_T

    def calculate_d_starch(self, starch):
        """Rate of starch degradation from triosesP (µmol C starch s-1 g-1 MS * DELTA_T).
        """
        return max(0, PhotosyntheticOrgan.PARAMETERS.DELTA_DSTARCH * (starch/(self.mstruct*self.__class__.PARAMETERS.ALPHA))) * Organ.PARAMETERS.DELTA_T

    def calculate_s_sucrose(self, triosesP):
        """Rate of sucrose synthesis from triosesP (µmol C sucrose s-1 g-1 MS * DELTA_T).
        """
        return (((max(0,triosesP)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) * PhotosyntheticOrgan.PARAMETERS.VMAX_SUCROSE) / ((max(0, triosesP)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) + PhotosyntheticOrgan.PARAMETERS.K_SUCROSE)) * Organ.PARAMETERS.DELTA_T

    def calculate_loading_sucrose(self, sucrose, sucrose_phloem):
        """Rate of sucrose loading to phloem (µmol C sucrose s-1 g-1 MS * DELTA_T).
        """
        driving_sucrose_compartment = max(sucrose / (self.mstruct*self.__class__.PARAMETERS.ALPHA), sucrose_phloem/(Organ.PARAMETERS.MSTRUCT_AXIS*self.__class__.PARAMETERS.ALPHA_AXIS))
        diff_sucrose = sucrose/(self.mstruct*self.__class__.PARAMETERS.ALPHA) - sucrose_phloem/(Organ.PARAMETERS.MSTRUCT_AXIS*self.__class__.PARAMETERS.ALPHA_AXIS)
        conductance = PhotosyntheticOrgan.PARAMETERS.SIGMA * PhotosyntheticOrgan.PARAMETERS.BETA * self.mstruct**(2/3)
        return driving_sucrose_compartment * diff_sucrose * conductance * Organ.PARAMETERS.DELTA_T

    def calculate_s_fructan(self, sucrose, regul_s_fructan):
        """Rate of fructan synthesis (µmol C fructan s-1 g-1 MS * DELTA_T)
        """
        return (((max(0, sucrose)/(self.mstruct*self.__class__.PARAMETERS.ALPHA))**(PhotosyntheticOrgan.PARAMETERS.N_SFRUCTAN) * PhotosyntheticOrgan.PARAMETERS.VMAX_SFRUCTAN) / ((max(0, sucrose)/(self.mstruct*self.__class__.PARAMETERS.ALPHA))**(PhotosyntheticOrgan.PARAMETERS.N_SFRUCTAN) + PhotosyntheticOrgan.PARAMETERS.K_SFRUCTAN**(PhotosyntheticOrgan.PARAMETERS.N_SFRUCTAN))) * regul_s_fructan * Organ.PARAMETERS.DELTA_T

    def calculate_d_fructan(self, sucrose, fructan):
        """Rate of fructan degradation (µmol C fructan s-1 g-1 MS)
        """
        return min((PhotosyntheticOrgan.PARAMETERS.K_DFRUCTAN * PhotosyntheticOrgan.PARAMETERS.VMAX_DFRUCTAN) / ((max(0, sucrose)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) + PhotosyntheticOrgan.PARAMETERS.K_DFRUCTAN) , max(0, fructan)) * Organ.PARAMETERS.DELTA_T

    def calculate_nitrates_import(self, roots_uptake_nitrate, organ_transpiration, total_transpiration):
        """Total nitrates imported from roots (through xylem) distributed relatively to organ transpiration (µmol N nitrates integrated over delta_t [already accounted in transpiration])
        """
        if total_transpiration>0:
            nitrates_import = roots_uptake_nitrate * (organ_transpiration/total_transpiration)* Organ.PARAMETERS.RATIO_EXPORT_NITRATES_ROOTS # Proportion of uptaked nitrates exported from roots to shoot
        else:
            nitrates_import = 0
        #print self.name, 'Ratio transpiration', nitrates_import
        return nitrates_import

    def calculate_amino_acids_import(self, roots_exported_amino_acids, organ_transpiration, total_transpiration):
        """Total Amino acids imported from roots (through xylem) distributed relatively to organ transpiration (µmol N Amino acids integrated over delta_t [already accounted in transpiration])
        """
        if total_transpiration>0:
            amino_acids_import = roots_exported_amino_acids * (organ_transpiration/total_transpiration)
        else:
            amino_acids_import = 0
        return amino_acids_import

    def calculate_s_amino_acids(self, nitrates, triosesP):
        """Rate of amino acid synthesis (µmol N amino acids s-1 g-1 MS * DELTA_T)
        """
        calculate_s_amino_acids = PhotosyntheticOrgan.PARAMETERS.VMAX_AMINO_ACIDS / ((1 + PhotosyntheticOrgan.PARAMETERS.K_AMINO_ACIDS_NITRATES/(nitrates/(self.mstruct*self.__class__.PARAMETERS.ALPHA))) * (1 + PhotosyntheticOrgan.PARAMETERS.K_AMINO_ACIDS_TRIOSESP/(triosesP/(self.mstruct*self.__class__.PARAMETERS.ALPHA)))) * Organ.PARAMETERS.DELTA_T
        #print calculate_s_amino_acids
        return calculate_s_amino_acids

    def calculate_s_proteins(self, amino_acids):
        """Rate of protein synthesis (µmol N proteins s-1 g-1 MS * DELTA_T)
        """
        calculate_s_proteins = (((max(0,amino_acids)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) * PhotosyntheticOrgan.PARAMETERS.VMAX_SPROTEINS) / ((max(0, amino_acids)/(self.mstruct*self.__class__.PARAMETERS.ALPHA)) + PhotosyntheticOrgan.PARAMETERS.K_SPROTEINS)) * Organ.PARAMETERS.DELTA_T
        #print calculate_s_proteins
        return calculate_s_proteins

    def calculate_d_proteins(self, proteins):
        """Rate of protein degradation (µmol N proteins s-1 g-1 MS * DELTA_T)
        """
        return max(0, PhotosyntheticOrgan.PARAMETERS.DELTA_DPROTEINS * (proteins/(self.mstruct*self.__class__.PARAMETERS.ALPHA))) * Organ.PARAMETERS.DELTA_T

    def calculate_loading_amino_acids(self, amino_acids, amino_acids_phloem):
        """Rate of amino acids loading to phloem (µmol N amino acids s-1 g-1 MS * DELTA_T) # TODO: formalism to be tested
        """
        driving_amino_acids_compartment = max(amino_acids / (self.mstruct*self.__class__.PARAMETERS.ALPHA), amino_acids_phloem/(Organ.PARAMETERS.MSTRUCT_AXIS*self.__class__.PARAMETERS.ALPHA_AXIS))
        diff_amino_acids = amino_acids/(self.mstruct*self.__class__.PARAMETERS.ALPHA) - amino_acids_phloem/(Organ.PARAMETERS.MSTRUCT_AXIS*self.__class__.PARAMETERS.ALPHA_AXIS)
        conductance = PhotosyntheticOrgan.PARAMETERS.SIGMA * PhotosyntheticOrgan.PARAMETERS.BETA * self.mstruct**(2/3)
        return driving_amino_acids_compartment * diff_amino_acids * conductance * Organ.PARAMETERS.DELTA_T

    # COMPARTMENTS

    def calculate_triosesP_derivative(self, photosynthesis, s_sucrose, s_starch, s_amino_acids):
        """ delta triosesP of organ integrated over delat_t (µmol C triosesP).
        """
        triosesP_consumption_AA = (s_amino_acids / Organ.PARAMETERS.AMINO_ACIDS_N_RATIO) * Organ.PARAMETERS.AMINO_ACIDS_C_RATIO #: Contribution of triosesP to the synthesis of amino_acids
        return max(0, photosynthesis) - (s_sucrose + s_starch + triosesP_consumption_AA) * (self.mstruct*self.__class__.PARAMETERS.ALPHA)

    def calculate_starch_derivative(self, s_starch, d_starch):
        """delta starch of organ integrated over delat_t (µmol C starch).
        """
        return (s_starch - d_starch) * (self.mstruct*self.__class__.PARAMETERS.ALPHA)

    def calculate_sucrose_derivative(self, s_sucrose, d_starch, loading_sucrose, s_fructan, d_fructan):
        """delta sucrose of organ integrated over delat_t (µmol C sucrose)
        """
        return (s_sucrose + d_starch + d_fructan - s_fructan - loading_sucrose) * (self.mstruct*self.__class__.PARAMETERS.ALPHA)

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
##        if self. name == 'lamina1':
##            print 'AA deriv+', amino_acids_import, s_amino_acids, d_proteins,'AA deriv-',s_proteins, loading_amino_acids
        return amino_acids_import + (s_amino_acids + d_proteins - s_proteins - loading_amino_acids) * (self.mstruct*self.__class__.PARAMETERS.ALPHA)

    def calculate_proteins_derivative(self, s_proteins, d_proteins):
        """delta proteins integrated over delat_t (µmol N proteins)
        """
        return (s_proteins - d_proteins) * (self.mstruct*self.__class__.PARAMETERS.ALPHA)


class Chaff(PhotosyntheticOrgan):
    """
    The Chaff class defines the CN exchanges in a chaff.
    """

    PARAMETERS = parameters.ChaffParameters #: the parameters of the organ


class Lamina(PhotosyntheticOrgan):
    """
    The Lamina class defines the CN exchanges in a lamina.
    """

    PARAMETERS = parameters.LaminaParameters #: the parameters of the organ

    # VARIABLES

    def _calculate_green_area(self, t):
        """Compute green area of the organ.
        """
        t_inflexion, value_inflexion = Lamina.PARAMETERS.INFLEXION_POINTS.get(self.name, (float("inf"), None))
        if t <= t_inflexion: # Non-senescent lamina
            green_area = self.area
        else: # Senescent lamina
            green_area = ((-0.0721*t + value_inflexion)/10000)
        return green_area


class Internode(PhotosyntheticOrgan):
    """
    The Internode class defines the CN exchanges in an internode.
    """

    PARAMETERS = parameters.InternodeParameters #: the parameters of the organ


class Peduncle(PhotosyntheticOrgan):
    """
    The Peduncle class defines the CN exchanges in a peduncle.
    """

    PARAMETERS = parameters.PeduncleParameters #: the parameters of the organ


class Sheath(PhotosyntheticOrgan):
    """
    The Sheath class defines the CN exchanges in a sheath.
    """

    PARAMETERS = parameters.SheathParameters #: the parameters of the organ


class Phloem(Organ):
    """
    The Phloem class defines the CN exchanges in a phloem.
    """

    PARAMETERS = parameters.PhloemParameters #: the parameters of the organ

    def __init__(self, t, sucrose, amino_acids, name=''):
        super(Phloem, self).__init__(t, name)

        # variables
        self.sucrose = sucrose
        self.amino_acids = amino_acids
        

    # VARIABLES

    def calculate_conc_sucrose(self, sucrose):
        """sucrose concentration (µmol sucrose g-1 MS)
        """
        return (sucrose/Organ.PARAMETERS.MSTRUCT_AXIS)/12

    def calculate_conc_c_sucrose(self, sucrose):
        """sucrose concentration expressed in C (µmol C sucrose g-1 MS)
        """
        return sucrose/(Organ.PARAMETERS.MSTRUCT_AXIS)

    def calculate_conc_amino_acids(self, amino_acids):
        """Amino_acid concentration (µmol amino_acids g-1 MS)
        """
        return (amino_acids/Organ.PARAMETERS.AMINO_ACIDS_N_RATIO) / Organ.PARAMETERS.MSTRUCT_AXIS

    # COMPARTMENTS

    def calculate_sucrose_derivative(self, organs):
        """delta sucrose of phloem integrated over delat_t (µmol C sucrose)
        """
        sucrose_derivative = 0
        # compute the sucrose contribution of each organ
        logger = logging.getLogger(__name__)
        for organ in organs:
            if isinstance(organ, Phloem):
                continue
            if isinstance(organ, PhotosyntheticOrgan):
                flow = organ.loading_sucrose * organ.mstruct * organ.__class__.PARAMETERS.ALPHA
            elif isinstance(organ, Grains):
                flow = - (organ.s_grain_structure + (organ.s_grain_starch * ((organ.structure*1E-6) * Organ.PARAMETERS.C_MOLAR_MASS))) #: Conversion of structure from umol of C to g of C
                logger.debug('organ name = {} ; contrib to phloem sucrose struct = {}'.format(organ.name, organ.s_grain_structure))
                logger.debug('organ name = {} ; contrib to phloem sucrose starch = {}; mstruct={}'.format(organ.name, organ.s_grain_starch * ((organ.structure*1E-6) * Organ.PARAMETERS.C_MOLAR_MASS), organ.structure))
            elif isinstance(organ, Roots):
                flow = -(organ.unloading_sucrose * organ.mstruct * organ.__class__.PARAMETERS.ALPHA)

            logger.debug('organ name = {} ; contrib to phloem sucrose = {}'.format(organ.name, flow))
            sucrose_derivative += flow
        logger.debug('total phloem sucrose_derivative = {}'.format(sucrose_derivative))
        return sucrose_derivative

    def calculate_amino_acids_derivative(self, organs):
        """delta amino acids of phloem integrated over delat_t (µmol N amino acids)
        """
        amino_acids_derivative = 0
        # compute the amino acids contribution of each organ
        for organ in organs:
            if isinstance(organ, PhotosyntheticOrgan):
                amino_acids_derivative += organ.loading_amino_acids * organ.mstruct * organ.__class__.PARAMETERS.ALPHA
            elif isinstance(organ, Grains):
                amino_acids_derivative -= organ.s_proteins * ((organ.structure/1E6)*12) #: Conversion of structure from umol of C to g of C
            elif isinstance(organ, Roots):
                amino_acids_derivative -= organ.unloading_amino_acids * organ.mstruct * organ.__class__.PARAMETERS.ALPHA
        return amino_acids_derivative

class Grains(Organ):
    """
    The Grains class defines the CN exchanges in grains.
    """

    PARAMETERS = parameters.GrainsParameters #: the parameters of the organ

    def __init__(self, t, starch, structure, proteins, name=''):
        super(Grains, self).__init__(t, name)

        # variables
        self.starch = starch
        self.structure = structure
        self.proteins = proteins
        
        # fluxes from phloem
        self.s_grain_structure = None            #: current rate of grain structure synthesis
        self.s_grain_starch = None             #: current rate of grain starch C synthesis
        self.s_proteins = None                   #: current rate of grain protein synthesis
        

    # VARIABLES

    def calculate_dry_mass(self, structure, starch):
        """Grain total dry mass (g) # TODO: ajouter la masse des prot?
        """
        return ((structure + starch)*1E-6) * Organ.PARAMETERS.C_MOLAR_MASS

    def calculate_protein_mass(self, proteins):
        """Grain total protein mass                                                            # TODO trouver stoechiometrie proteines grains
        """
        mass_N_proteins = proteins*1E-6 * Organ.PARAMETERS.N_MOLAR_MASS                        #: Mass of nitrogen in proteins (g)
        #mass_proteins = mass_N_proteins / Organ.PARAMETERS.AMINO_ACIDS_MOLAR_MASS_N_RATIO     #: Total mass of proteins (g)
        return mass_N_proteins

    def calculate_RGR_structure(self, sucrose_phloem):
        """Relative Growth Rate of grain structure, regulated by phloem concentrations
        """
        return ((max(0, sucrose_phloem)/(Organ.PARAMETERS.MSTRUCT_AXIS*Organ.PARAMETERS.ALPHA_AXIS)) * Grains.PARAMETERS.VMAX_RGR) / ((max(0, sucrose_phloem)/(Organ.PARAMETERS.MSTRUCT_AXIS*Organ.PARAMETERS.ALPHA_AXIS)) + Grains.PARAMETERS.K_RGR)

    def calculate_unloading_sucrose(self, s_grain_structure, s_grain_starch, structure):
        """Unloading of sucrose from phloem to grains integrated over delta_t (µmol sucrose)
        """
        return (s_grain_structure + s_grain_starch * (structure*1E-6) * Organ.PARAMETERS.C_MOLAR_MASS)/12


    # FLUXES

    def calculate_s_grain_structure(self, t, prec_structure, RGR_structure):
        """Synthesis of grain structure integrated over delta_t (µmol C structure s-1 * DELTA_T). Rate regulated by phloem concentrations
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

    def calculate_s_proteins(self, s_grain_structure, s_grain_starch, amino_acids_phloem, sucrose_phloem, structure):
        """Synthesis of grain proteins over delta_t (µmol N proteins). Rate regulated by phloem concentrations and unloading. Co-transported with sucrose relatively to the ratio amino acids:sucrose in phloem
        """
        if sucrose_phloem >0:
            s_proteins = (s_grain_structure + s_grain_starch*((structure*1E-6) * Organ.PARAMETERS.C_MOLAR_MASS)) * (amino_acids_phloem / sucrose_phloem)
        else:
            s_proteins = 0
        return s_proteins

    # COMPARTMENTS

    def calculate_structure_derivative(self, s_grain_structure):
        """delta grain structure integrated over delat_t (µmol C structure)
        """
        return s_grain_structure * Grains.PARAMETERS.Y_GRAINS

    def calculate_starch_derivative(self, s_grain_starch, structure):
        """delta grain starch integrated over delat_t (µmol C starch)
        """
        return s_grain_starch * Grains.PARAMETERS.Y_GRAINS * ((structure*1E-6)*Organ.PARAMETERS.C_MOLAR_MASS) #: Conversion of grain structure from µmol of C to g of C

    def calculate_proteins_derivative(self, s_proteins):
        """delta grain proteins integrated over delat_t (µmol N proteins)
        """
        return s_proteins


class Roots(Organ):
    """
    The Roots class defines the CN exchanges in roots.
    """

    PARAMETERS = parameters.RootsParameters #: the parameters of the organ

    def __init__(self, t, mstruct, sucrose, nitrates, amino_acids, name=''):
        super(Roots, self).__init__(t, name)

        # variables
        self.mstruct = mstruct  #: Structural mass (g)
        self.sucrose = sucrose
        self.nitrates = nitrates
        self.amino_acids = amino_acids

        # fluxes from phloem
        self.unloading_sucrose = None          #: current unloading of sucrose from phloem to roots
        self.unloading_amino_acids = None      #: current unloading of amino acids from phloem to roots
        

    # VARIABLES

    def calculate_conc_sucrose(self, sucrose):
        """Sucrose concentration (µmol sucrose g-1 MS).
        This is a concentration output (expressed in amount of substance g-1 MS).
        """
        return (sucrose/self.mstruct)/12

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
        return (amino_acids/Organ.PARAMETERS.AMINO_ACIDS_N_RATIO)/self.mstruct

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

    def calculate_sucrose_derivative(self, unloading_sucrose, s_amino_acids):
        """delta root sucrose integrated over delat_t (µmol C sucrose)
        """
        sucrose_consumption_AA = (s_amino_acids / Organ.PARAMETERS.AMINO_ACIDS_N_RATIO) * Organ.PARAMETERS.AMINO_ACIDS_C_RATIO      #: Contribution of sucrose to the synthesis of amino_acids

        return (unloading_sucrose - sucrose_consumption_AA) * self.mstruct

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