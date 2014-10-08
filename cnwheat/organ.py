# -*- coding: latin-1 -*-
"""
    cnwheat.organ
    ~~~~~~~~~~~~~

    The classes of the organs.

    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.
"""

from __future__ import division # use "//" to do integer division


class Organ(object):

    MSTRUCT_AXIS = 2.08             #: Structural mass  of a plant (g) (Bertheloot, 2011)
    ALPHA_AXIS = 1                  #: Proportion of the structural mass containing the substrates
    DELTA_T = 3600                  #: Timestep of the model (s)

    def __init__(self, name):
        if name is None:
            name = self.__class__.__name__
        self.name = name


class PhotosyntheticOrgan(Organ):

    # sucrose
    VMAX_SUCROSE = 1                #: Maximal rate of sucrose synthesis (µmol C s-1 g-1 MS)
    K_SUCROSE = 0.66                #: Affinity coefficient of sucrose synthesis (µmol C g-1 MS)

    # storage
    VMAX_STORAGE = 2                #: Maximal rate of storage synthesis (µmol C s-1 g-1 MS)
    K_STORAGE = 20                  #: Affinity coefficient of storage synthesis (µmol C g-1 MS)
    DELTA_DSTORAGE = 0.0001         #: Rate of storage degradation (µmol C storage s-1 g-1 MS)

    # fructans
    VMAX_SFRUCTAN = 0.2             #: Maximal rate of fructan synthesis (µmol C s-1 g-1 MS)
    K_SFRUCTAN = 20000              #: Affinity coefficient of fructan synthesis (µmol C g-1 MS)
    N_SFRUCTAN = 3                  #: Number of "substrates" for fructan synthesis (dimensionless)
    VMAX_REGUL_SFRUCTAN = 1         #: Maximal value of the regulation function of fructan synthesis (dimensionless)
    K_REGUL_SFRUCTAN = 8            #: Affinity coefficient of the regulation function of fructan synthesis (dimensionless)
    N_REGUL_SFRUCTAN = 15           #: Parameter of the regulation function of fructan synthesis (dimensionless)
    VMAX_DFRUCTAN = 0.035           #: Maximal rate of fructan degradation (µmol C s-1 g-1 MS)
    K_DFRUCTAN = 100                #: Affinity coefficient of fructan degradation (µmol C g-1 MS)

    SIGMA = 1.85e-07                #: Conductivity of an organ-phloem pathway (g mol-1 m-2 s-1) ; used to compute the sucrose loaded to the phloem

    def __init__(self, area, mstruct, PAR, storage_0,
                 sucrose_0, triosesP_0, fructan_0, name=None):

        super(PhotosyntheticOrgan, self).__init__(name)

        # parameters
        self.area = area                    #: area (m-2)
        self.mstruct = mstruct              #: Structural mass (g)
        self.PAR = PAR                      #: the PAR. Must be a :class:`pandas.Series` which index is time in hours 

        self.PAR_linear_interpolation = None #: linear interpolation of PAR

        self.photosynthesis_mapping = {} #: mapping to store the computed photosynthesis

        self.loading_sucrose = 0 #: current rate of sucrose loading to phloem

        # initialize the compartments
        self.triosesP_0 = triosesP_0 #: initial value of compartment triosesP
        self.storage_0 = storage_0 #: initial value of compartment storage
        self.sucrose_0 = sucrose_0 #: initial value of compartment sucrose
        self.fructan_0 = fructan_0 #: initial value of compartment fructan

    def get_initial_conditions(self):
        return [self.storage_0, self.sucrose_0, self.triosesP_0, self.fructan_0] # keep this order ! TODO: pas le plus "logique" mais ok si trop complique a modifier

    # VARIABLES

    def calculate_photosynthesis(self, t, An):
        """Total photosynthesis of an organ integrated over DELTA_T (µmol CO2 on organ area integrated over delat_t)
        """
        return An * self._calculate_green_area(t) * Organ.DELTA_T

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
        """sucrose concentration (µmol sucrose g-1 MS).
        This is a concentration output (expressed in amount of substance g-1 MS).
        """
        return (sucrose/self.mstruct)/12

    def calculate_conc_storage(self, storage):
        """storage concentration (µmol storage g-1 MS (eq glucose)).
        This is a concentration output (expressed in amount of substance g-1 MS).
        """
        return (storage/self.mstruct)/6

    def calculate_conc_fructan(self, fructan):
        """fructan concentration (µmol fructan g-1 MS (eq glucose))
        """
        return (fructan/self.mstruct)/6

    def calculate_regul_s_fructan(self, loading_sucrose):
        """Inhibition of fructan synthesis by the loading of sucrose to phloem
        """
        return ((PhotosyntheticOrgan.VMAX_REGUL_SFRUCTAN * PhotosyntheticOrgan.K_REGUL_SFRUCTAN**(PhotosyntheticOrgan.N_REGUL_SFRUCTAN)) / (max(0, loading_sucrose**(PhotosyntheticOrgan.N_REGUL_SFRUCTAN)) + PhotosyntheticOrgan.K_REGUL_SFRUCTAN**(PhotosyntheticOrgan.N_REGUL_SFRUCTAN)))

    # FLOWS

    def calculate_s_sucrose(self, triosesP):
        """Rate of sucrose synthesis from triosesP (µmol C sucrose s-1 g-1 MS * DELTA_T).
        This is a flow (expressed in amount of C substance g-1 MS integrated over DELTA_T).
        """
        return (((max(0,triosesP)/(self.mstruct*self.ALPHA)) * PhotosyntheticOrgan.VMAX_SUCROSE) / ((max(0, triosesP)/(self.mstruct*self.ALPHA)) + PhotosyntheticOrgan.K_SUCROSE)) * Organ.DELTA_T

    def calculate_s_storage(self, triosesP):
        """Rate of storage synthesis from triosesP (µmol C storage s-1 g-1 MS * DELTA_T).
        This is a flow (expressed in amount of C substance g-1 MS integrated over DELTA_T).
        """
        return (((max(0, triosesP)/(self.mstruct*self.ALPHA)) * PhotosyntheticOrgan.VMAX_STORAGE) / ((max(0, triosesP)/(self.mstruct*self.ALPHA)) + PhotosyntheticOrgan.K_STORAGE)) * Organ.DELTA_T

    def calculate_d_storage(self, storage):
        """Rate of storage degradation from triosesP (µmol C storage s-1 g-1 MS * DELTA_T).
        This is a flow (expressed in amount of C substance g-1 MS integrated over DELTA_T).
        """
        return max(0, PhotosyntheticOrgan.DELTA_DSTORAGE * (storage/(self.mstruct*self.ALPHA))) * Organ.DELTA_T

    def calculate_loading_sucrose(self, sucrose, sucrose_phloem):
        """Rate of sucrose loading to phloem (µmol C sucrose s-1 g-1 MS * DELTA_T).
        This is a flow (expressed in amount of C substance g-1 MS integrated over DELTA_T).
        """
        driving_sucrose_compartment = max(sucrose / (self.mstruct*self.ALPHA), sucrose_phloem/(Organ.MSTRUCT_AXIS*self.ALPHA_AXIS))
        diff_sucrose = sucrose/(self.mstruct*self.ALPHA) - sucrose_phloem/(Organ.MSTRUCT_AXIS*self.ALPHA_AXIS)
        conductance = PhotosyntheticOrgan.SIGMA * self.mstruct**(2/3)
        return driving_sucrose_compartment * diff_sucrose * conductance * Organ.DELTA_T

    def calculate_s_fructan(self, sucrose, regul_s_fructan):
        """Rate of fructan synthesis (µmol C fructan s-1 g-1 MS * DELTA_T)
        """
        return (((max(0, sucrose)/(self.mstruct*self.ALPHA))**(PhotosyntheticOrgan.N_SFRUCTAN) * PhotosyntheticOrgan.VMAX_SFRUCTAN) / ((max(0, sucrose)/(self.mstruct*self.ALPHA))**(PhotosyntheticOrgan.N_SFRUCTAN) + PhotosyntheticOrgan.K_SFRUCTAN**(PhotosyntheticOrgan.N_SFRUCTAN))) * regul_s_fructan * Organ.DELTA_T

    def calculate_d_fructan(self, sucrose, fructan):
        """Rate of fructan degradation (µmol C fructan s-1 g-1 MS)
        """
        return min((PhotosyntheticOrgan.K_DFRUCTAN * PhotosyntheticOrgan.VMAX_DFRUCTAN) / ((max(0, sucrose)/(Organ.MSTRUCT_AXIS*self.ALPHA)) + PhotosyntheticOrgan.K_DFRUCTAN) , max(0, fructan)) * Organ.DELTA_T

    # COMPARTMENTS

    def calculate_storage_derivative(self, s_storage, d_storage):
        """delta storage of organ integrated over delta-1 (µmol C storage).
        This is a differential equation of compartment expressed as a variation of the total amount of C substance in an organ per DELTA_T.
        """
        return (s_storage - d_storage) * (self.mstruct*self.ALPHA)

    def calculate_sucrose_derivative(self, s_sucrose, d_storage, loading_sucrose, s_fructan, d_fructan):
        """delta sucrose of organ integrated over delta-1 (µmol C sucrose)
        """
        return (s_sucrose + d_storage + d_fructan - s_fructan - loading_sucrose) * (self.mstruct*self.ALPHA)

    def calculate_triosesP_derivative(self, photosynthesis, s_sucrose, s_storage):
        """ delta triosesP of organ integrated over delta-1 (µmol C triosesP).
        This is a differential equation of compartment expressed as a variation of the total amount of C substance in an organ per DELTA_T.
        """
        return max(0, photosynthesis) - (s_sucrose + s_storage) * (self.mstruct*self.ALPHA)

    def calculate_fructan_derivative(self, s_fructan, d_fructan):
        """delta fructan integrated over delta-1 (µmol C fructan)
        """
        return (s_fructan - d_fructan)* (self.mstruct*self.ALPHA)


class Chaff(PhotosyntheticOrgan): 

    ALPHA = 1 #: Proportion of leaf structural mass containing substrate


class Lamina(PhotosyntheticOrgan):
    
    ALPHA = 1 #: Proportion of leaf structural mass containing substrate
    
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

    ALPHA = 1 #: Proportion of leaf structural mass containing substrate


class Peduncle(PhotosyntheticOrgan): 

    ALPHA = 1 #: Proportion of leaf structural mass containing substrate


class Sheath(PhotosyntheticOrgan): 

    ALPHA = 1 #: Proportion of leaf structural mass containing substrate


class Phloem(Organ):
    
    ALPHA = 1 #: Proportion of leaf structural mass containing substrate

    def __init__(self, sucrose_0, name=None):
        super(Phloem, self).__init__(name)

        # intialize the compartment
        self.sucrose_0 = sucrose_0  #: initial value of compartment sucrose

    def get_initial_conditions(self):
        return [self.sucrose_0]

    # VARIABLES

    def calculate_conc_sucrose(self, sucrose):
        """sucrose concentration (µmol sucrose g-1 MS)
        """
        return (sucrose/Organ.MSTRUCT_AXIS)/12

    def calculate_conc_c_sucrose(self, sucrose):
        """sucrose concentration expressed in C (µmol C sucrose g-1 MS)
        """
        return sucrose/(Organ.MSTRUCT_AXIS*Organ.ALPHA_AXIS)

    # COMPARTMENTS

    def calculate_sucrose_derivative(self, organs):
        """delta sucrose of phloem integrated over delta-1 (µmol C sucrose)
        """
        sucrose_derivative = 0
        for organ_ in organs:
            if isinstance(organ_, PhotosyntheticOrgan):
                sucrose_derivative += organ_.loading_sucrose*organ_.mstruct*self.ALPHA
            elif isinstance(organ_, Grains):
                sucrose_derivative -= (organ_.unloading_sucrose_structure + (organ_.unloading_sucrose_storage * ((organ_.structure/1E6)*12)))
            elif isinstance(organ_, Roots):
                sucrose_derivative -= (organ_.unloading_sucrose * organ_.mstruct)
        return sucrose_derivative


class Grains(Organ):
    
    ALPHA = 1 #: Proportion of leaf structural mass containing substrate

    # Structure
    GRAIN_STRUCTURE = 0     #: Initial value of structural mass of grains (µmol of C sucrose) TODO: initial value given in ModelMaker?
    VMAX_RGR = 1.9e-06      #: Maximal value of the Relative Growth Rate of grain structure (dimensionless)
    K_RGR = 300             #: Affinity coefficient of the Relative Growth Rate of grain structure (dimensionless)

    # storage
    GRAIN_STORAGE = 0       #: Initial value of grain storage (µmol of C) TODO: needed?
    VMAX_STORAGE = 0.5      #: Maximal rate of grain filling (µmol C s-1 g-1 MS)
    K_STORAGE = 100         #: Affinity coefficient of grain filling (µmol C g-1 MS)

    Y_GRAINS = 0.75         #: Proportion of C loaded from phloem actually used for grain structure and storage (1 - Y_GRAINS is a kind of growth respiration)
    FILLING_INIT = 360      #: Time (h) at which phloem loading switch from grain structure to grain storage

    def __init__(self, storage_0, structure_0, name=None):
        super(Grains, self).__init__(name)

        # flow to phloem
        self.unloading_sucrose_structure = 0 #: current unloading of sucrose from phloem to grain structure
        self.unloading_sucrose_storage = 0 #: current unloading of sucrose from phloem to grain storage
        self.structure = 0# TODO: utlisation valeur init??

        # initialize the compartments
        self.storage_0 = storage_0 #: initial value of compartment storage
        self.structure_0 = structure_0 #: initial value of compartment structure

    def get_initial_conditions(self):
        return [self.storage_0, self.structure_0]

    # VARIABLES

    def calculate_dry_mass(self, structure, storage):
        """Grain total dry mass (g)
        """
        return ((structure + storage)/1000000)*12

    def calculate_RGR_structure(self, sucrose_phloem):
        """Relative Growth Rate of grain structure
        """
        return ((max(0, sucrose_phloem)/(Organ.MSTRUCT_AXIS*Organ.ALPHA_AXIS)) * Grains.VMAX_RGR) / ((max(0, sucrose_phloem)/(Organ.MSTRUCT_AXIS*Organ.ALPHA_AXIS)) + Grains.K_RGR)

    # FLOWS

    def calculate_unloading_sucrose_structure(self, t, structure, RGR_structure):
        """Unloading of sucrose from phloem to grain structure (µmol C unloaded sucrose s-1 g-1 MS * DELTA_T)
        """
        if t<=Grains.FILLING_INIT:
            unloading_sucrose_structure = structure * RGR_structure * Organ.DELTA_T
        else:
            unloading_sucrose_structure = 0
        return unloading_sucrose_structure

    def calculate_unloading_sucrose_storage(self, t, sucrose_phloem):
        """Unloading of sucrose from phloem to grain storage (µmol C unloaded sucrose s-1 g-1 MS * DELTA_T)
        """
        if t<=Grains.FILLING_INIT:
            unloading_sucrose_storage = 0
        else:
            unloading_sucrose_storage = (((max(0, sucrose_phloem)/(Organ.MSTRUCT_AXIS*Organ.ALPHA_AXIS)) * Grains.VMAX_STORAGE) / ((max(0, sucrose_phloem)/(Organ.MSTRUCT_AXIS*Organ.ALPHA_AXIS)) + Grains.K_STORAGE)) * Organ.DELTA_T
        return unloading_sucrose_storage

    # COMPARTMENTS

    def calculate_structure_derivative(self, unloading_sucrose_structure):
        """delta grain structure integrated over delta-1 (µmol C structure)
        """
        return unloading_sucrose_structure*Grains.Y_GRAINS

    def calculate_storage_derivative(self, unloading_sucrose_storage, structure):
        """delta grain storage integrated over delta-1 (µmol C storage)
        """
        return unloading_sucrose_storage* Grains.Y_GRAINS * ((structure/1E6)*12) # Conversion of grain structure from µmol of C to g of C


class Roots(Organ):
    
    ALPHA = 1 #: Proportion of leaf structural mass containing substrate
    
    VMAX_ROOTS = 0.015  #: Maximal rate of sucrose unloading from phloem to roots (µmol C unloaded sucrose s-1 g-1 MS)
    K_ROOTS = 100       #: Affinity coefficient of sucrose unloading from phloem to roots (µmol C unloaded sucrose g-1 MS)

    def __init__(self, mstruct, sucrose_0, name=None):
        super(Roots, self).__init__(name)

        # parameters
        self.mstruct = mstruct #: Structural mass (g)

        self.unloading_sucrose = 0 #: current unloading of sucrose from phloem to roots

        # initialize the compartment
        self.sucrose_0 = sucrose_0 #: initial value of compartment sucrose

    def get_initial_conditions(self):
        return [self.sucrose_0] # keep this order !

    # VARIABLES

    def calculate_dry_mass(self, sucrose):
        """Dry mass of roots (g)
        """
        return (sucrose*12)/1000000

    # FLOWS

    def calculate_unloading_sucrose(self, sucrose_phloem):
        """Unloading of sucrose from phloem to roots (µmol C unloaded sucrose s-1 g-1 MS * DELTA_T)
        """
        return (((max(0, sucrose_phloem)/(Organ.MSTRUCT_AXIS*Organ.ALPHA_AXIS)) * Roots.VMAX_ROOTS) / ((max(0, sucrose_phloem)/(Organ.MSTRUCT_AXIS*Organ.ALPHA_AXIS)) + Roots.K_ROOTS))*Organ.DELTA_T

    # COMPARTMENTS

    def calculate_sucrose_derivative(self, unloading_sucrose):
        """delta root sucrose integrated over delta-1 (µmol C sucrose)
        """
        return unloading_sucrose*self.mstruct
