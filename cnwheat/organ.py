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

    Mstruct_axis = 2.08             #: Structural mass  of a plant (g) (Bertheloot, 2011)
    alpha_axis = 1                  #: Proportion of the structural mass containing the substrates
    delta_t = 3600                  #: Timestep of the model (s)

    # Sucrose
    Vmax_sucrose = 1                #: Maximal rate of sucrose synthesis (µmol C s-1 g-1 MS)
    K_sucrose = 0.66                #: Affinity coefficient of sucrose synthesis (µmol C g-1 MS)

    # Storage
    Vmax_storage = 2                #: Maximal rate of storage synthesis (µmol C s-1 g-1 MS)
    K_storage = 20                  #: Affinity coefficient of storage synthesis (µmol C g-1 MS)
    delta_Dstorage = 0.0001         #: Rate of storage degradation (µmol C storage s-1 g-1 MS)

    # Fructans
    Vmax_Sfructan = 0.2             #: Maximal rate of fructan synthesis (µmol C s-1 g-1 MS)
    K_Sfructan = 20000              #: Affinity coefficient of fructan synthesis (µmol C g-1 MS)
    n_Sfructan = 3                  #: Number of "substrates" for fructan synthesis (dimensionless)
    Vmax_regul_Sfructan = 1         #: Maximal value of the regulation function of fructan synthesis (dimensionless)
    K_regul_Sfructan = 8            #: Affinity coefficient of the regulation function of fructan synthesis (dimensionless)
    n_regul_Sfructan = 15           #: Parameter of the regulation function of fructan synthesis (dimensionless)
    Vmax_Dfructan = 0.035           #: Maximal rate of fructan degradation (µmol C s-1 g-1 MS)
    K_Dfructan = 100                #: Affinity coefficient of fructan degradation (µmol C g-1 MS)

    sigma = 1.85e-07                #: Conductivity of an organ-phloem pathway (g mol-1 m-2 s-1) ; used to compute the sucrose loaded to the phloem


    alpha = 1 #: Proportion of leaf structural mass containing substrate # TODO: is it really common to all organs? NOTE: I think this parameter is unneeded

    def __init__(self, name):
        self.name = name

    # VARIABLES

    def calculate_Photosynthesis(self, t, An):
        """Total photosynthesis of an organ integrated over delta_t (µmol CO2 on organ area integrated over delat_t)
        """
        return An * self._calculate_green_area(t) * Organ.delta_t
    
    def calculate_Conc_TriosesP(self, TRIOSESP):
        """Trioses Phosphate concentration (µmol TriosesP g-1 MS).
        This is a concentration output (expressed in amount of substance g-1 MS).
        """
        return (TRIOSESP/self.Mstruct)/3

    def calculate_Conc_Sucrose(self, SUCROSE):
        """Sucrose concentration (µmol sucrose g-1 MS).
        This is a concentration output (expressed in amount of substance g-1 MS).
        """
        return (SUCROSE/self.Mstruct)/12

    def calculate_Conc_Storage(self, STORAGE):
        """Storage concentration (µmol storage g-1 MS (eq glucose)).
        This is a concentration output (expressed in amount of substance g-1 MS).
        """
        return (STORAGE/self.Mstruct)/6
    
    def _calculate_green_area(self, t):
        """Compute green area of the organ.
        """
        return self.Area
    
    def calculate_Conc_Fructan(self, FRUCTAN):
        """Fructan concentration (µmol fructan g-1 MS (eq glucose))
        """
        return (FRUCTAN/self.Mstruct)/6

    def calculate_Regul_Sfructan(self, Loading_sucrose):
        """Inhibition of fructan synthesis by the loading of sucrose to phloem
        """
        return ((Organ.Vmax_regul_Sfructan * Organ.K_regul_Sfructan**(Organ.n_regul_Sfructan)) / ((max(0, Loading_sucrose)/(self.Mstruct*Organ.alpha))**(Organ.n_regul_Sfructan) + Organ.K_regul_Sfructan**(Organ.n_regul_Sfructan)))

    # FLOWS

    def calculate_S_sucrose(self, TRIOSESP):
        """Rate of SUCROSE synthesis from TRIOSESP (µmol C sucrose s-1 g-1 MS * delta_t).
        This is a flow (expressed in amount of C substance g-1 MS integrated over delta_t).
        """
        return (((max(0,TRIOSESP)/(self.Mstruct*Organ.alpha)) * Organ.Vmax_sucrose) / ((max(0, TRIOSESP)/(self.Mstruct*Organ.alpha)) + Organ.K_sucrose)) * Organ.delta_t

    def calculate_S_storage(self, TRIOSESP):
        """Rate of STORAGE synthesis from TRIOSESP (µmol C storage s-1 g-1 MS * delta_t).
        This is a flow (expressed in amount of C substance g-1 MS integrated over delta_t).
        """
        return (((max(0, TRIOSESP)/(self.Mstruct*Organ.alpha)) * Organ.Vmax_storage) / ((max(0, TRIOSESP)/(self.Mstruct*Organ.alpha)) + Organ.K_storage)) * Organ.delta_t

    def calculate_D_storage(self, STORAGE):
        """Rate of STORAGE degradation from TRIOSESP (µmol C storage s-1 g-1 MS * delta_t).
        This is a flow (expressed in amount of C substance g-1 MS integrated over delta_t).
        """
        return max(0, Organ.delta_Dstorage * (STORAGE/(self.Mstruct*Organ.alpha))) * Organ.delta_t

    def calculate_Loading_sucrose(self, SUCROSE, SUCROSE_phloem):
        """Rate of SUCROSE loading to phloem (µmol C sucrose s-1 g-1 MS * delta_t).
        This is a flow (expressed in amount of C substance g-1 MS integrated over delta_t).
        """
        return ((max(0.1, SUCROSE)/(self.Mstruct*Organ.alpha)) * ((max(0.1, SUCROSE)/(self.Mstruct*Organ.alpha)) - (max(0, SUCROSE_phloem)/(Organ.Mstruct_axis*Organ.alpha_axis))) * (Organ.sigma * self.Mstruct**(2/3))) * Organ.delta_t
    
    def calculate_S_fructan(self, SUCROSE, Regul_Sfructan):
        """Rate of fructan synthesis (µmol C fructan s-1 g-1 MS * delta_t)
        """
        return (((max(0, SUCROSE)/(self.Mstruct*Organ.alpha))**(Organ.n_Sfructan) * Organ.Vmax_Sfructan) / ((max(0, SUCROSE)/(self.Mstruct*Organ.alpha))**(Organ.n_Sfructan) + Organ.K_Sfructan**(Organ.n_Sfructan))) * Regul_Sfructan * Organ.delta_t

    def calculate_D_fructan(self, SUCROSE, FRUCTAN):
        """Rate of fructan degradation (µmol C fructan s-1 g-1 MS)
        """
        return min((Organ.K_Dfructan * Organ.Vmax_Dfructan) / ((max(0, SUCROSE)/(self.Mstruct_axis*Organ.alpha)) + Organ.K_Dfructan) , max(0, FRUCTAN)) * Organ.delta_t

    # COMPARTMENTS

    def calculate_TRIOSESP_derivative(self, Photosynthesis, S_sucrose, S_storage):
        """ delta TRIOSESP of organ integrated over delta-1 (µmol C TRIOSESP).
        This is a differential equation of compartment expressed as a variation of the total amount of C substance in an organ per delta_t.
        """
        return Photosynthesis - (S_sucrose + S_storage) * (self.Mstruct*Organ.alpha)

    def calculate_STORAGE_derivative(self, S_storage, D_storage):
        """delta STORAGE of organ integrated over delta-1 (µmol C STORAGE).
        This is a differential equation of compartment expressed as a variation of the total amount of C substance in an organ per delta_t.
        """
        return (S_storage - D_storage) * (self.Mstruct*Organ.alpha)
    
    def calculate_FRUCTAN_derivative(self, S_fructan, D_fructan):
        """delta FRUCTAN of internode integrated over delta-1 (µmol C FRUCTAN)
        """
        return (S_fructan - D_fructan)* (self.Mstruct*Organ.alpha)

    def calculate_SUCROSE_derivative(self, S_sucrose, D_storage, Loading_sucrose, S_fructan=0, D_fructan=0):
        """delta SUCROSE of organ integrated over delta-1 (µmol C SUCROSE)
        """
        return (S_sucrose + D_storage + D_fructan - S_fructan - Loading_sucrose) * (self.Mstruct*Organ.alpha)
    

class Lamina(Organ):
    #: Temporary estimation of lamina senescence ({'lamina_order': (time of senescence beginning (h), offset of the linear regression)})
    laminae_inflexion_points = {'lamina1': (600, 78.75),
                                'lamina2': (480, 68.61),
                                'lamina3': (360, 48.76)}

    def __init__(self, Area, Mstruct, PAR, STORAGE_0,
                 SUCROSE_0, TRIOSESP_0, name='lamina'):
        super(Lamina, self).__init__(name)
        # parameters
        self.Area = Area                    #: Area (m-2)
        self.Mstruct = Mstruct              #: Structural mass (g)
        self.PAR = PAR    #: PAR estimated from photosynthesis model, TODO: homogeneiser les termes pr la photosynthese (unite?)

        self.PAR_linear_interpolation = None #: linear interpolation of PAR
        
        self.photosynthesis_mapping = {} #: mapping to store the computed photosynthesis

        self.Loading_sucrose = 0 #: current rate of SUCROSE loading to phloem
        
        # initialize the compartments
        self.TRIOSESP_0 = TRIOSESP_0 #: initial value of compartment TRIOSESP
        self.STORAGE_0 = STORAGE_0 #: initial value of compartment STORAGE
        self.SUCROSE_0 = SUCROSE_0 #: initial value of compartment SUCROSE


    def get_initial_conditions(self):
        return [self.STORAGE_0, self.SUCROSE_0, self.TRIOSESP_0] # keep this order ! TODO: pas le plus "logique" mais ok si trop complique a modifier

    # VARIABLES

    def _calculate_green_area(self, t):
        """Compute green area of the organ.
        """
        t_inflexion, value_inflexion = Lamina.laminae_inflexion_points.get(self.name, (float("inf"), None))
        if t <= t_inflexion: # Non-senescent lamina
            green_area = self.Area
        else: # Senescent lamina
            green_area = ((-0.0721*t + value_inflexion)/10000)
        return green_area


class Phloem(Organ):

    def __init__(self, SUCROSE_0, name='phloem'):
        super(Phloem, self).__init__(name)
        
        # intialize the compartment
        self.SUCROSE_0 = SUCROSE_0  #: initial value of compartment SUCROSE

    def get_initial_conditions(self):
        return [self.SUCROSE_0]

    # VARIABLES

    def calculate_Conc_Sucrose(self, SUCROSE):
        """Sucrose concentration (µmol sucrose g-1 MS)
        """
        return (SUCROSE/Organ.Mstruct_axis)/12

    def calculate_Conc_C_Sucrose(self, SUCROSE):
        """Sucrose concentration expressed in C (µmol C sucrose g-1 MS)
        """
        return SUCROSE/(Organ.Mstruct_axis*Organ.alpha_axis)

    # COMPARTMENTS

    def calculate_SUCROSE_derivative(self, organs):
        """delta SUCROSE of phloem integrated over delta-1 (µmol C SUCROSE)
        """
        SUCROSE_derivative = 0
        for organ_ in organs:
            if isinstance(organ_, (Lamina, Sheath, Internode, Peduncle, Chaff)):
                SUCROSE_derivative += organ_.Loading_sucrose*organ_.Mstruct*Organ.alpha
            elif isinstance(organ_, Grains):
                SUCROSE_derivative -= (organ_.Unloading_sucrose_structure + (organ_.Unloading_sucrose_storage * ((organ_.STRUCTURE/1E6)*12)))
            elif isinstance(organ_, Roots):
                SUCROSE_derivative -= (organ_.Unloading_sucrose * organ_.Mstruct)
        return SUCROSE_derivative


class Chaff(Organ):

    def __init__(self, Area, Mstruct, PAR, STORAGE_0, SUCROSE_0,
                 TRIOSESP_0, name='chaff'):
        super(Chaff, self).__init__(name)
        # parameters
        self.Area = Area                    #: Area (m-2)
        self.Mstruct = Mstruct              #: Structural mass (g)
        self.PAR = PAR

        self.PAR_linear_interpolation = None #: linear interpolation of PAR
        
        self.photosynthesis_mapping = {} #: mapping to store the computed photosynthesis

        self.Loading_sucrose = 0 #: current rate of SUCROSE loading to phloem

        # initialize the compartments
        self.TRIOSESP_0 = TRIOSESP_0 #: initial value of compartment TRIOSESP
        self.STORAGE_0 = STORAGE_0 #: initial value of compartment STORAGE
        self.SUCROSE_0 = SUCROSE_0 #: initial value of compartment SUCROSE

    def get_initial_conditions(self):
        return [self.STORAGE_0, self.SUCROSE_0, self.TRIOSESP_0] # keep this order !


class Internode(Organ):

    def __init__(self, Area, Mstruct, PAR, FRUCTAN_0, STORAGE_0,
                 SUCROSE_0, TRIOSESP_0, name='internode'):
        super(Internode, self).__init__(name)
        # parameters
        self.Area = Area                    #: Area (m-2)
        self.Mstruct = Mstruct              #: Structural mass (g)
        self.PAR = PAR
        
        self.PAR_linear_interpolation = None #: linear interpolation of PAR
        
        self.photosynthesis_mapping = {} #: mapping to store the computed photosynthesis

        self.Loading_sucrose = 0 #: current rate of SUCROSE loading to phloem

        # initialize the compartments
        self.TRIOSESP_0 = TRIOSESP_0 #: initial value of compartment TRIOSESP
        self.STORAGE_0 = STORAGE_0 #: initial value of compartment STORAGE
        self.SUCROSE_0 = SUCROSE_0 #: initial value of compartment SUCROSE
        self.FRUCTAN_0 = FRUCTAN_0 #: initial value of compartment FRUCTAN

    def get_initial_conditions(self):
        return [self.FRUCTAN_0, self.STORAGE_0, self.SUCROSE_0, self.TRIOSESP_0]


class Peduncle(Organ):

    def __init__(self, Area, Mstruct, PAR, FRUCTAN_0, STORAGE_0,
                 SUCROSE_0, TRIOSESP_0, name='peduncle'):
        super(Peduncle, self).__init__(name)
        # parameters
        self.Area = Area                    #: Area (m-2)
        self.Mstruct = Mstruct              #: Structural mass (g)
        self.PAR = PAR

        self.PAR_linear_interpolation = None #: linear interpolation of PAR
        
        self.photosynthesis_mapping = {} #: mapping to store the computed photosynthesis

        self.Loading_sucrose = 0 #: current rate of SUCROSE loading to phloem

        # initialize the compartments
        self.TRIOSESP_0 = TRIOSESP_0 #: initial value of compartment TRIOSESP
        self.STORAGE_0 = STORAGE_0 #: initial value of compartment STORAGE
        self.SUCROSE_0 = SUCROSE_0 #: initial value of compartment SUCROSE
        self.FRUCTAN_0 = FRUCTAN_0 #: initial value of compartment FRUCTAN

    def get_initial_conditions(self):
        return [self.FRUCTAN_0, self.STORAGE_0, self.SUCROSE_0, self.TRIOSESP_0]


class Sheath(Organ):

    def __init__(self, Area, Mstruct, PAR, FRUCTAN_0, STORAGE_0,
                 SUCROSE_0, TRIOSESP_0, name='sheath'):
        super(Sheath, self).__init__(name)

        # parameters
        self.Area = Area                    #: Area (m-2)
        self.Mstruct = Mstruct              #: Structural mass (g)
        self.PAR = PAR

        self.PAR_linear_interpolation = None #: linear interpolation of PAR
        
        self.photosynthesis_mapping = {} #: mapping to store the computed photosynthesis

        self.Loading_sucrose = 0 #: current rate of SUCROSE loading to phloem

        # initialize the compartments
        self.TRIOSESP_0 = TRIOSESP_0 #: initial value of compartment TRIOSESP
        self.STORAGE_0 = STORAGE_0 #: initial value of compartment STORAGE
        self.SUCROSE_0 = SUCROSE_0 #: initial value of compartment SUCROSE
        self.FRUCTAN_0 = FRUCTAN_0 #: initial value of compartment FRUCTAN

    def get_initial_conditions(self):
        return [self.FRUCTAN_0, self.STORAGE_0, self.SUCROSE_0, self.TRIOSESP_0]


class Grains(Organ):

    # Structure
    Grain_structure = 0     #: Initial value of structural mass of grains (µmol of C sucrose) TODO: initial value given in ModelMaker?
    Vmax_RGR = 1.9e-06      #: Maximal value of the Relative Growth Rate of grain structure (dimensionless)
    K_RGR = 300             #: Affinity coefficient of the Relative Growth Rate of grain structure (dimensionless)

    # Storage
    Grain_storage = 0       #: Initial value of grain storage (µmol of C) TODO: needed?
    Vmax_storage = 0.5      #: Maximal rate of grain filling (µmol C s-1 g-1 MS)
    K_storage = 100         #: Affinity coefficient of grain filling (µmol C g-1 MS)

    Y_grains = 0.75         #: Proportion of C loaded from phloem actually used for grain structure and storage (1 - Y_grains is a kind of growth respiration)
    filling_init = 360      #: Time (h) at which phloem loading switch from grain structure to grain storage

    def __init__(self, STORAGE_0, STRUCTURE_0, name='grains'):
        super(Grains, self).__init__(name)

        # flow to phloem
        self.Unloading_sucrose_structure = 0 #: current unloading of sucrose from phloem to grain structure
        self.Unloading_sucrose_storage = 0 #: current unloading of sucrose from phloem to grain storage
        self.STRUCTURE = 0# TODO: utlisation valeur init??
        
        # initialize the compartments
        self.STORAGE_0 = STORAGE_0 #: initial value of compartment STORAGE
        self.STRUCTURE_0 = STRUCTURE_0 #: initial value of compartment STRUCTURE

    def get_initial_conditions(self):
        return [self.STORAGE_0, self.STRUCTURE_0]

    # VARIABLES

    def calculate_Dry_mass(self, STRUCTURE, STORAGE):
        """Grain total dry mass (g)
        """
        return ((STRUCTURE + STORAGE)/1000000)*12

    def calculate_RGR_structure(self, SUCROSE_phloem):
        """Relative Growth Rate of grain structure
        """
        return ((max(0, SUCROSE_phloem)/(Organ.Mstruct_axis*Organ.alpha_axis)) * Grains.Vmax_RGR) / ((max(0, SUCROSE_phloem)/(Organ.Mstruct_axis*Organ.alpha_axis)) + Grains.K_RGR)

    # FLOWS

    def calculate_Unloading_sucrose_structure(self, t, STRUCTURE, RGR_structure):
        """Unloading of sucrose from phloem to grain structure (µmol C unloaded sucrose s-1 g-1 MS * delta_t)
        """
        if t<=Grains.filling_init:
            Unloading_sucrose_structure = STRUCTURE * RGR_structure * Organ.delta_t
        else:
            Unloading_sucrose_structure = 0
        return Unloading_sucrose_structure

    def calculate_Unloading_sucrose_storage(self, t, SUCROSE_phloem):
        """Unloading of sucrose from phloem to grain storage (µmol C unloaded sucrose s-1 g-1 MS * delta_t)
        """
        if t<=Grains.filling_init:
            Unloading_sucrose_storage = 0
        else:
            Unloading_sucrose_storage = (((max(0, SUCROSE_phloem)/(Organ.Mstruct_axis*Organ.alpha_axis)) * Grains.Vmax_storage) / ((max(0, SUCROSE_phloem)/(Organ.Mstruct_axis*Organ.alpha_axis)) + Grains.K_storage)) * Organ.delta_t
        return Unloading_sucrose_storage

    # COMPARTMENTS

    def calculate_STRUCTURE_derivative(self, Unloading_sucrose_structure):
        """delta grain STRUCTURE integrated over delta-1 (µmol C STRUCTURE)
        """
        return Unloading_sucrose_structure*Grains.Y_grains

    def calculate_STORAGE_derivative(self, Unloading_sucrose_storage, STRUCTURE):
        """delta grain STORAGE integrated over delta-1 (µmol C STORAGE)
        """
        return Unloading_sucrose_storage* Grains.Y_grains * ((STRUCTURE/1E6)*12) # Conversion of grain structure from µmol of C to g of C


class Roots(Organ):
    Vmax_roots = 0.015  #: Maximal rate of sucrose unloading from phloem to roots (µmol C unloaded sucrose s-1 g-1 MS)
    K_roots = 100       #: Affinity coefficient of sucrose unloading from phloem to roots (µmol C unloaded sucrose g-1 MS)

    def __init__(self, Mstruct, Sucrose_0, name='roots'):
        super(Roots, self).__init__(name)

        # parameters
        self.Mstruct = Mstruct #: Structural mass (g)

        self.Unloading_sucrose = 0 #: current unloading of sucrose from phloem to roots

        # initialize the compartment
        self.Sucrose_0 = Sucrose_0 #: initial value of compartment Sucrose

    def get_initial_conditions(self):
        return [self.Sucrose_0] # keep this order !

    # VARIABLES

    def calculate_Dry_mass(self, Sucrose):
        """Dry mass of roots (g)
        """
        return (Sucrose*12)/1000000

    # FLOWS

    def calculate_Unloading_sucrose(self, SUCROSE_phloem):
        """Unloading of sucrose from phloem to roots (µmol C unloaded sucrose s-1 g-1 MS * delta_t)
        """
        return (((max(0, SUCROSE_phloem)/(Organ.Mstruct_axis*Organ.alpha_axis)) * Roots.Vmax_roots) / ((max(0, SUCROSE_phloem)/(Organ.Mstruct_axis*Organ.alpha_axis)) + Roots.K_roots))*Organ.delta_t

    # COMPARTMENTS

    def calculate_Sucrose_derivative(self, Unloading_sucrose):
        """delta root Sucrose integrated over delta-1 (µmol C Sucrose)
        """
        return Unloading_sucrose*self.Mstruct