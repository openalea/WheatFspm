#!/usr/bin/env python
# -*- coding: latin-1 -*-

from __future__ import division # use '//' to do integer division

"""
    farquharwheat.model
    ~~~~~~~~~~~~~~~~~~~

    Model of photosynthesis based on Farquhar's approach. The model includes the dependance of photosynthesis to organ temperature and nitrogen content. Internal CO2 and organ temperature are found numerically.

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

from math import sqrt, log,  exp

class PhotosynthesisModel(object):

    O = 21000      #: Photosynthetic parameter: Intercellular O2 concentration, umol mol(air)-1 or Pa, from Bernacchi et al. (2001)
    KC25 = 404     #: Photosynthetic parameter: Affinity constant of RuBisCO for C, umol mol-1 or Pa, from Bernacchi et al. (2001) (estimation in Braune et al. (2009) not enough accurate)
    KO25 = 278.4E3 #: Photosynthetic parameter: Affinity constant of RuBisCO for O, umol mol-1 or Pa, from Bernacchi et al. (2001) (estimation in Braune et al. (2009) not enough accurate)
    GAMMA25 = 39   #: Photosynthetic parameter: CO2 compensation point, umol(CO2) mol-1 (air), from Braune et al. (2009)
    THETA = 0.72    #: Photosynthetic parameter: curvature parameter of J, dimensionless

    #: Nitrogen dependance of photosynthetic parameters
    #: Derived from Braune et al. (2009):
    #:
    #:     * S_Na: slope of the relation between Na and the parameters (umol g-1 s-1)
    #:     * Na_min: minimum amount of leaf nitrogen below which photosynthesis rate is zero (g (N) m-2 leaf)
    #:     * Gamma_Na1 and Gamma_Na2: parameters of ALPHA dependance to Na (mol mol-1 and m2 g-1 respectively)
    #:     * delta1 and delta2: parameters of m (scaling factor of gs) dependance to Na (m2 g-1 and dimensionless respectively)
    PARAM_N = {'S_Na': {'Vc_max25': 63.2, 'Jmax25': 151, 'TPU25': 9.25, 'Rdark25': 0.493}, 'Na_min': {'Vc_max25': 0.198, 'Jmax25': 0.225, 'TPU25': 0.229, 'Rdark25': 0.118},
                'Gamma_Na1': 0.437, 'Gamma_Na2': 2.29, 'delta1': 14.7, 'delta2': -0.548}

    GSMIN = 0.05            #: Stomatal conductance parameter: Minimum gs, measured in the dark (mol m-2 s-1). Braune et al. (2009).
    GB = 3.5                #: Stomatal conductance parameter: Boundary layer conductance (mol m-2 s-1). Muller et al., (2005)

    SIGMA = 5.6704E-8       #: Physical parameter: Stefan-Bolzmann constant (W-2 K-4)
    I0 = 1370               #: Physical parameter: Extraterrestrial solar radiation (W m-2)
    LAMBDA = 2260E3         #: Physical parameter: Latent heat for vaporisation of water (J kg-1)
    RHOCP = 1256            #: Physical parameter: Volumetric heat capacity of air (J m-3 K-1)
    GAMMA = 66E-3           #: Physical parameter: Psychrometric constant (KPa K-1). Mean value
    A = 2.5                 #: Physical parameter: Attenuation coefficient of wind within a wheat canopy. From Campbell and Norman (1998), second edition. Can also be estimaed by: A = sqrt((0.2*LAI*h)/sqrt((4*width*h)/(pi*LAI))
    R = 8.3144              #: Physical parameter: Gas constant (J mol-1 K-1)
    PATM = 1.01325E5        #: Physical parameter: Atmospheric pressure (Pa)

    FR = 0.1                #: Physical parameter: Leaf radiation reflectance
    FT = 0.1                #: Physical parameter: Leaf radiation transmittance

    #: Temperature dependance of photosynthetic parameters
    #:
    #: Parameter values derived from Braune et al. (2009) except for Kc, Ko, and Rdark (Bernacchi et al., 2001)
    #:     * deltaHa, deltaHd: enthalpie of activation and deactivation respectively (kJ mol-1)
    #:     * deltaS: entropy term (kJ mol-1 K-1)
    #:     * Tref: reference temperature (K)
    #:     * R: universal gas constant (kJ mol-1 K-1)
    PARAM_TEMP = {'deltaHa': {'Vc_max': 89.7, 'Jmax': 48.9, 'TPU': 47., 'Kc': 79.43, 'Ko': 36.38, 'Gamma': 35., 'Rdark': 46.39},
                  'deltaHd': {'Vc_max': 149.3, 'Jmax': 152.3, 'TPU': 152.3},
                  'deltaS': {'Vc_max': 0.486, 'Jmax': 0.495, 'TPU': 0.495},
                  'Tref': 298.15, 'R': 8.3145E-03}

    @classmethod
    def _leaf_temperature(cls, leaf_width, z, H, wind0, PAR, gs, Ta, Tleaf, RH):
        """
        Energy balance for the estimation of leaf temperature
            - leaf_width (m)
            - z: height of leaf from soil (m)
            - H: canopy height (m)
            - Wind0: wind at the top of the canopy (m s-1)
            - PAR (umol m-2 s-1)
            - gs: stomatal conductance (mol m-2 s-1)
            - Ta: air temperature (degree C)
            - Tleaf: leaf temperature (degree C). Tleaf = Ta at the first iteration of the numeric resolution
            - RH: Relative humidity (decimal fraction)

        :class:`_leaf_temperature` should NOT be INSTANTIATE directly.
        """

        # Wind speed (m s-1)
        wind = wind0 * exp(cls.A*(z/H -1))                      # From Campbell and Norman (1998), second edition. z: organ height, H: height of the canopy

        # Boundary layer restistance to heat (s m-1)
        rbh = 100*sqrt(leaf_width/wind)

        # Turbulence resistance to heat (s m-1)
        rt = (0.74 * log(((2 - 0.7*H) / (0.1*H))**2)) / (0.16*wind)

        # Net absorbed radiation Rn (PAR and NIR, J m-2 s-1)
        Iabs = (PAR*(1-cls.FR-cls.FT))/(0.55*4.55)          # Global absorbed radiation by leaf (J m-2 s-1). TODO: relation a verifier
        es_Ta = 0.611 * exp((17.4*Ta)/(239+Ta))             # Saturated vapour pressure of the air (kPa), Ta in degree Celsius
        V = RH * es_Ta                                      # Vapour pressure of the air (kPa)
        fvap = 0.56 - 0.079*sqrt(10*V)                      # Fraction of vapour pressure

        tau = Iabs/cls.I0                                   # Atmospheric transmissivity (dimensionless)
        fclear = 0.1 + 0.9*max(0, min(1, (tau-0.2)/0.5))    # Fraction sky clearness

        Rn = Iabs - cls.SIGMA * (Tleaf+273)**4*fvap*fclear

        # Transpiration (mm s-1), Penman-Monteith
        if Tleaf == Ta:
            Ta_K = Ta + 273.15                              # Ta in kelvin
            s = ((17.4*239)/(Ta_K + 239)**2)*es_Ta          # Slope of the curve relating saturation vapour pressure to temperature (kPa K-1)
        else:
            es_Tl = 0.611 * exp((17.4*Tleaf)/(239+Tleaf))   # Saturated vapour pressure at leaf (kPa), Tleaf in degree Celsius
            Tleaf_K, Ta_K = Tleaf + 273.15, Ta + 273.15     # Temperatures in kelvin
            s = (es_Tl - es_Ta)/(Tleaf_K - Ta_K)            # Slope of the curve relating saturation vapour pressure to temperature (kPa K-1)

        VPDa = es_Ta - V
        rbw = 0.96 * rbh                                    # Boundary layer resistance for water (s m-1)
        gsw = (1.6*gs * cls.R * (Tleaf+273.15)) / cls.PATM  # Stomatal conductance to water (m s-1). 1.6 convert gs_CO2 in gs_water. Relation given by A. Tuzet (2003)
        rswp = 1/gsw                                        # Stomatal resistance for water (s m-1)

        Tr = max(0, (s * Rn + (cls.RHOCP * VPDa)/(rbh + rt)) / (cls.LAMBDA * (s + cls.GAMMA*((rbw + rt + rswp)/(rbh + rt)))))

        # Leaf temperature
        Tleaf = Ta + ((rbh + rt) * (Rn - cls.LAMBDA*Tr)) / cls.RHOCP
        return Tleaf, Tr

    @classmethod
    def _stomatal_conductance(cls, Ag, An, Na, ambient_CO2, RH):
        """
        BWB model of stomatal conductance
            - Ag: global assimilation (umol m-2 s-1)
            - An: net assimilation (umol m-2 s-1)
            - Na: nitrogen content of leaf (g m-2)
            - ambient_CO2: Air CO2 (umol mol-1)
            - RH: Relative humidity (decimal fraction)
        :class:`_stomatal_conductance` should NOT be INSTANTIATE directly.
        """

        Cs = ambient_CO2 - An *(1.37/(cls.GB))                   # CO2 concentration at leaf surface (umol mol-1 or Pa). From Prieto et al. (2012). GB in mol m-2 s-1
        m = cls.PARAM_N['delta1'] * Na**cls.PARAM_N['delta2']    # Scaling factor dependance to Na (dimensionless). This focntion is maintained although I'm not sure that it should be taken into account
        gs = (cls.GSMIN + m*((Ag*RH)/(Cs)))                      # Stomatal conductance (mol m-2 s-1), from Braune et al. (2009), Muller et al. (2005): using Ag rather than An. Would be better with a function of VPD and with (Ci-GAMMA) instead of Cs.
        return gs

    @classmethod
    def _f_temperature(cls, pname, p25, T):
        """
        Photosynthetic parameters relation to temperature
            - pname: name of parameter
            - p25: parameter value at 25 degree C
            - T: leaf temperature (degree C)
        :class:`_f_temperature` should NOT be INSTANTIATE directly.
        """
        Tk = T + 273.15
        deltaHa = cls.PARAM_TEMP['deltaHa'][pname]
        Tref = cls.PARAM_TEMP['Tref']
        R = cls.PARAM_TEMP['R']

        f_activation = exp((deltaHa * (Tk - Tref))/(R * Tref * Tk))

        if pname in ('Vc_max', 'Jmax', 'TPU'):
            deltaS = cls.PARAM_TEMP['deltaS'][pname]
            deltaHd = cls.PARAM_TEMP['deltaHd'][pname]
            f_deactivation = (1 + exp((Tref*deltaS - deltaHd) / (Tref*R))) / (1 + exp((Tk*deltaS - deltaHd) / (Tk*R)))
        else:
            f_deactivation = 1

        p = p25 * f_activation * f_deactivation

        return p

    @classmethod
    def _photosynthesis(cls, PAR, Na, Tleaf, Ci):
        """
        Estimation of C02 assimilation. In this version, most of parameters are derived from Braune et al. (2009) on barley
            - PAR: PAR intercepted by leaf (umol m-2 s-1)
            - Na: nitrogen content of leaf (g m-2)
            - Tleaf: leaf temperature (degree C)
            - Ci: internal CO2 (umol mol-1), by default = 0.7*CO2air
        :class:`_photosynthesis` should NOT be INSTANTIATE directly.
        """

        ### # RuBisCO-limited carboxylation rate ###
        # RuBisCO parameters dependance to temperature
        Kc = cls._f_temperature('Kc', cls.KC25, Tleaf)
        Ko = cls._f_temperature('Ko', cls.KO25, Tleaf)
        Gamma = cls._f_temperature('Gamma', cls.GAMMA25, Tleaf)

        # Vcmax
        Vc_max25 = 84.965 * (Na - 0.17)                                                     # Relation between Vc_max25 and Na
        Vc_max = cls._f_temperature ('Vc_max', Vc_max25, Tleaf)                             # Relation between Vc_max and temperature
        Ac = (Vc_max * (Ci-Gamma)) / (Ci + Kc * (1 + cls.O/Ko))                             # Rate of assimilation under Vc_max limitation
        ### RuBP regeneration-limited carboxylation rate via electron transport ###
        ALPHA = 0.0413 * Na + 0.2101                                                        # Relation between ALPHA and Na
        Jmax25 = 117.6 * (Na - 0.17)                                                        # Relation between Jmax25 and Na
        Jmax = cls._f_temperature('Jmax', Jmax25, Tleaf)                                    # Relation between Jmax and temperature

        # Electron transport rate
        J = ((Jmax+ALPHA*PAR) - sqrt((Jmax+ALPHA*PAR)**2 - 4*cls.THETA*ALPHA*PAR*Jmax))/(2*cls.THETA) # Muller et al. (2005), Evers et al. (2010)
        Aj = (J * (Ci-Gamma)) / (4*Ci + 8*Gamma)                                            # Rate of assimilation under RuBP regeneration limitation
        ### Triose phosphate utilisation-limited carboxylation rate ###
        TPU25 = cls.PARAM_N['S_Na']['TPU25'] * (Na - cls.PARAM_N['Na_min']['TPU25'])        # Relation between TPU25 and Na
        TPU = cls._f_temperature('TPU', TPU25, Tleaf)                                       # Relation between TPU and temperature
        Vomax = (Vc_max*Ko*Gamma)/(0.5*Kc*cls.O)                                            # Maximum rate of Vo (umol m-2 s-1)
        Vo = (Vomax * cls.O) / (cls.O + Ko*(1+Ci/Kc))                                       # Rate of oxygenation of RuBP (umol m-2 s-1)
        Ap = (1-Gamma/Ci)*(3*TPU) + Vo                                                      # Rate of assimilation under TPU limitation

        # Gross assimilation rate
        Ag = min(Ac, Aj, Ap)

        # Mitochondrial respiration rate of leaf in light Rd (processes other than photorespiration)
        Rdark25 = cls.PARAM_N['S_Na']['Rdark25'] * (Na - cls.PARAM_N['Na_min']['Rdark25'])  # Relation between Rdark25 (respiration in obscurity at 25 degree C) and Na
        Rdark = cls._f_temperature('Rdark', Rdark25, Tleaf)                                 # Relation between Rdark and temperature
        Rd = Rdark * (0.33 + (1-0.33)*(0.5)**(PAR/15))                                      # Found in Muller et al. (2005), eq. 19
        # Net C assimilation
        An = Ag - Rd
        return An, Ag, Rd


    @classmethod
    def calculate_An(cls, organ_width, organ_height, PAR, Ta, ambient_CO2, RH, wind0):
        """
        Compute An (µmol m-2 s-1) and Tr (mm s-1).

        :Parameters:

            - `organ_width` (:class:`float`) - width of the organ (or diameter for stem organ) (m)

            - `organ_height` (:class:`float`) - height of the organ from soil (m)

            - `PAR` (:class:`float`) - PAR (µmol m-2 s-1)

            - `Ta` (:class:`float`) - air temperature (degree Celsius)

            - `ambient_CO2` (:class:`float`) - air CO2 (umol mol-1)

            - `RH` (:class:`float`) - relative humidity (decimal fraction)

            - `wind0` (:class:`float`) - wind speed at the top of the canopy (m s-1)

        :Returns:
            An (µmol m-2 s-1) and Tr (mm s-1).

        :Returns Type:
            :class:`float`

        """
        #: Organ parameters
        LEAF_WIDTH = organ_width                    #: m
        H_CANOPY = 0.8                              #: m
        H_ORGAN = organ_height                      #: m
        NA_INIT = 2.5                               #: g m-2

        ### Iterations to find leaf temperature and Ci ###
        Ci, Tleaf = 0.7*ambient_CO2, Ta # Initial values
        prec_Ci, prec_Tleaf = 0.1, 0.1
        count = 0
        while abs((Ci - prec_Ci)/prec_Ci) >= 0.01 or abs((Tleaf - prec_Tleaf)/prec_Tleaf) >= 0.01:
            if count >=30: # TODO: test a faire? Semble prendre du tps de calcul
                if abs((Ci - prec_Ci)/prec_Ci) >= 0.01:
                    print 'Ci cannot converge, prec_Ci= {}, Ci= {}'.format(prec_Ci, Ci)
                else:
                    print 'Tleaf cannot converge, prec_Tleaf= {}, Tleaf= {}'.format(prec_Tleaf, Tleaf)
                break
            else:
                prec_Ci, prec_Tleaf = Ci, Tleaf
                An, Ag, Rd = cls._photosynthesis(PAR, NA_INIT, Tleaf, Ci)
                # Stomatal conductance
                gs = cls._stomatal_conductance(Ag, An, NA_INIT, ambient_CO2, RH)
                # New value of Ci
                Ci = ambient_CO2 - An * ((1.6/gs) + (1.37/cls.GB)) # gs and GB in mol m-2 s-1
                # New value of Tleaf
                Tleaf, Tr = cls._leaf_temperature(LEAF_WIDTH, H_ORGAN, H_CANOPY, wind0, PAR, gs, Ta, Tleaf, RH)
                count +=1

        return An, Tr

