#!/usr/bin/env python
# -*- coding: latin-1 -*-

from __future__ import division # use '//' to do integer division

"""
    farquharwheat.model
    ~~~~~~~~~~~~~~~~~~~

    Model of photosynthesis based on Farquhar's approach.
    The model includes the dependence of photosynthesis to organ temperature and nitrogen content.
    Internal CO2 and organ temperature are found numerically.

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

class Model(object):

    O = 21000      #: Photosynthetic parameter: Intercellular O2 concentration, µmol mol(air)-1 or Pa, from Bernacchi et al. (2001)
    KC25 = 404     #: Photosynthetic parameter: Affinity constant of RuBisCO for C, µmol mol-1 or Pa, from Bernacchi et al. (2001) (estimation in Braune et al. (2009) not enough accurate)
    KO25 = 278.4E3 #: Photosynthetic parameter: Affinity constant of RuBisCO for O, µmol mol-1 or Pa, from Bernacchi et al. (2001) (estimation in Braune et al. (2009) not enough accurate)
    GAMMA25 = 39   #: Photosynthetic parameter: CO2 compensation point, µmol(CO2) mol-1 (air), from Braune et al. (2009)
    THETA = 0.72   #: Photosynthetic parameter: curvature parameter of J, dimensionless

    MM_WATER = 18  #: Molar mass of water (g mol-1)

    #: Nitrogen dependance of photosynthetic parameters (derived from Braune et al. (2009) and Evers et al. (2010):
    #:     * S_surfacic_nitrogen: slope of the relation between surfacic_nitrogen and the parameter
    #:         * alpha: mol e- m2 mol-1 photon g-1 N
    #:         * Vc_max25: µmol CO2 g-1 N s-1
    #:         * Jmax25: µmol e- g-1 N s-1
    #:         * TPU25: µmol CO2 g-1 N s-1
    #:         * Rdark25: µmol CO2 g-1 N s-1
    #:     * surfacic_nitrogen_min: minimum amount of nitrogen below which photosynthesis rate is zero (g (N) m-2)
    #:     * beta: intercept parameter of the relation between alpha and surfacic_nitrogen (mol e- mol-1 photons)
    #:     * delta1 and delta2: parameters of m (scaling factor of gs) dependance to surfacic_nitrogen (m2 g-1 and dimensionless respectively)

    PARAM_N = {'S_surfacic_nitrogen': {'Vc_max25': 84.965, 'Jmax25': 117.6, 'alpha': 0.0413, 'TPU25': 9.25, 'Rdark25': 0.493}, 'surfacic_nitrogen_min': {'Vc_max25': 0.17, 'Jmax25': 0.17, 'TPU25': 0.229, 'Rdark25': 0.118},
               'beta': 0.2101, 'delta1': 14.7, 'delta2': -0.548}
    NA_0 = 2                #: Initial value of surfacic_nitrogen (g m-2), used if no surfacic_nitrogen is provided by user

    GSMIN = 0.05            #: Stomatal conductance parameter: Minimum gsw, measured in the dark (mol m-2 s-1). Braune et al. (2009).
    GB = 3.5                #: Stomatal conductance parameter: Boundary layer conductance to water vapour (mol m-2 s-1). Muller et al., (2005)

    A = 2.5                 #: Physical parameter: Attenuation coefficient of wind within a wheat canopy. From Campbell and Norman (1998), second edition. Can also be estimated by: A = sqrt((0.2*LAI*h)/sqrt((4*width*h)/(pi*LAI))
    GAMMA = 66E-3           #: Physical parameter: Psychrometric constant (KPa K-1). Mean value
    I0 = 1370               #: Physical parameter: Extraterrestrial solar radiation (W m-2)
    K =0.40                 #: Physical parameter: Von Kármán's constant (dimensionless)
    LAMBDA = 2260E3         #: Physical parameter: Latent heat for vaporisation of water (J kg-1)
    RHOCP = 1256            #: Physical parameter: Volumetric heat capacity of air (J m-3 K-1)
    SIGMA = 5.6704E-8       #: Physical parameter: Stefan-Bolzmann constant (W-2 K-4)
    ZR = 2                  #: Physical parameter: Height above canopy at which reference wind (Ur) is measured (m)

    R = 8.3144              #: Physical parameter: Gas constant (J mol-1 K-1)
    PATM = 1.01325E5        #: Physical parameter: Atmospheric pressure (Pa)

    FR = 0.1                #: Physical parameter: radiation reflectance
    FT = 0.1                #: Physical parameter: radiation transmittance

    #: Temperature dependance of photosynthetic parameters (parameter values derived from Braune et al. (2009) except for Kc, Ko, and Rdark (Bernacchi et al., 2001))
    #:     * deltaHa, deltaHd: enthalpie of activation and deactivation respectively (kJ mol-1)
    #:     * deltaS: entropy term (kJ mol-1 K-1)
    #:     * Tref: reference temperature (K)
    #:     * R: universal gas constant (kJ mol-1 K-1)

    PARAM_TEMP = {'deltaHa': {'Vc_max': 89.7, 'Jmax': 48.9, 'TPU': 47., 'Kc': 79.43, 'Ko': 36.38, 'Gamma': 35., 'Rdark': 46.39},
                  'deltaHd': {'Vc_max': 149.3, 'Jmax': 152.3, 'TPU': 152.3},
                  'deltaS': {'Vc_max': 0.486, 'Jmax': 0.495, 'TPU': 0.495},
                  'Tref': 298.15, 'R': 8.3145E-03}
    KELVIN_DEGREE = 273.15                #: Conversion factor from degree C to Kelvin

    DELTA_CONVERGENCE = 0.01 #: The relative delta for Ci and Ts convergence.

    @classmethod
    def _organ_temperature(cls, w, z, Zh, Ur, PAR, gsw, Ta, Ts, RH, organ_name):
        """
        Energy balance for the estimation of organ temperature

        :Parameters:
            - `w` (:class:`float`) - organ characteristic dimension (m) to be considered for heat transfer through forced convection (by wind).
                 For a leaf: its width (more related to wind direction than length), for cylindric stem elements: diameter.
            - `z` (:class:`float`) - organ height from soil (m)
            - `Zh` (:class:`float`) - canopy height (m)
            - `Ur` (:class:`float`) - wind speed (m s-1) at the reference height (zr), e.g. top of the canopy + 2m (in the case of wheat, Ur can be approximated as the wind speed at 2m from soil)
            - `PAR` (:class:`float`) - absorbed PAR(µmol m-2 s-1)
            - `gsw` (:class:`float`) - stomatal conductance to water vapour (mol m-2 s-1)
            - `Ta` (:class:`float`) - air temperature (degree C)
            - `Ts` (:class:`float`) - organ temperature (degree C). Ts = Ta at the first iteration of the numeric resolution
            - `RH` (:class:`float`) - Relative humidity (decimal fraction)
            - `organ_name` (:class:`string`) - organ name (used to distinguish lamina from cylindric organs)

        :Returns:
            Ts (organ temperature, degree C), Tr (organ transpiration rate, mm s-1)
        :Returns Type:
            :class:`float`
        """

        d = 0.7*Zh                                          #: Zero plane displacement height (m)
        Zo = 0.1*Zh                                         #: Roughness length (m)

        #: Wind speed
        u_star = (Ur * cls.K) / log((cls.ZR - d)/Zo)        #: Friction velocity (m s-1)
        Uh = (u_star/cls.K) * log((Zh-d)/Zo)                #: Wind speed at the top of canopy (m s-1)
        u = Uh * exp(cls.A*(z/Zh -1))                       #: Wind speed at organ height (m s-1), from Campbell and Norman (1998), second edition.

        #: Boundary layer resistance to heat (s m-1). See Finnigan J, Raupach M. 1987 and Monteith JL. 1973 for basic equations.
        if organ_name == 'Lamina':
            rbh = 154*sqrt(w/u)                             #: Case of horizontal planes submitted to forced convection
        else:
            rbh = w / (1.2E-5 * ((u*w)/1.5E-5)**0.47)       #: Case of vertical cylinders submitted to forced convection

        #: Turbulence resistance to heat (s m-1)
        ra = 1/(cls.K**2 * Ur) * (log((cls.ZR - d)/Zo))**2  #: Aerodynamic resistance integrated from zr to z0 + d

        #: Net absorbed radiation Rn (PAR and NIR, J m-2 s-1)
        Iabs = (PAR*(1-cls.FR-cls.FT))/(0.55*4.55)          #: Global absorbed radiation by organ (J m-2 s-1). TODO: relation a verifier
        es_Ta = 0.611 * exp((17.4*Ta)/(239+Ta))             #: Saturated vapour pressure of the air (kPa), Ta in degree Celsius
        V = RH * es_Ta                                      #: Vapour pressure of the air (kPa)
        fvap = 0.56 - 0.079*sqrt(10*V)                      #: Fraction of vapour pressure

        tau = Iabs/cls.I0                                   #: Atmospheric transmissivity (dimensionless)
        fclear = 0.1 + 0.9*max(0, min(1, (tau-0.2)/0.5))    #: Fraction sky clearness

        Rn = Iabs - cls.SIGMA * (Ts + cls.KELVIN_DEGREE)**4*fvap*fclear

        #: Transpiration (mm s-1), Penman-Monteith
        if Ts == Ta:
            Ta_K = Ta + cls.KELVIN_DEGREE
            s = ((17.4*239)/(Ta_K + 239)**2)*es_Ta          #: Slope of the curve relating saturation vapour pressure to temperature (kPa K-1)
        else:
            es_Tl = 0.611 * exp((17.4*Ts)/(239+Ts))     #: Saturated vapour pressure at organ level (kPa), Ts in degree Celsius
            Ts_K, Ta_K = Ts + cls.KELVIN_DEGREE, Ta + cls.KELVIN_DEGREE
            s = (es_Tl - es_Ta)/(Ts - Ta_K)               #: Slope of the curve relating saturation vapour pressure to temperature (kPa K-1)

        VPDa = es_Ta - V
        rbw = 0.96 * rbh                                                   #: Boundary layer resistance for water (s m-1)
        gsw_physic = (gsw * cls.R * (Ts+cls.KELVIN_DEGREE)) / cls.PATM   #: Stomatal conductance to water in physical units (m s-1). Relation given by A. Tuzet (2003)
        rswp = 1/gsw_physic                                                #: Stomatal resistance for water (s m-1)

        Tr = max(0, (s * Rn + (cls.RHOCP * VPDa)/(rbh + ra)) / (cls.LAMBDA * (s + cls.GAMMA*((rbw + ra + rswp)/(rbh + ra))))) # mm s-1

        #: Organ temperature
        Ts = Ta + ((rbh + ra) * (Rn - cls.LAMBDA*Tr)) / cls.RHOCP

        return Ts, Tr

    @classmethod
    def _stomatal_conductance(cls, Ag, An, surfacic_nitrogen, ambient_CO2, RH):
        """
        Ball, Woodrow, and Berry model of stomatal conductance (1987)

        :Parameters:
            - `Ag` (:class:`float`) - gross assimilation rate (µmol m-2 s-1)
            - `An` (:class:`float`) - net assimilation rate (µmol m-2 s-1)
            - `surfacic_nitrogen` (:class:`float`) - surfacic nitrogen content(g m-2)
            - `ambient_CO2` (:class:`float`) - Air CO2 (µmol mol-1)
            - `RH` (:class:`float`) - Relative humidity (decimal fraction)
        :Returns:
            gsw (mol m-2 s-1)
        :Returns Type:
            :class:`float`
        """

        Cs = ambient_CO2 - An *(1.37/(cls.GB))                   #: CO2 concentration at organ surface (µmol mol-1 or Pa). From Prieto et al. (2012). GB in mol m-2 s-1
        m = cls.PARAM_N['delta1'] * surfacic_nitrogen**cls.PARAM_N['delta2']    #: Scaling factor dependance to surfacic_nitrogen (dimensionless). This focntion is maintained although I'm not conviced that it should be taken into account
        gsw = (cls.GSMIN + m*((Ag*RH)/(Cs)))                     #: Stomatal conductance to water vapour (mol m-2 s-1), from Braune et al. (2009), Muller et al. (2005): using Ag rather than An. Would be better with a function of VPD and with (Ci-GAMMA) instead of Cs.
        return gsw

    @classmethod
    def _calculate_Ci(cls, ambient_CO2, An, gsw):
        """
        Calculates the internal CO2 concentration (Ci)

        :Parameters:
            - `ambient_CO2` (:class:`float`) - air CO2 (µmol mol-1)
            - `An` (:class:`float`) -net assimilation rate of CO2 (µmol m-2 s-1)
            - `gsw` (:class:`float`) - stomatal conductance to water vapour (mol m-2 s-1)
        :Returns:
            Ci (µmol mol-1)
        :Returns Type:
            :class:`float`
        """
        Ci = ambient_CO2 - An * ((1.6/gsw) + (1.37/cls.GB)) #: Intercellular concentration of CO2 (µmol mol-1)
        # gsw and GB in mol m-2 s-1 so that  (An * ((1.6/gs) + (1.37/cls.GB)) is thus in µmol mol-1 as ambient_CO2
        # 1.6 converts gsw to gs_CO2, and 1.37 comes from (1.6)^(2/3)
        return Ci

    @classmethod
    def _f_temperature(cls, pname, p25, T):
        """
        Photosynthetic parameters relation to temperature

        :Parameters:
            - `pname` (:class:`string`) - name of parameter
            - `p25` (:class:`float`) - parameter value at 25 degree C
            - `T` (:class:`float`) - organ temperature (degree C)
        :Returns:
            p (parameter value at organ temperature)
        :Returns Type:
            :class:`float`
        """
        Tk = T + cls.KELVIN_DEGREE
        deltaHa = cls.PARAM_TEMP['deltaHa'][pname]                  #: Enthalpie of activation of parameter pname (kJ mol-1)
        Tref = cls.PARAM_TEMP['Tref']
        R = cls.PARAM_TEMP['R']

        f_activation = exp((deltaHa * (Tk - Tref))/(R * Tref * Tk)) #: Energy of activation (normalized to unity)

        if pname in ('Vc_max', 'Jmax', 'TPU'):
            deltaS = cls.PARAM_TEMP['deltaS'][pname]                #: entropy term of parameter pname (kJ mol-1 K-1)
            deltaHd = cls.PARAM_TEMP['deltaHd'][pname]              #: Enthalpie of deactivation of parameter pname (kJ mol-1)
            f_deactivation = (1 + exp((Tref*deltaS - deltaHd) / (Tref*R))) / (1 + exp((Tk*deltaS - deltaHd) / (Tk*R))) #: Energy of deactivation (normalized to unity)
        else:
            f_deactivation = 1

        p = p25 * f_activation * f_deactivation

        return p

    @classmethod
    def calculate_photosynthesis(cls, PAR, surfacic_nitrogen, Ts, Ci):
        """
        Computes photosynthesis rate following Farquhar's model with regulation by organ temperature and nitrogen content.
        In this version, most of parameters are derived from Braune et al. (2009) on barley and Evers et al. (2010) for N dependencies.

        :Parameters:
            - `PAR` (:class:`float`) - PAR absorbed (µmol m-2 s-1)
            - `surfacic_nitrogen` (:class:`float`) - surfacic nitrogen content(g m-2)
            - `Ts` (:class:`float`) - organ temperature (degree C)
            - `Ci` (:class:`float`) - internal CO2 (µmol mol-1), Ci = 0.7*CO2air for the first iteration
        :Returns:
            Ag (µmol m-2 s-1), An (µmol m-2 s-1), Rd (µmol m-2 s-1)
        :Returns Type:
            :class:`float`
        """

        #: RuBisCO parameters dependance to temperature
        Kc = cls._f_temperature('Kc', cls.KC25, Ts)
        Ko = cls._f_temperature('Ko', cls.KO25, Ts)
        Gamma = cls._f_temperature('Gamma', cls.GAMMA25, Ts)

        #: RuBisCO-limited carboxylation rate
        Sna_Vcmax25 = cls.PARAM_N['S_surfacic_nitrogen']['Vc_max25']
        surfacic_nitrogen_min_Vcmax25 = cls.PARAM_N['surfacic_nitrogen_min']['Vc_max25']
        Vc_max25 = Sna_Vcmax25 * (surfacic_nitrogen - surfacic_nitrogen_min_Vcmax25)                                      #: Relation between Vc_max25 and surfacic_nitrogen (µmol m-2 s-1)
        Vc_max = cls._f_temperature ('Vc_max', Vc_max25, Ts)                              #: Relation between Vc_max and temperature (µmol m-2 s-1)
        Ac = (Vc_max * (Ci-Gamma)) / (Ci + Kc * (1 + cls.O/Ko))                             #: Rate of assimilation under Vc_max limitation (µmol m-2 s-1)

        #: RuBP regeneration-limited carboxylation rate via electron transport
        ALPHA = cls.PARAM_N['S_surfacic_nitrogen']['alpha'] * surfacic_nitrogen + cls.PARAM_N['beta']                     #: Relation between ALPHA and surfacic_nitrogen (mol e- mol-1 photon)
        Sna_Jmax25 = cls.PARAM_N['S_surfacic_nitrogen']['Jmax25']
        surfacic_nitrogen_min_Jmax25 = cls.PARAM_N['surfacic_nitrogen_min']['Jmax25']
        Jmax25 = Sna_Jmax25 * (surfacic_nitrogen - surfacic_nitrogen_min_Jmax25)                                          #: Relation between Jmax25 and surfacic_nitrogen (µmol m-2 s-1)
        Jmax = cls._f_temperature('Jmax', Jmax25, Ts)                                     #: Relation between Jmax and temperature (µmol m-2 s-1)

        J = ((Jmax+ALPHA*PAR) - sqrt((Jmax+ALPHA*PAR)**2 - 4*cls.THETA*ALPHA*PAR*Jmax))/(2*cls.THETA) #: Electron transport rate (Muller et al. (2005), Evers et al. (2010)) (µmol m-2 s-1)
        Aj = (J * (Ci-Gamma)) / (4*Ci + 8*Gamma)                                            #: Rate of assimilation under RuBP regeneration limitation (µmol m-2 s-1)

        #: Triose phosphate utilisation-limited carboxylation rate
        Sna_TPU25 = cls.PARAM_N['S_surfacic_nitrogen']['TPU25']
        surfacic_nitrogen_min_TPU25 = cls.PARAM_N['surfacic_nitrogen_min']['TPU25']
        TPU25 = Sna_TPU25 * (surfacic_nitrogen - surfacic_nitrogen_min_TPU25)                                             #: Relation between TPU25 and surfacic_nitrogen (µmol m-2 s-1)
        TPU = cls._f_temperature('TPU', TPU25, Ts)                                        #: Relation between TPU and temperature (µmol m-2 s-1)
        Vomax = (Vc_max*Ko*Gamma)/(0.5*Kc*cls.O)                                            #: Maximum rate of Vo (µmol m-2 s-1) (µmol m-2 s-1)
        Vo = (Vomax * cls.O) / (cls.O + Ko*(1+Ci/Kc))                                       #: Rate of oxygenation of RuBP (µmol m-2 s-1)
        Ap = (1-Gamma/Ci)*(3*TPU + Vo)                                                      #: Rate of assimilation under TPU limitation (µmol m-2 s-1). I think there was a mistake in the paper of Braune t al. (2009) where they wrote Ap = (1-Gamma/Ci)*(3*TPU) + Vo
        # A more recent expression of Ap was given by S. v Caemmerer in her book (2000): AP = (3TPU * (Ci-Gamma))/(Ci-(1+3alpha)*Gamma),
        # where 0 < alpha > 1 is the fraction of glycolate carbon not returned to the chloroplast, but I couldn't find any estimation of alpha for wheat

        #: Gross assimilation rate (µmol m-2 s-1)
        Ag = min(Ac, Aj, Ap)

        #: Mitochondrial respiration rate of organ in light Rd (processes other than photorespiration)
        Rdark25 = cls.PARAM_N['S_surfacic_nitrogen']['Rdark25'] * (surfacic_nitrogen - cls.PARAM_N['surfacic_nitrogen_min']['Rdark25'])  #: Relation between Rdark25 (respiration in obscurity at 25 degree C) and surfacic_nitrogen (µmol m-2 s-1)
        Rdark = cls._f_temperature('Rdark', Rdark25, Ts)                                  #: Relation between Rdark and temperature (µmol m-2 s-1)
        Rd = Rdark * (0.33 + (1-0.33)*(0.5)**(PAR/15))                                      # Found in Muller et al. (2005), eq. 19 (µmol m-2 s-1)

        #: Net C assimilation (µmol m-2 s-1)
        if Ag <= 0: # Occurs when Ci is lower than Gamma, in these cases there is no net assimilation (Farquhar, 1980; Caemmerer, 2000), so I set Ag = Rd
            Ag = Rd
        An = Ag - Rd

        return Ag, An, Rd


    @classmethod
    def calculate_An(cls, surfacic_nitrogen, width, height, PAR, Ta, ambient_CO2, RH, Ur, organ_name):
        """
        quelle doc?

        :Parameters:
            - `surfacic_nitrogen` (:class:`float`) - total surfacic nitrogen content of organs (g m-2), obtained by the sum of nitrogen, amino acids, proteins and structural N.
               Properly speaking, photosynthesis should be related to proteins (RubisCO), but parameters of most Farquhar models are calibrated on total N measurements (DUMAS method).
               If None, surfacic_nitrogen = :attr:`NA_0`

            - `width` (:class:`float`) - width of the organ (or diameter for stem organ) (m),
               characteristic dimension to be considered for heat transfer through forced convection (by wind).

            - `height` (:class:`float`) - height of the organ from soil (m)

            - `PAR` (:class:`float`) - absorbed PAR (µmol m-2 s-1)

            - `Ta` (:class:`float`) - air temperature (°C)

            - `ambient_CO2` (:class:`float`) - air CO2 (µmol mol-1)

            - `RH` (:class:`float`) - relative humidity (decimal fraction)

            - `Ur` (:class:`float`) - Ur: wind at the reference height (zr) (m s-1), e.g. top of the canopy + 2m
               (in the case of wheat, Ur can be approximated as the wind speed at 2m from soil)

            - `organ_name` (:class:`string`) - type of the organ (e.g. "Lamina", "Sheath", etc.)

        :Returns:
            Ag (µmol m-2 s-1), An (µmol m-2 s-1), Rd (µmol m-2 s-1), Tr (mmol m-2 s-1), Ts (°C) and  gsw (mol m-2 s-1)

        :Returns Type:
            :class:`float`
        """

        if surfacic_nitrogen is None:
            surfacic_nitrogen = cls.NA_0

        #: Organ parameters
        H_CANOPY = 0.8                              #: m, temporary

        ### Iterations to find organ temperature and Ci ###
        Ci, Ts = 0.7*ambient_CO2, Ta # Initial values
        count = 0

        while True:
            prec_Ci, prec_Ts = Ci, Ts
            Ag, An, Rd = cls.calculate_photosynthesis(PAR, surfacic_nitrogen, Ts, Ci)
            # Stomatal conductance to water
            gsw = cls._stomatal_conductance(Ag, An, surfacic_nitrogen, ambient_CO2, RH)

            # New value of Ci
            Ci = cls._calculate_Ci(ambient_CO2, An, gsw)

            # New value of Ts
            Ts, Tr = cls._organ_temperature(width, height, H_CANOPY, Ur, PAR, gsw, Ta, Ts, RH, organ_name)
            count +=1

            if count >=30: # TODO: test a faire? Semble prendre du tps de calcul
                if abs((Ci - prec_Ci)/prec_Ci) >= cls.DELTA_CONVERGENCE:
                    print '{}, Ci cannot converge, prec_Ci= {}, Ci= {}'.format(organ_name, prec_Ci, Ci)
                if abs((Ts - prec_Ts)/prec_Ts) >= cls.DELTA_CONVERGENCE:
                    print '{}, Ts cannot converge, prec_Ts= {}, Ts= {}'.format(organ_name, prec_Ts, Ts)
                break
                #Ci = cls._calculate_Ci(ambient_CO2, An, gsw)

            if abs((Ci - prec_Ci)/prec_Ci) < cls.DELTA_CONVERGENCE and abs((Ts - prec_Ts)/prec_Ts) < cls.DELTA_CONVERGENCE:
                break

        #: Conversion of Tr from mm s-1 to mmol m-2 s-1 (more suitable for further use of Tr)
        Tr = (Tr * 1E6) / cls.MM_WATER # Using 1 mm = 1kg m-2

        if organ_name != 'Lamina':
            Ag = Ag*0.75
        return Ag, An, Rd, Tr, Ts, gsw