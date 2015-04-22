#!/usr/bin/env python
# -*- coding: latin-1 -*-

from __future__ import division # use '//' to do integer division

"""
    respiwheat.model
    ~~~~~~~~~~~~~~~~~~~

    Model of respiration based on Thornley and Cannell, 2000.
    The model computes the respiration associated with the main biological processes.

    R_total = sum(R_growth) + R_phloem + R_Namm_upt + R_Nnit_upt + R_Nnit_red(shoot + root) + R_N2fix + R_min_upt + sum(R_residual)

    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.
"""

class RespirationModel(object):

    ### R_growth###
    YG = 0.80             # Growth yield (units of C appearing in new biomass per unit of C substrate utilized for growth)

    ### R_phloem###
    CPHLOEM = 0.006       # Units C respired per unit C substrate loaded into the phloem

    ### R_Namm_upt###
    C_AMM_UPT = 0.198     # µmol of C substrate respired per µmol of N ammonium taken up

    ### R_Nnit_upt###
    C_NIT_UPT = 0.397     # µmol of C substrate respired per µmol of N nitrates taken up

    ### R_Nnit_red###
    F_NIT_RED_SH_CS = 0.5 # fraction of nitrate reduced in the shoot using C substrate rather than using excess ATP and reducing power obtained directly from photosynthesis
    C_NIT_RED = 1.98      # µmol of C substrate per µmol of N nitrates reduced

    ### R_N2fix ###
    C_NFIX = 6            # kg substrate C respired (kg N fixed)-1 (in the range  5 to 12)

    ### R_min_upt ###
    CMIN_UPT = 5000       # µmol of C substrate respired per g of minerals taken up
    ASHE_CONTENT = 0.05   # g minerals per g of structural dry mass

    ### R_residual ###
    KM_MAX = 4.1E-6       # Maximum value of the maintenance constant when C is much greater than KM (µmol of C substrate respired per µmol N s-1)
    KM = 1.67E-3          # The Michaelis-Menten constant affinity i.e. the C substrate concentration at half the value of KM_MAX (µmol of C substrate per g of structural mass)

    @classmethod
    def R_growth(cls, G, mstruct):
        """
        Local growth respiration

        : Parameters:
            - `G` (:class:`float`) - gross growth of plant tissue (µmol C added in new plant tissue g-1 mstruct)
            - `mstruct` (:class:`float`) -  structural dry mass of organ (g)

        : Returns:
            _R_growth (µmol C respired)

        :Returns Type:
            :class:`float`
        """
        _R_growth = ((1 - cls.YG)/cls.YG) * (G*mstruct)
        return _R_growth

    @classmethod
    def R_grain_growth(cls, mstruct_growth, starch_filling, mstruct):
        """
        Grain growth respiration

        : Parameters:
            - `mstruct_growth` (:class:`float`) - gross growth of grain structure (µmol C added in grain structure)
            - `starch_filling` (:class:`float`) - gross growth of grain starch (µmol C added in grain starch g-1 mstruct)
            - `mstruct` (:class:`float`) -  structural dry mass of organ (g)

        : Returns:
            R_grain_growth (µmol C respired)

        :Returns Type:
            :class:`float`
        """
        R_grain_growth_struct = ((1 - cls.YG)/cls.YG) * mstruct_growth
        R_grain_growth_starch = ((1 - cls.YG)/cls.YG) * (starch_filling*mstruct)
        return R_grain_growth_struct, R_grain_growth_starch

    @classmethod
    def R_phloem(cls, sucrose_loading, sucrose, mstruct):
        """
        Phloem loading respiration

        : Parameters:
            - `sucrose_loading` (:class:`float`) -  Loading flux from the C substrate pool to phloem (µmol C g-1 mstruct)
            - `sucrose` (:class:`float`) -  amount of C sucrose in organ (µmol C)
            - `mstruct` (:class:`float`) -  structural dry mass of organ (g)

        : Returns:
            _R_phloem (µmol C respired)

        :Returns Type:
            :class:`float`
        """

        if sucrose_loading >0 and sucrose >0:
            _R_phloem = cls.CPHLOEM * sucrose_loading * mstruct
        else:
            _R_phloem = 0
        return _R_phloem

    @classmethod
    def R_Namm_upt(cls, U_Namm):
        """
        Ammonium uptake respiration

        : Parameters:
            - `U_Namm` (:class:`float`) -  uptake of N ammonium (µmol N)

        : Returns:
            _R_Namm (µmol C respired)

        :Returns Type:
            :class:`float`
        """
        R_Namm_upt = cls.C_AMM_UPT * U_Namm
        return R_Namm_upt

    @classmethod
    def R_Nnit_upt(cls, U_Nnit, sucrose):
        """
        Nitrate uptake respiration

        : Parameters:
            - `U_Nnit` (:class:`float`) - uptake of N nitrates (µmol N)
            - `sucrose` (:class:`float`) -  amount of C sucrose in organ (µmol C)

        : Returns:
            _R_Nnit_upt (µmol C respired)

        :Returns Type:
            :class:`float`
        """
        if sucrose > 0:
            _R_Nnit_upt = cls.C_NIT_UPT * U_Nnit
        else:
            _R_Nnit_upt = 0
        return _R_Nnit_upt

    @classmethod
    def R_Nnit_red(cls, s_amino_acids, sucrose, mstruct, root=False):
        """
        Nitrate reduction-linked respiration
        Distinction is made between nitrate realised in roots or in shoots where a part of the energy required is derived from ATP
        and reducing power obtained directly from photosynthesis (rather than C substrate)

        : Parameters:
            - `s_amino_acids` (:class:`float`) - consumption of N for the synthesis of amino acids (µmol N g-1 mstruct)
              (in the present version, this is used to approximate nitrate reduction needed in the original model of Thornley and Cannell, 2000)
            - `sucrose` (:class:`float`) -  amount of C sucrose in organ (µmol C)
            - `mstruct` (:class:`float`) -  structural dry mass of organ (g)
            - `root` (:class:`bool`) - specifies if the nitrate reduction-linked respiration is computed for shoot (False) or root (True) tissues.

        : Returns:
            _R_Nnit_upt (µmol C respired)

        :Returns Type:
            :class:`float`
        """
        if sucrose > 0:
            if not root:
                _R_Nnit_red = cls.F_NIT_RED_SH_CS * cls.C_NIT_RED *  s_amino_acids * mstruct # Respiration in shoot tissues
            else:
                _R_Nnit_red = cls.C_NIT_RED * s_amino_acids * mstruct                        # Respiration in root tissues
        else:
            _R_Nnit_red = 0
        return _R_Nnit_red

    @classmethod
    def R_N2fix(cls, I_Nfix):
        """
        N2-fixation respiration

        : Parameters:
            - `I_Nfix` (:class:`float`) - flux of fixed N into the root substrate N pool (kg fixed N)

        : Returns:
            _R_N2fix (µmol C respired)

        :Returns Type:
            :class:`float`
        """

        _R_N2fix = cls.C_NFIX * I_Nfix
        return _R_N2fix

    @classmethod
    def R_min_upt(cls, delta_BMstruct):
        """
        Mineral ion (other than N) uptake-linked respiration

        : Parameters:
            - `delta_BMstruct` (:class:`float`) - gross plant structural growth rate (g)

        : Returns:
            _R_min_upt (µmol C respired)

        :Returns Type:
            :class:`float`
        """
        # Uptake of minerals (g)
        Umin = (cls.ASHE_CONTENT * delta_BMstruct) * 0.5 # 0.5: minerals internally recycled from senescing leaves
        # Respiratory cost
        _R_min_upt = cls.CMIN_UPT * Umin
        return _R_min_upt

    @classmethod
    def R_residual(cls, sucrose, mstruct, Ntot, delta_t):
        """
        Residual maintenance respiration (cost from protein turn-over, cell ion gradients, futile cycles...)

        : Parameters:
            - `sucrose` (:class:`float`) - amount of C sucrose (µmol C)
            - `mstruct` (:class:`float`) - structural dry mass of organ (g)
            - `Ntot` (:class:`float`) - total N in plant (µmol N)
            - `delta_t` (:class:`float`) - timestep (s)

        : Returns:
            _R_residual (µmol C respired)

        :Returns Type:
            :class:`float`
        """
        conc_sucrose = sucrose / mstruct
        if conc_sucrose >0 :
            _R_residual = ((cls.KM_MAX * conc_sucrose)/(cls.KM + conc_sucrose)) * Ntot * delta_t
        else:
            _R_residual = 0
        return _R_residual