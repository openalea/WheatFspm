#!/usr/bin/env python
# -*- coding: latin-1 -*-

from __future__ import division # use '//' to do integer division

"""
    respiwheat.model
    ~~~~~~~~~~~~~~~~~~~

    Model of respiration based on on Thornley and Cannell, 2000.
    The model computes the respiration associated with the main biological processes.

    R_total = sum(R_growth) + R_phloem + R_Namm_upt + R_Nnit_upt + R_Nnit_red(shoot + root) + R_N2fix + R_min_upt + sum(R_residual)

    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.
"""

class RespirationModel(object):

    ### R_growth###
    Yg = 0.75             # Growth yield (units of C appearing in new biomass per unit of C substrate utilized for growth)

    ### R_phloem###
    Cphloem = 0.006       # Units C respired per unit C substrate loaded into the phloem

    ### R_Namm_upt###
    C_Namm_upt = 0.198    # µmol of C substrate respired per µmol of N ammonium taken up

    ### R_Nnit_upt###
    C_Nnit_upt = 0.397    # µmol of C substrate respired per µmol of N nitrates taken up

    ### R_Nnit_red###
    fNnit_red_sh_CS = 0.5 # fraction of nitrate reduced in the shoot using C substrate rather than using excess ATP and reducing power obtained directly from photosynthesis
    C_Nnit_red = 1.98     # µmol of C substrate per µmol of N nitrates reduced

    ### R_N2fix ###
    C_Nfix = 6            # kg substrate C respired (kg N fixed)-1 (in the range  5 to 12)

    ### R_min_upt ###
    Cmin_upt = 5000       # µmol of C substrate respired per g of minerals taken up
    ashe_content = 0.05   # g minerals per g of structural dry mass

    ### R_residual ###
    km_max = 4.1E-6       # Maximum value of the maintenance constant when C is much greater than Km (µmol of C substrate respired per µmol N s-1)
    Km = 1.67E3           # The Michaelis-Menten constant affinity i.e. the C substrate concentration at half the value of km_max (µmol of C substrate per g of structural mass)

    @classmethod
    def R_growth(cls, G):
        """
        Local growth respiration
        - G: gross growth of plant tissue (µmol C added in new plant tissue)
        """
        Rg = ((1 - cls.Yg)/cls.Yg)*G
        return Rg

    @classmethod
    def R_phloem(cls, L, sucrose):
        """
        Phloem loading respiration
        - L: Loading flux from the C substrate pool to phloem (µmol C)
        """
        if sucrose > 0:
            Rp = cls.Cphloem * L
        else:
            Rp = 0
        return Rp

    @classmethod
    def R_Namm_upt(cls, U_Namm):
        """
        Ammonium uptake respiration
        - U_Namm: uptake of N ammonium (µmol N)
        Neglected for wheat model
        """
        R_Namm = cls.C_Namm_upt * U_Namm
        return R_Namm

    @classmethod
    def R_Nnit_upt(cls, U_Nnit, sucrose):
        """
        Nitrate uptake respiration
        - U_Nnit: uptake of N nitrates (µmol N)
        """
        if sucrose > 0:
            R_Nnit_upt_ = cls.C_Nnit_upt * U_Nnit
        else:
            R_Nnit_upt_ = 0
        return R_Nnit_upt_

    @classmethod
    def R_Nnit_red(cls, Nit_red, sucrose, root=False):
        """
        Nitrate reduction-linked respiration
        - Nit_red: µmol of N nitrates reduced (supposed to equal consumption of N for the synthesis of amino acids)
        - root: specifies if the nitrate reduction-linked respiration is computed for shoot (False) or root (True) tissues.
        Distinction is made between nitrate realised in roots or in shoots where a part of the energy required is derived from ATP
        and reducing power obtained directly from photosynthesis (rather than C substrate)
        """
        if sucrose > 0:
            if not root:
                R_Nnit_red_ = cls.fNnit_red_sh_CS * cls.C_Nnit_red *  Nit_red #: Respiration in shoot tissues
            else:
                R_Nnit_red_ = cls.C_Nnit_red * Nit_red                         #: Respiration in root tissues
        else:
            R_Nnit_red_ = 0

        return R_Nnit_red_

    @classmethod
    def R_N2fix(cls, I_Nfix):
        """
        N2-fixation respiration
        - I_Nfix: flux of fixed N into the root substrate N pool (kg fixed N)
        Not used in wheat model
        """

        R_N2 = cls.C_Nfix * I_Nfix
        return R_N2

    @classmethod
    def R_min_upt(cls, delta_BMstruct):
        """
        Mineral ion (other than N) uptake-linked respiration
        - delta_BMstruct: gross plant structural growth rate  (g)
        """
        # Uptake of minerals (g)
        Umin = (cls.ashe_content * delta_BMstruct) * 0.5 #: 0.5: minerals internally recycled from senescing leaves
        # Respiratory cost
        Rmin = cls.Cmin_upt * Umin
        return Rmin
    # c bien delta_BMstruct? (cls.ashe_content * created_tissues_BMstruct) - (cls.ashe_content * senesced_tissues_Mstruct)*0.5
    # senesced_part_Mstruct = Mstruct des tissus qui deviennet senescents (progression d'un 'front de senscence' sur les organes)
    # a chq pas de tps, calcul delta_BMstruct, puis C respiré depuis sucrose phloem?: OK car uptake + transport
    @classmethod
    def R_residual(cls, C, Ntot):
        """
        Residual maintenance respiration (cost from protein turn-over, cell ion gradients, futile cycles...)
        - C: C substrate concentration (µmol C g-1 structural mass)
        - Ntot: total N in plant (µmol N)= structural N + meristematic N + substrate N + photosynthetic N
        """
        if C > 0:
            Rres = ((cls.km_max * C)/(cls.Km + C)) * Ntot
        else:
            Rres = 0
        return Rres
# A chq pas de tps, faire ce calcul pr chq org_ph, roots et grains?

#Prise en compte température:
#    - pr la croissance: dépendance temp pr RER et beta (=ttes les vitesses pr les dimensions) type Arrhenius modifiee pr diminuer a 1 temp. Va ainsi affecter les flux de masse
#    - pr respi residuelle: mettre Arrhenius classique (chercher ds biblio)
#    - respi phloeme: a reguler + tard a travers les flux de C chargé