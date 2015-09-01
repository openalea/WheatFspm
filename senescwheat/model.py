#!/usr/bin/env python
# -*- coding: latin-1 -*-

from __future__ import division # use '//' to do integer division

"""
    senescwheat.model
    ~~~~~~~~~~~~~~~~~~~

    Model of senescenece.

    :copyright: Copyright 2015 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.
"""

import pandas as pd


class SenescenceModel(object):

    N_MOLAR_MASS = 14             #: Molar mass of nitrogen (g mol-1)
    SENESCENCE_ROOTS = 3.5E-7     #: Rate of root turnover at 20°C (s-1). Value coming from Johnson and Thornley (1985), see also Asseng et al. (1997).

    @classmethod
    def calculate_forced_relative_delta_green_area(cls, green_area_df, group_id, prev_green_area):
        """relative green_area variation due to senescence

        : Parameters:
            - `green_area_df` (:class:`pandas.core.frame.DataFrame`) - a pandas DataFrame containing the green area values for each photosynthetic element at each time
            - `group_id` (:class:`tuple`) - the group id to be used to select data in the DataFrame
            - `prev_green_area` (:class:`float`) - previous value of an organ green area (m-2)

        : Returns:
            new_green_area (m-2), relative_delta_green_area (dimensionless)

        :Returns Type:
            :class:`float`
        """
        new_green_area = green_area_df.get_group(group_id).green_area.values[0]
        relative_delta_green_area = (prev_green_area - new_green_area) / prev_green_area
        return new_green_area, relative_delta_green_area


    @classmethod
    def calculate_relative_delta_green_area(cls, group_id, prev_green_area, proteins, max_proteins, DELTA_T):
        """relative green_area variation due to senescence

        : Parameters:
            - `group_id` (:class:`tuple`) - the group id to be used to select data in the DataFrame
            - `prev_green_area` (:class:`float`) - previous value of an organ green area (m-2)

        : Returns:
            new_green_area (m-2), relative_delta_green_area (dimensionless)

        :Returns Type:
            :class:`float`
        """
        senescence_max_rate = 0.5E-8 # maximal senescence m² s-1
        fraction_N_max = 0.55

        # First run
        if max_proteins.has_key(group_id)==False:
            max_proteins[group_id] = proteins
            new_green_area = prev_green_area
            relative_delta_green_area = 0
        # Overwrite max proteins
        elif max_proteins[group_id] < proteins:
            max_proteins[group_id] = proteins
            new_green_area = prev_green_area
            relative_delta_green_area = 0
        # Senescence if (actual proteins/max_proteins) < fraction_N_max
        elif (proteins / max_proteins[group_id]) < fraction_N_max :
            senesced_area = min(prev_green_area, senescence_max_rate * DELTA_T)
            new_green_area = max(0, prev_green_area - senesced_area)
            relative_delta_green_area = senesced_area / prev_green_area
        else:
            new_green_area = prev_green_area
            relative_delta_green_area = 0

        return new_green_area, relative_delta_green_area, max_proteins

    @classmethod
    def calculate_delta_mstruct_shoot(cls, relative_delta_green_area, prev_mstruct, prev_Nstruct):
        """delta of structural mass due to senescence of photosynthetic elements

        : Parameters:
            - `relative_delta_green_area` (:class:`float`) - relative variation of a photosynthetic element green area
            - `prev_mstruct` (:class:`float`) - previous value of an organ structural mass (g)
            - `prev_Nstruct` (:class:`float`) - previous value of an organ structural N (g)

        : Returns:
            new_mstruct (g), new_Nstruct (g)

        :Returns Type:
            :class:`float`
        """
        new_mstruct = prev_mstruct - prev_mstruct*relative_delta_green_area
        new_Nstruct = prev_Nstruct - prev_Nstruct*relative_delta_green_area
        return new_mstruct, new_Nstruct

    @classmethod
    def calculate_remobilisation(cls, metabolite, relative_delta_structure):
        """Metabolite remobilisation due to senescence over delta_t (µmol).
        : Parameters:
            - `relative_delta_structure` (:class:`float`) - could be relative variation of a photosynthetic element green area or relative variation of mstruct
        """
        return metabolite * relative_delta_structure

    @classmethod
    def calculate_surfacic_nitrogen(cls, nitrates, amino_acids, proteins, Nstruct, green_area):
        # TODO: to be moved in a particular model (shouldn't be included in SenescenceModel)

        """Surfacic content of nitrogen

        : Parameters:
            - `nitrates` (:class:`float`) - amount of nitrates (µmol N)
            - `amino_acids` (:class:`float`) - amount of amino_acids (µmol N)
            - `proteins` (:class:`float`) - amount of proteins (µmol N)
            - `Nstruct` (:class:`float`) - structural N (g)
            - `green_area` (:class:`float`) - green area (m-2)

        : Returns:
            Surfacic nitrogen (g m-2)

        :Returns Type:
            :class:`float`
        """
        mass_N_tot = (nitrates + amino_acids + proteins)*1E-6 * cls.N_MOLAR_MASS + Nstruct
        return (mass_N_tot / green_area)

    # Roots
    @classmethod
    def calculate_roots_mstruct_growth(cls, sucrose, amino_acids, mstruct, DELTA_T):
        # TODO: to be moved in a particular model (shouldn't be included in SenescenceModel)

        """Root structural dry matter growth

        : Parameters:
            - `sucrose` (:class:`float`) - amount of sucrose (µmol C)
            - `mstruct` (:class:`float`) - structural mass (g)
            - `DELTA_T` (:class:`float`) - value of the timestep (s)

        : Returns:
            mstruct_C_growth (µmol C), mstruct_growth (g), Nstruct_growth (g), Nstruct_N_growth (µmol N)

        :Returns Type:
            :class:`float`
        """
        ALPHA = 1
        VMAX_GROWTH = 0.015             #: Maximal rate of root structural dry matter growth (µmol C s-1 g-1 MS)
        K_GROWTH = 1250                 #: Affinity coefficient of root structural dry matter growth (µmol C g-1 MS)
        C_MOLAR_MASS = 12               #: Carbon molar mass (g mol-1)
        RATIO_C_MSTRUCT = 0.384         #: Mean contribution of carbon to structural dry mass (g C g-1 Mstruct)
        RATIO_N_MSTRUCT = 0.02          #: Mean contribution of nitrogen to structural dry mass (g N g-1 Mstruct)

        mstruct_C_growth = (((max(0, sucrose)/(mstruct*ALPHA)) * VMAX_GROWTH) / ((max(0, sucrose)/(mstruct*ALPHA)) + K_GROWTH)) * DELTA_T * mstruct     #: root growth in C (µmol of C)
        mstruct_growth = (mstruct_C_growth*1E-6 * C_MOLAR_MASS) / RATIO_C_MSTRUCT                                                                       #: root growth (g of structural dry mass)
        Nstruct_growth = mstruct_growth*RATIO_N_MSTRUCT                                                                                                 #: root growth in nitrogen (g)
        Nstruct_N_growth = min(amino_acids, (Nstruct_growth/cls.N_MOLAR_MASS)*1E6)                                                                      #: root growth in nitrogen (µmol N)

        return mstruct_C_growth, mstruct_growth, Nstruct_growth, Nstruct_N_growth

    @classmethod
    def calculate_roots_senescence(cls, mstruct, Nstruct, DELTA_T):
        """Root senescence

        : Parameters:
            - `mstruct` (:class:`float`) - structural mass (g)
            - `Nstruct` (:class:`float`) - structural N (g)
            - `DELTA_T` (:class:`float`) - value of the timestep (s)

        : Returns:
            Amount of C lost by root senescence over DELTA_T (g mstruct), amount of N lost by root senescence over DELTA_T (g Nstruct)

        :Returns Type:
            :class:`float`
        """
        return mstruct * cls.SENESCENCE_ROOTS * DELTA_T, Nstruct * cls.SENESCENCE_ROOTS * DELTA_T         #: Structral mass lost by turn_over (g)


    @classmethod
    def calculate_delta_mstruct_roots(cls, mstruct_growth, Nstruct_growth, mstruct_senescence, Nstruct_senescence, root_mstruct):
        """delta of root structural dry matter (g)

        : Parameters:
            - `mstruct_growth` (:class:`float`) - root structural mass growth (g)
            - `Nstruct_growth` (:class:`float`) - root structural N growth (g)
            - `mstruct_senescence` (:class:`float`) - amount of C lost by root senescence (g mstruct)
            - `Nstruct_senescence` (:class:`float`) - amount of N lost by root senescence (g Nstruct)
            - `root_mstruct` (:class:`float`) - actual mstruct of roots (g)

        : Returns:
            delta_mstruct (g), delta_mstruct (g), relative_delta_mstruct

        :Returns Type:
            :class:`float`
        """
        delta_mstruct, delta_Nstruct = (mstruct_growth - mstruct_senescence), (Nstruct_growth - Nstruct_senescence)
        relative_delta_mstruct = mstruct_senescence/root_mstruct
        return delta_mstruct, delta_Nstruct, relative_delta_mstruct