# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    cnwheat.postprocessings
    ~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`cnwheat.postprocessings` defines post-processings to apply 
    on CN-Wheat outputs, and provides a front-end to automatize the generation of graphs 
    for validation of the outputs.
    
    Please use preferably the front-ends :func:`postprocessings` and :func:`generate_graphs`.

    :copyright: Copyright 2014-2017 INRA-ECOSYS, see AUTHORS.
    :license: CeCILL-C, see LICENSE for details.
    
    **Acknowledgments**: The research leading these results has received funding through the 
    Investment for the Future programme managed by the Research National Agency 
    (BreedWheat project ANR-10-BTBR-03).
    
    .. seealso:: Barillot et al. 2016.
"""

"""
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
"""

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from cnwheat import simulation as cnwheat_simulation, model as cnwheat_model, parameters as cnwheat_parameters, tools as cnwheat_tools
from respiwheat import model as respiwheat_model

#: the time index
T_INDEX = cnwheat_simulation.Simulation.T_INDEX

#: the index to locate the plants in the modeled system
PLANTS_INDEXES = cnwheat_simulation.Simulation.PLANTS_INDEXES
#: concatenation of :attr:`T_INDEX` and :attr:`PLANTS_INDEXES`
PLANTS_T_INDEXES = cnwheat_simulation.Simulation.PLANTS_T_INDEXES
#: plants post-processing variables 
PLANTS_POSTPROCESSING_VARIABLES = []
#: concatenation of :attr:`PLANTS_T_INDEXES`, :attr:`PLANTS_RUN_VARIABLES <cnwheat_simulation.Simulation.PLANTS_RUN_VARIABLES>` and :attr:`PLANTS_POSTPROCESSING_VARIABLES`
PLANTS_RUN_POSTPROCESSING_VARIABLES = PLANTS_T_INDEXES + cnwheat_simulation.Simulation.PLANTS_RUN_VARIABLES + PLANTS_POSTPROCESSING_VARIABLES

#: the indexes to locate the axes in the modeled system
AXES_INDEXES = cnwheat_simulation.Simulation.AXES_INDEXES
#: concatenation of :attr:`T_INDEX` and :attr:`AXES_INDEXES`
AXES_T_INDEXES = cnwheat_simulation.Simulation.AXES_T_INDEXES
#: axes post-processing variables 
AXES_POSTPROCESSING_VARIABLES = []
#: concatenation of :attr:`AXES_T_INDEXES`, :attr:`AXES_RUN_VARIABLES <cnwheat_simulation.Simulation.AXES_RUN_VARIABLES>` and :attr:`AXES_POSTPROCESSING_VARIABLES`
AXES_RUN_POSTPROCESSING_VARIABLES = AXES_T_INDEXES + cnwheat_simulation.Simulation.AXES_RUN_VARIABLES + AXES_POSTPROCESSING_VARIABLES

#: the indexes to locate the phytomers in the modeled system
PHYTOMERS_INDEXES = cnwheat_simulation.Simulation.PHYTOMERS_INDEXES
#: concatenation of :attr:`T_INDEX` and :attr:`PHYTOMERS_INDEXES`
PHYTOMERS_T_INDEXES = cnwheat_simulation.Simulation.PHYTOMERS_T_INDEXES
#: phytomers post-processing variables 
PHYTOMERS_POSTPROCESSING_VARIABLES = []
#: concatenation of :attr:`PHYTOMERS_T_INDEXES`, :attr:`PHYTOMERS_RUN_VARIABLES <cnwheat_simulation.Simulation.PHYTOMERS_RUN_VARIABLES>` and :attr:`PHYTOMERS_POSTPROCESSING_VARIABLES`
PHYTOMERS_RUN_POSTPROCESSING_VARIABLES = PHYTOMERS_T_INDEXES + cnwheat_simulation.Simulation.PHYTOMERS_RUN_VARIABLES + PHYTOMERS_POSTPROCESSING_VARIABLES

#: the indexes to locate the organs in the modeled system
ORGANS_INDEXES = cnwheat_simulation.Simulation.ORGANS_INDEXES
#: concatenation of :attr:`T_INDEX` and :attr:`ORGANS_INDEXES`
ORGANS_T_INDEXES = cnwheat_simulation.Simulation.ORGANS_T_INDEXES
#: organs post-processing variables 
ORGANS_POSTPROCESSING_VARIABLES = ['Conc_Amino_Acids', 'Conc_Nitrates', 'Conc_Sucrose', 'Conc_cytokinins', 'Dry_Mass', 'Proteins_N_Mass', 'R_maintenance']
#: concatenation of :attr:`ORGANS_T_INDEXES`, :attr:`ORGANS_RUN_VARIABLES <cnwheat_simulation.Simulation.ORGANS_RUN_VARIABLES>` and :attr:`ORGANS_POSTPROCESSING_VARIABLES`
ORGANS_RUN_POSTPROCESSING_VARIABLES = ORGANS_T_INDEXES + cnwheat_simulation.Simulation.ORGANS_RUN_VARIABLES + ORGANS_POSTPROCESSING_VARIABLES

#: the indexes to locate the hidden zones in the modeled system
HIDDENZONE_INDEXES = cnwheat_simulation.Simulation.HIDDENZONE_INDEXES
#: concatenation of :attr:`T_INDEX` and :attr:`HIDDENZONE_INDEXES`
HIDDENZONE_T_INDEXES = cnwheat_simulation.Simulation.HIDDENZONE_T_INDEXES
#: hidden zones post-processing variables 
HIDDENZONE_POSTPROCESSING_VARIABLES = ['Conc_Amino_Acids', 'Conc_Fructan', 'Conc_Proteins', 'Conc_Sucrose']
#: concatenation of :attr:`HIDDENZONE_T_INDEXES`, :attr:`HIDDENZONE_RUN_VARIABLES <cnwheat_simulation.Simulation.HIDDENZONE_RUN_VARIABLES>` and :attr:`HIDDENZONE_POSTPROCESSING_VARIABLES`
HIDDENZONE_RUN_POSTPROCESSING_VARIABLES = HIDDENZONE_T_INDEXES + cnwheat_simulation.Simulation.HIDDENZONE_RUN_VARIABLES + HIDDENZONE_POSTPROCESSING_VARIABLES

#: the indexes to locate the elements in the modeled system
ELEMENTS_INDEXES = cnwheat_simulation.Simulation.ELEMENTS_INDEXES
#: concatenation of :attr:`T_INDEX` and :attr:`ELEMENTS_INDEXES`
ELEMENTS_T_INDEXES = cnwheat_simulation.Simulation.ELEMENTS_T_INDEXES
#: elements post-processing variables 
ELEMENTS_POSTPROCESSING_VARIABLES = ['Conc_Amino_Acids', 'Conc_Fructan', 'Conc_Nitrates', 'Conc_Proteins', 'Conc_Starch', 'Conc_Sucrose', 'Conc_TriosesP', 
                                     'Conc_cytokinins', 'R_maintenance']
#: concatenation of :attr:`ELEMENTS_T_INDEXES`, :attr:`ELEMENTS_RUN_VARIABLES <cnwheat_simulation.Simulation.ELEMENTS_RUN_VARIABLES>` and :attr:`ELEMENTS_POSTPROCESSING_VARIABLES`
ELEMENTS_RUN_POSTPROCESSING_VARIABLES = ELEMENTS_T_INDEXES + cnwheat_simulation.Simulation.ELEMENTS_RUN_VARIABLES + ELEMENTS_POSTPROCESSING_VARIABLES

#: the indexes to locate the soils in the modeled system
SOILS_INDEXES = cnwheat_simulation.Simulation.SOILS_INDEXES
#: concatenation of :attr:`T_INDEX` and :attr:`SOILS_INDEXES`
SOILS_T_INDEXES = cnwheat_simulation.Simulation.SOILS_T_INDEXES
#: soils post-processing variables 
SOILS_POSTPROCESSING_VARIABLES = []
#: concatenation of :attr:`SOILS_T_INDEXES`, :attr:`SOILS_RUN_VARIABLES <cnwheat_simulation.Simulation.SOILS_RUN_VARIABLES>` and :attr:`SOILS_POSTPROCESSING_VARIABLES`
SOILS_RUN_POSTPROCESSING_VARIABLES = SOILS_T_INDEXES + cnwheat_simulation.Simulation.SOILS_RUN_VARIABLES + SOILS_POSTPROCESSING_VARIABLES


###################################################
############ POST-PROCESSINGS FUNCTIONS ###########
############# DO NOT USE THEM DIRECTLY ############
###################################################  

class Roots:
    """
    Post-processings to apply on Roots outputs.
    """
    @staticmethod
    def calculate_conc_nitrates(nitrates, mstruct):
        """Nitrate concentration.

        :Parameters:
            - `nitrates` (:class:`float`) - Amount of nitrates (µmol N)
        :Returns:
            Nitrate concentration (µmol nitrates g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return (nitrates/mstruct)
    
    @staticmethod
    def calculate_conc_amino_acids(amino_acids, mstruct):
        """Amino_acid concentration.

        :Parameters:
            - `amino_acids` (:class:`float`) - Amount of amino acids (µmol N)
        :Returns:
            Amino_acid concentration (µmol amino_acids g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return (amino_acids/cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_N_RATIO)/mstruct
    
    @staticmethod
    def calculate_conc_sucrose(sucrose, mstruct):
        """Sucrose concentration.

        :Parameters:
            - `sucrose` (:class:`float`) - Amount of sucrose (µmol C)
        :Returns:
            Sucrose concentration (µmol sucrose g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return (sucrose/mstruct)/cnwheat_model.EcophysiologicalConstants.NB_C_SUCROSE
    
    @staticmethod
    def calculate_conc_cytokinins(cytokinins, mstruct):
        """Cytokinin concentration.

        :Parameters:
            - `cytokinins` (:class:`float`) - Amount of cytokinins (AU)
        :Returns:
            cytokinins concentration (AU cytokinins g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return cytokinins/mstruct
    
    
class Phloem:
    """
    Post-processings to apply on Phloem outputs.
    """
    @staticmethod
    def calculate_conc_amino_acids(amino_acids, mstruct_axis):
        """Amino_acids concentration. Related to the structural dry mass of the culm.

        :Parameters:
            - `amino_acids` (:class:`float`) - Amount of amino_acids in phloem (µmol N)
            - `mstruct_axis` (:class:`float`) -The structural dry mass of the axis (g)
        :Returns:
            Amino_acids concentration (µmol amino_acids g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return (amino_acids / cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_N_RATIO) / mstruct_axis
    
    @staticmethod
    def calculate_conc_sucrose(sucrose, mstruct_axis):
        """Sucrose concentration. Related to the structural dry mass of the culm

        :Parameters:
            - `sucrose` (:class:`float`) - Amount of sucrose in phloem (µmol C)
            - `mstruct_axis` (:class:`float`) -The structural dry mass of the axis (g)
        :Returns:
            Sucrose concentration (µmol sucrose g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return (sucrose / cnwheat_model.EcophysiologicalConstants.NB_C_SUCROSE) / mstruct_axis


class Grains:
    """
    Post-processings to apply on Grains outputs.
    """
    @staticmethod
    def calculate_dry_mass(structure, starch, proteins):
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
        C_mass = ((structure + starch)*1E-6*cnwheat_model.EcophysiologicalConstants.C_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.HEXOSE_MOLAR_MASS_C_RATIO

        #: N mass, grain proteins were supposed to be gluten mainly composed of Glu, Gln and Pro
        N_mass = (proteins*1E-6*cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS) / cnwheat_model.Grains.AMINO_ACIDS_MOLAR_MASS_N_RATIO

        return  C_mass + N_mass
    
    @staticmethod
    def calculate_protein_N_mass(proteins):
        """Grain total protein mass.

        :Parameters:
            - `proteins` (:class:`float`) - Grain protein content (µmol N)
        :Returns:
            Grain total protein mass (g)
        :Returns Type:
            :class:`float`
        """
        mass_N_proteins = proteins*1E-6 * cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS                        #: Mass of nitrogen in proteins (g)
        #masS_proteins = mass_N_proteins / EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_N_RATIO     #: Total mass of proteins (g)
        return mass_N_proteins
    
    
class HiddenZone:
    """
    Post-processings to apply on HiddenZone outputs.
    """
    @staticmethod
    def calculate_conc_amino_acids(amino_acids, mstruct):
        """Amino acid concentration.

        :Parameters:
            - `amino_acids` (:class:`float`) - N amino acids (µmol N)
        :Returns:
            Amino_acid concentration (µmol amino acids g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return (amino_acids/cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_N_RATIO) / mstruct
    
    @staticmethod
    def calculate_conc_sucrose(sucrose, mstruct):
        """Sucrose concentration.

        :Parameters:
            - `sucrose` (:class:`float`) - C sucrose (µmol C)
        :Returns:
            Sucrose concentration (µmol sucrose g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return (sucrose/mstruct) / cnwheat_model.EcophysiologicalConstants.NB_C_SUCROSE
    
    @staticmethod
    def calculate_conc_fructan(fructan, mstruct):
        """Fructan concentration.

        :Parameters:
            - `fructan` (:class:`float`) - C fructan (µmol C)
        :Returns:
            Fructan concentration (µmol fructan g-1 mstruct, eq. glucose).
        :Returns Type:
            :class:`float`
        """
        return (fructan/mstruct)/cnwheat_model.EcophysiologicalConstants.NB_C_HEXOSES

    @staticmethod
    def calculate_conc_protein(proteins, mstruct):
        """Proteins concentration.

        :Parameters:
            - `proteins` (:class:`float`) - N proteins (µmol N)
            - `mstruct` (:class:`float`) - Structural mass (g)
        :Returns:
            Protein concentration (g proteins g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        mass_N_proteins = proteins*1E-6 * cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS                        #: Mass of N in proteins (g)
        masS_proteins = mass_N_proteins / cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_N_RATIO      #: Total mass of proteins (g)
        return (masS_proteins / mstruct)


class Element:
    """
    Post-processings to apply on Element outputs.
    """
    @staticmethod
    def calculate_conc_triosesP(triosesP, mstruct):
        """Triose Phosphates concentration.

        :Parameters:
            - `triosesP` (:class:`float`) - Amount of triose phosphates (µmol C)
        :Returns:
            Triose phosphates concentration (µmol triosesP g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return (triosesP/mstruct)/cnwheat_model.EcophysiologicalConstants.NB_C_TRIOSEP

    @staticmethod
    def calculate_conc_sucrose(sucrose, mstruct):
        """Sucrose concentration.

        :Parameters:
            - `sucrose` (:class:`float`) - Amount of sucrose (µmol C)
        :Returns:
            Sucrose concentration (µmol sucrose g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return (sucrose/mstruct)/cnwheat_model.EcophysiologicalConstants.NB_C_SUCROSE

    @staticmethod
    def calculate_conc_starch(starch, mstruct):
        """Starch concentration.

        :Parameters:
            - `starch` (:class:`float`) - Amount of sucrose (µmol C)
        :Returns:
            Starch concentration (µmol starch g-1 mstruct, eq. glucose).
        :Returns Type:
            :class:`float`
        """
        return (starch/mstruct)/cnwheat_model.EcophysiologicalConstants.NB_C_HEXOSES

    @staticmethod
    def calculate_conc_fructan(fructan, mstruct):
        """Fructan concentration.

        :Parameters:
            - `fructan` (:class:`float`) - Amount of fructan (µmol C)
        :Returns:
            Fructan concentration (µmol fructan g-1 mstruct, eq. glucose).
        :Returns Type:
            :class:`float`
        """
        return (fructan/mstruct)/cnwheat_model.EcophysiologicalConstants.NB_C_HEXOSES
    
    @staticmethod
    def calculate_conc_nitrates(nitrates, mstruct):
        """Nitrate concentration.

        :Parameters:
            - `nitrates` (:class:`float`) - Amount of nitrates (µmol N)
        :Returns:
            Nitrate concentration (µmol nitrates g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return (nitrates/mstruct)
    
    @staticmethod
    def calculate_conc_amino_acids(amino_acids, mstruct):
        """Amino_acid concentration.

        :Parameters:
            - `amino_acids` (:class:`float`) - Amount of amino acids (µmol N)
        :Returns:
            Amino_acid concentration (µmol amino acids g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return (amino_acids/cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_N_RATIO) / mstruct
    
    @staticmethod
    def calculate_conc_proteins(proteins, mstruct):
        """Protein concentration.

        :Parameters:
            - `proteins` (:class:`float`) - Amount of proteins (µmol N)
        :Returns:
            Protein concentration (g proteins g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        mass_N_proteins = proteins*1E-6 * cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS                        #: Mass of N in proteins (g)
        masS_proteins = mass_N_proteins / cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_N_RATIO      #: Total mass of proteins (g)
        return (masS_proteins / mstruct)
    
    @staticmethod
    def calculate_conc_cytokinins(cytokinins, mstruct):
        """Cytokinin concentration.

        :Parameters:
            - `cytokinins` (:class:`float`) - Amount of cytokinins (AU)
        :Returns:
            Cytokinin concentration (AU g-1 mstruct)
        :Returns Type:
            :class:`float`
        """
        return cytokinins/mstruct


###################################################
############ POST-PROCESSINGS FRONT-END ###########
###################################################

def postprocessings(plants_df=None, axes_df=None, metamers_df=None, hiddenzones_df=None, 
                    organs_df=None, elements_df=None, soils_df=None, delta_t=1):
    """
    Compute post-processings from CN-Wheat outputs, and format the post-processings to :class:`dataframes <pandas.DataFrame>`.
    
    For each post-processing output dataframe:
    
        * compute post-processings from CN-Wheat outputs, 
        * concatenate CN-Wheat outputs and post-processings and place the results in a jointed dataframe,
        * reorder the columns of the dataframes according to :attr:`PLANTS_RUN_POSTPROCESSING_VARIABLES`, 
          :attr:`AXES_RUN_POSTPROCESSING_VARIABLES`, :attr:`PHYTOMERS_RUN_POSTPROCESSING_VARIABLES`, 
          :attr:`ORGANS_RUN_POSTPROCESSING_VARIABLES`, :attr:`HIDDENZONE_RUN_POSTPROCESSING_VARIABLES`,
          :attr:`ELEMENTS_RUN_POSTPROCESSING_VARIABLES` and :attr:`SOILS_RUN_POSTPROCESSING_VARIABLES`, 
        * and convert the indexes of plants and metamers to integers (if relevant).

    :Parameters:
            - `plants_df` (:class:`pandas.DataFrame`) - CN-Wheat outputs at plant scale (see :attr:`simulation.Simulation.PLANTS_RUN_VARIABLES`)
            - `axes_df` (:class:`pandas.DataFrame`) - CN-Wheat outputs at axis scale (see :attr:`simulation.Simulation.AXES_RUN_VARIABLES`)
            - `metamers_df` (:class:`pandas.DataFrame`) - CN-Wheat outputs at phytomer scale (see :attr:`simulation.Simulation.PHYTOMERS_RUN_VARIABLES`)
            - `hiddenzones_df` (:class:`pandas.DataFrame`) - CN-Wheat outputs at hidden zone scale (see :attr:`simulation.Simulation.HIDDENZONE_RUN_VARIABLES`)
            - `organs_df` (:class:`pandas.DataFrame`) - CN-Wheat outputs at organ scale (see :attr:`simulation.Simulation.ORGANS_RUN_VARIABLES`)
            - `elements_df` (:class:`pandas.DataFrame`) - CN-Wheat outputs at element scale (see :attr:`simulation.Simulation.ELEMENTS_RUN_VARIABLES`)
            - `soils_df` (:class:`pandas.DataFrame`) - CN-Wheat outputs at soil scale (see :attr:`simulation.Simulation.SOILS_RUN_VARIABLES`)
            - `delta_t` (:class:`pandas.DataFrame`) - the delta t between 2 outputs (in seconds).

    :Returns:
        :class:`dataframes <pandas.DataFrame>` of post-processings for each scale:

            * plant (see :attr:`PLANTS_RUN_POSTPROCESSING_VARIABLES`)
            * axis (see :attr:`AXES_RUN_POSTPROCESSING_VARIABLES`)
            * metamer (see :attr:`PHYTOMERS_RUN_POSTPROCESSING_VARIABLES`)
            * organ (see :attr:`ORGANS_RUN_POSTPROCESSING_VARIABLES`)
            * hidden zone (see :attr:`HIDDENZONE_RUN_POSTPROCESSING_VARIABLES`)
            * element (see :attr:`ELEMENTS_RUN_POSTPROCESSING_VARIABLES`)
            * and soil (see :attr:`SOILS_RUN_POSTPROCESSING_VARIABLES`)
        
        depending of the dataframes given as argument.
        For example, if user passes only dataframes `plants_df`, `axes_df` and `metamers_df`, 
        then only post-processings dataframes of plants, axes and metamers are returned.  
            
    :Returns Type:
        :class:`tuple` of :class:`pandas.DataFrame`
        
    """
    
    returned_dataframes = []
    
    # plants
    if plants_df is not None:
        pp_plants_df = pd.concat([plants_df, pd.DataFrame(columns=PLANTS_POSTPROCESSING_VARIABLES)])
        pp_plants_df = pp_plants_df.reindex_axis(PLANTS_RUN_POSTPROCESSING_VARIABLES, axis=1, copy=False)
        pp_plants_df['plant'] = pp_plants_df['plant'].astype(int)
        returned_dataframes.append(pp_plants_df)
                
    # axes
    if axes_df is not None:
        pp_axes_df = pd.concat([axes_df, pd.DataFrame(columns=AXES_POSTPROCESSING_VARIABLES)])
        pp_axes_df = pp_axes_df.reindex_axis(AXES_RUN_POSTPROCESSING_VARIABLES, axis=1, copy=False)
        pp_axes_df['plant'] = pp_axes_df['plant'].astype(int)
        returned_dataframes.append(pp_axes_df)
    
    # metamers
    if metamers_df is not None:
        pp_metamers_df = pd.concat([metamers_df, pd.DataFrame(columns=PHYTOMERS_POSTPROCESSING_VARIABLES)])
        pp_metamers_df = pp_metamers_df.reindex_axis(PHYTOMERS_RUN_POSTPROCESSING_VARIABLES, axis=1, copy=False)
        pp_metamers_df[['plant', 'metamer']] = pp_metamers_df[['plant', 'metamer']].astype(int)
        returned_dataframes.append(pp_metamers_df)
        
    # hidden zones
    if hiddenzones_df is not None:
        pp_hiddenzones_df = pd.concat([hiddenzones_df, pd.DataFrame(columns=HIDDENZONE_POSTPROCESSING_VARIABLES)])
        pp_hiddenzones_df.loc[:, 'Conc_Amino_Acids'] = HiddenZone.calculate_conc_amino_acids(hiddenzones_df['amino_acids'], hiddenzones_df['mstruct'])
        pp_hiddenzones_df.loc[:, 'Conc_Fructan'] = HiddenZone.calculate_conc_fructan(hiddenzones_df['fructan'], hiddenzones_df['mstruct'])
        pp_hiddenzones_df.loc[:, 'Conc_Proteins'] = HiddenZone.calculate_conc_protein(hiddenzones_df['proteins'], hiddenzones_df['mstruct'])
        pp_hiddenzones_df.loc[:, 'Conc_Sucrose'] = HiddenZone.calculate_conc_sucrose(hiddenzones_df['sucrose'], hiddenzones_df['mstruct'])
        pp_hiddenzones_df = pp_hiddenzones_df.reindex_axis(HIDDENZONE_RUN_POSTPROCESSING_VARIABLES, axis=1, copy=False)
        pp_hiddenzones_df[['plant', 'metamer']] = pp_hiddenzones_df[['plant', 'metamer']].astype(int)
        returned_dataframes.append(pp_hiddenzones_df)
    
    # organs
    if organs_df is not None:
        pp_organs_df = pd.concat([organs_df, pd.DataFrame(columns=ORGANS_POSTPROCESSING_VARIABLES)])
        # roots
        roots_df = organs_df.loc[organs_df.organ == 'roots']
        pp_organs_df.loc[pp_organs_df.organ == 'roots', 'Conc_Nitrates'] = Roots.calculate_conc_nitrates(roots_df['nitrates'], roots_df['mstruct'])
        pp_organs_df.loc[pp_organs_df.organ == 'roots', 'Conc_Amino_Acids'] = Roots.calculate_conc_amino_acids(roots_df['amino_acids'], roots_df['mstruct'])
        pp_organs_df.loc[pp_organs_df.organ == 'roots', 'Conc_Sucrose'] = Roots.calculate_conc_sucrose(roots_df['sucrose'], roots_df['mstruct'])
        pp_organs_df.loc[pp_organs_df.organ == 'roots', 'Conc_cytokinins'] = Roots.calculate_conc_cytokinins(roots_df['cytokinins'], roots_df['mstruct'])
        pp_organs_df.loc[pp_organs_df.organ == 'roots', 'Conc_cytokinins'] = Roots.calculate_conc_cytokinins(roots_df['cytokinins'], roots_df['mstruct'])
        R_residual = np.array(map(respiwheat_model.RespirationModel.R_residual, roots_df['sucrose'], roots_df['mstruct']*cnwheat_model.Roots.PARAMETERS.ALPHA, roots_df['Total_Organic_Nitrogen'], [delta_t]*len(roots_df), soils_df['Tsoil']))
        pp_organs_df.loc[pp_organs_df.organ == 'roots', 'R_maintenance'] = R_residual[:,1]    
        # phloem
        phloems_df = organs_df.loc[organs_df.organ == 'phloem']
        pp_organs_df.loc[pp_organs_df.organ == 'phloem', 'Conc_Amino_Acids'] = Phloem.calculate_conc_amino_acids(phloems_df['amino_acids'], axes_df['mstruct'])
        pp_organs_df.loc[pp_organs_df.organ == 'phloem', 'Conc_Sucrose'] = Phloem.calculate_conc_sucrose(phloems_df['sucrose'], axes_df['mstruct'])
        # grains
        grains_df = organs_df.loc[organs_df.organ == 'grain']
        pp_organs_df.loc[pp_organs_df.organ == 'grain', 'Dry_Mass'] = Grains.calculate_dry_mass(grains_df['structure'], grains_df['starch'], grains_df['proteins'])
        pp_organs_df.loc[pp_organs_df.organ == 'grain', 'Proteins_N_Mass'] = Grains.calculate_protein_N_mass(grains_df['proteins'])
        pp_organs_df = pp_organs_df.reindex_axis(ORGANS_RUN_POSTPROCESSING_VARIABLES, axis=1, copy=False)
        pp_organs_df['plant'] = pp_organs_df['plant'].astype(int)
        returned_dataframes.append(pp_organs_df)
    
    # elements
    if elements_df is not None:
        pp_elements_df = pd.concat([elements_df, pd.DataFrame(columns=ELEMENTS_POSTPROCESSING_VARIABLES)])
        pp_elements_df.loc[:, 'Conc_TriosesP'] = Element.calculate_conc_triosesP(elements_df['triosesP'], elements_df['mstruct'])
        pp_elements_df.loc[:, 'Conc_Starch'] = Element.calculate_conc_starch(elements_df['starch'], elements_df['mstruct'])
        pp_elements_df.loc[:, 'Conc_Sucrose'] = Element.calculate_conc_sucrose(elements_df['sucrose'], elements_df['mstruct'])
        pp_elements_df.loc[:, 'Conc_Fructan'] = Element.calculate_conc_fructan(elements_df['fructan'], elements_df['mstruct'])
        pp_elements_df.loc[:, 'Conc_Nitrates'] = Element.calculate_conc_nitrates(elements_df['nitrates'], elements_df['mstruct'])
        pp_elements_df.loc[:, 'Conc_Amino_Acids'] = Element.calculate_conc_amino_acids(elements_df['amino_acids'], elements_df['mstruct'])
        pp_elements_df.loc[:, 'Conc_Proteins'] = Element.calculate_conc_proteins(elements_df['proteins'], elements_df['mstruct'])
        pp_elements_df.loc[:, 'Conc_cytokinins'] = Element.calculate_conc_cytokinins(elements_df['cytokinins'], elements_df['mstruct'])
        grouped = elements_df.groupby('organ')
        for organ_type, parameters_class in \
            (('ear', cnwheat_parameters.CHAFF_ELEMENT_PARAMETERS),
             ('blade', cnwheat_parameters.LAMINA_ELEMENT_PARAMETERS),
             ('internode', cnwheat_parameters.INTERNODE_ELEMENT_PARAMETERS),
             ('peduncle', cnwheat_parameters.PEDUNCLE_ELEMENT_PARAMETERS),
             ('sheath', cnwheat_parameters.SHEATH_ELEMENT_PARAMETERS)):
            if organ_type not in grouped.groups:
                continue  
            group = grouped.get_group(organ_type)
            if len(group) == 0:
                continue
            curr_organ_elements_df = elements_df.loc[group.index]
            pp_curr_organ_elements_df = pp_elements_df.loc[group.index]
            R_residual = np.array(map(respiwheat_model.RespirationModel.R_residual, curr_organ_elements_df['sucrose'], curr_organ_elements_df['mstruct'] * parameters_class.ALPHA, curr_organ_elements_df['Total_Organic_Nitrogen'], [delta_t] * len(curr_organ_elements_df), curr_organ_elements_df['Ts']))
            pp_curr_organ_elements_df.loc[:, 'R_maintenance'] = R_residual[:,1]
        pp_elements_df = pp_elements_df.reindex_axis(ELEMENTS_RUN_POSTPROCESSING_VARIABLES, axis=1, copy=False)
        pp_elements_df[['plant', 'metamer']] = pp_elements_df[['plant', 'metamer']].astype(int)
        returned_dataframes.append(pp_elements_df)
    
    # soils
    if soils_df is not None:
        pp_soils_df = pd.concat([soils_df, pd.DataFrame(columns=SOILS_POSTPROCESSING_VARIABLES)])
        pp_soils_df = pp_soils_df.reindex_axis(SOILS_RUN_POSTPROCESSING_VARIABLES, axis=1, copy=False)
        pp_soils_df[['plant']] = pp_soils_df[['plant']].astype(int)
        returned_dataframes.append(pp_soils_df)
    
    return tuple(returned_dataframes)


###################################################
############ GRAPHS GENERATION FRONT-END ##########
###################################################

def generate_graphs(axes_df=None, hiddenzones_df=None, organs_df=None, elements_df=None, soils_df=None, graphs_dirpath='.'):
    """
    Generate graphs to validate the outputs of CN-Wheat, and save them in directory `graphs_dirpath`.
    
    :Parameters:
        - `axes_df` (:class:`pandas.DataFrame`) - CN-Wheat outputs and post-processings at axis scale (see :attr:`PLANTS_RUN_POSTPROCESSING_VARIABLES`)
        - `hiddenzones_df` (:class:`pandas.DataFrame`) - CN-Wheat outputs at hidden zone scale (see :attr:`HIDDENZONE_RUN_POSTPROCESSING_VARIABLES`)
        - `organs_df` (:class:`pandas.DataFrame`) - CN-Wheat outputs at organ scale (see :attr:`ORGANS_RUN_POSTPROCESSING_VARIABLES`)
        - `elements_df` (:class:`pandas.DataFrame`) - CN-Wheat outputs at element scale (see :attr:`ELEMENTS_RUN_POSTPROCESSING_VARIABLES`)
        - `soils_df` (:class:`pandas.DataFrame`) - CN-Wheat outputs at soil scale (see :attr:`SOILS_RUN_POSTPROCESSING_VARIABLES`)
        - `graphs_dirpath` (:class:`pandas.DataFrame`) - the path of the directory to save the generated graphs in
    
    """
    
    x_name = 't'
    x_label='Time (Hour)'
    
    # 1) Photosynthetic organs
    # Obsolete variables? 
    # 'PARa': u'Absorbed PAR (µmol m$^{-2}$ s$^{-1}$)',
    # 'An': u'Net photosynthesis (µmol m$^{-2}$ s$^{-1}$)', 
    # 'Rd': u'Mitochondrial respiration rate of organ in light (µmol C h$^{-1}$)', 
    # 'gs': u'Conductance stomatique (mol m$^{-2}$ s$^{-1}$)',
    graph_variables_ph_elements = {'Ag': u'Gross photosynthesis (µmol m$^{-2}$ s$^{-1}$)', 'Tr':u'Organ surfacic transpiration rate (mmol H$_{2}$0 m$^{-2}$ s$^{-1}$)', 'Transpiration':u'Organ transpiration rate (mmol H$_{2}$0 s$^{-1}$)', 'Ts': u'Temperature surface (°C)', 
                                   'Conc_TriosesP': u'[TriosesP] (µmol g$^{-1}$ mstruct)', 'Conc_Starch':u'[Starch] (µmol g$^{-1}$ mstruct)', 'Conc_Sucrose':u'[Sucrose] (µmol g$^{-1}$ mstruct)', 'Conc_Fructan':u'[Fructan] (µmol g$^{-1}$ mstruct)',
                                   'Conc_Nitrates': u'[Nitrates] (µmol g$^{-1}$ mstruct)', 'Conc_Amino_Acids': u'[Amino_Acids] (µmol g$^{-1}$ mstruct)', 'Conc_Proteins': u'[Proteins] (g g$^{-1}$ mstruct)',
                                   'Nitrates_import': u'Total nitrates imported (µmol h$^{-1}$)', 'Amino_Acids_import': u'Total amino acids imported (µmol N h$^{-1}$)',
                                   'S_Amino_Acids': u'[Rate of amino acids synthesis] (µmol N g$^{-1}$ mstruct h$^{-1}$)', 'S_Proteins': u'Rate of protein synthesis (µmol N g$^{-1}$ mstruct h$^{-1}$)', 'D_Proteins': u'Rate of protein degradation (µmol N g$^{-1}$ mstruct h$^{-1}$)', 'k_proteins': u'Relative rate of protein degradation (s$^{-1}$)',
                                   'Loading_Sucrose': u'Loading Sucrose (µmol C sucrose h$^{-1}$)', 'Loading_Amino_Acids': u'Loading Amino acids (µmol N amino acids h$^{-1}$)',
                                   'green_area': u'Green area (m$^{2}$)', 'R_phloem_loading': u'Respiration phloem loading (µmol C h$^{-1}$)', 'R_Nnit_red': u'Respiration nitrate reduction (µmol C h$^{-1}$)', 'R_residual': u'Respiration residual (µmol C h$^{-1}$)', 'R_maintenance': u'Respiration residual (µmol C h$^{-1}$)',
                                   'mstruct': u'Structural mass (g)', 'Nstruct': u'Structural N mass (g)',
                                   'Conc_cytokinins':u'[cytokinins] (UA g$^{-1}$ mstruct)', 'D_cytokinins':u'Cytokinin degradation (UA g$^{-1}$ mstruct)', 'cytokinins_import':u'Cytokinin import (UA)'}
    
    for org_ph in (['blade'], ['sheath'], ['internode'], ['peduncle', 'ear']):
        for variable_name, variable_label in graph_variables_ph_elements.iteritems():
            graph_name = variable_name + '_' + '_'.join(org_ph) + '.PNG'
            cnwheat_tools.plot_cnwheat_ouputs(elements_df,
                          x_name = x_name,
                          y_name = variable_name,
                          x_label=x_label,
                          y_label=variable_label,
                          filters={'organ': org_ph},
                          plot_filepath=os.path.join(graphs_dirpath, graph_name),
                          explicit_label=False)
    
    # 2) Roots, grains and phloem
    # Obsolete variables? 
    # 'R_growth': u'Growth respiration of roots (µmol C h$^{-1}$)', 
    graph_variables_organs = {'Conc_Sucrose':u'[Sucrose] (µmol g$^{-1}$ mstruct)', 'Dry_Mass':'Dry mass (g)',
                        'Conc_Nitrates': u'[Nitrates] (µmol g$^{-1}$ mstruct)', 'Conc_Amino_Acids':u'[Amino Acids] (µmol g$^{-1}$ mstruct)', 'Proteins_N_Mass': u'[N Proteins] (g)',
                        'Uptake_Nitrates':u'Nitrates uptake (µmol h$^{-1}$)', 'Unloading_Sucrose':u'Unloaded sucrose (µmol C g$^{-1}$ mstruct h$^{-1}$)', 'Unloading_Amino_Acids':u'Unloaded Amino Acids (µmol N AA g$^{-1}$ mstruct h$^{-1}$)',
                        'S_Amino_Acids': u'Rate of amino acids synthesis (µmol N g$^{-1}$ mstruct h$^{-1}$)', 'S_Proteins': u'Rate of protein synthesis (µmol N h$^{-1}$)', 'Export_Nitrates': u'Total export of nitrates (µmol N h$^{-1}$)', 'Export_Amino_Acids': u'Total export of Amino acids (µmol N h$^{-1}$)',
                        'R_Nnit_upt': u'Respiration nitrates uptake (µmol C h$^{-1}$)', 'R_Nnit_red': u'Respiration nitrate reduction (µmol C h$^{-1}$)', 'R_residual': u'Respiration residual (µmol C h$^{-1}$)', 'R_maintenance': u'Respiration residual (µmol C h$^{-1}$)',
                        'R_grain_growth_struct': u'Respiration grain structural growth (µmol C h$^{-1}$)', 'R_grain_growth_starch': u'Respiration grain starch growth (µmol C h$^{-1}$)',
                        'mstruct': u'Structural mass (g)',
                        'C_exudation': u'Carbon lost by root exudation (µmol C g$^{-1}$ mstruct h$^{-1}$', 'N_exudation': u'Nitrogen lost by root exudation (µmol N g$^{-1}$ mstruct h$^{-1}$',
                        'Conc_cytokinins':u'[cytokinins] (UA g$^{-1}$ mstruct)', 'S_cytokinins':u'Rate of cytokinins synthesis (UA g$^{-1}$ mstruct)', 'Export_cytokinins': 'Export of cytokinins from roots (UA h$^{-1}$)',
                        'HATS_LATS': u'Potential uptake (µmol h$^{-1}$)' , 'regul_transpiration':'Regulating transpiration function'}
    
    for org in (['roots'], ['grains'], ['phloem']):
        for variable_name, variable_label in graph_variables_organs.iteritems():
            graph_name = variable_name + '_' + '_'.join(org) + '.PNG'
            cnwheat_tools.plot_cnwheat_ouputs(organs_df,
                          x_name = x_name,
                          y_name = variable_name,
                          x_label=x_label,
                          y_label=variable_label,
                          filters={'organ': org},
                          plot_filepath=os.path.join(graphs_dirpath, graph_name),
                          explicit_label=False)
    
    # 3) Soil
    _, (ax1) = plt.subplots(1)
    conc_nitrates_soil = soils_df['Conc_Nitrates_Soil']*14E-6
    ax1.plot(soils_df['t'], conc_nitrates_soil)
    ax1.set_ylabel(u'[Nitrates] (g m$^{-3}$)')
    ax1.set_xlabel('Time from flowering (hour)')
    ax1.set_title = 'Conc Nitrates Soil'
    plt.savefig(os.path.join(graphs_dirpath, 'Conc_Nitrates_Soil.PNG'), format='PNG', bbox_inches='tight')
    plt.close()
    
    # 4) Hidden zones
    # Obsolete variables? 
    # 'Respi_growth': u'Growth respiration (µmol C)',
    # 'delta_leaf_L':u'Delta leaf length (m)',
    # 'leaf_L': u'Leaf length (m)',
    # 'leaf_dist_to_emerge': u'Length for leaf emergence (m)',
    # 'sucrose_consumption_mstruct': u'Consumption of sucrose for growth (µmol C)' 
    graph_variables_hiddenzones = {'Conc_Sucrose':u'[Sucrose] (µmol g$^{-1}$ mstruct)', 'Conc_Amino_Acids':u'[Amino Acids] (µmol g$^{-1}$ mstruct)', 'Conc_Proteins': u'[Proteins] (g g$^{-1}$ mstruct)', 'Conc_Fructan':u'[Fructan] (µmol g$^{-1}$ mstruct)',
                                    'Unloading_Sucrose':u'Sucrose unloading (µmol C)', 'Unloading_Amino_Acids':u'Amino_acids unloading (µmol N)', 'mstruct': u'Structural mass (g)'}
    
    for variable_name, variable_label in graph_variables_hiddenzones.iteritems():
        graph_name = variable_name + '_hz' + '.PNG'
        cnwheat_tools.plot_cnwheat_ouputs(hiddenzones_df,
                      x_name = x_name,
                      y_name = variable_name,
                      x_label = x_label,
                      y_label = variable_label,
                      filters={'plant': 1, 'axis': 'MS'},
                      plot_filepath=os.path.join(graphs_dirpath, graph_name),
                      explicit_label=False)
    
    # 5) Organs
    # Obsolete variables? 
    # 'visible_length': u'Length (m)'
    graph_variables_organs = {}
    for variable_name, variable_label in graph_variables_organs.iteritems():
        graph_name = variable_name + '.PNG'
        cnwheat_tools.plot_cnwheat_ouputs(organs_df,
                      x_name = x_name,
                      y_name = variable_name,
                      x_label = x_label,
                      y_label = variable_label,
                      filters={'plant': 1, 'axis': 'MS'},
                      plot_filepath=os.path.join(graphs_dirpath, graph_name),
                      explicit_label=False)
    
    
    

