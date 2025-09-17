
.. _farquharwheat_user:

Farquhar-Wheat User Guide
#########################

.. contents::

Introduction
============

This is the documentation for Farquhar-Wheat, a model of leaf photosynthesis based on Farquhar's approach.
The models computes C assimilation, stomatal conductance and organ temperature
for individual shoot organs. The model also provides an option to calculate the photosynthesis
at an infra-leaf scale (scene primitive).

- Calculation of C assimilation was adapted from the Farquhar model proposed by Müller et al. (2005) \
  and Braune et al. (2009) for wheat. The model was adapted to account
  for the metabolic contents provided by CN-Wheat model :
  (i) the actual surfacic N amount (stimulation of the assimilation rates) and
  (ii) the actual non-structural carbohydrate content (TPU limitation).

- Calculation of the stomatal conductance is obtained from the BWB model

- Leaf erergy balance is computed from the weather data inputs, C assimilation and
  stomatal conductance in order to simulate organ temperature.

The 3 processes described above are therefore coupled and iteratively computed until
internal CO2 (Ci) and leaf temperature converge at each time step.


Inputs of Farquhar-Wheat
========================

- ambient_CO2: atmospheric CO2 concentration (ppm)
- RH : relative air humidity (decimal fraction)
- Ur : wind speed at 2m height above canopy (m s-1)
- height_canopy (m)
- height : organ height (m)
- width : lamina width or sheath diamter (m)
- PARa : PAR absorbed (µmol m-2 s-1)
- nitrates : amount of nitrates (µmol N)
- amino_acids : amount of amino acids (µmol N)
- proteins : amount of proteins (µmol N)
- Nstruct : structural N (g)
- sucrose : amount of sucrose (µmol C)
- starch : amount of starch (µmol C)
- fructan : amount of fructan (µmol C)
- PARa_prim: PAR absorbed by the primitive (µmol m-2 s-1)
- area_prim: area of the primitive (m²)
- organ_label

Outputs of Farquhar-Wheat
=========================
- Ag: gross C assimilation (µmol m-2 s-1)
- An : Net C assimilation (µmol m-2 s-1)
- Rd : respiration (µmol m-2 s-1)
- Tr : transpiration (mmol m-2 s-1)
- Ts : organ surface temperature (°C)
- gs : stomatal conductance for water (mol m-2 s-1)


Package architecture
=====================

Farquhar-Wheat is a Python package which consists of several Python modules:

* :mod:`openalea.farquharwheat.model`: the state and the equations of the model,
* :mod:`openalea.farquharwheat.parameters`: the parameters of the model,
* :mod:`openalea.farquharwheat.simulation`: the simulator (front-end) to run the model,
* and :mod:`openalea.farquharwheat.converter`: functions to convert Farquhar-Wheat inputs/outputs to/from Pandas dataframes.