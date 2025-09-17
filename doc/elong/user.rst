
.. _elongwheat_user:

elong-Wheat User Guide
#########################

.. contents::

Introduction
============

Each growing leaves is represented by a set of compartments:

- a hiddenzone: encompasses growing and mature tissues of the growing leaf that are
  enclosed in the pseudo stem.
- any exposed (emerged) part of the lamina (photosynthetically active)
- any exposed part of the sheath (photosynthetically active).

Tissue emergence are calculated according to the length of the pseudostem which depends on
the above sheath and internodes of each growing leaf. Leaf elongation rate is regulated
by the amino acids and sucrose concentration inside the hiddenzone. The dynamics of leaf elongation
is coordinated with the emergence of the previous leaf:
- before leaf n-1 emergence, leaf n follows an exponential-like function with a strong metabolic regulation
- after leaf n-1 emergence, leaf n follows a predefined sigmoidal kinetics with low metabolic control

At leaf n emergence, its maximal width and surfacic mass are defined according to the sucrose concentration
of the hiddenzone averaged during to phyllochrons.

Elong-Wheat also simulates some functions of the shoot apical meristem : leaf primordia emission and floral transition
Elong-Wheat has on option to force the final dimensions of each leaves to those defined in :mod:`openalea.elongwheat.parameters` module

Inputs of elong-Wheat
========================

- Initial dimensions (m) and structural masses (g) of each hiddenzone, laminae, sheaths and internodes
- Initial content of the metabolites calculated by CN-Wheat (µmol)
- air and soil temperature (°C) : used to calculate the temperature of the hiddenzones and shoot apical meristem

Details on each inputs are given in the docstring.

Outputs of elong-Wheat
=========================
Updated leaf dimensions and shoot apical meristem status

Details on each outputs are given in the docstring.


Package architecture
=====================

Elong-Wheat is a Python package which consists of several Python modules:

* :mod:`openalea.elongwheat.model`: the state and the equations of the model,
* :mod:`openalea.elongwheat.parameters`: the parameters of the model,
* :mod:`openalea.elongwheat.simulation`: the simulator (front-end) to run the model,
* and :mod:`openalea.elongwheat.converter`: functions to convert Elong-Wheat inputs/outputs to/from Pandas dataframes.