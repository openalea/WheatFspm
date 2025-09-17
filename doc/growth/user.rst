
.. _growthwheat_user:

Growth-Wheat User Guide
#########################

.. contents::

Introduction
============

Growth-Wheat simulates leaf growth in mass according to leaf elongation calculated in Elong-Wheat.
Leaf dimension increase of first interpreted by Adel-wheat model, which provides a variation in leaf area.
For non-emerged leaves, a empiric and constant relation between leaf length and its strucral mass is used.
For emerged leaves, the structural mass growth is computed from their green area and
a surfacic mass defined by Elong-Wheat. According to the computed increase in structral mass, Growth-Wheat
calculates the related consumption of sucrose and amino acids. Growth-Wheat also handles the transfer of CN metabolites between the
hiddenzone and the newly emerged part of the leaf if any.
Root growth and the related metabolite consumption are regulated by the local concentration of
sucrose. The maximal rate of root growth also depends on the reproductive status of the plant.


Inputs of Growth-Wheat
========================

- Initial structral masses (g), dimensions (m) of each hiddenzone + the green area for laminae, sheaths and internodes
- Initial content of the metabolites calculated by CN-Wheat (µmol)

Details on each inputs are given in the docstring.

Outputs of Growth-Wheat
=========================

Updated structural masses for hiddenzones, shoot organs and roots
Updated CN metabolites contents

Details on each outputs are given in the docstring.


Package architecture
=====================

Growth-Wheat is a Python package which consists of several Python modules:

* :mod:`openalea.growthwheat.model`: the state and the equations of the model,
* :mod:`openalea.growthwheat.parameters`: the parameters of the model,
* :mod:`openalea.growthwheat.simulation`: the simulator (front-end) to run the model,
* and :mod:`openalea.growthwheat.converter`: functions to convert Growth-Wheat inputs/outputs to/from Pandas dataframes.