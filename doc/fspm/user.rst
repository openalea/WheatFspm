
.. _fspmwheat_user:

fspm-Wheat User Guide
#########################

.. contents::

Introduction
============

Fspm-Wheat is a subpackage used to couple the models cn-wheat, farquhar-wheat, elong-wheat,
growth-wheat, respi-wheat, senesc-wheat, adel-wheat and caribu. Each model has a specific facade used to read
its inputs from a MTG object and write its outputs after running.
The MTG is provided by adel-wheat model.


Package architecture
=====================

fspm-Wheat is a Python package which consists of several Python modules:

* :mod:`openalea.fspmwheat.caribu_facade`: the interface between Caribu and the MTG
* :mod:`openalea.fspmwheat.cnwheat_facade`: the interface between CN-Wheat and the MTG
* :mod:`openalea.fspmwheat.elong_facade`: the interface between Elong-Wheat and the MTG
* :mod:`openalea.fspmwheat.farquhar_facade`: the interface between Farquhar-Wheat and the MTG
* :mod:`openalea.fspmwheat.fspmwheat_facade`: the interface between FSPM-Wheat and the MTG
* :mod:`openalea.fspmwheat.fspmwheat_postprocessing`: the postprocessings of Fspm-Wheat
* :mod:`openalea.fspmwheat.growth_facade`: the interface between Growth-Wheat and the MTG
* :mod:`openalea.fspmwheat.senesc_facade`: the interface between Senesc-Wheat and the MTG
* :mod:`openalea.fspmwheat.tools`: some tools
