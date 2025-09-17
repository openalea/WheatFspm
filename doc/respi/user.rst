
.. _respiwheat_user:

Respi-Wheat User Guide
#########################

.. contents::

Introduction
============

Respi-Wheat simulates the different respiration fluxes related to the main physiological processes
described in Thornley and Cannell (2000).
The respiration fluxes described are related:
- the growth of shoot and root organs and grains
- phloem loading of sucrose by source organs
- nitrates, ammonium and other ions uptake by roots
- nitrates reduction
- N fixation for legumes only
- residual processes ((cost from protein turn-over, cell ion gradients, futile cycles...)


Inputs of Respi-Wheat
========================

Details on each inputs are given in the docstring.


Outputs of Respi-Wheat
=========================

Amount of C lost by respiration (µmol C) for each described process.

Details on each outputs are given in the docstring.


Package architecture
=====================

Respi-Wheat is a Python package which consists of several Python modules:

* :mod:`openalea.respiwheat.model`: the state and the equations of the model,