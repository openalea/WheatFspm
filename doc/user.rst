
.. _farquharwheat_user:

Farquhar-Wheat User Guide
#########################

.. contents::

Introduction
============

This is the documentation for Farquhar-Wheat, a model of photosynthesis based on Farquhar's approach. 

Prerequisites
-------------

Farquhar-Wheat needs at least Python_ (2.7 or newer) to run and Pandas_ (0.14.0 or newer) to format outputs.
Farquhar-Wheat has not been tested with **Python 3**. 

.. _Python: http://www.python.org/
.. _Pandas: http://pandas.pydata.org/


Usage
-----

See :ref:`getting_started` for an introduction. 


Installing
==========

First get the sources using ``svn``:: 

  svn checkout https://subversion.renater.fr/farquhar-wheat
  
This creates the directory ``farquhar-wheat``.

Then, in the directory ``farquhar-wheat``, run::

  python setup.py install
  
Or, to install in develop mode, run::

  python setup.py develop
  

.. _getting_started:


Getting started
===============

TODO


Inputs of Farquhar-Wheat
========================

TODO


Outputs of Farquhar-Wheat
=========================

TODO


Farquhar inputs from MTG properties 
===================================

===========================  ==========================================================  
   Farquhar inputs									MTG properties     											
===========================  ==========================================================  
surfacic_nitrogen		 	 surfacic_nitrogen											
width					 	 blade						
width 					 	 other organs diameter
height of hidden element	 height of the precedent sheath					
height of visible element	 mean(plantgl_utils.get_height(geometry))
organ_type					 label of the vertex									
PAR							 exposed_area / area * PARi 									
===========================  ==========================================================  


MTG properties from Farquhar outputs 
====================================

======================  ====================
   MTG properties		  Farquhar outputs     	
======================  ====================
Ag						Ag
An						An
Tr						Tr
Ts						Ts
gs						gs
======================  ====================
