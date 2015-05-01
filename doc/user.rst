
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


To couple Farquhar-Wheat with other models, you also need NumPy_ (1.7.2 or newer), Alinea.Astk_ and 
OpenAlea.MTG_ libraries.

.. _NumPy: http://www.numpy.org/
.. _Alinea.Astk: https://scm.gforge.inria.fr/svn/openaleapkg/trunk/astk
.. _OpenAlea.MTG: https://scm.gforge.inria.fr/svn/vplants/vplants/trunk/newmtg/


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

