
.. _cnwheat_user:

CN-Wheat User Guide
###################

.. contents::

Introduction
============

TODO


Installation
============

Prerequisites
-------------

CN-Wheat requires the following software installed for your platform:

1) Python__ 2.7.x or newer

__ http://www.python.org/

2) NumPy__ 1.7.2 or newer

__ http://www.numpy.org/

3) SciPy__ 0.12.1 or newer

__ http://docs.scipy.org/

4) Pandas__ 0.14.0 or newer

__ http://pandas.pydata.org/

5) If you want to plot graph using :mod:`cnwheat.post_processing`: Matplotlib__ 1.2.0 or newer

__ http://matplotlib.sourceforge.net/

6) If you want to run the tests: Nose__ 1.3.0 or newer

__ http://nose.readthedocs.org/

7) If you want to build the documentation: Sphinx__ 1.1.3 or newer

__ http://sphinx-doc.org/


Getting the sources
===================
Use the command::

  svn checkout https://subversion.renater.fr/cn-wheat
  
This creates the directory ``cn-wheat``.


Install in user mode
--------------------
::

  cd cn-wheat
  python setup.py install
  
  
Install in developer mode
-------------------------
::

  cd cn-wheat
  python setup.py develop


Run the tests
-------------
Install CN-Wheat and run::

  cd cn-wheat
  nosetests


Build the documentation
-----------------------

To build the documentation, just run::

    python setup.py build_sphinx

Then open the file :download:`doc/_build/html/index.html <_build/html/index.html>` 
with your web browser.


Getting started
===============

TODO


Inputs of CN-Wheat
==================

TODO


Outputs of CN-Wheat
===================

TODO

