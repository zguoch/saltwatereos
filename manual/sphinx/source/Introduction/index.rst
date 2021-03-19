.. _introduction:

*******************
Introduction
*******************

.. include:: ../include.rst_

Realistic simulations of fluid flow in natural hydrothermal systems require accurate formulation of fluid properties
especially for seawater convection at mid-ocean ridges. 
In order to explore circulation patterns and realistic phase seperation phenomenon, 
we have to calculate equation of state of binary salt-water fluids over pressure-temperature-salinity ranges encountered in the Earth's crust.
Fortunately, pure water can be described by the IAPS-84 equation of state and :cite:`Driesner2007Part1,Driesner2007Part2` have developed a set of correction formulations of phase relations and thermodynamic properties for NaCl-H2O system.
Further, we have developed a set of multi-language ( `C++ <https://www>`_ , `Swift <https://www>`_ , `Python <https://www>`_ , `Matlab <https://www>`_ ) and multi-platform ( `Windows <https://www>`_ , `MacOS <https://www>`_ , `Linux <https://www>`_ , `IOS <https://www>`_ ) tools, including callable `C++ library <https://www>`_ , `desktop application <../Application/desktop.html>`_ with graphical user interface (GUI), `command line tools <https://www>`_ (just like `gmt <http://gmt.soest.hawaii.edu>`_ style), and Mobile apps for `iphone <https://www>`_ and `ipad <https://www>`_ . In addition, parallel computing is available for desktop application and command line tool.

Platforms
=============

   +----------------+----------------+----------------+-----------+
   | |mac|          | |linux|        | |windows|      | |ios|     |
   +----------------+----------------+----------------+-----------+

Applications
=================

   +----------------------------+------------------------------+---------------------+
   | |desktop|                  | |cmd|                        | |mobile|            |
   +----------------------------+------------------------------+---------------------+

Supported Programing Languages
==================================

   +----------------+----------------+----------------+----------------+----------------+----------------+
   | |c++|          | |matlab|       | |python|       |    |swift|     |    |JS|        |      |tcl|     |
   +----------------+----------------+----------------+----------------+----------------+----------------+

Authors and Developers
============================


Cite
=======

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4603878.svg
   :target: https://doi.org/10.5281/zenodo.4603878

.. code:: bib

   @software{zhikui_guo_2021_4603878,
      author       = {Zhikui Guo and
                     Lars Rüpke},
      title        = {{swEOS: multi-platform multi-language package of 
                        salt-water equation of state}},
      month        = mar,
      year         = 2021,
      publisher    = {Zenodo},
      version      = {1.7.0},
      doi          = {10.5281/zenodo.4603878},
      url          = {https://doi.org/10.5281/zenodo.4603878}
   }

.. |c++| image:: /_static/logo/cpp.*
   :width: 50 px
   :align: middle
   :target: #

.. |matlab| image:: /_static/logo/matlab.*
   :width: 50 px
   :align: middle
   :target: #

.. |python| image:: /_static/logo/python.*
   :width: 50 px
   :align: middle
   :target: #

.. |swift| image:: /_static/logo/swift.*
   :width: 50 px
   :align: middle
   :target: #

.. |JS| image:: /_static/logo/js.*
   :width: 50 px
   :align: middle
   :target: #

.. |tcl| image:: /_static/logo/tcl.*
   :width: 30 px
   :align: middle
   :target: #

.. |mac| image:: /_static/logo/mac.*
   :width: 80 px
   :align: middle
   :target: #

.. |linux| image:: /_static/logo/linux.*
   :width: 80 px
   :align: middle
   :target: #

.. |windows| image:: /_static/logo/windows.*
   :width: 80 px
   :align: middle
   :target: #

.. |ios| image:: /_static/logo/ios.*
   :width: 80 px
   :align: middle
   :target: #

.. |desktop| image:: /_static/logo/apps_desktop.png
   :width: 300 px
   :align: middle
   :target: #

.. |cmd| image:: /_static/logo/apps_cmd.png
   :width: 300 px
   :align: middle
   :target: #

.. |mobile| image:: /_static/logo/apps_mobile.png
   :width: 300 px
   :align: middle
   :target: #
