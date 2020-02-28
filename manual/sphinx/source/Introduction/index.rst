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
Further, we have developed a set of multi-language ( `C++ <https://www>`_ , `Swift <https://www>`_ , `Python <https://www>`_ , `Matlab <https://www>`_ ) and multi-platform ( `Windows <https://www>`_ , `MacOS <https://www>`_ , `Linux <https://www>`_ , `IOS <https://www>`_ ) tools, including callable `C++ library <https://www>`_ , `desktop application <https://www>`_ with graphical user interface (GUI), `command line tools <https://www>`_ (just like `gmt <http://gmt.soest.hawaii.edu>`_ style), and Mobile apps for `iphone <https://www>`_ and `ipad <https://www>`_ . In addition, parallel computing is available for desktop application and command line tool.

Language 
=============

+-----------+----------------+-----------------+
| |c++|     | |matlab|       | |python|        |
+-----------+----------------+-----------------+

.. |c++| image:: /_images/c++.png
   :height: 150 px
   :target: https://stackoverflow.com/questions/14087784/linked-image-in-restructuredtext

.. |matlab| image:: /_images/matlab.png
   :height: 150 px

.. |python| image:: /_images/python.png
   :height: 150 px

Platform 
=============

Application 
=============

.. figure:: /_figures/filetree_main.*
   :width: 500 px
   :align: center

   Folder/file layout of HydrothermalFoam.

The second part of the script shows how to make the color version of

Features of the HydrothermalFoam are summarized as following:

- Original characteristics of OpenFOAM: HydrothermalFoam keeps all the original characteristics of Open- Foam, for example, file structure of case, syntax of all the input files and output files, mesh, utilities, even part of varable names in the source code are kept the same. This principle has two advantages, one is that it is easy to understand and to use HydrothermalFoam if you are a OpenFoam user. The other is that it is easy to compar and understand HydrothermalFoam solver and the other standard solvers in OpenFoam.

- Mesh: HydrothermalFoam supports both structured regular mesh and unsructured mesh. The internal structured regular mesh tool, blockMesh, is recommended for new users. While Gmsh_ is also an excellent open source unstuctured mesh generator, and there is a utility named gmshToFoam can transfoam gmsh to OpenFoam mesh.

- Boundary conditions: even though OpenFoam has a lot of build-in boundary conditions, we also developed some specific boundary conditions, e.g. fixedHeatFlux, fixedMassFlux for the specific problem — Hydrothermal system.



Acknowledgments
-----------------

Cite
=======

