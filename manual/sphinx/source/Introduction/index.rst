.. _introduction:

*******************
Introduction
*******************

.. include:: ../include.rst_

HydrothermalFoam —combination of Hydrothermal and OpenFOAM — is a series of programs or toolbox tended to solve the equations that describe natural hydrothermal convection and minerial reactions in porous media, fractured-porous media. It is primarily developed by Zhikui Guo and Lars Rüpke based on OpenFOAM_, which is a open source C++ library for CFD (Computational Fluid Dynamics). 
See :cite:`kawada2010formation` for an introduction to non-standard analysis :cite:`kawada2010formation` .



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

