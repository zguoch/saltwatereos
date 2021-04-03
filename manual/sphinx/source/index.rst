.. SaltWater EOS documentation front page

.. only:: html 

   .. raw:: html 
   
      <div class="bs-header">
         <div class="container">
            <h2 font face="Times">
            <b>swEOS:</b>
            <b style="color: #50B0F0">s</b>alt-
            <b style="color: #50B0F0">w</b>ater (NaCl-H2O)
            <b style="color: #50B0F0">E</b>quation
            <b style="color: #50B0F0">o</b>f
            <b style="color: #50B0F0">S</b>tate
            </h2>

            <table>
               <tbody><tr valign="top">
               <td width="10%"><b>What it is:</b></td>
               <td>
                  An extensible code written in C++ to support research
                  in using EOS and thermodynamic properties of H2O-NaCl system in both p-T-X and p-H-X coordinate space.
               </td>
               </tr>
               <tr valign="top">
               <td><b>Motivation:</b></td>
               <td>
                  To provide easy-to-use application, library and API for cross-platforms and programing languages. The users can implement their own research calculation using the library and APIs.
               </td>
               </tr>
            </tbody></table>
         </div>
      </div>


.. download page 

.. only:: html

   .. tab:: Phase changes animation

      .. raw:: html

         <video width=100% autoplay muted controls loop>
            <source src="./_static/video/PhaseChanges.mp4" type="video/mp4">
            Your browser does not support HTML video.
         </video>

   .. tab:: Phase diagram in P-T-X space

      .. raw:: html 

         <iframe src="_static/glance/swPhaseDiagram.html" width="100%" height="700px"></iframe>

   .. tab:: Phase diagram in X-H-P space

      .. raw:: html 

         <iframe src="_static/glance/swPhaseDiagram_XHP.html" width="100%" height="700px"></iframe>

Welcome to the SaltWater EOS Docs! Here you'll find resources for using SaltWater EOS and examples of what
it can do.

swEOS User Manual
=======================================


.. toctree::
   :maxdepth: 3
   :caption: Contents

   Introduction/index.rst
   Install/apps/index.rst
   Compile/index.rst
   Tutorial/index.rst
   Cookbooks/index.rst

.. bibliography:: manual.bib
    :cited: 
    :style: apa 