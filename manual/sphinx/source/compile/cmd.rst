.. _build_cmd:

Command line tool
=============================

.. include:: ../include.rst_

It's pretty easy to build command line version of swEOS from source code.

1. Download or clone source code 

.. code-block:: bash 

    git clone https://gitlab.com/hydrothermal-openfoam/saltwatereos.git

2. Change directory to the source code 

.. code-block:: bash 

    cd saltwatereos

3. Compile c++ library

.. code-block:: bash

    cd Library
    mkdir build 
    cd build 
    cmake ..
    make install

Then you will get :code:`libeosH2ONaCl.a` static library in :code:`saltwatereos/Library/lib` folder.

4. Build command line tool

.. code-block:: bash 

    cd ../commandline
    mkdir build
    cd build
    cmake ..
    cmake ..
    make 

.. attention::

    If there are some error information about openmp during cmake, please just type :code:`cmake ..` again.

Then you will get the executable file :code:`swEOS` in the current path. In tutorial of :ref:`cmd` you can find detail usage and example.


