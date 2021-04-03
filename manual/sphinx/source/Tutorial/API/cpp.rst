.. _API_tutorial_cpp:

.. include:: /include.rst_

C++
=========

The C++ API is provided by head files and a library file named :code:`libeosH2ONaCl.a` for macOS and linux,  :code:`libeosH2ONaCl.lib` for windows system, respectively. The head files and library file are contained in the app installer (in the :code:`include` and :code:`lib` folder), see also :ref:`install` section.

Before starting to use c++ api of swEOS, the CMake_ and c++ compiler have to be installed.

How to start ?
---------------------

Assuming the library file and head files has been downloaded and saved in :code:`~/Download/swEOS` (please use your own path) path, it means that the :code:`lib` folder and :code:`include` folder are in :code:`~/Download/swEOS` path.

**Step 1.** Create a source code folder and check directory to this folder.

**Step 2.** Create c++ source code, e.g. :download:`main.cpp <cpp/main.cpp>`, see :numref:`lst:tutorial:api:cpp`. Include the head file(line 1) and instantiate a object (line 5) of :code:`cH2ONaCl` class, then all the properties and member functions can be accessed through object :code:`eos`, e.g. density (line 9).

.. literalinclude:: cpp/main.cpp
    :language: c++
    :linenos:
    :emphasize-lines: 1,5,9
    :caption: Source code of C++ api for calculating density of H2ONaCl.
    :name: lst:tutorial:api:cpp


**Step 3.** Create :download:`CMakeLists.txt <cpp/CMakeLists.txt>`, see :numref:`lst:tutorial:api:cpp:cmake`.

.. literalinclude:: cpp/CMakeLists.txt
    :language: cmake
    :linenos:
    :emphasize-lines: 21-22
    :caption: CMake file
    :name: lst:tutorial:api:cpp:cmake

**Step 4.** Configure and generate project using CMake_ command line tool or GUI app.

Please remember set cmake cache variable of :code:`SWEOS_DIR` to specify the swEOS library path which contains :code:`lib` and :code:`include` folders.

.. code-block:: bash

    mkdir build
    cd build
    cmake cmake -DSWEOS_DIR=~/Download/swEOS ..

.. note:: 

    If :code:`SWEOS_DIR` is not set or set incorrectly, you will get the following error information. Then you just make sure set a correct path for :code:`SWEOS_DIR` to fix the problem.

    .. code:: console 

        CMake Error at CMakeLists.txt:14 (message):
        Please specify path of SWEOS library which contains lib and include paths
    
            cmake -DSWEOS_DIR=path_of_SWEOS .. 
    
        -- Configuring incomplete, errors occurred!


**Step 5.** Compile and build the program.

In macOS or Linux system, just run :code:`make` to compile and build the program. While, in Windows system and using |vs2017|, you will get a :code:`.sln` file. Now, every VS user should know how to do!

**Step 6.** Run the program.

The executable program name is :code:`test_swEOS`, which is configured in the CMakeLists file (see line 4 of :numref:`lst:tutorial:api:cpp:cmake`). Run this program in the terminal, you will get the output like this,

.. code-block:: console 

    ./test_swEOS
    
    Pressure(bar): 200
    Temperature(deg.C): 400
    Salinity (wt.% NaCl): 0.2
    Density(kg/m3): 188.056

More features
---------------------

All right, now you must know how to use c++ api of swEOS package in your own program, see :ref:`cookbooks` for usage of more functions.