.. _build:

*******************
Build from source
*******************

.. include:: ../include.rst_

.. _build_lib:

Build c++ library 
=====================

One can download the source code and save it to wherever you like, assuming save the source to :code:`Download` download folder.

2. How to do ?

.. only:: html

   .. tab:: Animation on macOS
      
      .. raw:: html

         <script id="asciicast-9cN1baxW50xXDPBDKKjmU5mH0" src="https://asciinema.org/a/9cN1baxW50xXDPBDKKjmU5mH0.js" async></script>

.. tab:: macOS

   .. admonition:: Requirements: basic development environment

      * CMake_ : >=3.3
      * C++ compiler
      * make 
      * git (optional)

   .. code-block:: bash

      # 1. clone source code from github, or just download the source code and skip this step.
      git clone https://github.com/zguoch/saltwatereos.git

      # 2. check directory in to Library folder of the source code
      cd Library
      
      # 3. create a build folder
      mkdir build
      cd build
      
      # 4. cmake without building other APIs
      cmake -DBuild_API_MultiLanguage=OFF ..
      
      # 5. build the lib: you will get libeosH2ONaCl.a in the build folder
      make 
      
      # 6. install the lib to ../lib path, so the ../lib path will be created and libeosH2ONaCl.a is copied in this path
      make install

.. tab:: Linux

   .. admonition:: Requirements: basic development environment

      * CMake_ : >=3.3
      * C++ compiler
      * make 
      * git (optional)

   .. code-block:: bash

      # 1. clone source code from github, or just download the source code and skip this step.
      git clone https://github.com/zguoch/saltwatereos.git

      # 2. check directory in to Library folder of the source code
      cd Library
      
      # 3. create a build folder
      mkdir build
      cd build
      
      # 4. cmake without building other APIs
      cmake -DBuild_API_MultiLanguage=OFF ..
      
      # 5. build the lib: you will get libeosH2ONaCl.a in the build folder
      make 
      
      # 6. install the lib to ../lib path, so the ../lib path will be created and libeosH2ONaCl.a is copied in this path
      make install

.. tab:: Windows 10: Visual Studio

   .. admonition:: Requirements: basic development environment

      * CMake_ : >=3.3
      * |vs2017|
      * git (optional)

      The |vs2017| have to be installed and the path of `MSBUILD.exe` is added in the system environment PATH variable. All the following steps are performed in PowerShell.

   .. code-block:: bash

      # 1. clone source code from github, or just download the source code and skip this step.
      git clone https://github.com/zguoch/saltwatereos.git

      # 2. check directory in to Library folder of the source code
      cd Library
      
      # 3. create a build folder
      mkdir build
      cd build
      
      # 4. cmake without building other APIs
      cmake -DBuild_API_MultiLanguage=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_GENERATOR_PLATFORM=x64 ..

      # 5. build using msbuild.exe: then you will find eosH2ONaCl.lib is generated in the Release folder
      msbuild /m /p:Configuration=Release eosH2ONaCl.vcxproj

      # 6. install (optional): the eosH2ONaCl.lib is copied to ../lib folder
      msbuild /m /p:Configuration=Release INSTALL.vcxproj


3. What will you get ?

.. tab:: macOS

   .. image:: images/build_lib_macOS.*

.. tab:: Linux

   need to make snapshot

.. tab:: Windows 10: Visual Studio

   .. image:: images/build_lib_windows.*


.. _build_cmd:

Build standalone command line tool
=======================================

.. important:: 

   If you want to compile the standalone command line tool by yourself, you have to finish the previous step of :ref:`build_lib` firstly!

.. only:: html

   .. tab:: Animation on macOS
      
      .. raw:: html

         <script id="asciicast-AS1QAbLASVACpYpUHbMZaC4s0" src="https://asciinema.org/a/AS1QAbLASVACpYpUHbMZaC4s0.js" async></script>

.. tab:: macOS

   .. admonition:: Requirements: basic development environment

      * All requirements in :ref:`build_lib` step
      * OpenMP: could use `brew install libomp` and check `/usr/local/Cellar/libomp/*/include` path. The cmake will automatically detect the OpenMP include files and library files.

   .. code-block:: bash

      # 1. go to the commandline folder in the source code path 
      cd commandline

      # 2. create build folder and check directory to the build folder
      mkdir build
      cd build

      # 3. cmake 
      cmake ..

      # 4. build 
      make 

      # 5. check if swEOS is generated 
      ls 

      # 6. test it 
      ./swEOS -h

.. tab:: Linux

   .. admonition:: Requirements: basic development environment

      * All requirements in :ref:`build_lib` step
      * OpenMP: could use `sudo apt-get install libomp`. The cmake will automatically detect the OpenMP include files and library files.

   .. code-block:: bash

      # 1. go to the commandline folder in the source code path 
      cd commandline

      # 2. create build folder and check directory to the build folder
      mkdir build
      cd build

      # 3. cmake 
      cmake ..

      # 4. build 
      make 

      # 5. check if swEOS is generated 
      ls 

      # 6. test it 
      ./swEOS -h

.. tab:: Windows 10: Visual Studio

   .. admonition:: Requirements: basic development environment

      * All requirements in :ref:`build_lib` step

      
   .. code-block:: bash

      # 1. go to the commandline folder in the source code path 
      cd commandline

      # 2. create build folder and check directory to the build folder
      mkdir build
      cd build

      # 4. cmake 
      cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_GENERATOR_PLATFORM=x64 ..

      # 5. build using msbuild.exe: the swEOS.exe will generated in the Release folder
      msbuild /m /p:Configuration=Release swEOS.vcxproj

      # 6. Test it
      ./Release/swEOS.exe -h

Build APIs
=====================

.. warning:: 

   APIs of other programing language, e.g., python, tcl, js, depends on a lot of tools (e.g. swig, npm, python) and need a lot of programing skills. 
   We don't recommend users to try out this unless they master cross-platform and multi-language programing skills. But users can use the APIs for research, for example the python API :code:`pyswEOS` python users.

   Of course, users can find details in the source code if they are interested in that.

Build Desktop App with GUI 
================================

.. warning:: 

   Again, compilation of the desktop app with GUI also not that easy, because the GUI version is depends on Qt_ and VTK_.
   Therefore, it also needs a lot of programing skills to reach that.
   If users are interested in building the GUI version by themselves, please read the source code.  
