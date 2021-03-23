.. _install_linux:

*******************
Linux
*******************

.. include:: /include.rst_


Download and install
==================================

For Linux system, the install process is pretty easy.
Similar to :ref:`install_mac`, download the :code:`.zip` installer from |app_linux_download| and unzip it to wherever you like, e.g. :code:`/home/swEOS`. 

Desktop app with GUI
-------------------------

Run the following command in terminal to launch the swEOS GUI version.

.. code:: bash 

    /home/swEOS/swEOS.sh

Command line tool
-------------------------
    
There are two ways to use the command line tool.

1. Set arguments and options directory after :code:`swEOS.sh`, for example :code:`/home/swEOS/swEOS.sh -h`, because the command line arguments and options is also integrated in the GUI version.
    
2. Use the standalone command line tool, for example :code:`/home/swEOS/CommandLineTool/swEOS`

.. tip:: 

    1. If there is a error information of :code:`permission denied` after run :code:`/home/swEOS/CommandLineTool/swEOS` command, please set the swEOS command line tool as executable mode. This can be done by running :code:`chmod 755 /home/swEOS/CommandLineTool/swEOS`
   
    2. If you want to access swEOS much easier, please add the command line tool path to the environment variable, just run command of :code:`echo "PATH="/home/swEOS/CommandLineTool:$PATH"" >>~/.bashrc`.