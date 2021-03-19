.. _install_windows:

*******************
Windows
*******************

.. include:: /include.rst_

Download 
=================

Please go to the |app_mac_download| to download the proper version (e.g. Catalina or newer) according to your system version.
If there is not a proper version, please look at the :ref:`build` section for help.

Install
===========

The downloaded installer is a :code:`.dmg` file.
Then you can get all the files (see :numref:`fig:installer:mac`) of swEOS app by simply double clicking the .dmg file and accepting the licence.

.. tab:: what's in the installer (snapshot)

    .. figure:: mac_installer_contents.png
        :align: center
        :name: fig:installer:win

        Files in the installer

.. tab:: what's in the installer (file tree)

    .. literalinclude:: mac_installer_contents.txt
        :language: terminfo
        :linenos:
        :caption: File trees of macOS installer of swEOS app
        :name: lst:installer:win
        :emphasize-lines: 42-48, 27-28

Desktop app with GUI
-------------------------

If you want to use the desktop app, just drag the :code:`swEOS` (see :numref:`fig:installer:mac`) to :code:`Applications` folder.
This "drag" install process is the same as any other macOS app.

.. tip:: 

    The app has not been notarized by Apple because there is no funding to support a Apple Developer ID. 
    Therefore the :code:`swEOS` app will be blocked by the Gatekeeper of macOS (see `Apple support <https://support.apple.com/en-us/HT202491>`_ for more details).
    In order to allow swEOS to run on your macOS, you have to run the following command with superuser permission(:code:`sudo`) in the terminal:

    .. code:: bash
    
        sudo xattr -r -d com.apple.quarantine /Applications/swEOS.app


Command line tool
---------------------

The standalone command line tool(cmd) is also included in the installer, 
if you want to use this cmd app, 
just need to copy :code:`commandline/swEOS` (see lines 27-28 in :numref:`lst:installer:mac`) file to some directory (e.g. /usr/local/bin) in your file system.
The you can use the cmd app in the for batch calculation purpose (see ).

.. tip::

    If there is a error information of :code:`permission denied: swEOS`, one can run command of :code:`chmod 755 /usr/local/bin/swEOS` to change its file model as an executable.



