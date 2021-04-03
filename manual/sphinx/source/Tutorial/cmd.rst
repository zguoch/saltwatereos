.. _cmd:

Command line tool
=============================

Synopsis
----------

.. include:: common_SYN_OPTs.rst_

.. include:: ../include.rst_


**swEOS** [ |-D|  [ *dimension* ]] 
[ |-V|  [ *variables* ] ]
[ |-P|  [:math:`p_{bar}`] ]
[ |-T|  [:math:`T_{^{\circ}C}`] ]
[ |-X|  [:math:`x_{wt. NaCl}`] ]
[ |-H|  [:math:`h_{kJ/kg}`] ]
[ |-R|  *min1/delta1/max1/min2/delta2/max2/min3/delta3/max3* ]
[ |-G|  [*inputfile*] ]
[ |-O|  [*outputfile*] ]
[ |-t|  [*threads*] ]
[ |-n|  ]
[ |-hh|  ]

|No-spaces|


Description
------------

**swEOS** can calculate phase relations and thermodynamic properties of salt water in single point, one-dimension, two-dimension and three-dimension, respectively.

Required Arguments
------------------

.. _-D:

|-D| [ *dimension* ]
    Sets dimension. This is the first key option, the available arguments are **0**, **1**, **2** and **3**.

    * **0** means **single point** calculation. If |-D| is set to **0**, the |-V| option only support **PTX** or **PHX**. In addition, the pressure, salinity, temperature or enthalpy must be specified by |-P|, |-X|, |-T| or |-H| option, respectively.

    * **1** means only one variable changes and the others are set to fixed value. If |-D| is set to **1**, the |-V| option only support **P**, **T**, **X**, or **H**. The range of variable specified by |-V| option is set by |-R| option. And the other variables are set to fixed value by |-P|, |-X|, |-T| or |-H| options. The result will be saved to file. If the output file name is not specified by |-O| option, swEOS will use the default file name and print the file path in terminal.

    * **2** means change two variables and fixed the third variable. Similar to **-D1** case.

    * **3** means no fixed variable. Similar to **-D1** case.

.. _-V:

|-V|  [ *variables* ]
    Sets variables accroding to |-D| option. This is the second key option, the available arguments are **PTX**, **PHX**, **P**, **T**, **X**, **H**, **PT**, **PX**, **TX**, **PH** and **HX**.

    * **PTX**, if argument of |-D| option is **0** (single point or multiple points case), then the pressure, temperature and salinity is set by |-P|, |-T| and |-X| option, respectively. While if argument of |-D| is **3** (three dimension case), the range of pressure, temperature and salinity must be specified by |-R| option **in the same order** of argument of |-V| option. Therefore, for the same calculation, |-V| option can be **PTX**, **PXT**, **TPX**, **TXP**, **XPT** and **XTP** unless you set |-R| option in the save order. For example, :code:`-VXPT -R0/0.1/1/5/1/400/0/1/100` means salinity in range of [0, 1] with interval of 0.1, pressure in range of [5, 400] bar with interval of 1 bar, temperature in range of [0, 100] :math:`^{\circ}\text{C}` with interval 1 :math:`^{\circ}\text{C}` . Alternately, you can also do the same thing using command of |-V|  *TPX* |-R|  *0/1/100/5/1/400/0/0.1/1* .

    * **PHX**, similar to **PTX**, it just replaced temperature with enthalpy.

    * **T**. This is only valid when argument of |-D| option is **1** (one dimension case). It means temperature is the independent variable, its range is specified by |-R| option, e.g. :code:`-R0/1/100` means temperature in range of [0, 100] |ssd| with interval of 1 |ssd|,  pressure and salinity are fixed by |-P| and |-X| option, respectively. In addition, the output file name has to be specified by |-O| option, will write as csv file format. **P**, **X**, **H** similar to **T**.

    * **PT**. This is only valid when argument of |-D| option is set to **2** (two dimension case). It means pressure and temperature are the independent variable, their range are specified by |-R| option, e.g. :code:`-R5/1/500/0/1/100` means temperature in range of [0, 100] |ssd| with interval of 1 |ssd|, pressure in range of [5, 500] bar with interval of 1 bar. Salinity is the fixed variable and the fixed value is specified by |-X| option. Note that it doesn't matter what the order of variable follow |-V| option, e.g. :code:`-VPT` and :code:`-VTP` are all valid, but the order of range follow |-R| option matters, it must be the same order with |-V| option. For example, :code:`-VPT -R5/1/500/0/100` and :code:`-VTP -R0/1/100/5/1/500` are equivalent. Of course, you also have to specify output file by |-O| option. **PX**, **TX**, **PH**, **HX** are similar with **PT**.

.. _-P:

|-P|  [ *pressure* ]
    Sets fixed pressure value, it should be a float or integer number. The unit is **bar**. 

.. _-T:

|-T|  [ *temperature* ]
    Sets fixed temperature value, it should be a float or integer number. The unit is |ssd|. 

.. _-X:

|-X|  [ *salinity* ]
    Sets fixed salinity value, it should be a float or integer number. The unit is :math:`wt. NaCl`. 

.. _-H:

|-H|  [ *enthalpy* ]
    Sets fixed enthalpy value, it should be a float or integer number. 

.. _-R:

|-R|  [ *min/delta/max* ]
    Sets range and interval of independent variable(s), |-R| option must correspond to |-V| and |-D| options. For example, :code:`-D1 -VT -R0/1/100` , :code:`-D2 -VTX -R0/1/100/0/0.1/0.8` , :code:`-D3 -VTXP -R0/1/100/0/0.1/0.8/5/1/500` .

.. _-G:

|-G|  [ *inputfile* ]
    Sets input file for multiple points calculation, the input file with three columns with delimiter of space or table( :code:`\t` ) correspondint to |-V| option. For example, :code:`-D0 -VPTX -Ginput.txt` means calculate EOS of some points, independent variables of each point are pressure, temperature and salinity which are listed in :code:`input.txt` as belows

    .. tab:: PTX

        .. code-block:: bash
            :name: lst:inputfile:ptx
            :caption: Input file for multi-points calculation, three columns are :code:`p(bar) T(deg.C) X(wt% NaCl)`

            316	10	0.032
            316	11	0.032
            316	12	0.032
            316	13	0.032
            316	14	0.032
            316	15	0.032
            316	16	0.032
            316	17	0.032
            316	18	0.032
            316	19	0.032
            316	20	0.032
            316	21	0.032
            316	22	0.032
            316	23	0.032
    
    .. tab:: PHX

            .. code-block:: bash
                :name: lst:inputfile:phx
                :caption: Input file for multi-points calculation, three columns are :code:`p(bar) H(kJ/kg) X(wt% NaCl)`
    
                316	74.3015	0.032
                316	78.2695	0.032
                316	82.2389	0.032
                316	86.2098	0.032
                316	90.1822	0.032
                316	94.1561	0.032
                316	98.1316	0.032
                316	102.109	0.032
                316	106.087	0.032
                316	110.067	0.032
                316	114.049	0.032
                316	118.031	0.032
                316	122.016	0.032
                316	126.001	0.032
                316	126.001	0.032

.. _-O:

|-O|  [ *outputfile* ]
    Sets output file name for one, two and three-dimensional calculation. The supported file format for 1D is csv, for 2D and 3D can be one of txt(delimiter is :code:`\t` ), csv(delimiter is :code:`,` ) and vtk. Because temperature, pressure and salinity have different scaling, if open the vtk file directlly by paraview, you can not see anything in salinity dimension. Therefore there is a python script generated by :code:`swEOS` can deal with the scaling issue automatically. You just need to type a command in terminal to visualize the result. e.g.  :code:`paraview --script=test3D.vtk.py` .

.. _-tt:

|-t|  [ *thread* ]
    Sets number of threads for parallel computing.

.. _-n:

|-n| 
    If need to normalize the coordinates of result in vtk file, it is only valid in 3D calculation. Because the coordinates of the results are (P,T,X), (P,H,X), ... with different scale, if use paraview to visualize the result, one have to deal with the coordinate scale. One option is set using |-n| to normalize the coordinate. Otherwise, a python script will be generated by swEOS automatically, one can use this python script to visualize results (see also :ref:`cmd_3D`).

    .. code-block:: console
        :emphasize-lines: 5,6

        3D calculation using 8 threads, H in [100, 600] kJ/kg, P in [10, 500] bar, X in [0.001, 1]
        [################################]100%
        Writing results to file ...
        Results have been saved to file: PTX.vtk
        Paraview-python script is generated as : PTX.vtk.py
        You can use command of  paraview --script=PTX.vtk.py to visualize result in paraview

.. _-hh:

|-hh| 
    Print a short message about the syntax of the command. 


Examples
------------------------

Single point calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash 

    swEOS -D 0 -V PXT -P 316 -T 100 -X 0.032
    swEOS -D 0 -V PXH -H 438 -P 316 -X 0.032

.. raw:: html

    <script id="asciicast-z13UZmFuRN7U39tVguaBkXF12" src="https://asciinema.org/a/zCb09bq5gI4UKZYGPHMOxYdrK.js" async data-autoplay="false"></script>

Multi-points calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash 

    swEOS -D 0 -V PHX -G phx.txt -O PHX_0D.csv
    swEOS -D 0 -V PTX -G ptx.txt -O PTX_0D.csv

.. only:: html

    .. tab:: Usage animation 

        .. raw:: html

            <script id="asciicast-vevYoK7t4jfuODEfBKfQfGBZ6" src="https://asciinema.org/a/zCb09bq5gI4UKZYGPHMOxYdrK.js" async data-autoplay="false"></script>

    .. tab:: Input: PTX

        .. code-block:: bash
    
            316	10	0.032
            316	11	0.032
            316	12	0.032
            316	13	0.032
            316	14	0.032
            316	15	0.032
            316	16	0.032
            316	17	0.032
            316	18	0.032
            316	19	0.032
            316	20	0.032
            316	21	0.032
            316	22	0.032
            316	23	0.032
        
    .. tab:: Result: PTX
    
        .. csv-table:: Calculation result with `PTX` input file
            :file: res_0D_ptx.csv

    .. tab:: Input: PHX

        .. code-block:: bash

            316	74.3015	0.032
            316	78.2695	0.032
            316	82.2389	0.032
            316	86.2098	0.032
            316	90.1822	0.032
            316	94.1561	0.032
            316	98.1316	0.032
            316	102.109	0.032
            316	106.087	0.032
            316	110.067	0.032
            316	114.049	0.032
            316	118.031	0.032
            316	122.016	0.032
            316	126.001	0.032
            316	126.001	0.032
        
    .. tab:: Result: PHX

        .. csv-table:: Calculation result with `PHX` input file
            :file: res_0D_phx.csv

One-dimensional calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash 

    swEOS -D 1 -V H -X 0.032 -P 399 -R 43/1/100 -O H_1D.csv
    swEOS -D 1 -V T -X 0.032 -P 399 -R 0/1/100 -O T_1D.csv
    swEOS -D 1 -V P -X 0.032 -T 100 -R 5/1/100 -O P_1D.csv
    swEOS -D 1 -V X -T 100 -P 399 -R 0/0.001/1 -O X_1D.csv

.. only:: html

    .. tab:: Usage animation 

        .. raw:: html

            <script id="asciicast-iIdDC7uC5lvoHVvsBVgWyi64Y" src="https://asciinema.org/a/zCb09bq5gI4UKZYGPHMOxYdrK.js" async data-autoplay="false"></script>

    .. tab:: Results: change T

        .. csv-table:: Fix P and X, change T
            :file: T_1D.csv

    .. tab:: Results: change P

        .. csv-table:: Fix T and X, change P
            :file: P_1D.csv

    .. tab:: Results: change X

        .. csv-table:: Fix P and T, change X
            :file: X_1D.csv

    .. tab:: Results: change H

        .. csv-table:: Fix P and X, change H
            :file: H_1D.csv

Two-dimensional calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash 

    swEOS -D 2 -V PT -R 10/0.1/100/1/1/500 -X 0.032 -O PT_2D.vtk
    swEOS -D 2 -V PX -R 100/0.1/800/0/0.01/1 -T 100  -O PX_2D.vtk
    swEOS -D 2 -V TX -R 1/1/800/0/0.01/1 -P 100  -O TX_2D.vtk
    swEOS -D 2 -V PH -R 100/1/800/100/1/700 -X 0.032  -O PH_2D.vtk
    swEOS -D 2 -V XH -R 0.001/0.001/1/100/1/700 -P 200  -O XH_2D.vtk

.. raw:: html

    <script id="asciicast-7p6v2o8fux1lvW2T48o3oBsev" src="https://asciinema.org/a/zCb09bq5gI4UKZYGPHMOxYdrK.js" async data-autoplay="false"></script>

.. _cmd_3D:

Three-dimensional calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash 

    swEOS -D 3 -V PTX -R 10/10/500/1/10/600/0/0.01/1 -O PTX.vtk
    swEOS -D 3 -V PHX -R 10/10/500/100/10/600/0.001/0.01/1 -O PTX.vtk

.. raw:: html

    <script id="asciicast-t2PDdWGa13XRh76tlB4OYGpFu" src="https://asciinema.org/a/zCb09bq5gI4UKZYGPHMOxYdrK.js" async data-autoplay="false"></script>


.. admonition:: Phase region index table

    +-------------------+--------------------------------------+
    |Phase region index | Phase region                         |
    +-------------------+--------------------------------------+
    |0                  | Single phase(Liquid)                 |
    +-------------------+--------------------------------------+
    |1                  | Liquid + Vapor at X=0                |
    +-------------------+--------------------------------------+
    |2                  | Pure vapour phase                    |
    +-------------------+--------------------------------------+
    |3                  | Liquid + Halite                      |
    +-------------------+--------------------------------------+
    |4                  | Vapour + Halite                      |
    +-------------------+--------------------------------------+
    |5                  | Vapour + Liquid + Halite             |
    +-------------------+--------------------------------------+
    |6                  | Vapour + Liquid on the liquid side   |
    +-------------------+--------------------------------------+
    |7                  | Vapour + Liquid on the vapour side   |
    +-------------------+--------------------------------------+