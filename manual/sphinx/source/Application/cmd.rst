.. _cmd:

Command line tool
=============================

Synopsis
----------

.. include:: common_SYN_OPTs.rst_

**swEOS** [ |-D|\ [ *dimension* ]] 
[ |-V|\ [ *variables* ] ]
[ |-P|\ [:math:`p_{bar}`] ]
[ |-T|\ [:math:`T_{^{\circ}C}`] ]
[ |-X|\ [:math:`x_{wt. NaCl}`] ]
[ |-H|\ [:math:`h_{kJ/kg}`] ]
[ |-R|\ *min1/delta1/max1/min2/delta2/max2/min3/delta3/max3* ]
[ |-G|\ [*inputfile*] ]
[ |-O|\ [*outputfile*] ]

|No-spaces|


Description
------------

**swEOS** can calculate phase relations and thermodynamic properties of salt water in single point, one-dimension, two-dimension and three-dimension, respectively.

Required Arguments
------------------

.. _-D:

|-D|\ [ *dimension* ]
    Sets dimension. This is the first key option, the available arguments are **0**, **1**, **2** and **3**.

    * **0** means **single point** calculation. If |-D| is set to **0**, the |-V| option only support **PTX** or **PHX**. In addition, the pressure, salinity, temperature or enthalpy must be specified by |-P|, |-X|, |-T| or |-H| option, respectively.

    * **1** means only one variable changes and the others are set to fixed value. If |-D| is set to **1**, the |-V| option only support **P**, **T**, **X**, or **H**. The range of variable specified by |-V| option is set by |-R| option. And the other variables are set to fixed value by |-P|, |-X|, |-T| or |-H| options. The result will be saved to file. If the output file name is not specified by |-O| option, swEOS will use the default file name and print the file path in terminal.

    * **2** means change two variables and fixed the third variable. Similar to **-D1** case.

    * **3** means no fixed variable. Similar to **-D1** case.

.. _-V:

|-V|\ [ *variables* ]
    Sets variables accroding to |-D| option. This is the second key option, the available arguments are **PTX**, **PHX**, **P**, **T**, **X**, **H**, **PT**, **PX**, **TX**, **PH** and **HX**.

    * **PTX**, if argument of |-D| option is **1**, then the pressure, temperature and salinity is set by |-P|, |-T| and |-X| option, respectively. While if argument of |-D| is **3**, the range of pressure, temperature and salinity must be specified by |-R| option **in the same order** of argument of |-V| option. 