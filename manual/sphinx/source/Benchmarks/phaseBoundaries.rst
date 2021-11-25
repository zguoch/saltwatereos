.. _phaseboundaries:

**************************************
Phase boundaries
**************************************

.. include:: /include.rst_

Critical curve of H2O-NaCl system
=======================================

.. Note that the minimum valid temperature of the critical curve formula is the critical point of H2O, i.e. 373.946 |ssd|. 


**P_X_Critical**
    |P_X_Critical| calculate the critical pressure (P) and salinity (X) given temperature.

.. tab:: Chart

    .. include:: echarts/criticalCruve_H2ONaCl.rst

.. tab:: python

    .. code:: python

        import pyswEOS
        from pyswEOS import H2ONaCl
        from pyswEOS import H2O
        sw=H2ONaCl.cH2ONaCl()
        # calculate critical curve by giving temperature
        T=np.linspace(H2O.T_Critic,H2ONaCl.TMAX_C,100)
        p,x = sw.P_X_Critical(T)
        x,p = np.array(x)*100,np.array(p)

.. tab:: c++

    .. code:: c++

        #include "H2ONaCl.H"
        H2ONaCl::cH2ONaCl eos;
        double dT = (H2ONaCl::TMAX_C - H2O::T_Critic)/100;
        for (double T = H2O::T_Critic; T <= H2ONaCl::TMAX_C; T=T+dT)
        {
            double P,X;
            eos.P_X_Critical(T,P,X);
            cout<<T<<" "<<P<<" "<<X<<endl;
        }

.. tab:: Below water critical point

    The critical salinity at critical point of water is zero, 
    so what's the relationship between critical pressure given by Equation 5a of Driesner(2007a) and boiling pressure given by IAPWS formula (e.g. IAPWS-IF97 |TSat_P|) ? 
    The calculation results are quite different (up to 9 bar).

    .. include:: echarts/criticalPressure_belowCriticalPoint_H2O.rst