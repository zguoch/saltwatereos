from libcpp cimport bool

# dependence header file
cdef extern from "Fluid.H":
    pass

# Declare the class with cdef
cdef extern from "H2O.H" namespace "H2O":
    cdef cppclass cH2O:
        cH2O() except +
        double P_Boiling(double T); 
        double Rho_Liquid_Saturated(double T);
        double Rho_Vapor_Saturated(double T);
        double Rho(double T, double P);
        double Phi_r_delta(double delta, double tau);
        double Pressure_T_Rho(double T, double Rho);
        double SublimationCurve(double T);
        double BoilingCurve(double T);
        double MeltingCurve(double T, bool isIceI);