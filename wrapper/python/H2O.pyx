# distutils: language = c++
# distutils: sources ="../../Library/src/H2O.C"

# 此处的名字取决于相应的.pxd文件的名字
cimport H2O

cdef class IAPWS95:
    cdef cH2O c_H2O  
    def __cinit__(self):
        self.c_H2O = cH2O()
    def P_Boiling(self, P):
        return self.c_H2O.P_Boiling(P)
    def Rho_Liquid_Saturated(self, T):
        return self.c_H2O.Rho_Liquid_Saturated(T)
    def Rho_Vapor_Saturated(self, T):
        return self.c_H2O.Rho_Vapor_Saturated(T)
    def Rho(self, T, P):
        return self.c_H2O.Rho(T,P)
    def Pressure_T_Rho(self, T, Rho):
        return self.c_H2O.Pressure_T_Rho(T, Rho)
    def SublimationCurve(self, T):
        return self.c_H2O.SublimationCurve(T)
    def BoilingCurve(self, T):
        return self.c_H2O.BoilingCurve(T)
    def MeltingCurve(self, T, isIceI=False):
        return self.c_H2O.MeltingCurve(T, isIceI)