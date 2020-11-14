' Akinfief and Diamond 2009
Public Function Qtz_solubility(xnacli#, Ti#, rhoi#) As Double
Dim a#(0 To 3), At#, b#(0 To 3), Bt#, mH2O#, mNaCl#, VH2O#, T!, LogSolub#, xNaCl!

Dim i%

mH2O = 18.015268
mNaCl = 58.4428
xNaCl = xnacli

T = Ti + 273.15
a(0) = 4.262:  a(1) = -5764.2: a(2) = 1751300#: a(3) = -228690000#
b(0) = 2.8454: b(1) = -1006.9: b(2) = 356890#: b(3) = 0

For i = 0 To 3
    At = At + a(i) * T ^ -i
    Bt = Bt + b(i) * T ^ -i
Next i
VH2O = 1 / rhoi * 1000 * (mNaCl * xNaCl + mH2O * (1 - xNaCl))
VH2O = (VH2O - 30.8 * xNaCl) / (1 - xNaCl)

LogSolub = At + Bt * Log10(18.0152 / VH2O) + 2 * Log10(1 - xNaCl)
Qtz_solubility = 10 ^ LogSolub

End Function
' Qtz solubility derivative at constant temperature
Public Function dQtzDP_const_T(xnacli#, Ti#, Pi#, rhoi#) As Double
Dim Rho2#, dP#, tmpPPPPP#
'step in pressure for (dQ/dP)T
dP = 0.5

Rho2 = Rho_Brine(xnacli, Ti, Pi + dP)
tmpPPPPP = Qtz_solubility(xnacli, Ti, rhoi)
dQtzDP_const_T = (tmpPPPPP - Qtz_solubility(xnacli, Ti, Rho2)) / (-dP) / tmpPPPPP

End Function
' Qtz solubility derivative at constant pressure
Public Function dQtzDT_const_P(xnacli#, Ti#, Pi#, rhoi#) As Double
Dim Rho2#, dT#, tmpPPPPP#
'step in pressure for (dQ/dT)P
dT = -0.5

Rho2 = Rho_Brine(xnacli, Ti + dT, Pi)
tmpPPPPP = Qtz_solubility(xnacli, Ti, rhoi)
dQtzDT_const_P = (tmpPPPPP - Qtz_solubility(xnacli, Ti + dT, Rho2)) / (dT) / tmpPPPPP

End Function


Public Function mnslb#(Min_ID$, x_mol_frac#, T_in_C#, Rho_br#)
'Brooks and Steele-MacInnis 2019 here and in following functions
'Min_ID must be one of these: Ap Calc Cor Fl Qtz Ru
Dim K1#, K2#, RhoSt#, r#, T_K#, K_coef1#(), K_coef2#(), Tmp#(), q#, j#
Dim V_Mix#, m_NaCl#

T_K = T_in_C + 273.15
r = 8.3144598
K_coef1 = Coef_per_Min(Min_ID, True)
K_coef2 = Coef_per_Min(Min_ID, False)

Tmp = Qu_Jey(Min_ID)
j = Tmp(0)
q = Tmp(1)

V_Mix = 1 / Rho_br * 1000 * (58.4428 * x_mol_frac + 18.015268 * (1 - x_mol_frac))
RhoSt = 18.0152 / ((1 - x_mol_frac) * V_Mix)

K1 = -1 / r * (K_coef1(0) / T_K + K_coef1(1) + K_coef1(2) * LogExp(T_K) + K_coef1(3) * T_K + K_coef1(4) * LogExp(RhoSt))
K1 = Exp(K1)

If Min_ID <> "qtz" Then
    K2 = -1 / r * (K_coef2(0) / T_K + K_coef2(1) + K_coef2(4) * LogExp(RhoSt))
    K2 = Exp(K2)
    m_NaCl = S_Unit_Converter(x_mol_frac * 100, "MolPer", "Molal")
Else
    K2 = 0
End If

mnslb = K1 * (1 - x_mol_frac) ^ j + K1 * K2 * (1 - x_mol_frac) ^ j * m_NaCl ^ q

End Function

' first value represents J, second - Q
Private Function Qu_Jey(Min_ID$) As Double()
Dim Q_J(1) As Double

Select Case Min_ID
Case "ap"
    Q_J(0) = 6: Q_J(1) = 1
Case "calc"
    Q_J(0) = 1: Q_J(1) = 1
Case "cor"
    Q_J(0) = 0.5: Q_J(1) = 0.5
Case "fl"
    Q_J(0) = 0: Q_J(1) = 1
Case "qtz"
    Q_J(0) = 4: Q_J(1) = 0
Case "ru"
    Q_J(0) = 2: Q_J(1) = 1
End Select

Qu_Jey = Q_J

End Function


Private Function Coef_per_Min(Min_ID$, First_Eq As Boolean) As Double()
'a-e coefficients from returned as an array indexed from 0 to 4.
'Boolean switcher activates either K1 or K2 coefficients


Dim coeff#(4)

If First_Eq Then
    Select Case Min_ID
    Case "ap"
        coeff(0) = 63400:  coeff(1) = 3.9:      coeff(2) = 0:       coeff(3) = 0:       coeff(4) = -89.32
    Case "calc"
        coeff(0) = 57400:  coeff(1) = -35.71:   coeff(2) = 0:       coeff(3) = 0:       coeff(4) = -72.98
    Case "cor"
        coeff(0) = 80300:  coeff(1) = -29.31:   coeff(2) = 0:       coeff(3) = 0:       coeff(4) = -37.01
    Case "fl"
        coeff(0) = 56400:  coeff(1) = -24.89:   coeff(2) = 0:       coeff(3) = 0:       coeff(4) = -59.73
    Case "qtz"
        coeff(0) = 23600:  coeff(1) = -52.92:   coeff(2) = 10.93:   coeff(3) = -0.0463: coeff(4) = -18.52
    Case "ru"
        coeff(0) = 104000: coeff(1) = -33.72:   coeff(2) = 0:       coeff(3) = 0:       coeff(4) = -13.32
    End Select
Else
    Select Case Min_ID
    Case "ap"
        coeff(0) = 41000:   coeff(1) = -45.95: coeff(4) = 25.24
    Case "calc"
        coeff(0) = 13300:   coeff(1) = -6.07: coeff(4) = 57.29
    Case "cor"
        coeff(0) = 0:       coeff(1) = -3.31: coeff(4) = 54.61
    Case "fl"
        coeff(0) = 0:       coeff(1) = 2.62: coeff(4) = 56.16
    Case "qtz"
        coeff(0) = 0:       coeff(1) = 0:      coeff(4) = 0
    Case "ru"
        coeff(0) = -54600:  coeff(1) = 37.08: coeff(4) = 31.09
    End Select
End If

Coef_per_Min = coeff

End Function



