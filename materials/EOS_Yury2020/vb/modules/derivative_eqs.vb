' reduced susceptibility parameter to estimate critical region
Public Function XsiCapital(bbeta#, xnacli#, Ti#, Pi#, rhoi#) As Double
Dim CP_array(0, 0 To 5) As Double, PCrit#, RhoCrit#, xNaCl#, mH2O#, mNaCl#, SltWtFrac#

mH2O = 18.015268
mNaCl = 58.4428

xNaCl = xnacli
SltWtFrac = MolToWt(xNaCl * 100) / 100

Crit_Array_Filling SltWtFrac, 0, CP_array
PCrit = CP_array(0, 2)
RhoCrit = CP_array(0, 3)

XsiCapital = rhoi ^ 2 * bbeta * PCrit / RhoCrit ^ 2

End Function

'isotermal compressibility beta
Public Function Beta(xnacli#, Ti#, Pi#, rhoi#) As Double
Dim Delta#, p2#, Rho2#, point2#

Delta = 0.1  '0.5 '
p2 = Pi + Delta

If xnacli = 0 Then
    point2 = Water_Boiling_Curve(Ti)
    If Sgn(Pi - point2) <> Sgn(p2 - point2) And Sgn(rhoi - 322) = -1 Then p2 = Pi - Delta
End If

Rho2 = Rho_Brine(xnacli, Ti, p2)
Beta = 1 / rhoi * (rhoi - Rho2) / (Pi - p2)

End Function
'isobaric expansivity alpha
Public Function Alpha(xnacli#, Ti#, Pi#, rhoi#) As Double
Dim Delta#, p2#, t2#, point1#, point2#

Delta = -0.5
t2 = Ti + Delta
If xnacli = 0 Then
    point1 = Water_Boiling_Curve(Ti)
    point2 = Water_Boiling_Curve(t2)
    If Sgn(Pi - point1) <> Sgn(Pi - point2) And Sgn(Pi - p1) <> 0 Then t2 = Ti - Delta
End If

Rho2 = Rho_Brine(xnacli, t2, Pi)
Alpha = -1 / rhoi * (rhoi - Rho2) / (Ti - t2)

End Function
