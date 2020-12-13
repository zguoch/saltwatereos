
Public Function New_LV_P#(Ti#, Pi#, PressureUpperLimit#, XNaCl_Moli#, FindWhereItIs)
'salinity in mole fractions, temperature in C, pressure in bar

' This whole function designed to verify that PTx point of interest is in one phase field
' OR to identify (boolean FindWhereItIs) at which pressure Tx line is in one phase
' At very low salinities may misbehave as it not designed/tested to operate at such conditions

' if the returned results from this function is 5000 - then there is no single phase fluid for a given TX

Dim n%, a#, b#, C#, Toler#, chckr1 As Boolean, chckr2 As Boolean

If PressureUpperLimit = 0 Then PressureUpperLimit = 5000
If High_Accuracy Then Toler = 1E-10 Else Toler = 0.001

If XNaCl_Moli = 0 Then
    If Ti <= 373.946 Then New_LV_P = P_H2O_Boiling_Curve(Ti) Else New_LV_P = 0
ElseIf XNaCl_Moli = 1 Then
    If Ti <= 800.7 Then
        New_LV_P = P_Subl(Ti)
        Exit Function
    ElseIf Ti > 800.7 Then
        a = P_Boil(Ti)
        b = 5000
        If Ti > T_hm(5000) Then
            New_LV_P = a
            Exit Function
        Else
            ' This equation is solved for pressure halite melting eq.1 from Driesner & Heinrich 2007
            ' described here in Driesners_eqs, T_hm function
            a = (Ti - 800.7) / (0.024726) - 0.0005
            
            If FindWhereItIs Then
                New_LV_P = a
            Else
                If Pi < b And Pi > a Then New_LV_P = a Else New_LV_P = 5000
            End If
            Exit Function
        End If
    End If
Else
'defining upper P limit along L+H curve
    b = X_L_Sat(Ti, 5000)
    a = X_L_Sat(Ti, P_VLH(Ti))
    If a > 1 Then a = 1
    If a > b Then
        aa = a
        a = b
        b = aa
    End If

    If XNaCl_Moli >= 0.0989944153489242 And a <= XNaCl_Moli And b >= XNaCl_Moli Then
        If Ti < 800.7 Then a = P_VLH(Ti) Else a = Toler
        b = PressureUpperLimit
        n = 1
        While n <= 1000
            C = (a + b) / 2
            If (b - a) / 2 < Toler Then
                n = 1001
            End If

            n = n + 1
            aa = X_L_Sat(Ti, C) - XNaCl_Moli
            bb = X_L_Sat(Ti, b) - XNaCl_Moli
            If Sgn(aa) = Sgn(bb) Then b = C Else a = C
        Wend

    Else
        b = 5000
    End If
    'defining lower P limit
    If FindWhereItIs Then a = Toler Else a = Pi
    
    n = 1
    While n <= 1000
        C = (a + b) / 2
        If (b - a) / 2 < Toler Then
            New_LV_P = b
            n = 1001
        End If
        n = n + 1
        aa = Single_Phase_Checker(Ti, C, XNaCl_Moli, chckr1)
        bb = Single_Phase_Checker(Ti, a, XNaCl_Moli, chckr2)
        
        If chckr1 = chckr2 And chckr1 = True Then
            If aa = True Then b = C Else a = C
        ElseIf chckr1 <> chckr2 And chckr1 = False Then
            b = C
        End If
    Wend
    
    If Not High_Accuracy Then New_LV_P = Round_Down(New_LV_P, 3)
End If

End Function

Public Function LH_pressure_finder#(LowP_boundary_BAR#, T_in_C#, xNaCl_frac#)
Dim n%, Toler#, aa As Boolean, bb As Boolean, chckr1 As Boolean, chckr2 As Boolean
Dim a#, b#, C#, tmp1#, Tmp2#, FlipVal#
a = LowP_boundary_BAR
n = 1
b = 5000

tmp1 = X_L_Sat(T_in_C, a)
Tmp2 = X_L_Sat(T_in_C, b)

If tmp1 > Tmp2 Then
    FlipVal = tmp1
    tmp1 = Tmp2
    Tmp2 = FlipVal
End If

If High_Accuracy Then Toler = 1E-10 Else Toler = 0.001

While n <= 1000
    C = (a + b) / 2
    
    tmp1 = X_L_Sat(T_in_C, C) - xNaCl_frac
    Tmp2 = X_L_Sat(T_in_C, a) - xNaCl_frac
    
    If Sgn(tmp1) = Sgn(Tmp2) Then a = C Else b = C
    
    If (b - a) / 2 < Toler Then
        LH_pressure_finder = b
        n = 1001
    End If
    n = n + 1

Wend
LH_pressure_finder = Round(LH_pressure_finder, 3)
End Function
Public Function Single_Phase_Checker(Ti#, Pi#, XNaCl_Moli#, NoErrors As Boolean) As Boolean
'Role of this function to check are there any phase boundaries above given PTx or no
'T- C, P- bar, x- mol frac. NoErrors to check if screw up or no

Dim n%, j%, T#, P#, xNaCl#, XNaCl_Crit#, PCrit#, a#, b#, C#, Toler#, Tmp#, Tmp2#
Dim TLSat#, TVSat#, TVLliq#, TVLvap#
T = Ti
If Pi < 0 Then P = 1E-12 Else P = Pi
xNaCl = XNaCl_Moli
X_and_P_crit T, XNaCl_Crit, PCrit
If High_Accuracy Then Toler = 1E-10 Else Toler = 0.001
If xNaCl = 0 Then
    If P <> P_H2O_Boiling_Curve(T) Then Single_Phase_Checker = True: NoErrors = True
ElseIf xNaCl = 1 Then
    'modified to check for the phase boundary conditions only
    If T > 800.7 Then Tmp = P_Boil(T): Tmp2 = (T - 800.7) / 2.4726 / 10 ^ -2 + 0.0005 Else Tmp = P_Subl(T)
    If P <> Tmp Or P <> Tmp2 Then
        Single_Phase_Checker = True
        If P <> 0 Then NoErrors = True Else Exit Function
    End If
Else
    Tmp = X_L_Sat(Ti, Pi)
    If Tmp >= xNaCl Then
        NoErrors = True
        Single_Phase_Checker = True
    Else
        NoErrors = True
        Single_Phase_Checker = False
        If P >= P_VLH(T) Then Exit Function
    End If

    If P >= PCrit Then
        NoErrors = True
        Single_Phase_Checker = True
    Else
        If T <= 373.946 Then
            Tmp = Water_Boiling_Curve(T)
            If P > Tmp Then
                NoErrors = True
                Single_Phase_Checker = True
                Exit Function
            Else
                Tmp2 = P_VLH(T)
                If P < Tmp2 Then
                    If xNaCl <= X_V_Sat(T, P) Then
                        NoErrors = True
                        Single_Phase_Checker = True
                    Else
                        NoErrors = True
                        Single_Phase_Checker = False
                    End If
                ElseIf xNaCl <= Abs(X_VL_Vap(T, P)) Or xNaCl >= Abs(X_VL_Liq(T, P)) Then
                    NoErrors = True
                    Single_Phase_Checker = True
                Else
                    NoErrors = True
                    Single_Phase_Checker = False
                End If
            End If
        Else
            If T <= 800.7 Then Tmp = P_VLH(T) Else Tmp = -1
            If P <= Tmp Then
                If xNaCl <= X_V_Sat(T, P) Then
                    NoErrors = True
                    Single_Phase_Checker = True
                Else
                    NoErrors = True
                    Single_Phase_Checker = False
                End If
            Else
                If xNaCl <= X_VL_Vap(T, P) Or xNaCl >= X_VL_Liq(T, P) Then
                    NoErrors = True
                    Single_Phase_Checker = True
                Else
                    NoErrors = True
                    Single_Phase_Checker = False
                End If
            End If
        End If
        
    End If
    
End If
End Function



Public Function Rho_Brine#(xNaCl_frac#, T_in_C#, P_in_Bar#)
Dim T#, P#, T_Star#, V_water#
SetupExcelForCalc True
Dim mH2O#, mNaCl#
mH2O = 18.015268
mNaCl = 58.4428

T = T_in_C
P = P_in_Bar
T_Star = T_Star_V(xNaCl_frac, T, P)
V_water = V_Extrapol(xNaCl_frac, T, P)
'V_water = 0
If V_water = 0 Then V_water = mH2O / Rho_Water(T_Star, P) * 1000#

Rho_Brine = (mH2O * (1 - xNaCl_frac) + mNaCl * xNaCl_frac) / V_water * 1000#
' FOR CRITICAL POINT OF PURE WATER 373.946 C, 220.64 B RHO RETURNED AS 316.29
' Opposed to expected 322
' Problem here is that approaching critical density from PT  is very sensitive
' to decimals, and MS VBA is not precise enough
End Function


Public Function Isob_Heat_cap#(xNaCl#, T_in_C#, P_in_Bar#)
Dim ThSt#, q2#, RhoW#
ThSt = T_star_H(xNaCl, q2, T_in_C, P_in_Bar)
RhoW = Rho_Water(ThSt - 273.15, P_in_Bar)
'/(1+xnacl) correction factor doesn't mention in the original GCA Driesner's paper
'but it exist, according to calcs completed by SoWat and personal communications
Isob_Heat_cap = (Water_Isobaric_Heat_capacity_calc(ThSt, RhoW) * q2) / (1 + xNaCl)

End Function

Public Function WtToMol#(wtPecent#)
    ' this and next functions are designed to convert wt to mol % and back.
    ' these simple equations works faster compared to complicated S_Unit_Converter function
    Dim mH2O#, mNaCl#, Tmp#
    mH2O = 18.015268
    mNaCl = 58.4428
    Tmp = wtPecent / 100
    WtToMol = Tmp / mNaCl / (Tmp / mNaCl + (1 - Tmp) / mH2O) * 100
End Function
Public Function MolToWt#(MolPercent#)
    Dim mH2O#, mNaCl#, Tmp#
    mH2O = 18.015268
    mNaCl = 58.4428
    Tmp = MolPercent / 100
    MolToWt = mNaCl * Tmp / (mNaCl * Tmp + (1 - Tmp) * mH2O) * 100
End Function

Public Function SuplFuncs_LV_P_ActualFinder#(T#, x_wt#, VaporSide As Boolean, LowerEnd As Boolean, UnderBeak As Boolean)

'This function designed to identify low and high end of pressure for a given Tx
Dim EndLoop As Boolean, i%
Dim p1#, pTmp#, p2#, x#, xTmp#, x1#, x2#, P_CP#, x_CP#, vlh#, CaseID%
x = WtToMol(x_wt) / 100
If VaporSide Then
    If T <= 800.7 Then vlh = P_VLH(T)
    If LowerEnd Then
        CaseID = 1
        p1 = SuplFuncs_LV_PMinFinder(T)
        If T > 800.7 Then p2 = P_Boil(T) Else p2 = P_Subl(T)
    Else
        If T > 373.946 Then
            If T > 376.17893 Then
                UnderBeak = False
                X_and_P_crit T, x_CP, p1
                p2 = SuplFuncs_LV_PMinFinder(T)
                CaseID = 2
            Else
                X_and_P_crit T, x_CP, p1
                p2 = SuplFuncs_LV_PMaxFinder(T)
                CaseID = 3
            End If
        Else
            If UnderBeak Then
                p1 = Water_Boiling_Curve(T)
                p2 = SuplFuncs_LV_PMaxFinder(T)
                CaseID = 3
            Else
                p1 = SuplFuncs_LV_PMaxFinder(T)
                p2 = SuplFuncs_LV_PMinFinder(T)
                CaseID = 4
            End If
        End If
    End If
Else
    CaseID = 5
    If T > 373.946 Then X_and_P_crit T, pTmp, p1 Else p1 = Water_Boiling_Curve(T)
    If T < 800.7 Then p2 = vlh Else p2 = P_Boil(T)
End If

Select Case CaseID
Case 1, 3
    If p1 < vlh Then x1 = X_V_Sat(T, p1) Else x1 = X_VL_Vap(T, p1)
    If p2 < vlh Then x2 = X_V_Sat(T, p2) Else x2 = X_VL_Vap(T, p2)
    
    If x < x1 Or x > x2 Then
        Exit Function
    Else
        pTmp = (p1 + p2) / 2
        While Not EndLoop
            If pTmp < vlh Then xTmp = X_V_Sat(T, pTmp) - x Else xTmp = X_VL_Vap(T, pTmp) - x
            If p2 < vlh Then x2 = X_V_Sat(T, p2) - x Else x2 = X_VL_Vap(T, p2) - x
            If Sgn(xTmp) <> Sgn(x2) Then p1 = pTmp Else p2 = pTmp
            pTmp = (p1 + p2) / 2
            If Abs(xTmp / x * 100) < 1E-05 Or i = 200 Then EndLoop = True
            i = i + 1
        Wend
    End If
    
Case 2, 4
    If p1 < vlh Then x1 = X_V_Sat(T, p1) Else x1 = X_VL_Vap(T, p1)
    If p2 < vlh Then x2 = X_V_Sat(T, p2) Else x2 = X_VL_Vap(T, p2)
    
    If Not VaporSide Then
        If x > x_CP Or x < x2 Then
            Exit Function
        End If
    ElseIf x < x2 Or x > x1 Then
        Exit Function
    Else
        pTmp = (p1 + p2) / 2
        While Not EndLoop
            If pTmp < vlh Then xTmp = X_V_Sat(T, pTmp) - x Else xTmp = X_VL_Vap(T, pTmp) - x
            If p2 < vlh Then x2 = X_V_Sat(T, p2) - x Else x2 = X_VL_Vap(T, p2) - x
            
            If Sgn(xTmp) <> Sgn(x2) Then p1 = pTmp Else p2 = pTmp
            pTmp = (p1 + p2) / 2
            If Abs(xTmp / x * 100) < 1E-05 Or i = 200 Then EndLoop = True
            i = i + 1
        Wend
    End If
Case 5
     x1 = X_VL_Liq(T, p1)
     x2 = X_VL_Liq(T, p2)
    
    If x < x1 Or x > x2 Then
        Exit Function
    Else
        pTmp = (p1 + p2) / 2
        While Not EndLoop
            xTmp = X_VL_Liq(T, pTmp) - x
            x2 = X_VL_Liq(T, p2) - x
            If Sgn(xTmp) <> Sgn(x2) Then p1 = pTmp Else p2 = pTmp
            pTmp = (p1 + p2) / 2
            If Abs(xTmp / x * 100) < 1E-05 Or i = 200 Then EndLoop = True
            i = i + 1
        Wend
    End If
Case Else
    SuplFuncs_LV_P_ActualFinder = 0
    Exit Function
End Select

SuplFuncs_LV_P_ActualFinder = pTmp

End Function

Public Function SuplFuncs_LV_PMinFinder#(T#)
Dim EndLoop As Boolean, i%
Dim p1#, p2#, x1#, x2#, p_i#, x_i#, p_lvh#, slp1#, slp2#

If T <= 800.7 Then
    p1 = P_VLH(T)
    p_lvh = p1
Else
    X_and_P_crit T, x1, p1
    p1 = p1 / 2
End If
If T > 800.7 Then p2 = P_Boil(T) Else p2 = P_Subl(T)
p_i = (p1 + p2) / 2

If p1 > p_lvh Then x1 = X_VL_Vap(T, p1) Else x1 = X_V_Sat(T, p1)
If p2 > p_lvh Then x2 = X_VL_Vap(T, p2) Else x2 = X_V_Sat(T, p2)
If p_i > p_lvh Then x_i = X_VL_Vap(T, p_i) Else x_i = X_V_Sat(T, p_i)

While Not EndLoop
    slp1 = (p1 - p_i) / (x1 - x_i)
    slp2 = (p2 - p_i) / (x2 - x_i)
    If Sgn(slp1) <> Sgn(slp2) Then
        If slp1 >= 0 Then
            p1 = p_i
        Else
            If x1 < x2 Then p2 = p_i Else p1 = p_i
        End If
    ElseIf slp1 >= 0 Then
        p1 = p_i
    Else
        p2 = p_i
    End If
    p_i = (p1 + p2) / 2
    
    If p1 > p_lvh Then x1 = X_VL_Vap(T, p1) Else x1 = X_V_Sat(T, p1)
    If p2 > p_lvh Then x2 = X_VL_Vap(T, p2) Else x2 = X_V_Sat(T, p2)
    If p_i > p_lvh Then x_i = X_VL_Vap(T, p_i) Else x_i = X_V_Sat(T, p_i)

    
    If Abs(p1 - p2) / p1 < 1E-06 Or i > 200 Then EndLoop = True
    i = i + 1
Wend

SuplFuncs_LV_PMinFinder = (p2 + p1) / 2
End Function

Public Function SuplFuncs_LV_PMaxFinder#(T#)
Dim EndLoop As Boolean, i%
Dim p1#, p2#, x1#, x2#, p_i#, x_i#, p_lvh#, slp1#, slp2#

If T > 376.17893 Then   'Critical point here replaced by this empirical number,
                        'as this T appears to be an end of the bird's beak
    X_and_P_crit T, Xi, p_i
    SuplFuncs_LV_PMaxFinder = p_i
Else
    If T >= 373.946 Then X_and_P_crit T, Xi, p1 Else p1 = Water_Boiling_Curve(T)
    p2 = SuplFuncs_LV_PMinFinder(T)
    p_i = (p1 + p2) / 2
    
    If p1 > p_lvh Then x1 = X_VL_Vap(T, p1) Else x1 = X_V_Sat(T, p1)
    If p2 > p_lvh Then x2 = X_VL_Vap(T, p2) Else x2 = X_V_Sat(T, p2)
    If p_i > p_lvh Then x_i = X_VL_Vap(T, p_i) Else x_i = X_V_Sat(T, p_i)
    
    While Not EndLoop
        slp1 = (p1 - p_i) / (x1 - x_i)
        slp2 = (p2 - p_i) / (x2 - x_i)
        If Sgn(slp2) <= 0 And x1 > x2 Then
            p2 = p_i
        ElseIf Sgn(slp2) <= 0 Then
            p1 = p_i
        Else
            p2 = p_i
        End If
        
        p_i = (p1 + p2) / 2
        
        If p1 > p_lvh Then x1 = X_VL_Vap(T, p1) Else x1 = X_V_Sat(T, p1)
        If p2 > p_lvh Then x2 = X_VL_Vap(T, p2) Else x2 = X_V_Sat(T, p2)
        If p_i > p_lvh Then x_i = X_VL_Vap(T, p_i) Else x_i = X_V_Sat(T, p_i)
    
        
        If Abs(p1 - p2) < 0.0001 Or i > 200 Then EndLoop = True
        i = i + 1
    Wend
    
    SuplFuncs_LV_PMaxFinder = (p2 + p1) / 2
End If
End Function
