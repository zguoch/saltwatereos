''' functions in this module can be used from Excel cells directly (see sheet Microthermometry as an example)

Public Function Phases_PTx_Surveyor(Optional WtPercent As Variant = "none", Optional Temp_in_C As Variant = "none", Optional Pres_in_Bar As Variant = "none")
' this function focused on identifying phase equilibria at any given PTX/PT/PX/Tx coordinates
Dim v_name(), LwLim(), ULim(), Vl(), msg$
Dim x#, T#, P#, i%, j%, PropInUse As Byte
Dim x_test#, T_Test#, P_Test#, Temp1#, Temp2#

If IsEmpty(WtPercent) Then WtPercent = "text"
If IsEmpty(Temp_in_C) Then Temp_in_C = "text"
If IsEmpty(Pres_in_Bar) Then Pres_in_Bar = "text"

v_name = Array("x", "T", "P")
LwLim = Array(0, 100, 0)

ULim = Array(100, 1000, 5000)

'PropInUse used as indicator which case will be followed:
    '0 for single PTx point
    '1 variable salinity. PT case.
If IsNumeric(WtPercent) And Not IsEmpty(WtPercent) Then i = i + 1: x = WtPercent Else x = 0: PropInUse = 1
    '2 variable temperature, Px case.
If IsNumeric(Temp_in_C) And Not IsEmpty(Temp_in_C) Then i = i + 1: T = Temp_in_C Else T = 100: PropInUse = 2
    '3 variable pressure, Tx case.
If IsNumeric(Pres_in_Bar) And Not IsEmpty(Pres_in_Bar) Then i = i + 1: P = Pres_in_Bar Else P = 5000: PropInUse = 3


If i < 2 Then
    Phases_PTx_Surveyor = "Function requires at least 2 input arguments"
    Exit Function
End If

Vl = Array(x, T, P)
msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    Phases_PTx_Surveyor = msg
    Exit Function
End If

Select Case PropInUse
Case 0
    'PTx single point
    ' Only pressure for corresponding TX will be returned
    
    'Pure H2O
    If WtPercent = 0 Then
        If T < 373.946 Then
            P_Test = Round(P_H2O_Boiling_Curve(T), 3)
            If Round(P, 3) = P_Test Then
                Phases_PTx_Surveyor = "Water boiling curve"
            ElseIf Round(P, 3) > P_Test Then
                Phases_PTx_Surveyor = "Single phase Liquid state, boiling pressure " & SuplFuncs_Report_Number(P_Test) & "."
            Else
                Phases_PTx_Surveyor = "Single phase Vapor state, boiling pressure " & SuplFuncs_Report_Number(P_Test) & "."
            End If
        Else
            Temp1 = Rho_Water(T, P)
            If Temp1 > 322 Then Phases_PTx_Surveyor = "Single phase Liquid state." Else Phases_PTx_Surveyor = "Single phase Vapor state."
        End If
        Exit Function
        
    End If
    'pure NaCl
    If WtPercent = 100 Then
        If T < 800.7 Then
            P_Test = P_Subl(T)
            If P > P_Test Then Phases_PTx_Surveyor = "NaCl is in Solid phase. " Else Phases_PTx_Surveyor = "NaCl in Vapor phase. "
            Phases_PTx_Surveyor = Phases_PTx_Surveyor & SuplFuncs_Report_Number(P_Test) & " bar is the sublimation curve."
            Exit Function
        Else
            P_Test = P_Boil(T)
            Temp1 = (T - 800.7) / 2.4726 / 10 ^ -2 + 0.0005
            
            If P > Temp1 Then
                Phases_PTx_Surveyor = "NaCl is in Solid phase. "
            ElseIf P < P_Test Then
                Phases_PTx_Surveyor = "NaCl in Vapor phase. "
            Else
                Phases_PTx_Surveyor = "NaCl is in Liquid phase. "
            End If
            Phases_PTx_Surveyor = Phases_PTx_Surveyor & SuplFuncs_Report_Number(P_Test) & " bar NaCl boiling, " & SuplFuncs_Report_Number(Temp1) & " bar NaCl melting."
            Exit Function
        End If
    End If
    
    If T < 800.7 Then
        P_Test = P_VLH(T)
        If P > P_Test Then
            x_test = MolToWt(100 * X_L_Sat(T, P))
            X_and_P_crit T, Temp1, Temp2
            Temp1 = MolToWt(100 * Temp1)
            If x > x_test Then
                Phases_PTx_Surveyor = "The fluid is in L+H state."
                Exit Function
            ElseIf P < Temp2 Then
                Temp1 = MolToWt(100 * X_VL_Vap(T, P))
                Temp2 = MolToWt(100 * X_VL_Liq(T, P))
                If x < Temp1 Or x > Temp2 Then
                    Phases_PTx_Surveyor = "The fluid is in single phase state."
                    Exit Function
                Else
                    Phases_PTx_Surveyor = "The fluid is in L+V state."
                    Exit Function
                End If
            Else
                Phases_PTx_Surveyor = "The fluid is in single phase state."
                Exit Function
            End If
        Else
            Temp1 = MolToWt(100 * X_VL_Vap(T, P))
            If x < Temp1 Then
                Phases_PTx_Surveyor = "The fluid is in single phase state."
                Exit Function
            Else
                Phases_PTx_Surveyor = "The fluid is in V+H state."
                Exit Function
            End If
        End If
    Else
        x_test = MolToWt(100 * X_L_Sat(T, P))
        If x > x_test Then
            Phases_PTx_Surveyor = "The fluid is in L+H state."
            Exit Function
        End If
        
        P_Test = SuplFuncs_LV_P_ActualFinder(T, x, True, True, False)
        Temp1 = MolToWt(100 * X_VL_Vap(T, P))
        Temp2 = MolToWt(100 * X_VL_Liq(T, P))
        If P < P_Test Then
            Phases_PTx_Surveyor = "The fluid is in single phase state."
            Exit Function
        End If

        If x >= Temp1 And x <= Temp2 Then
            Phases_PTx_Surveyor = "The fluid is in L+V state."
            Exit Function
        Else
            Phases_PTx_Surveyor = "The fluid is in single phase state."
            Exit Function
        End If
    End If
    
Case 1
    '1 Variable salinity
    X_and_P_crit T, x_test, P_Test
    x_test = MolToWt(100 * x_test)
    Temp1 = P_VLH(T)
    Temp2 = MolToWt(100 * X_L_Sat(T, P))
    
    If P > P_Test Then
        If Temp2 >= 100 Then
            Phases_PTx_Surveyor = "The fluid is in single phase between 0-100 NaCl."
        Else
            Phases_PTx_Surveyor = "The fluid is in single phase up to " & SuplFuncs_Report_Number(Temp2) & " wt% NaCl. Above that is in two phase L+H."
        End If
    Else
        If T < 800.7 Then
            If P > Temp1 Then
                If Temp2 < 100 Then Phases_PTx_Surveyor = SuplFuncs_Report_Number(Temp2)
                Temp1 = MolToWt(100 * X_VL_Vap(T, P))
                Temp2 = MolToWt(100 * X_VL_Liq(T, P))
                If Phases_PTx_Surveyor <> "" Then
                    Phases_PTx_Surveyor = "Below " & SuplFuncs_Report_Number(Temp1) & " and between " & SuplFuncs_Report_Number(Temp2) & " to " & Phases_PTx_Surveyor & _
                        " wt. % NaCl the fluid is in single phase. From " & SuplFuncs_Report_Number(Temp1) & " to " & SuplFuncs_Report_Number(Temp2) & " is L+V, above " & _
                        Phases_PTx_Surveyor & " is L+H."
                Else
                    Phases_PTx_Surveyor = "Below " & SuplFuncs_Report_Number(Temp1) & " and above " & SuplFuncs_Report_Number(Temp2) & " wt. % NaCl the fluid is in single phase, between - L+V."
                End If
            Else
                Temp1 = MolToWt(100 * X_VL_Vap(T, P))
                Phases_PTx_Surveyor = "Below " & SuplFuncs_Report_Number(Temp1) & " wt. % NaCl the fluid is in single (vapor) phase, above that - V+H."
            End If
        Else
            If Temp2 < 100 Then Phases_PTx_Surveyor = CStr(SuplFuncs_Report_Number(Temp2))
            Temp1 = MolToWt(100 * X_VL_Vap(T, P))
            Temp2 = MolToWt(100 * X_VL_Liq(T, P))
            If Phases_PTx_Surveyor <> "" Then
                Phases_PTx_Surveyor = "Below " & SuplFuncs_Report_Number(Temp1) & " and between " & SuplFuncs_Report_Number(Temp2) & " to " & Phases_PTx_Surveyor & _
                    " wt. % NaCl the fluid is in single phase. From " & SuplFuncs_Report_Number(Temp1) & " to " & SuplFuncs_Report_Number(Temp2) & _
                    " is L+V, above " & Phases_PTx_Surveyor & " is L+H."
            Else
                Phases_PTx_Surveyor = "Below " & SuplFuncs_Report_Number(Temp1) & " and above " & SuplFuncs_Report_Number(Temp2) & _
                " wt. % NaCl the fluid is in single phase, between - L+V."
            End If
        End If
    End If

Case 2
    '2 variable T
    If x > 26.2 Then
        'l+h
        Temp2 = SuplFuncs_LH_TFinder(P, x)
        Phases_PTx_Surveyor = "L+H surface crossed at " & SuplFuncs_Report_Number(Temp2) & " C. "
    End If

    Temp1 = SuplFuncs_LV_TFinder(P, x)
    If Temp1 > Temp2 Then
        If Round(Temp1, 4) <> 1000 Then
            Phases_PTx_Surveyor = Phases_PTx_Surveyor & "L+V surface entered at " & SuplFuncs_Report_Number(Temp1) & " C, "
            X_and_P_crit SuplFuncs_CP_TFinder(P), Temp1, Temp2
            Temp1 = MolToWt(100 * Temp1)
            If x > Temp1 Then
                Phases_PTx_Surveyor = Phases_PTx_Surveyor & "below fluid is in Liquid state. "
            Else
                Phases_PTx_Surveyor = Phases_PTx_Surveyor & "below fluid is in Vapor state. "
            End If
        ElseIf Temp2 = 0 And P > 390.1474443 Then
            Phases_PTx_Surveyor = "The fluid is in single phase between 0-1000 C."
            Exit Function
        End If
    End If
    
    'T 594.6324458 at max P_LVH=390.1474443
    '26.278 is a lowest LH salinity, at 0 C
    If P <= 390.1474443 Then
        Temp1 = SuplFuncs_VH_TFinder(P, False)
        Temp2 = SuplFuncs_VH_TFinder(P, True)
        Phases_PTx_Surveyor = Phases_PTx_Surveyor & " V+H crossed at " & SuplFuncs_Report_Number(Temp1) & " and " & SuplFuncs_Report_Number(Temp2) & " C."
    End If

Case 3
    '3 variable P
    
    X_and_P_crit T, x_test, P_Test
    x_test = MolToWt(100 * x_test)
    If x < x_test Then
    ''Case for vapor branch, it will work only at T above CP of H2O
        P_Test = SuplFuncs_LV_PMinFinder(T)
        Temp1 = MolToWt(100 * X_VL_Vap(T, P_Test))
        'check for salinity too low to have 2 phase
        If Temp1 > x Then
            Phases_PTx_Surveyor = "The fluid is in single phase state at any pressure."
            Exit Function
        Else
            Temp2 = SuplFuncs_LV_P_ActualFinder(T, x, True, False, False)
            If T < 800.7 Then
                Temp1 = P_VLH(T)
                If Temp1 > Temp2 Then
                    Phases_PTx_Surveyor = "Between " & SuplFuncs_Report_Number(Temp2) & " and " & SuplFuncs_Report_Number(SuplFuncs_LV_P_ActualFinder(T, x, True, True, False)) _
                        & " bar fluid is V+H, above fluid is in single phase."
                Else
                    Phases_PTx_Surveyor = "Between " & SuplFuncs_Report_Number(Temp2) & " and " & SuplFuncs_Report_Number(Temp1) & _
                        " bar fluid is in L+V, between " & Round(P_VLH(T), 1) & " and " & SuplFuncs_Report_Number(SuplFuncs_LV_P_ActualFinder(T, x, True, True, False)) _
                        & " it is in V+H."
                End If
            Else
                Temp1 = P_Boil(T)
                Phases_PTx_Surveyor = "Between " & SuplFuncs_Report_Number(Temp2) & " and " & SuplFuncs_Report_Number(Temp1) & " bar fluid is in L+V."
            End If
            Exit Function
        End If
    Else
        Temp1 = SuplFuncs_LV_PMinFinder(T)
        Temp1 = MolToWt(100 * X_VL_Vap(T, Temp1))
        
        If x < Temp1 Then
            Phases_PTx_Surveyor = "The fluid is in the single phase at any pressure, crossing water boining curve at " _
                & SuplFuncs_Report_Number(Water_Boiling_Curve(T)) & " bar."
            Exit Function
        End If
        
        Temp2 = SuplFuncs_LV_P_ActualFinder(T, x, False, False, False)
        
        If T < 800.7 Then
            Temp1 = P_VLH(T)
            If Round(Temp2, 4) >= Round(Temp1, 4) Then
                Phases_PTx_Surveyor = "Between " & SuplFuncs_Report_Number(SuplFuncs_LV_P_ActualFinder(T, x, True, True, False)) & _
                    " and " & SuplFuncs_Report_Number(Temp1) & " bar fluid is in V+H field, up to " & SuplFuncs_Report_Number(SuplFuncs_LV_P_ActualFinder(T, x, False, True, False)) & _
                    " is in L+V. "
            Else
                Phases_PTx_Surveyor = "Between " & SuplFuncs_Report_Number(SuplFuncs_LV_P_ActualFinder(T, x, True, True, False)) & _
                    " and " & SuplFuncs_Report_Number(SuplFuncs_LV_P_ActualFinder(T, x, True, False, False)) & " bar fluid is in V+H field. Up to "
            End If
        Else
            Phases_PTx_Surveyor = "Between " & SuplFuncs_Report_Number(SuplFuncs_LV_P_ActualFinder(T, x, False, False, False)) & _
                " and " & SuplFuncs_Report_Number(SuplFuncs_LV_P_ActualFinder(T, x, False, True, False)) & " bar fluid is in V+H field. "
        End If
    End If
    
    If T < 925 Then
        Temp1 = MolToWt(100 * X_L_Sat(T, P_VLH(T)))
        Temp2 = MolToWt(100 * X_L_Sat(T, 5000))
        x_test = 0
        If Temp1 > Temp2 Then
            x_test = Temp1
            Temp1 = Temp2
            Temp2 = x_test
        End If
        If x >= Temp1 And x <= Temp2 Then
            P_Test = LH_pressure_finder(P_VLH(T), T, WtToMol(x) / 100)
            If x_test = 0 Then
                If P_Test > 4999 Then
                    Phases_PTx_Surveyor = Phases_PTx_Surveyor & "Above it, the fluid is in single phase state."
                Else
                    Phases_PTx_Surveyor = Phases_PTx_Surveyor & "Up to " & SuplFuncs_Report_Number(P_Test) & " bar fluid is in L+H field, above it is in single phase state."
                End If
            Else
                If P_Test > 4999 Then
                    Phases_PTx_Surveyor = Phases_PTx_Surveyor & "Above it, the fluid is in L+H."
                Else
                    Phases_PTx_Surveyor = Phases_PTx_Surveyor & "Up to " & SuplFuncs_Report_Number(P_Test) & " bar fluid is in single phase state field, above it is L+H."
                End If
            End If
        ElseIf x > Temp2 Then
            Phases_PTx_Surveyor = Phases_PTx_Surveyor & " above it, the fluid is in L+H."
        End If
    End If

End Select

End Function
Public Function FlInc_Salinity(T_in_C_melt#, Optional ByVal T_in_C_Homog# = 0#, Optional ByVal Dissolv_Phase$ = "ice")
'governing funciton which identifies specific model needed for salinity estimation
If Dissolv_Phase = "" Then Dissolv_Phase = "ice"
Dissolv_Phase = LCase(Left(Dissolv_Phase, 2))

If Left(Dissolv_Phase, 1) = "i" Then T_in_C_Homog = 0#

If T_in_C_Homog <> 0# Then
    If T_in_C_Homog < T_in_C_melt Then
        FlInc_Salinity = FlInc_Salinity_LH(T_in_C_melt, T_in_C_Homog)
    Else
        FlInc_Salinity = FlInc_Salinity_Vap_Sat(T_in_C_melt, Dissolv_Phase)
    End If
Else
    FlInc_Salinity = FlInc_Salinity_Vap_Sat(T_in_C_melt, Dissolv_Phase)
End If

If IsEmpty(FlInc_Salinity) Then FlInc_Salinity = "Error in incoming data"


End Function
Private Function FlInc_Salinity_Vap_Sat(Temp_in_C#, Dissolv_Phase$)

Select Case Dissolv_Phase
Case "ic", "i"
    'ice melting - Bodnar 1993
    If Temp_in_C < -21.2 Then
        FlInc_Salinity_Vap_Sat = "ice melting temperature cannot be lower than H2O-NaCl eutectic point (-21.2 C)"
    ElseIf Temp_in_C > 0 Then
        FlInc_Salinity_Vap_Sat = "ice melting temperature should be negative"
    Else
        FlInc_Salinity_Vap_Sat = Round(1.78 * Abs(Temp_in_C) - 0.0442 * Abs(Temp_in_C) ^ 2 + 0.000557 * Abs(Temp_in_C) ^ 3, 2)
    End If
Case "ha", "h"
'This and following - Sterner 1988
    If Temp_in_C < 0.1 Then
        FlInc_Salinity_Vap_Sat = "halite melting temperature should be greater than 0.1"
    ElseIf Temp_in_C > 925 Then
        FlInc_Salinity_Vap_Sat = "halite melting temperature cannot exceed halite liquidus (925 C at 5000 bar)"
    Else
        Temp_in_C = Temp_in_C / 100
        FlInc_Salinity_Vap_Sat = Round(26.242 + 0.4928 * Temp_in_C + 1.42 * Temp_in_C ^ 2 - 0.223 * Temp_in_C ^ 3 + 0.04129 * Temp_in_C ^ 4 + 0.006295 * Temp_in_C ^ 5 - 0.001967 * Temp_in_C ^ 6 + 0.0001112 * Temp_in_C ^ 7, 2)
    End If
Case "hh", "hy"
    If Temp_in_C < -21.2 Then
        FlInc_Salinity_Vap_Sat = "hydrohalite melting temperature cannot be lower than H2O-NaCl eutectic point (-21.2)"
    ElseIf Temp_in_C > 0 Then
        FlInc_Salinity_Vap_Sat = "hydrohalite melting temperature should be negative"
    Else
        Temp_in_C = Temp_in_C / 100
        FlInc_Salinity_Vap_Sat = Round(40.36947594 + 14.80771966 * Temp_in_C - 14.08238722, 2)
    End If
Case Else
    FlInc_Salinity_Vap_Sat = "Dissolv_phase argument accepts only ice, halite, hydrohalite (i, h, hh) values"
End Select

End Function
Private Function FlInc_Salinity_LH(Temp_in_C_halite#, Temp_in_C_homog#)
'Lecumberri-Sanchez et al. 2012
Dim a#(0 To 1), i%
If Temp_in_C_homog > Temp_in_C_halite Then
    FlInc_Salinity_LH = "specified Thom>T_halite indicates vapor-saturated halite melting, please use function FlInc_Salinity_Vap_Sat"
    Exit Function
ElseIf Temp_in_C_halite < 0.1 Or Temp_in_C_halite# > 925 Then
    FlInc_Salinity_LH = "specified T_halite is outside of 0.1-925 C temperature range"
    Exit Function
End If

a(0) = 26.4575 - 0.000361 * Temp_in_C_homog ^ 2 + 5.5302E-07 * Temp_in_C_homog ^ 3
a(1) = 0.010765 + 0.0003697 * Temp_in_C_homog - 1.544E-07 * Temp_in_C_homog ^ 2 - 3.79E-10 * Temp_in_C_homog ^ 3

For i = 0 To 1
    FlInc_Salinity_LH = FlInc_Salinity_LH + a(i) * Temp_in_C_halite ^ i
Next i

FlInc_Salinity_LH = Round(FlInc_Salinity_LH, 1)

End Function

Public Function Phases_Single_ph_pressure_at_TX(WtPercent#, Temp_in_C#, Optional Precision As Boolean = False)

Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("T", "x")
LwLim = Array(0, 0)
ULim = Array(1000, 100)
Vl = Array(Temp_in_C, WtPercent)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)
 
If msg <> "" Then
    Phases_Single_ph_pressure_at_TX = msg
Else
    High_Accuracy = Precision
    Phases_Single_ph_pressure_at_TX = New_LV_P(Temp_in_C, 1E-05, 0, WtToMol(WtPercent) / 100, True)
End If

If Phases_Single_ph_pressure_at_TX = 5000 Then Phases_Single_ph_pressure_at_TX = "2 phases in range from 1E-5 to 5000 bar"
End Function
Public Function Phases_x_VL_coex_Liquid_side(Temp_in_C#, Pres_in_Bar#)
Dim v_name(), LwLim(), ULim(), Vl(), msg$, tmp1#, Tmp2#
v_name = Array("T", "P")
LwLim = Array(0, 0)
ULim = Array(1000, 2300)
Vl = Array(Temp_in_C, Pres_in_Bar)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    If Pres_in_Bar >= 2300 Then Phases_x_VL_coex_Liquid_side = "The fluid is in single phase" Else Phases_x_VL_coex_Liquid_side = msg
Else
    X_and_P_crit Temp_in_C, tmp1, Tmp2
    If Pres_in_Bar > Tmp2 Then
        Phases_x_VL_coex_Liquid_side = "The fluid is in single phase"
    ElseIf Pres_in_Bar < P_VLH(Temp_in_C) Then
        Phases_x_VL_coex_Liquid_side = "the fluid is below V+H surface"
    Else
        Phases_x_VL_coex_Liquid_side = MolToWt(100 * X_VL_Liq(Temp_in_C, Pres_in_Bar))
    End If
End If
End Function

Public Function Phases_x_VL_coex_Vapor_side(Temp_in_C#, Pres_in_Bar#)
Dim v_name(), LwLim(), ULim(), Vl(), msg$, tmp1#, Tmp2#
v_name = Array("T", "P")
LwLim = Array(0, 0)
ULim = Array(1000, 2300)
Vl = Array(Temp_in_C, Pres_in_Bar)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    If Pres_in_Bar >= 2300 Then Phases_x_VL_coex_Vapor_side = "The fluid is in single phase" Else Phases_x_VL_coex_Vapor_side = msg
Else
    X_and_P_crit Temp_in_C, tmp1, Tmp2
    If Pres_in_Bar > Tmp2 Then
        Phases_x_VL_coex_Vapor_side = "The fluid is in single phase"
    Else
        Phases_x_VL_coex_Vapor_side = MolToWt(100 * X_VL_Vap(Temp_in_C, Pres_in_Bar))
    End If
    
End If
End Function

Public Function Phases_P_LVH_Coexistence(Temp_in_C#)
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("T")
LwLim = Array(0)
ULim = Array(800.7)
Vl = Array(Temp_in_C)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    Phases_P_LVH_Coexistence = msg
Else
    Phases_P_LVH_Coexistence = P_VLH(Temp_in_C)
End If
End Function

Public Function Phases_x_Halite_saturated_vapor(Temp_in_C#, Pres_in_Bar#)
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("T", "P")
LwLim = Array(0, 0)
ULim = Array(802, 400)
Vl = Array(Temp_in_C, Pres_in_Bar)


msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    Phases_x_Halite_saturated_vapor = msg
Else
    Phases_x_Halite_saturated_vapor = MolToWt(100 * X_V_Sat(Temp_in_C, Pres_in_Bar))
End If
End Function

Public Function Phases_x_Halite_Liquidus(Temp_in_C#, Pres_in_Bar#)
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("T", "P")
LwLim = Array(0, 0)
ULim = Array(1000, 5000)
Vl = Array(Temp_in_C, Pres_in_Bar)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    Phases_x_Halite_Liquidus = msg
Else
    If Pres_in_Bar < P_VLH(Temp_in_C) Then
        Phases_x_Halite_Liquidus = "The fluid is below L+H fluid field"
    Else
        Phases_x_Halite_Liquidus = MolToWt(100 * X_L_Sat(Temp_in_C, Pres_in_Bar))
    End If
End If
End Function

Public Function Phases_P_crit(Temp_in_C#)
Dim tmp1#, Tmp2#
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("T")
LwLim = Array(373.946)
ULim = Array(1000)
Vl = Array(Temp_in_C)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    Phases_P_crit = msg
Else
    X_and_P_crit Temp_in_C, tmp1, Tmp2
    Phases_P_crit = Tmp2
End If
End Function

Public Function Phases_x_crit(Temp_in_C#)
Dim tmp1#, Tmp2#
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("T")
LwLim = Array(373.946)
ULim = Array(1000)
Vl = Array(Temp_in_C)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    Phases_x_crit = msg
Else
    X_and_P_crit Temp_in_C, tmp1, Tmp2
    Phases_x_crit = MolToWt(tmp1 * 100)
End If
End Function

Public Function Phases_P_Halite_Boil(Temp_in_C#)
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("T")
LwLim = Array(800.7)
ULim = Array(1000)
Vl = Array(Temp_in_C)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    Phases_P_Halite_Boil = msg
Else
    Phases_P_Halite_Boil = P_Boil(Temp_in_C)
End If
End Function

Public Function Phases_P_Halite_Subl(Temp_in_C#)
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("T")
LwLim = Array(0)
ULim = Array(800.7)
Vl = Array(Temp_in_C)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    Phases_P_Halite_Subl = msg
Else
    Phases_P_Halite_Subl = P_Subl(Temp_in_C)
End If
End Function

Public Function Phases_T_Halite_Melt(Pres_in_Bar#)
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("P")
LwLim = Array(0)
ULim = Array(5000)
Vl = Array(Pres_in_Bar)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    Phases_T_Halite_Melt = msg
Else
    Phases_T_Halite_Melt = T_hm(Pres_in_Bar)
End If
End Function
Public Function Phases_x_of_coexisting_phases(WtPercent#, Temp_in_C#, Pres_in_Bar#)
Dim FstPhase#, ScndPhase#, tmp1#, Tmp2#, n%, k%
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("x", "T", "P")
LwLim = Array(0, 0, 0)
ULim = Array(100, 1000, 5000)
Vl = Array(WtPercent, Temp_in_C, Pres_in_Bar#)
n = 1
k = 1

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)
If msg <> "" Then
    Phases_x_of_coexisting_phases = msg
Else
    If Single_Phase_Checker(Temp_in_C, Pres_in_Bar, WtToMol(WtPercent) / 100, False) Then
        Phases_x_of_coexisting_phases = " The fluid is in the single state"
    Else
        If Temp_in_C < 800.7 Then tmp1 = P_VLH(Temp_in_C)
        If Pres_in_Bar < tmp1 Then
            FstPhase = MolToWt(X_VL_Vap(Temp_in_C, Pres_in_Bar) * 100)
            ScndPhase = 100
        Else
            Tmp = MolToWt(X_L_Sat(Temp_in_C, Pres_in_Bar) * 100)
            If Tmp > WtPercent Then
                ScndPhase = MolToWt(X_VL_Liq(Temp_in_C, Pres_in_Bar) * 100)
                FstPhase = MolToWt(X_VL_Vap(Temp_in_C, Pres_in_Bar) * 100)
            Else
                ScndPhase = 100
                FstPhase = Tmp
            End If
        End If
        
        If FstPhase > 1 Then
            n = 0
        Else
            While Round(FstPhase, n) = 0
                n = n + 1
            Wend
        End If
        
        If ScndPhase > 1 Then
            k = 0
        Else
            While Round(ScndPhase, k) = 0
                k = k + 1
            Wend
        End If
        
        Phases_x_of_coexisting_phases = "First phase have " & SuplFuncs_Report_Number(FstPhase) & ", second " & SuplFuncs_Report_Number(ScndPhase) & " of wt. % NaCl salinity"
        Exit Function
    End If
End If

End Function

Public Function Brine_dP_dT_for_Isochore(WtPercent#, Temp_in_C#, Pres_in_Bar#, Dens_in_kg#)
'Works properly only above 1 bar, and note that isochores at low PT have some curvature

Dim T_wa1#, T_wa2#, T_Inc#, xNaCl#, Dens_wa#
Dim mH2O#, mNaCl#, p1#, p2#, Rho1#, Rho2#, der#, EndLoop As Boolean, i%
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("x", "T", "P")
LwLim = Array(0, 0, 1)
ULim = Array(100, 1000, 5000)
Vl = Array(WtPercent, Temp_in_C, Pres_in_Bar)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    Brine_dP_dT_for_Isochore = msg
Else
    mH2O = 18.015268
    mNaCl = 58.4428
    T_Inc = 1
    
    xNaCl = WtToMol(WtPercent) / 100
    Dens_wa = Dens_in_kg * mH2O / (mH2O * (1 - xNaCl) + mNaCl * xNaCl)
    T_wa1 = T_Star_V(xNaCl, Temp_in_C + T_Inc, Pres_in_Bar)
    
    p1 = Water_Pressure_calc(T_wa1 + 273.15, Dens_wa) * 10
    p2 = p1 + 0.1
    
    While EndLoop = False
    
        T_wa1 = T_Star_V(xNaCl, Temp_in_C + T_Inc, p1)
        T_wa2 = T_Star_V(xNaCl, Temp_in_C + T_Inc, p2)
        Rho1 = Rho_Water(T_wa1, p1) - Dens_wa
        Rho2 = Rho_Water(T_wa2, p2) - Dens_wa
        p1 = Water_Pressure_calc(T_wa1 + 273.15, Rho1 + Dens_wa) * 10
        p2 = Water_Pressure_calc(T_wa2 + 273.15, Rho2 + Dens_wa) * 10
    
        der = (Rho1 - Rho2) / (p1 - p2)
        p1 = p1 - Rho1 / der
        p2 = p1 + 0.1
    
        If Abs(Rho1) < 1E-05 Or i > 250 Then EndLoop = True
        i = i + 1
    Wend
    Brine_dP_dT_for_Isochore = (p1 - Pres_in_Bar) / T_Inc
End If
End Function

Public Function Brine_Enthalpy(WtPercent#, Temp_in_C#, Pres_in_Bar#, Optional Precision As Boolean = False)
Dim IDoNno#, Tmp#
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("x", "T", "P")
LwLim = Array(0, 0, 1)
ULim = Array(100, 1000, 5000)
Vl = Array(WtPercent, Temp_in_C, Pres_in_Bar)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    Brine_Enthalpy = msg
Else
    IDoNno = 1
    Tmp = T_star_H(WtToMol(WtPercent) / 100, IDoNno, Temp_in_C, Pres_in_Bar)
    If Single_Phase_Checker(Temp_in_C, Pres_in_Bar, WtToMol(WtPercent) / 100, False) Then
        High_Accuracy = Precision
        Brine_Enthalpy = 1000 * Water_Enthalpy_calc(Tmp, Rho_Water(Tmp - 273.15, Pres_in_Bar))
    Else
        Brine_Enthalpy = "more than 1 phase"
    End If
End If
End Function

Public Function Brine_Isobaric_heat_caps(WtPercent#, Temp_in_C#, Pres_in_Bar#, Optional Precision As Boolean = False)
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("x", "T", "P")
LwLim = Array(0, 0, 1)
ULim = Array(100, 1000, 5000)
Vl = Array(WtPercent, Temp_in_C, Pres_in_Bar)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    Brine_Isobaric_heat_caps = msg
Else
    High_Accuracy = Precision
    If Single_Phase_Checker(Temp_in_C, Pres_in_Bar, WtToMol(WtPercent) / 100, False) Then
        Brine_Isobaric_heat_caps = 1000 * Isob_Heat_cap(WtToMol(WtPercent) / 100, Temp_in_C, Pres_in_Bar)
    Else
        Brine_Isobaric_heat_caps = "more than 1 phase"
    End If
End If
End Function

Public Function SiO2_isobaric_change_solubility(WtPercent#, Temp_in_C#, Pres_in_Bar#, Dens_in_kg#)
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("x", "T", "P")
LwLim = Array(0, 0, 1)
ULim = Array(100, 1000, 5000)
Vl = Array(WtPercent, Temp_in_C, Pres_in_Bar)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    SiO2_isobaric_change_solubility = msg
Else
    If Single_Phase_Checker(Temp_in_C, Pres_in_Bar, WtToMol(WtPercent) / 100, False) Then
        SiO2_isobaric_change_solubility = dQtzDT_const_P(WtToMol(WtPercent) / 100, Temp_in_C, Pres_in_Bar, Dens_in_kg)
    Else
        SiO2_isobaric_change_solubility = "more than 1 phase"
    End If
End If
End Function

Public Function SiO2_isothermal_change_solubility(WtPercent#, Temp_in_C#, Pres_in_Bar#, Dens_in_kg#)
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("x", "T", "P")
LwLim = Array(0, 0, 1)
ULim = Array(100, 1000, 5000)
Vl = Array(WtPercent, Temp_in_C, Pres_in_Bar)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    SiO2_isothermal_change_solubility = msg
Else
    If Single_Phase_Checker(Temp_in_C, Pres_in_Bar, WtToMol(WtPercent) / 100, False) Then
        SiO2_isothermal_change_solubility = dQtzDP_const_T(WtToMol(WtPercent) / 100, Temp_in_C, Pres_in_Bar, Dens_in_kg)
    Else
        SiO2_isothermal_change_solubility = "more than 1 phase"
    End If
End If
End Function

Public Function SiO2_solubility(WtPercent#, Temp_in_C#, Pres_in_Bar#, Dens_in_kg#)
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("x", "T", "P")
LwLim = Array(0, 0, 1)
ULim = Array(100, 1000, 5000)
Vl = Array(WtPercent, Temp_in_C, Pres_in_Bar)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    SiO2_solubility = msg
Else
    If Single_Phase_Checker(Temp_in_C, Pres_in_Bar, WtToMol(WtPercent) / 100, False) Then
        SiO2_solubility = Qtz_solubility(WtToMol(WtPercent) / 100, Temp_in_C, Dens_in_kg)
    Else
        SiO2_solubility = "more than 1 phase"
    End If
End If
End Function

Public Function Mineral_Solubility(Min_ID$, WtPercent#, Temp_in_C#, Pres_in_Bar#, Dens_in_kg#)
Dim v_name(), LwLim(), ULim(), Vl(), msg$, tmp1#
v_name = Array("x", "T", "P")
LwLim = Array(0, 0, 1)
ULim = Array(100, 1000, 5000)
Vl = Array(WtPercent, Temp_in_C, Pres_in_Bar)
Min_ID = LCase(Min_ID)
tmp1 = WtPercent
If Min_ID <> "ap" And Min_ID <> "calc" And Min_ID <> "cor" And Min_ID <> "fl" And Min_ID <> "qtz" And Min_ID <> "ru" Then
msg = "Please use these mineral abbreviations : Ap/Calc/Cor/Fl/Qtz/Ru. "
End If

msg = msg + SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    Mineral_Solubility = msg
Else
    If Single_Phase_Checker(Temp_in_C, Pres_in_Bar, WtToMol(WtPercent) / 100, False) Then
        Mineral_Solubility = mnslb(Min_ID, WtToMol(WtPercent) / 100, Temp_in_C#, Dens_in_kg)
    Else
        Mineral_Solubility = "more than 1 phase"
    End If
End If
End Function


Public Function Brine_Xsi_Critical_Region(WtPercent#, Temp_in_C#, Pres_in_Bar#, Dens_in_kg#, b_beta#, Optional Precision As Boolean = False)
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("x", "T", "P")
LwLim = Array(0, 0, 1)
ULim = Array(100, 1000, 5000)
Vl = Array(WtPercent, Temp_in_C, Pres_in_Bar)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    Brine_Xsi_Critical_Region = msg
Else
    High_Accuracy = Precision
    If Single_Phase_Checker(Temp_in_C, Pres_in_Bar, WtToMol(WtPercent) / 100, False) Then
        Brine_Xsi_Critical_Region = XsiCapital(b_beta, WtToMol(WtPercent) / 100, Temp_in_C, Pres_in_Bar, Dens_in_kg)
    Else
        Brine_Xsi_Critical_Region = "more than 1 phase"
    End If
End If
End Function

Public Function Brine_Isothermal_Compressibility_beta(WtPercent#, Temp_in_C#, Pres_in_Bar#, Dens_in_kg#, Optional Precision As Boolean = False)
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("x", "T", "P")
LwLim = Array(0, 0, 1)
ULim = Array(100, 1000, 5000)
Vl = Array(WtPercent, Temp_in_C, Pres_in_Bar)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    Brine_Isothermal_Compressibility_beta = msg
Else
    High_Accuracy = Precision
    If Single_Phase_Checker(Temp_in_C, Pres_in_Bar, WtToMol(WtPercent) / 100, False) Then
        Brine_Isothermal_Compressibility_beta = Beta(WtToMol(WtPercent) / 100, Temp_in_C, Pres_in_Bar, Dens_in_kg)
    Else
        Brine_Isothermal_Compressibility_beta = "more than 1 phase"
    End If
End If
End Function

Public Function Brine_Isobaric_Expansivity_alpha(WtPercent#, Temp_in_C#, Pres_in_Bar#, Dens_in_kg#, Optional Precision As Boolean = False)
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("x", "T", "P")
LwLim = Array(0, 0, 1)
ULim = Array(100, 1000, 5000)
Vl = Array(WtPercent, Temp_in_C, Pres_in_Bar)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    Brine_Isobaric_Expansivity_alpha = msg
Else
    High_Accuracy = Precision
    If Single_Phase_Checker(Temp_in_C, Pres_in_Bar, WtToMol(WtPercent) / 100, False) Then
        Brine_Isobaric_Expansivity_alpha = Alpha(WtToMol(WtPercent) / 100, Temp_in_C, Pres_in_Bar, Dens_in_kg)
    Else
        Brine_Isobaric_Expansivity_alpha = "more than 1 phase"
    End If
End If
End Function

Public Function Brine_Density(WtPercent#, Temp_in_C#, Pres_in_Bar#, Optional Precision As Boolean = False)
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("x", "T", "P")
LwLim = Array(0, 0, 1)
ULim = Array(100, 1000, 5000)
Vl = Array(WtPercent, Temp_in_C, Pres_in_Bar)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    Brine_Density = msg
Else
    High_Accuracy = Precision
    If Single_Phase_Checker(Temp_in_C, Pres_in_Bar, WtToMol(WtPercent) / 100, False) Then
        Brine_Density = Rho_Brine(WtToMol(WtPercent) / 100, Temp_in_C, Pres_in_Bar)
    Else
        Brine_Density = "more than 1 phase"
    End If
End If
End Function

Public Function Brine_Viscosity(WtPercent#, Temp_in_C#, Pres_in_Bar#)
Dim tmp1#, Tmp2#
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("x", "T", "P")
LwLim = Array(0, 0, 1)
ULim = Array(100, 1000, 5000)
Vl = Array(WtPercent, Temp_in_C, Pres_in_Bar)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    Brine_Viscosity = msg
Else
    If Single_Phase_Checker(Temp_in_C, Pres_in_Bar, WtToMol(WtPercent) / 100, False) Then
        If Visc_Inc_data_Fit(T_star_Mu(WtPercent, Temp_in_C) - 273.15, Pres_in_Bar) Then
            Tmp2 = Rho_Brine(WtToMol(WtPercent) / 100, Temp_in_C, Pres_in_Bar)
            Brine_Viscosity = Viscosity_H2O_NaCl(WtPercent / 100, Temp_in_C, Pres_in_Bar, Tmp2)
            Exit Function
        Else
            Brine_Viscosity = "PTx are outside of model limits"
        End If
    Else
        Brine_Viscosity = "more than 1 phase"
    End If
End If

End Function

Public Function IAPWS_Water_Density(Temp_in_C#, Pres_in_Bar#, Optional Precision As Boolean = False)
Dim v_name(), LwLim(), ULim(), Vl(), msg$

v_name = Array("T", "P")
LwLim = Array(0, 0)
ULim = Array(1000, 100000)
Vl = Array(Temp_in_C, Pres_in_Bar)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    IAPWS_Water_Density = msg
Else
    High_Accuracy = Precision
    IAPWS_Water_Density = Rho_Water(Temp_in_C, Pres_in_Bar)
End If
End Function

Public Function IAPWS_Water_Vap_Dens_L_sat(Temp_in_C#)
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("T")
LwLim = Array(0)
ULim = Array(373.946)
Vl = Array(Temp_in_C)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    IAPWS_Water_Vap_Dens_L_sat = msg
Else
    IAPWS_Water_Vap_Dens_L_sat = Rho_Water_Vap_sat(Temp_in_C)
End If
End Function

Public Function IAPWS_Water_Liq_Dens_V_sat(Temp_in_C#)
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("T")
LwLim = Array(0)
ULim = Array(373.946)
Vl = Array(Temp_in_C)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    IAPWS_Water_Liq_Dens_V_sat = msg
Else
    IAPWS_Water_Liq_Dens_V_sat = Rho_Water_Liq_sat(Temp_in_C)
End If
End Function

Public Function IAPWS_Water_Isochoric_heat_capacity(Temp_in_C#, H2O_Density#)
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("T")
LwLim = Array(0)
ULim = Array(1000)
Vl = Array(Temp_in_C)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    IAPWS_Water_Isochoric_heat_capacity = msg
ElseIf Water_Pressure_calc(Temp_in_C + 273.15, H2O_Density) * 10 > 100000# Then
    IAPWS_Water_Isochoric_heat_capacity = "density is too high so it's above pressure limit of the model for a given temperature"
Else
    IAPWS_Water_Isochoric_heat_capacity = 1000 * Water_Isochoric_Heat_capacity_calc(Temp_in_C + 273.15, H2O_Density)
End If
End Function

Public Function IAPWS_Water_Isob_heat_capacity(Temp_in_C#, H2O_Density#)
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("T")
LwLim = Array(0)
ULim = Array(1000)
Vl = Array(Temp_in_C)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    IAPWS_Water_Isob_heat_capacity = msg
ElseIf Water_Pressure_calc(Temp_in_C + 273.15, H2O_Density) * 10 > 100000# Then
    IAPWS_Water_Isob_heat_capacity = "density is too high so it's above pressure limit of the model for a given temperature"
Else
    IAPWS_Water_Isob_heat_capacity = 1000 * Water_Isobaric_Heat_capacity_calc(Temp_in_C + 273.15, H2O_Density)
End If
End Function

Public Function IAPWS_Water_Enthalpy(Temp_in_C#, H2O_Density#)
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("T")
LwLim = Array(0)
ULim = Array(1000)
Vl = Array(Temp_in_C)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    IAPWS_Water_Enthalpy = msg
ElseIf Water_Pressure_calc(Temp_in_C + 273.15, H2O_Density) * 10 > 100000# Then
    IAPWS_Water_Enthalpy = "density is too high so it's above pressure limit of the model for a given temperature"
Else
    IAPWS_Water_Enthalpy = Water_Enthalpy_calc(Temp_in_C + 273.15, H2O_Density) * 1000
End If
End Function

Public Function IAPWS_Water_Pressure(Temp_in_C#, H2O_Density#)
Dim v_name(), LwLim(), ULim(), Vl(), msg$
v_name = Array("T")
LwLim = Array(0)
ULim = Array(1000)
Vl = Array(Temp_in_C)

msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    IAPWS_Water_Pressure = msg
Else
    IAPWS_Water_Pressure = Water_Pressure_calc(Temp_in_C + 273.15, H2O_Density) * 10
End If
End Function
Public Function Halite_liquid_density(Temp_in_C#, P_in_Bar#)
Dim v_name(), LwLim(), ULim(), Vl(), msg$

v_name = Array("T", "P")
LwLim = Array(800.7, 6E-06)
ULim = Array(1000, 5000)
Vl = Array(Temp_in_C, P_in_Bar)
msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    Halite_liquid_density = msg
Else
    Halite_liquid_density = NaCl_Rho_Liq(Temp_in_C, P_in_Bar) / 1000
End If

End Function

Public Function Halite_solid_density(Temp_in_C#, P_in_Bar#)
Dim v_name(), LwLim(), ULim(), Vl(), msg$

v_name = Array("T", "P")
LwLim = Array(0, 0)
ULim = Array(924.4, 5000)
Vl = Array(Temp_in_C, P_in_Bar)
msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    Halite_solid_density = msg
Else
    Halite_solid_density = NaCl_Rho_Solid(Temp_in_C, P_in_Bar) / 1000
End If

End Function

Public Function Halite_spec_enthalpy(Temp_in_C#, P_in_Bar#)
Dim v_name(), LwLim(), ULim(), Vl(), msg$

v_name = Array("T", "P")
LwLim = Array(0, 0)
ULim = Array(1000, 5000)
Vl = Array(Temp_in_C, P_in_Bar)
msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    Halite_spec_enthalpy = msg
Else
    Halite_spec_enthalpy = NaCl_specific_enthalpy(Temp_in_C, P_in_Bar)
End If

End Function
Public Function Halite_isob_heat_capacity(Temp_in_C#, P_in_Bar#)
Dim v_name(), LwLim(), ULim(), Vl(), msg$

v_name = Array("T", "P")
LwLim = Array(0, 0)
ULim = Array(1000, 5000)
Vl = Array(Temp_in_C, P_in_Bar)
msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)

If msg <> "" Then
    Halite_isob_heat_capacity = msg
Else
    Halite_isob_heat_capacity = NaCl_isobaric_heat_capacity(Temp_in_C, P_in_Bar)
End If

End Function



Public Function IAPWS_Water_Viscosity(Temp_in_C#, H2O_Density#)
Dim v_name(), LwLim(), ULim(), Vl(), msg$

If Visc_Inc_data_Fit(Temp_in_C#, Water_Pressure_calc(Temp_in_C + 273.15, H2O_Density) * 10) Then
    v_name = Array("T")
    LwLim = Array(0)
    ULim = Array(900)
    Vl = Array(Temp_in_C)
    msg = SuplFuncs_FoolProofMsg(v_name, LwLim, ULim, Vl)
    If msg <> "" Then
        IAPWS_Water_Viscosity = msg
    Else
        IAPWS_Water_Viscosity = Water_Viscosity_calc(Temp_in_C + 273.15, H2O_Density)
    End If
Else
    IAPWS_Water_Viscosity = "density is too high, it is above pressure limit of the model for a given temperature"
End If
End Function



'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'BELOW ARE SUPPLEMENTARY FUNCTIONS, DESIGNED TO                 '
'INJECT FUNCTION DESCRIPTION, GENERATE SHEET WITH EXAMPLES      '
'TEST THAT INCOMING VALUES ARE MATCHING THE LIMITS              '
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''


Private Function SuplFuncs_FoolProofMsg$(VarName(), LowLim(), upLim(), Val())
' the function that check if P, T, X, Rho are matching the model limit,
' and generating the message describing the problem
' if this private function returns nothing, the public functions proceed successfully
' in most cases density (rho) is exempted from check, assuming it is generated usinb brine_density functionDim i%, str$(), NotPass As Boolean

SuplFuncs_FoolProofMsg = ""

ReDim str(1 To 2, 0 To UBound(VarName, 1)) As String
For i = 0 To UBound(VarName, 1)
NotPass = False
    If Val(i) < LowLim(i) Then
        str(1, i) = "below the "
        NotPass = True
    ElseIf Val(i) > upLim(i) Then
        str(1, i) = "above the "
        NotPass = True
    End If

    If NotPass Then
        Select Case VarName(i)
            Case "T"
                str(2, i) = "Temperature "
            Case "P"
                str(2, i) = "Pressure "
            Case "x"
                str(2, i) = "Salinity "
            Case "rh_h2o"
                str(2, i) = "density of H2O "
        End Select
        SuplFuncs_FoolProofMsg = SuplFuncs_FoolProofMsg & " " & str(2, i) & "is " & str(1, i) & "limit ( " & CStr(Round(LowLim(i), 1) & " - " & Round(upLim(i), 1)) & ")."
    End If

Next i

End Function


'This and following sub are generating "Microthermometry" sheet with examples of each important functions
Sub SupSubs_List_of_User_Functions_for_Mac()
Dim a%, rowCount%, DescCol%, ExampCol%, ArgCol%
SetupExcelForCalc True
If SheetExist("Function examples") Then ActiveWorkbook.Sheets("Function examples").Delete
Sheets.Add.Name = "Function examples"


With ActiveWorkbook.Sheets("Function examples")
    Rows("3:3").Select
    ActiveWindow.FreezePanes = True
    ActiveWindow.SmallScroll Down:=-3

    rowCount = 1
    DescCol = 12
    ExampCol = 4
    ArgCol = 6
    
    Columns("D:D").ColumnWidth = 60
    Range("D1:D65").NumberFormat = "0.000E+0"
    
    .Cells(rowCount, 1) = "Function name": .Cells(rowCount, DescCol) = "Description"
    .Cells(rowCount, ExampCol) = "Example": .Cells(rowCount, ArgCol) = "Arguments for functions"
    
    .Cells(rowCount + 1, ArgCol) = "NaCl, wt. %": .Cells(rowCount + 1, ArgCol + 1) = "T C, H2O-NaCl": .Cells(rowCount + 1, ArgCol + 2) = "P bar, H2O-NaCl"
    .Cells(rowCount + 1, ArgCol + 3) = "H2O-NaCl density, kg/m3": .Cells(rowCount + 1, ArgCol + 4) = "beta coef., H2O-NaCl": .Cells(rowCount + 1, ArgCol + 5) = "H2O density, kg/m3"
    
    Range(CStr(ConvertToLetter(ArgCol) & rowCount + 1 & ":" & ConvertToLetter(ArgCol + 5) & rowCount + 1)).Font.Size = 8
    Range(CStr(ConvertToLetter(ArgCol) & rowCount + 1 & ":" & ConvertToLetter(ArgCol + 5) & rowCount + 1)).WrapText = True
    Range(CStr(ConvertToLetter(ArgCol) & rowCount) & ":" & CStr(ConvertToLetter(ArgCol + 5) & rowCount)).Merge
    rowCount = rowCount + 2
    
    SupSubs_GenerateExample rowCount, DescCol, ExampCol, ArgCol

    .Cells(69, 1) = "The list of standard Excel functions is extended, and can be accessed directly from function builder. Abobe are examples for each of the functions."
    .Cells(70, 1) = "in order to see the the function and its arguments, start typing the function by typing " & Chr(34) & Chr(61) & Chr(34) & " and name of the function, then select it from drop list and push keys Ctrl+Shift+A (for Mac - Cmd+Shift+A)"
    .Cells(71, 1) = "The default units for function arguments are - weight (mass) percent of NaCl, temperature in C, pressure in bar."
    .Cells(74, 1) = "The function S_Unit_Converter converts between different salinity units, with % in (range 0-100), and fractions out (0-1), unless molality units involved."

    .Cells(75, 1) = "=A75": .Cells(75, 4) = "mass percent NaCl --> to mol fraction"
    .Cells(76, 1) = "15.0"
    .Cells(77, 1) = "WtPer": .Cells(77, 2) = "MolPer": .Cells(77, 3) = "VolPer": .Cells(77, 4) = "Molal"
    .Cells(78, 1) = "In cells above are shown arguments for S_Unit_Converter function"
    .Cells(75, 1) = "=S_Unit_Converter(A76, A77, B77)": .Cells(75, 1).NumberFormat = "0.00"
    .Cells(1, 1).Select
    
End With
SetupExcelForCalc False

End Sub
Private Sub SupSubs_GenerateExample(FncRow%, FncDescCol%, exCol%, FncArgCol%)
' this sub is loading function examples on "Microthermometry" sheet

Dim F_n(), F_D(), DefVal$(0 To 5, 0 To 62), F_nn$(), F_Dd(0 To 62), FncExmpRef$(0 To 6, 0 To 62), FncRefStr$(0 To 62), i%
ReDim F_nn$(0 To 62)
'array with default values for properties, x-0, T-1, P-2, RhoNaCl-3, Beta-4, RhoH2O-5
'(0-31)*2 is index for the property #

'DefVal is a default value to load on sheet "functions"
DefVal(0, 0) = 10:  DefVal(1, 0) = 200:     DefVal(2, 0) = 500:     DefVal(3, 0) = "":                  DefVal(4, 0) = "":                      DefVal(5, 0) = ""
DefVal(0, 2) = 10:  DefVal(1, 2) = 200:     DefVal(2, 2) = 500:     DefVal(3, 2) = "":                  DefVal(4, 2) = "":                      DefVal(5, 2) = ""
DefVal(0, 4) = 10:  DefVal(1, 4) = 200:     DefVal(2, 4) = 500:     DefVal(3, 4) = 971.667853903444:    DefVal(4, 4) = "":                      DefVal(5, 4) = ""
DefVal(0, 6) = 10:  DefVal(1, 6) = 200:     DefVal(2, 6) = 500:     DefVal(3, 6) = "":                  DefVal(4, 6) = "":                      DefVal(5, 6) = ""
DefVal(0, 8) = 10:  DefVal(1, 8) = 200:     DefVal(2, 8) = 500:     DefVal(3, 8) = 971.667853903444:    DefVal(4, 8) = "":                      DefVal(5, 8) = ""
DefVal(0, 10) = 10: DefVal(1, 10) = 200:    DefVal(2, 10) = 500:    DefVal(3, 10) = "":                 DefVal(4, 10) = "":                     DefVal(5, 10) = ""
DefVal(0, 12) = 10: DefVal(1, 12) = 200:    DefVal(2, 12) = 500:    DefVal(3, 12) = 971.667853903444:   DefVal(4, 12) = 5.38324831142102E-05:   DefVal(5, 12) = ""
DefVal(0, 14) = "": DefVal(1, 14) = 200:    DefVal(2, 14) = 500:    DefVal(3, 14) = "":                 DefVal(4, 14) = "":                     DefVal(5, 14) = ""
DefVal(0, 16) = "": DefVal(1, 16) = 200:    DefVal(2, 16) = "":     DefVal(3, 16) = "":                 DefVal(4, 16) = "":                     DefVal(5, 16) = 900
DefVal(0, 18) = "": DefVal(1, 18) = 200:    DefVal(2, 18) = "":     DefVal(3, 18) = "":                 DefVal(4, 18) = "":                     DefVal(5, 18) = 900
DefVal(0, 20) = "": DefVal(1, 20) = 200:    DefVal(2, 20) = "":     DefVal(3, 20) = "":                 DefVal(4, 20) = "":                     DefVal(5, 20) = 900
DefVal(0, 22) = "": DefVal(1, 22) = 200:    DefVal(2, 22) = "":     DefVal(3, 22) = "":                 DefVal(4, 22) = "":                     DefVal(5, 22) = 900
DefVal(0, 24) = "": DefVal(1, 24) = 200:    DefVal(2, 24) = "":     DefVal(3, 24) = "":                 DefVal(4, 24) = "":                     DefVal(5, 24) = ""
DefVal(0, 26) = "": DefVal(1, 26) = 200:    DefVal(2, 26) = "":     DefVal(3, 26) = "":                 DefVal(4, 26) = "":                     DefVal(5, 26) = ""
DefVal(0, 28) = "": DefVal(1, 28) = 200:    DefVal(2, 28) = "":     DefVal(3, 28) = "":                 DefVal(4, 28) = "":                     DefVal(5, 28) = 900
DefVal(0, 30) = "": DefVal(1, 30) = 200:    DefVal(2, 30) = 500:    DefVal(3, 30) = "":                 DefVal(4, 30) = "":                     DefVal(5, 30) = ""
DefVal(0, 32) = "": DefVal(1, 32) = 500:    DefVal(2, 32) = 300:    DefVal(3, 32) = "":                 DefVal(4, 32) = "":                     DefVal(5, 32) = ""
DefVal(0, 34) = "": DefVal(1, 34) = 400:    DefVal(2, 34) = "":     DefVal(3, 34) = "":                 DefVal(4, 34) = "":                     DefVal(5, 34) = ""
DefVal(0, 36) = "": DefVal(1, 36) = 850:    DefVal(2, 36) = "":     DefVal(3, 36) = "":                 DefVal(4, 36) = "":                     DefVal(5, 36) = ""
DefVal(0, 38) = "": DefVal(1, 38) = 750:    DefVal(2, 38) = "":     DefVal(3, 38) = "":                 DefVal(4, 38) = "":                     DefVal(5, 38) = ""
DefVal(0, 40) = "": DefVal(1, 40) = 400:    DefVal(2, 40) = "":     DefVal(3, 40) = "":                 DefVal(4, 40) = "":                     DefVal(5, 40) = ""
DefVal(0, 42) = "": DefVal(1, 42) = "":     DefVal(2, 42) = 2500:   DefVal(3, 42) = "":                 DefVal(4, 42) = "":                     DefVal(5, 42) = ""
DefVal(0, 44) = "": DefVal(1, 44) = 400:    DefVal(2, 44) = "":     DefVal(3, 44) = "":                 DefVal(4, 44) = "":                     DefVal(5, 44) = ""
DefVal(0, 46) = 20: DefVal(1, 46) = 600:    DefVal(2, 46) = 390:    DefVal(3, 46) = "":                 DefVal(4, 46) = "":                     DefVal(5, 46) = ""
DefVal(0, 48) = "": DefVal(1, 48) = 400:    DefVal(2, 48) = 200:    DefVal(3, 48) = "":                 DefVal(4, 48) = "":                     DefVal(5, 48) = ""
DefVal(0, 50) = "": DefVal(1, 50) = 400:    DefVal(2, 50) = 200:    DefVal(3, 50) = "":                 DefVal(4, 50) = "":                     DefVal(5, 50) = ""
DefVal(0, 52) = 30: DefVal(1, 52) = 350:    DefVal(2, 52) = "":     DefVal(3, 52) = "":                 DefVal(4, 52) = "":                     DefVal(5, 52) = ""
DefVal(0, 54) = 10: DefVal(1, 54) = 200:    DefVal(2, 54) = 500:    DefVal(3, 54) = 971.667853903444:   DefVal(4, 54) = "":                     DefVal(5, 54) = ""
DefVal(0, 56) = 10: DefVal(1, 56) = 200:    DefVal(2, 56) = 500:    DefVal(3, 56) = 971.667853903444:   DefVal(4, 56) = "":                     DefVal(5, 56) = ""
DefVal(0, 58) = 10: DefVal(1, 58) = 200:    DefVal(2, 58) = 500:    DefVal(3, 58) = 971.667853903444:   DefVal(4, 58) = "":                     DefVal(5, 58) = ""
DefVal(0, 60) = 10: DefVal(1, 60) = 200:    DefVal(2, 60) = 500:    DefVal(3, 60) = 971.667853903444:   DefVal(5, 60) = "Calc"
DefVal(0, 62) = 15: DefVal(1, 62) = 400:    DefVal(2, 62) = "None": DefVal(3, 62) = "":                 DefVal(4, 62) = ""

'FncExmpRef is the text for referencing to other cells for a given function
FncExmpRef(0, 0) = CStr(ConvertToLetter(FncArgCol + 0) & FncRow + 0):   FncExmpRef(1, 0) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 0):   FncExmpRef(2, 0) = CStr(ConvertToLetter(FncArgCol + 2) & FncRow + 0): FncExmpRef(3, 0) = "": FncExmpRef(4, 0) = "": FncExmpRef(5, 0) = ""
FncExmpRef(0, 2) = CStr(ConvertToLetter(FncArgCol + 0) & FncRow + 2):   FncExmpRef(1, 2) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 2):   FncExmpRef(2, 2) = CStr(ConvertToLetter(FncArgCol + 2) & FncRow + 2): FncExmpRef(3, 2) = "": FncExmpRef(4, 2) = "": FncExmpRef(5, 2) = ""
FncExmpRef(0, 4) = CStr(ConvertToLetter(FncArgCol + 0) & FncRow + 4):   FncExmpRef(1, 4) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 4):   FncExmpRef(2, 4) = CStr(ConvertToLetter(FncArgCol + 2) & FncRow + 4): FncExmpRef(3, 4) = CStr(ConvertToLetter(FncArgCol + 3) & FncRow + 4): FncExmpRef(4, 4) = "": FncExmpRef(5, 4) = ""
FncExmpRef(0, 6) = CStr(ConvertToLetter(FncArgCol + 0) & FncRow + 6):   FncExmpRef(1, 6) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 6):   FncExmpRef(2, 6) = CStr(ConvertToLetter(FncArgCol + 2) & FncRow + 6): FncExmpRef(3, 6) = "": FncExmpRef(4, 6) = "": FncExmpRef(5, 6) = ""
FncExmpRef(0, 8) = CStr(ConvertToLetter(FncArgCol + 0) & FncRow + 8):   FncExmpRef(1, 8) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 8):   FncExmpRef(2, 8) = CStr(ConvertToLetter(FncArgCol + 2) & FncRow + 8): FncExmpRef(3, 8) = CStr(ConvertToLetter(FncArgCol + 3) & FncRow + 8): FncExmpRef(4, 8) = "": FncExmpRef(5, 8) = ""
FncExmpRef(0, 10) = CStr(ConvertToLetter(FncArgCol + 0) & FncRow + 10): FncExmpRef(1, 10) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 10): FncExmpRef(2, 10) = CStr(ConvertToLetter(FncArgCol + 2) & FncRow + 10): FncExmpRef(3, 10) = "": FncExmpRef(4, 10) = "": FncExmpRef(5, 10) = ""
FncExmpRef(0, 12) = CStr(ConvertToLetter(FncArgCol + 0) & FncRow + 12): FncExmpRef(1, 12) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 12): FncExmpRef(2, 12) = CStr(ConvertToLetter(FncArgCol + 2) & FncRow + 12): FncExmpRef(3, 12) = CStr(ConvertToLetter(FncArgCol + 3) & FncRow + 12): FncExmpRef(4, 12) = CStr(ConvertToLetter(FncArgCol + 4) & FncRow + 12): FncExmpRef(5, 12) = ""
FncExmpRef(0, 14) = "":                                                 FncExmpRef(1, 14) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 14): FncExmpRef(2, 14) = CStr(ConvertToLetter(FncArgCol + 2) & FncRow + 14): FncExmpRef(3, 14) = "": FncExmpRef(4, 14) = "": FncExmpRef(5, 14) = ""
FncExmpRef(0, 16) = "":                                                 FncExmpRef(1, 16) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 16): FncExmpRef(2, 16) = "": FncExmpRef(3, 16) = "": FncExmpRef(4, 16) = "": FncExmpRef(5, 16) = CStr(ConvertToLetter(FncArgCol + 5) & FncRow + 16)
FncExmpRef(0, 18) = "":                                                 FncExmpRef(1, 18) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 18): FncExmpRef(2, 18) = "": FncExmpRef(3, 18) = "": FncExmpRef(4, 18) = "": FncExmpRef(5, 18) = CStr(ConvertToLetter(FncArgCol + 5) & FncRow + 18)
FncExmpRef(0, 20) = "":                                                 FncExmpRef(1, 20) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 20): FncExmpRef(2, 20) = "": FncExmpRef(3, 20) = "": FncExmpRef(4, 20) = "": FncExmpRef(5, 20) = CStr(ConvertToLetter(FncArgCol + 5) & FncRow + 20)
FncExmpRef(0, 22) = "":                                                 FncExmpRef(1, 22) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 22): FncExmpRef(2, 22) = "": FncExmpRef(3, 22) = "": FncExmpRef(4, 22) = "": FncExmpRef(5, 22) = CStr(ConvertToLetter(FncArgCol + 5) & FncRow + 22)
FncExmpRef(0, 24) = "":                                                 FncExmpRef(1, 24) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 24): FncExmpRef(2, 24) = "": FncExmpRef(3, 24) = "": FncExmpRef(4, 24) = "": FncExmpRef(5, 24) = ""
FncExmpRef(0, 26) = "":                                                 FncExmpRef(1, 26) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 26): FncExmpRef(2, 26) = "": FncExmpRef(3, 26) = "": FncExmpRef(4, 26) = "": FncExmpRef(5, 26) = ""
FncExmpRef(0, 28) = "":                                                 FncExmpRef(1, 28) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 28): FncExmpRef(2, 28) = "": FncExmpRef(3, 28) = "": FncExmpRef(4, 28) = "": FncExmpRef(5, 28) = CStr(ConvertToLetter(FncArgCol + 5) & FncRow + 28)
FncExmpRef(0, 30) = "":                                                 FncExmpRef(1, 30) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 30): FncExmpRef(2, 30) = CStr(ConvertToLetter(FncArgCol + 2) & FncRow + 30): FncExmpRef(3, 30) = "": FncExmpRef(4, 30) = "": FncExmpRef(5, 30) = ""
FncExmpRef(0, 32) = "":                                                 FncExmpRef(1, 32) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 32): FncExmpRef(2, 32) = CStr(ConvertToLetter(FncArgCol + 2) & FncRow + 32): FncExmpRef(3, 32) = "": FncExmpRef(4, 32) = "": FncExmpRef(5, 32) = ""
FncExmpRef(0, 34) = "":                                                 FncExmpRef(1, 34) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 34): FncExmpRef(2, 34) = "": FncExmpRef(3, 34) = "": FncExmpRef(4, 34) = "": FncExmpRef(5, 34) = ""
FncExmpRef(0, 36) = "":                                                 FncExmpRef(1, 36) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 36): FncExmpRef(2, 36) = "": FncExmpRef(3, 36) = "": FncExmpRef(4, 36) = "": FncExmpRef(5, 36) = ""
FncExmpRef(0, 38) = "":                                                 FncExmpRef(1, 38) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 38): FncExmpRef(2, 38) = "": FncExmpRef(3, 38) = "": FncExmpRef(4, 38) = "": FncExmpRef(5, 38) = ""
FncExmpRef(0, 40) = "":                                                 FncExmpRef(1, 40) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 40): FncExmpRef(2, 40) = "": FncExmpRef(3, 40) = "": FncExmpRef(4, 40) = "": FncExmpRef(5, 40) = ""
FncExmpRef(0, 42) = "":                                                 FncExmpRef(1, 42) = "":                                                 FncExmpRef(2, 42) = CStr(ConvertToLetter(FncArgCol + 2) & FncRow + 42): FncExmpRef(3, 42) = "": FncExmpRef(4, 42) = "": FncExmpRef(5, 42) = ""
FncExmpRef(0, 44) = "":                                                 FncExmpRef(1, 44) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 44): FncExmpRef(2, 44) = "": FncExmpRef(3, 44) = "": FncExmpRef(4, 44) = "": FncExmpRef(5, 44) = ""

FncExmpRef(0, 46) = CStr(ConvertToLetter(FncArgCol + 0) & FncRow + 46):  FncExmpRef(1, 46) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 46):   FncExmpRef(2, 46) = CStr(ConvertToLetter(FncArgCol + 2) & FncRow + 46): FncExmpRef(3, 0) = "": FncExmpRef(4, 0) = "": FncExmpRef(5, 0) = ""

FncExmpRef(0, 48) = "":                                                 FncExmpRef(1, 48) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 48): FncExmpRef(2, 48) = CStr(ConvertToLetter(FncArgCol + 2) & FncRow + 48): FncExmpRef(3, 48) = "": FncExmpRef(4, 48) = "": FncExmpRef(5, 48) = ""
FncExmpRef(0, 50) = "":                                                FncExmpRef(1, 50) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 50): FncExmpRef(2, 50) = CStr(ConvertToLetter(FncArgCol + 2) & FncRow + 50): FncExmpRef(3, 50) = "": FncExmpRef(4, 50) = "": FncExmpRef(5, 50) = ""
FncExmpRef(0, 52) = CStr(ConvertToLetter(FncArgCol + 0) & FncRow + 52): FncExmpRef(1, 52) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 52): FncExmpRef(2, 52) = "": FncExmpRef(3, 52) = "": FncExmpRef(4, 52) = "": FncExmpRef(5, 52) = ""
FncExmpRef(0, 54) = CStr(ConvertToLetter(FncArgCol + 0) & FncRow + 54): FncExmpRef(1, 54) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 54): FncExmpRef(2, 54) = CStr(ConvertToLetter(FncArgCol + 2) & FncRow + 54): FncExmpRef(3, 54) = CStr(ConvertToLetter(FncArgCol + 3) & FncRow + 54): FncExmpRef(4, 54) = "": FncExmpRef(5, 54) = ""
FncExmpRef(0, 56) = CStr(ConvertToLetter(FncArgCol + 0) & FncRow + 56): FncExmpRef(1, 56) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 56): FncExmpRef(2, 56) = CStr(ConvertToLetter(FncArgCol + 2) & FncRow + 56): FncExmpRef(3, 56) = CStr(ConvertToLetter(FncArgCol + 3) & FncRow + 56): FncExmpRef(4, 56) = "": FncExmpRef(5, 56) = ""
FncExmpRef(0, 58) = CStr(ConvertToLetter(FncArgCol + 0) & FncRow + 58): FncExmpRef(1, 58) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 58): FncExmpRef(2, 58) = CStr(ConvertToLetter(FncArgCol + 2) & FncRow + 58): FncExmpRef(3, 58) = CStr(ConvertToLetter(FncArgCol + 3) & FncRow + 58): FncExmpRef(4, 58) = "": FncExmpRef(5, 58) = ""
FncExmpRef(0, 60) = CStr(ConvertToLetter(FncArgCol + 5) & FncRow + 60): FncExmpRef(1, 60) = CStr(ConvertToLetter(FncArgCol + 0) & FncRow + 60): FncExmpRef(2, 60) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 60): FncExmpRef(3, 60) = CStr(ConvertToLetter(FncArgCol + 2) & FncRow + 60): FncExmpRef(4, 60) = CStr(ConvertToLetter(FncArgCol + 3) & FncRow + 60)
FncExmpRef(0, 62) = CStr(ConvertToLetter(FncArgCol + 0) & FncRow + 62): FncExmpRef(1, 62) = CStr(ConvertToLetter(FncArgCol + 1) & FncRow + 62): FncExmpRef(2, 62) = CStr(ConvertToLetter(FncArgCol + 2) & FncRow + 62): FncExmpRef(3, 62) = "": FncExmpRef(4, 62) = "": FncExmpRef(5, 62) = ""

'accuracy parameter default values
FncExmpRef(6, 0) = "False"
FncExmpRef(6, 2) = "False"
FncExmpRef(6, 4) = "False"
FncExmpRef(6, 6) = "False"
FncExmpRef(6, 8) = "False"
FncExmpRef(6, 12) = "False"
FncExmpRef(6, 14) = "False"
FncExmpRef(6, 52) = "False"

F_n = SuplFuncs_f_and_d_Arrays(True)
F_D = SuplFuncs_f_and_d_Arrays(False)

For i = 0 To 62
    'compiling the text for formula that will be loaded on the sheet
    If i = 0 Then
        F_nn(i) = F_n(i)
        F_Dd(i) = F_D(i)
        For j = 0 To 6
            
            If FncExmpRef(j, i) <> "" Then FncRefStr(i) = FncRefStr(i) + "," + FncExmpRef(j, i)
        Next j
        FncRefStr(i) = "=" & F_nn(i) & "(" & Right(FncRefStr(i), Len(FncRefStr(i)) - 1) & ")"
    Else
        If (i Mod 2) = 0 Then
            F_nn(i) = F_n(i / 2)
            F_Dd(i) = F_D(i / 2)
            For j = 0 To 6
                If FncExmpRef(j, i) <> "" Then FncRefStr(i) = FncRefStr(i) + "," + FncExmpRef(j, i)
            Next j
            FncRefStr(i) = "=" & F_nn(i) & "(" & Right(FncRefStr(i), Len(FncRefStr(i)) - 1) & ")"
        End If
    End If
Next i
SuplFuncs_CellMerge FncRow, 32, FncDescCol, FncDescCol + 6

Range(CStr(ConvertToLetter(1) & FncRow & ":" & ConvertToLetter(1) & FncRow + 62)) = WorksheetFunction.Transpose(F_nn)
Range(CStr(ConvertToLetter(FncDescCol) & FncRow & ":" & ConvertToLetter(FncDescCol) & FncRow + 62)) = WorksheetFunction.Transpose(F_Dd)

Range(CStr(ConvertToLetter(FncArgCol) & FncRow & ":" & ConvertToLetter(FncArgCol + 5) & FncRow + 62)) = WorksheetFunction.Transpose(DefVal)
Range(CStr(ConvertToLetter(exCol) & FncRow & ":" & ConvertToLetter(exCol) & FncRow + 62)) = WorksheetFunction.Transpose(FncRefStr)


End Sub

Sub SupSubs_Function_Description_loader()
'This sub injects additional category of functions in Windows Excel
On Error GoTo Handler
Dim FuncName(), FuncDesc(), TotalProps(), Category$, ArgTXT$(0 To 34, 1 To 6), ArgDesc$()
Dim t_b$, t_d_n$, t_d_h$, t_p$, t_t$, t_wt$, Prec$, m_id$, s_lim$, i%, j%, tmp1$, Tmp2$

Categry = "H2O-NaCl properties"

t_b = "Beta coefficient":     t_d_n = "Density of H2O-NaCl mixture":     t_d_h = "Density of pure H2O"
t_p = "Pressure in bar":     t_t = "Temperature in degree C":              t_wt = "Weight (mass) percent of NaCl in H2O-NaCl"
m_id = "mineral ID, acceptable abbreviations: Ap/Calc/Cor/Fl/Qtz/Ru"
s_lim = "(text)Must be one of: WtPer/MolPer/VolPer/Molal"
Prec = "High precision, either True or False"

FuncName = SuplFuncs_f_and_d_Arrays(True)
    
FuncDesc = SuplFuncs_f_and_d_Arrays(False)
TotalProps = Array(4, 4, 5, 4, 5, 4, 6, 3, 3, 3, 3, 3, 2, 2, 3, 3, 3, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 5, 5, 5, 5, 3, 3, 4)

ArgTXT(0, 1) = t_wt:    ArgTXT(0, 2) = t_t:     ArgTXT(0, 3) = t_p
ArgTXT(1, 1) = t_wt:    ArgTXT(1, 2) = t_t:     ArgTXT(1, 3) = t_p
ArgTXT(2, 1) = t_wt:    ArgTXT(2, 2) = t_t:     ArgTXT(2, 3) = t_p: ArgTXT(2, 4) = t_d_n
ArgTXT(3, 1) = t_wt:    ArgTXT(3, 2) = t_t:     ArgTXT(3, 3) = t_p
ArgTXT(4, 1) = t_wt:    ArgTXT(4, 2) = t_t:     ArgTXT(4, 3) = t_p: ArgTXT(4, 4) = t_d_n
ArgTXT(5, 1) = t_wt:    ArgTXT(5, 2) = t_t:     ArgTXT(5, 3) = t_p
ArgTXT(6, 1) = t_wt:    ArgTXT(6, 2) = t_t:     ArgTXT(6, 3) = t_p: ArgTXT(6, 4) = t_d_n: ArgTXT(6, 5) = t_b
ArgTXT(7, 1) = t_t:     ArgTXT(7, 2) = t_p
ArgTXT(8, 1) = t_t:     ArgTXT(8, 2) = t_d_h
ArgTXT(9, 1) = t_t:     ArgTXT(9, 2) = t_d_h
ArgTXT(10, 1) = t_t:    ArgTXT(10, 2) = t_d_h
ArgTXT(11, 1) = t_t:    ArgTXT(11, 2) = t_d_h
ArgTXT(12, 1) = t_t
ArgTXT(13, 1) = t_t
ArgTXT(14, 1) = t_t:    ArgTXT(14, 2) = t_d_h
ArgTXT(15, 1) = t_t:    ArgTXT(15, 2) = t_p
ArgTXT(16, 1) = t_t:    ArgTXT(16, 2) = t_p
ArgTXT(17, 1) = t_t
ArgTXT(18, 1) = t_t
ArgTXT(19, 1) = t_t
ArgTXT(20, 1) = t_t
ArgTXT(21, 1) = t_p
ArgTXT(22, 1) = t_t

ArgTXT(23, 1) = t_wt:    ArgTXT(23, 2) = t_t:    ArgTXT(23, 3) = t_p
ArgTXT(24, 1) = t_t:    ArgTXT(24, 2) = t_p
ArgTXT(25, 1) = t_t:    ArgTXT(25, 2) = t_p
ArgTXT(26, 1) = t_wt:   ArgTXT(26, 2) = t_t
ArgTXT(27, 1) = t_wt:   ArgTXT(27, 2) = t_t:    ArgTXT(27, 3) = t_p: ArgTXT(27, 4) = t_d_n
ArgTXT(28, 1) = t_wt:   ArgTXT(28, 2) = t_t:    ArgTXT(28, 3) = t_p: ArgTXT(28, 4) = t_d_n
ArgTXT(29, 1) = t_wt:   ArgTXT(29, 2) = t_t:    ArgTXT(29, 3) = t_p: ArgTXT(29, 4) = t_d_n
ArgTXT(30, 1) = m_id:   ArgTXT(30, 2) = t_wt:   ArgTXT(30, 3) = t_t: ArgTXT(30, 4) = t_p: ArgTXT(30, 5) = t_d_n
ArgTXT(31, 1) = t_wt:   ArgTXT(31, 2) = t_t:    ArgTXT(31, 3) = t_p:
ArgTXT(32, 1) = "numeric value for salinity":   ArgTXT(32, 2) = s_lim:   ArgTXT(32, 3) = s_lim
ArgTXT(33, 1) = "Phase melting temperature, in degree C": ArgTXT(33, 2) = "(Optional, unnecessary for ice melting) Bubble dissapearance temperature, in degree C": ArgTXT(33, 3) = "(Optional, text) which phase melted, ice (default value), hydrohalite, halite"
ArgTXT(34, 1) = t_wt:   ArgTXT(34, 2) = t_t:    ArgTXT(34, 3) = t_p: ArgTXT(34, 4) = t_d_n


ArgTXT(0, 4) = Prec: ArgTXT(1, 4) = Prec: ArgTXT(2, 5) = Prec: ArgTXT(3, 4) = Prec: ArgTXT(4, 5) = Prec:
ArgTXT(6, 6) = Prec: ArgTXT(7, 3) = Prec: ArgTXT(25, 3) = Prec


For i = 0 To 34
    ReDim ArgDesc(1 To TotalProps(i))
    For j = 1 To UBound(ArgDesc, 1)
        ArgDesc(j) = ArgTXT(i, j)
    Next j
    tmp1 = FuncName(i)
    Tmp2 = FuncDesc(i)
    Application.MacroOptions Macro:=tmp1, Description:=Tmp2, Category:=Categry, ArgumentDescriptions:=ArgDesc

Next i
Handler:
End Sub

Private Function SuplFuncs_f_and_d_Arrays(F_n As Boolean)
'returns either the function name or function description
If F_n = True Then
    SuplFuncs_f_and_d_Arrays = Array("Brine_Density", "Brine_Enthalpy", "Brine_Isobaric_Expansivity_alpha", _
    "Brine_Isobaric_heat_caps", "Brine_Isothermal_Compressibility_beta", "Brine_Viscosity", _
    "Brine_Xsi_Critical_Region", "IAPWS_Water_Density", "IAPWS_Water_Enthalpy", "IAPWS_Water_Isob_heat_capacity", _
    "IAPWS_Water_Isochoric_heat_capacity", "IAPWS_Water_Pressure", "IAPWS_Water_Liq_Dens_V_sat", "IAPWS_Water_Vap_Dens_L_sat", _
    "IAPWS_Water_Viscosity", "Phases_x_Halite_Liquidus", "Phases_x_Halite_saturated_vapor", "Phases_P_crit", _
    "Phases_P_Halite_Boil", "Phases_P_Halite_Subl", "Phases_P_LVH_Coexistence", "Phases_T_Halite_Melt", _
    "Phases_x_crit", "Phases_x_of_coexisting_phases", "Phases_x_VL_coex_Liquid_side", "Phases_x_VL_coex_Vapor_side", "Phases_Single_ph_pressure_at_TX", _
    "SiO2_isobaric_change_solubility", "SiO2_isothermal_change_solubility", "SiO2_solubility", "Mineral_Solubility", "Phases_PTx_Surveyor", _
    "S_Unit_Converter", "FlInc_Salinity", "Brine_dP_dT_for_Isochore")
Else
    SuplFuncs_f_and_d_Arrays = Array("Returns the density of H2O-NaCl, in kg/m3" & Chr(10) & "eq. 7, Driesner (2007)", "Returns the specific enthalpy of H2O-NaCl, in J/kg" & Chr(10) & "eq. 21, Driesner (2007)", _
    "Returns the isobaric expansivity coefficient (alpha), in 1/K" & Chr(10) & "eq. 8, Klyukin et al (2016)", "Returns the isobaric heat capacity of H2O-NaCl, in J /(kg K)" & Chr(10) & "eq. 27 (Driesner 2007), divided by (1+Xnacl) where Xnacl is a molar fraction. This term was omitted in original publication but should be used, according to personal communications with Driesner. ", _
    "Returns the isothermal compressibility (beta), in 1/bar" & Chr(10) & "eq. 7, Klyukin et al (2016)", "Returns the viscosity of H2O-NaCl single phase mixture, in micropascal/second" & Chr(10) & "eq. 3, Klyukin et al (2017)", _
    "Returns the reduced succeptibility parameter, indicating critical region limits at value of 0.5, dimensionless" & Chr(10) & "eq. in text of 2.2.3.1 section, Anisimov et al (2004)", "Returns the density of pure H2O according to IAPWS-95 formulation, in kg/m3" & Chr(10) & "Root finding algorithms applied to eq. for H2O pressure at given temperature and density listed in Table 6.3, Wagner and Pruss (2002)", _
    "Returns the specific enthalpy of H2O, calculated according to IAPWS-95 formulation, in J/kg" & Chr(10) & "eq. listed in Table 6.3, Wagner and Pruss (2002)", "Returns the isobaric heat capacity of H2O, calculated according to IAPWS-95 formulation, in in Joule/(kg K)" & Chr(10) & "eq. listed in Table 6.3, Wagner and Pruss (2002)", _
    "Returns the isochoric heat capacity of H2O, calculated according to IAPWS-95 formulation, in in Joule/(kg K)" & Chr(10) & "eq. listed in Table 6.3, Wagner and Pruss (2002)", "Returns the pressure of H2O at given conditions, according to IAPWS-95 formulation, in bar" & Chr(10) & "eq. listed in Table 6.3, Wagner and Pruss (2002)", _
    "Returns the vapor-saturated liquid density of pure H2O according to IAPWS-95 formulation, in kg/m3" & Chr(10) & "eq.2.6, Wagner and Pruss (2002)", "Returns the liquid-saturated vapor density of pure H2O according to IAPWS-95 formulation, in kg/m3" & Chr(10) & "eq.2.7, Wagner and Pruss (2002)", _
    "Returns the dynamic viscosity of pure H2O according to the IAPWS Formulation 2008, extending IAPWS-95, in micropascal/sec" & Chr(10) & "eq.2, Huber et al (2009)", "Returns the concentration for halite liquidus at given pressure and temperature, in wt. % NaCl" & Chr(10) & "eq. 8, Driesner and Heinrich (2007)", _
    "Returns the concentration for halite-saturated vapor at given pressure and temperature, in wt. % NaCl" & Chr(10) & "eq. 9, Driesner and Heinrich (2007)", "Returns the critical pressure at given temperature, in bar" & Chr(10) & "eq. 5 and eq. 7, Driesner and Heinrich (2007)", _
    "Returns the pressure of the halite boiling at given temperature, in bar" & Chr(10) & "eq. 3, Driesner and Heinrich (2007)", "Returns the pressure of the halite sublimation at given temperature, in bar" & Chr(10) & "eq. 2, Driesner and Heinrich (2007)", _
    "Returns the pressure (bar) of liquid-vapor-halite coexisting curve for given temperature, in bar" & Chr(10) & "eq. 10, Driesner and Heinrich (2007)", "Returns the temperature of halite melting, in C" & Chr(10) & "eq. 1, Driesner and Heinrich (2007)", _
    "Returns the concentration for critical point at given temperature, in wt. % NaCl" & Chr(10) & "eq. 5 and eq. 7, Driesner and Heinrich (2007)", _
    "Returns the salinity for coexisting phases, if the fluid is in two phase fluid state" & Chr(10) & "eq. 8, 11, 12, Driesner and Heinrich (2007)", _
    "Returns the concentration of liquid-vapor coexistence surface at liquid branch, for given temperature and pressure, in wt. % NaCl" & Chr(10) & "eq. 12, Driesner and Heinrich (2007)", _
    "Returns the concentration of liquid-vapor coexistence surface at vapor branch, for given temperature and pressure, in wt. % NaCl" & Chr(10) & "eq. 11, Driesner and Heinrich (2007)", _
    "Returns the LV pressure at the point of interest, in bar", "Returns the coefficient of isobaric change in SiO2 solubility, in 1/K" & Chr(10) & "eq. 12, Klyukin et al (2016)", _
    "Returns the coefficient of isothermal change in SiO2 solubility, in 1/bar" & Chr(10) & "eq. 11, Klyukin et al (2016)", "Returns the SiO2 solubility in H2O-NaCl, in mol per kg H2O" & Chr(10) & "eq. 10, Akinfief and Diamond (2009)", _
    "Returns the mineral solubility, in mol per kg of H2O, for quartz, calcite, corundum, fluorapatite, fluorite and rutile" & Chr(10) & "from Brooks and Steele-MacInnis (in press)", _
    "Checks the state of the fluid for a given PTx point, or estimates coexisting phase along the lines of provided PT, Px, Tx coordinates (at least 2 parameters are required)" & Chr(10) & "eqs of Driesner and Heinrich (2007)", _
    "Converts salinity, works in any combintaions between weight(mass), molar, volume percents and molality", _
    "Estimates salinity based on melting of daughter phase, bubble homogenization, and phase ID" & Chr(10) & "from (Bodnar 1993), or (Sterner et al. 1988) or (Lecumberri-Sanchez et al. 2012)", _
    "Returns the dP/dT slope of the isochore at given parameters," & Chr(10) & "eq. 7, Driesner (2007)")
End If
End Function


Private Sub SuplFuncs_CellMerge(Starting_Row%, PropsCount%, StartCol%, EndCol%)

For i = 0 To PropsCount - 1
    Range(CStr(ConvertToLetter(StartCol) & Starting_Row + i * 2 & ":" & ConvertToLetter(EndCol) & Starting_Row + 1 + i * 2)).Select
    With Selection
        .HorizontalAlignment = xlLeft
        .VerticalAlignment = xlTop
        .WrapText = True
        .Orientation = 0
        .AddIndent = False
        .IndentLevel = 0
        .ShrinkToFit = False
        .ReadingOrder = xlContext
        .MergeCells = True
    End With
    
Next i

End Sub
Private Function SuplFuncs_Report_Number$(number#)
'rounding number to a the first meaningful digit, avoiding to return 0 for low values
Dim i%, toexport#
i = 1
While Round(number, i) = 0
    i = i + 1
Wend
toexport = Round(number, i)
SuplFuncs_Report_Number = CStr(toexport)
End Function

Private Function SuplFuncs_LV_TFinder#(P#, x#)
Dim EndLoop As Boolean, i%
Dim t1#, t2#, t_new#, x1#, x2#
Dim xCP#, pCP#
tempTemperature = SuplFuncs_CP_TFinder(P)
t1 = tempTemperature
t2 = 1000

If t1 > 1000 Then SuplFuncs_LV_TFinder = 1000: Exit Function

X_and_P_crit t1, xCP, pCP
xCP = MolToWt(100 * xCP)
While Not EndLoop
    t_new = (t1 + t2) / 2
    If x >= xCP Then
        x1 = MolToWt(100 * X_VL_Liq(t1, P)) - x
        x2 = MolToWt(100 * X_VL_Liq(t_new, P)) - x
    Else
        x1 = MolToWt(100 * X_VL_Vap(t1, P)) - x
        x2 = MolToWt(100 * X_VL_Vap(t_new, P)) - x
    End If
    If x1 > 100 Then x1 = -x1
    If x2 > 100 Then x2 = -x2
    If Sgn(x1) = Sgn(x2) Then t1 = t_new Else t2 = t_new
    If Abs(t1 - t2) < 1E-05 Or i > 100 Then EndLoop = True
    i = i + 1
Wend
SuplFuncs_LV_TFinder = t1

End Function
Private Function SuplFuncs_CP_TFinder#(P#)
Dim EndLoop As Boolean, i%
Dim t1#, t2#, p1#, p2#, der#

t1 = 400
t2 = 401
While Not EndLoop
    X_and_P_crit t1, der, p1
    X_and_P_crit t2, der, p2
    p1 = p1 - P
    p2 = p2 - P
    der = (p1 - p2) / (t1 - t2)
    t1 = t1 - p1 / der
    t2 = t1 + 0.1
    If Abs(p1) < 1E-05 Or i > 100 Then EndLoop = True
    i = i + 1
Wend
SuplFuncs_CP_TFinder = t1

End Function
Private Function SuplFuncs_VH_TFinder#(P#, LowEnd As Boolean)
Dim EndLoop
Dim t1#, t2#, p1#, p2#, T_Tmp#, der#
If LowEnd Then
    t1 = 590
    t2 = 801
Else
    t1 = 0
    t2 = 600
End If

While Not EndLoop
    T_Tmp = (t1 + t2) / 2
    p1 = P_VLH(t1) - P
    p2 = P_VLH(T_Tmp) - P
    
    If Sgn(p1) = Sgn(p2) Then t1 = T_Tmp Else t2 = T_Tmp

    If Abs(p1 - p2) < 0.0001 Or i > 250 Then EndLoop = True
    i = i + 1
Wend
    SuplFuncs_VH_TFinder = t1
End Function
Private Function SuplFuncs_LH_TFinder#(P#, x#)
Dim EndLoop
Dim t1#, t2#, x1#, x2#, VoidIt#, der#, Tmp#
t1 = 0
t2 = 930
Tmp = WtToMol(x) / 100
While Not EndLoop
    x1 = X_L_Sat(t1, P) - Tmp
    x2 = X_L_Sat(t2, P) - Tmp
    der = (x1 - x2) / (t1 - t2)
    If Abs(der) < 1E-09 Then GoTo WeReGood
    t1 = t1 - x1 / der
    t2 = t1 + 0.1

    If Abs(x1) < 1E-08 Or i > 250 Then EndLoop = True
    i = i + 1
Wend
WeReGood:
SuplFuncs_LH_TFinder = t1
End Function

Public Function S_Unit_Converter#(Salinity_Value#, Salinity_Unit_IN$, Salinity_Unit_OUT$)
'salinity unit out options should be defined from these 4 variants
' "WtPer"  "MolPer"  "VolPer"  "Molal"
'Salinity_Value should be between 0 to 100% no matter what units are (except molality)
'output values are shown in fractions, i.e. from 0 to 1 (except molality)

Dim mH2O#, mNaCl#, RhoH2O#, RhoNaCl#

mH2O = 18.015268
mNaCl = 58.4428
RhoH2O = 1
RhoNaCl = 2.16
Salinity_Value = Salinity_Value / 100
Salinity_Unit_IN = LCase(Salinity_Unit_IN)
Salinity_Unit_OUT = LCase(Salinity_Unit_OUT)

If Salinity_Unit_OUT = Salinity_Unit_IN Or Salinity_Value = 0 Then
    S_Unit_Converter = Salinity_Value
Else
    'converting anything to Weight %%
    If Salinity_Unit_OUT = "wtper" Then
        Select Case Salinity_Unit_IN
        Case "molper"
            S_Unit_Converter = Salinity_Value * mNaCl / (Salinity_Value * mNaCl + (1 - Salinity_Value) * mH2O)
        Case "volper"
             S_Unit_Converter = Salinity_Value * RhoNaCl / (Salinity_Value * RhoNaCl + (1 - Salinity_Value) * RhoH2O)
        Case Else
             S_Unit_Converter = Salinity_Value * mNaCl * 100 / (1000 + Salinity_Value * mNaCl * 100)
        End Select
    'converting anything to Mole %%
    ElseIf Salinity_Unit_OUT = "molper" Then
        Select Case Salinity_Unit_IN
        Case "wtper"
            S_Unit_Converter = Salinity_Value / mNaCl / (Salinity_Value / mNaCl + (1 - Salinity_Value) / mH2O)
        Case "volper"
            S_Unit_Converter = Salinity_Value * RhoNaCl / mNaCl / (Salinity_Value * RhoNaCl / mNaCl + (1 - Salinity_Value) * (RhoH2O / mH2O))
        Case Else
            S_Unit_Converter = Salinity_Value * mNaCl * 100 / (1000 + Salinity_Value * mNaCl * 100)
            S_Unit_Converter = S_Unit_Converter / mNaCl / (S_Unit_Converter / mNaCl + (1 - S_Unit_Converter) / mH2O)
        End Select
    'converting anything to volume %%
    ElseIf Salinity_Unit_OUT = "volper" Then
        Select Case Salinity_Unit_IN
        Case "molper"
            S_Unit_Converter = mNaCl * Salinity_Value / RhoNaCl / (mNaCl * Salinity_Value / RhoNaCl + (1 - Salinity_Value) * mH2O / RhoH2O)
        Case "wtper"
            S_Unit_Converter = Salinity_Value / RhoNaCl / (Salinity_Value / RhoNaCl + (1 - Salinity_Value) / RhoH2O)
        Case Else
            S_Unit_Converter = Salinity_Value * mNaCl * 100 / (1000 + Salinity_Value * mNaCl * 100)
            S_Unit_Converter = S_Unit_Converter / RhoNaCl / (S_Unit_Converter / RhoNaCl + (1 - S_Unit_Converter) / RhoH2O)
        End Select
    'converting anything to Mole/kg water
    Else
        Select Case Salinity_Unit_IN
        Case "wtper"
            S_Unit_Converter = 1000 * Salinity_Value / mNaCl / (1 - Salinity_Value)
        Case "molper"
            S_Unit_Converter = Salinity_Value * mNaCl / (Salinity_Value * mNaCl + (1 - Salinity_Value) * mH2O)
            S_Unit_Converter = 1000 * S_Unit_Converter / mNaCl / (1 - S_Unit_Converter)
        Case Else
            S_Unit_Converter = Salinity_Value * RhoNaCl / (Salinity_Value * RhoNaCl + (1 - Salinity_Value) * RhoH2O)
            S_Unit_Converter = 1000 * S_Unit_Converter / mNaCl / (1 - S_Unit_Converter)
        End Select
    End If
End If
S_Unit_Converter = Round(S_Unit_Converter, 9)
End Function
