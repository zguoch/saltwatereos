'Viscosity calculations are based on work published by Klyukin Lowell Bodnar 2017, Phase Equilibria

'this function verifies PT limits of IAPWS-95, which are imposed on H2O-NaCl viscosity, as well as PT limits for H2O-NaCl
'T in C, P in bar
Public Function Visc_Inc_data_Fit(T, P) As Boolean
If T > -22.15 And P >= 1 And P > -94.765 * T + 0.947 And P < 161.22 * T + 5671.1 And T <= 900 And P <= 10000 Then
            
    If T < 100 Then 'And T >= 0
        Visc_Inc_data_Fit = True
        Exit Function
    ElseIf T >= 100 And T < 160 And P < 5000 Then
        Visc_Inc_data_Fit = True
        Exit Function
    ElseIf T >= 160 And T < 600 And P < 3500 Then
        Visc_Inc_data_Fit = True
        Exit Function
    ElseIf T >= 600 And P < 3000 Then
        Visc_Inc_data_Fit = True
        Exit Function
    Else
        Visc_Inc_data_Fit = False
        Exit Function
    End If
Else
    Visc_Inc_data_Fit = False
End If

End Function

Public Function T_star_Mu(wt_Percent, T_in_C) As Double
'The equation to calculate T_star appears different compared to Klyukin et al 2017
'but it is simply rearranged equation
Dim z0#, b#, C#, d#, e#, f#, Wt_frac#
On Error GoTo Ffail
Wt_frac = wt_Percent / 100

z0 = 0
b = -35.9858
C = 0.80017
d = 1E-06
e = -0.05239
f = 1.32936

T_star_Mu = b * Wt_frac ^ C + T_in_C * (1 - d * T_in_C ^ e - f * Wt_frac ^ C * T_in_C ^ e)
T_star_Mu = Round(T_star_Mu + 273.15, 3)
Exit Function

Ffail:
T_star_Mu = -25
End Function

'MODEL LIMITATION: T_STAR HAVE TO BE <=900 AND ~>40 C
'PRESSURE LIMITS CANCELED BY ASSUMPTION THAT LOW-P REGION WITH P BELOW P H2O BOILING CURVE
'CAN BE EXTRAPOLATED UP TO H2O BOILING BECAUSE
'1) IT'S A LOW PRESSURE DIFFERENCE AND
'2) INCREASE IN PRESSURE DOESN'T HAVE ANY AFFECT ON VISCOSITY

'THIS EXTRAPOLATION PLANNED TO ACHIEVE BY DP/DT SLOPE FOR ISOVISCOSITY LINES
'EXTRAPOLATED LINEARY DOWN TO P_VLH
Public Function Viscosity_H2O_NaCl(S_in_Wt_frac, T_in_C, P_in_Bar, Rho_in_KGM) As Double
Dim T_Star#, Rho_W#, Pres#, Tmp#, S_in_Wt_Per#

Pres = P_in_Bar
S_in_Wt_Per = S_in_Wt_frac * 100

T_Star = T_star_Mu(S_in_Wt_Per, T_in_C)
If T_Star > 1173.15 Then
    Viscosity_H2O_NaCl = 0
Else
    Rho_W = Rho_Water(T_Star - 273.15, Pres)
    Viscosity_H2O_NaCl = Water_Viscosity_calc(T_Star, Rho_W)

End If
End Function
