'This module contains equations from:
''Driesner, T. and C. A. Heinrich (2007). "The system H2O-NaCl. I. Correlation
''     formulae for phase relations in temperature-pressure-composition space from
''     0 to 1000 °C, 0 to 5000 bar, and 0 to 1 XNaCl " Geochimica et Cosmochimica Acta 71(20): 4880-4901.
''Driesner, T. (2007). "The system H2O-NaCl. II. Correlations for molar volume, enthalpy,
''     and isobaric heat capacity from 0 to 1000 °C, 1 to 5000 bar,
''     and 0 to 1 X-NaCl." Geochimica et Cosmochimica Acta 71(20): 4902-4919.
    

Public Function T_Star_V(xNaCl, T, P) As Double
'Tv*  for density
Dim n10 As Double, n11 As Double, n12 As Double
Dim n20 As Double, n23 As Double, n21 As Double, n22 As Double, n300 As Double, n301 As Double
Dim n302 As Double, n310 As Double, n311 As Double, n312 As Double, n30 As Double, n31 As Double
Dim n1 As Double, n2 As Double, d As Double, n_oneNaCl As Double, n_twoNaCl As Double

xH2O = 1 - xNaCl

n11 = -54.2958 - 45.7623 * Exp(-0.000944785 * P)
n21 = -2.6142 - 0.000239092 * P
n22 = 0.0356828 + 4.37235E-06 * P + 2.0566E-09 * P ^ 2

n300 = 7606640# / ((P + 472.051) ^ 2)
n301 = -50 - 86.1446 * Exp(-0.000621128 * P)
n302 = 294.318 * Exp(-0.00566735 * P)
n310 = -0.0732761 * Exp(-0.0023772 * P) - 5.2948E-05 * P
n311 = -47.2747 + 24.3653 * Exp(-0.00125533 * P)
n312 = -0.278529 - 0.00081381 * P
n30 = n300 * (Exp(n301 * xNaCl) - 1) + n302 * xNaCl
n31 = n310 * Exp(n311 * xNaCl) + n312 * xNaCl

n_oneNaCl = 330.47 + 0.942876 * P ^ 0.5 + 0.0817193 * P - 2.47556E-08 * P ^ 2 + 3.45052E-10 * P ^ 3
n10 = n_oneNaCl
n12 = -n11 - n10
n20 = 1 - n21 * n22 ^ 0.5
n_twoNaCl = -0.0370751 + 0.00237723 * P ^ 0.5 + 5.42049E-05 * P + 5.84709E-09 * P ^ 2 - 5.99373E-13 * P ^ 3
n23 = n_twoNaCl - n20 - n21 * (1 + n22) ^ 0.5
n1 = n10 + n11 * xH2O + n12 * xH2O ^ 2
n2 = n20 + n21 * (xNaCl + n22) ^ 0.5 + n23 * xNaCl
d = n30 * Exp(n31 * T)

T_Star_V = n1 + n2 * T + d

End Function

Public Function T_star_H(xNaCl, q2h, T, P) As Double
'Th* for enthalpy
Dim q10#, q11#, q12#
Dim q20#, q21#, q22#, q23#
Dim q1#, q2#, q_oneNaCl#, q_twoNaCl#, xH2O#

xH2O = 1 - xNaCl

q11 = -32.1724 + 0.0621255 * P
q21 = -1.69513 - 0.000452781 * P - 6.04279E-08 * P ^ 2
q22 = 0.0612567 + 1.88082E-05 * P
q_oneNaCl = 47.9048 - 0.00936994 * P + 6.51059E-06 * P ^ 2
q_twoNaCl = 0.241022 + 3.45087E-05 * P - 4.28356E-09 * P ^ 2

q10 = q_oneNaCl
q12 = -q11 - q10

q20 = 1 - q21 * q22 ^ 0.5
q23 = q_twoNaCl - q20 - q21 * (1 + q22) ^ 0.5

q1 = q10 + q11 * xH2O + q12 * xH2O ^ 2
q2 = q20 + q21 * (xNaCl + q22) ^ 0.5 + q23 * xNaCl

T_star_H = q1 + q2 * T + 273.15
q2h = q2
End Function

'halite melting curve
'Eq 1 from Driesner and Heinrich 2007
Public Function T_hm(P) As Double
Dim a As Double, T_tr_NaCl As Single, P_tr_NaCl As Single
a = 0.024726
T_tr_NaCl = 800.7
P_tr_NaCl = 0.0005

T_hm = T_tr_NaCl + a * (P - P_tr_NaCl)

End Function

'Halite sublimation curve,
'eq 2 Driesner and Heinrich 2007
Public Function P_Subl(T) As Double
Dim B_subl As Single, P_Triple_NaCl As Single, T_Triple_NaCl As Single
B_subl = 11806.1
T_Triple_NaCl = 800.7
P_Triple_NaCl = 0.0005

P_Subl = 10 ^ (Log10(P_Triple_NaCl) + B_subl * (1 / (T_Triple_NaCl + 273.15) - 1 / (T + 273.15)))

End Function
'Halite boiling curve,
'eq 3 Driesner and Heinrich 2007
Public Function P_Boil(T) As Double
Dim B_boil As Single, P_Triple_NaCl As Single, T_Triple_NaCl As Single
B_boil = 9418.12
T_Triple_NaCl = 800.7
P_Triple_NaCl = 0.0005

P_Boil = 10 ^ (Log10(P_Triple_NaCl) + B_boil * (1 / (T_Triple_NaCl + 273.15) - 1 / (T + 273.15)))

End Function

'eq 5 and eq 7
Public Function X_and_P_crit(T_in, Xcrit_Out, Pcrit_Out)
Dim i As Integer, j As Integer, GotIt As Boolean
Dim Sum1 As Double, PH2O_Crit As Double, TH2O_Crit As Double
Dim RhoUpUp#
Dim C(1 To 14) As Double, CA(1 To 11) As Single, d(1 To 11) As Double
Dim T As Double, P_Crit As Double, x_crit As Double, X_Crit_2 As Double
''Comparison of IAPS-84 and IAPWS-95 for the code didn't reveal any noticable difference,
''and more recent, IAPWS-95 values then were choosen for use.

' Here are values from IAPS-84, excluded from the code:
' PH2O_Crit = 220.54915: TH2O_Crit = 373.976

'IAPWS-95 are in use
PH2O_Crit = 220.64: TH2O_Crit = 373.946

C(1) = -2.36:        C(2) = 0.128534:       C(3) = -0.023707:           C(4) = 0.00320089
C(5) = -0.000138917: C(6) = 1.02789E-07: C(7) = -4.8376E-11:            C(8) = 2.36
C(9) = -0.0131417:   C(10) = 0.00298491:    C(11) = -0.000130114:       C(14) = -0.000488336

CA(1) = 1: CA(2) = 1.5:  CA(3) = 2: CA(4) = 2.5
CA(5) = 3: CA(6) = 4:    CA(7) = 5: CA(8) = 1
CA(9) = 2: CA(10) = 2.5: CA(11) = 3

For i = 8 To 11
    Sum1 = Sum1 + C(i) * (500 - TH2O_Crit) ^ CA(i)
    C(13) = C(13) + C(i) * CA(i) * (500 - TH2O_Crit) ^ (CA(i) - 1)
Next i
C(12) = PH2O_Crit + Sum1

d(1) = 8E-05:        d(2) = 1E-05:              d(3) = -1.37125E-07:     d(4) = 9.46822E-10
d(5) = -3.50549E-12: d(6) = 6.57369E-15:        d(7) = -4.89423E-18:     d(8) = 0.0777761
d(9) = 0.00027042:   d(10) = -4.244821E-07:     d(11) = 2.580872E-10

Sum1 = 0
T = T_in
If T < TH2O_Crit Then
    'eq. 5a
    For j = 1 To 7
       Sum1 = Sum1 + C(j) * (TH2O_Crit - T) ^ CA(j)
    Next j
    P_Crit = PH2O_Crit + Sum1

Else
    If T >= TH2O_Crit And T <= 500 Then
        'eq 5b
        For j = 8 To 11
            Sum1 = Sum1 + C(j) * (T - TH2O_Crit) ^ CA(j)
        Next j
        P_Crit = PH2O_Crit + Sum1
    Else
        'eq. 5c
        For j = 12 To 14
            Sum1 = Sum1 + C(j) * (T - 500) ^ (j - 12)
        Next j
        P_Crit = Sum1
    End If
End If
Sum1 = 0

If T >= TH2O_Crit And T <= 600 Then
    'eq. 7a
    For j = 1 To 7
        Sum1 = Sum1 + d(j) * (T - TH2O_Crit) ^ j
    Next j
    x_crit = Sum1
ElseIf T > 600 Then ' And T <= 1000
    'eq. 7b
    For j = 8 To 11
        Sum1 = Sum1 + d(j) * (T - 600) ^ (j - 8)
    Next j
    x_crit = Sum1
End If

Xcrit_Out = x_crit
Pcrit_Out = P_Crit
X_and_P_crit = Pcrit_Out
End Function

'Halite Liquidus
'Eq 8 from Driesner and Heinrich 2007
Public Function X_L_Sat(T, P) As Double
Dim e(0 To 5) As Double, i As Integer, TmpUnt As Double
e(0) = 0.0989944 + 3.30796E-06 * P - 4.71759E-10 * P ^ 2
e(1) = 0.00947257 - 8.6646E-06 * P + 1.69417E-09 * P ^ 2
e(2) = 0.610863 - 1.51716E-05 * P + 1.1929E-08 * P ^ 2
e(3) = -1.64994 + 0.000203441 * P - 6.46015E-08 * P ^ 2
e(4) = 3.36474 - 0.000154023 * P + 8.17048E-08 * P ^ 2
e(5) = 1

For i = 0 To 4
    e(5) = e(5) - e(i)
Next i

TmpUnt = T / T_hm(P)
For i = 0 To 5
    X_L_Sat = X_L_Sat + e(i) * TmpUnt ^ i
Next i

If X_L_Sat > 1 Then X_L_Sat = 1

End Function

'eq 9  from Driesner and Heinrich 2007
Public Function X_V_Sat(T, P) As Double
Dim j(0 To 3) As Double, k(0 To 15) As Double, Log_K_supScr As Double, Log_K_Line As Double, P_Line As Double, P_NaCl As Double
Dim ii%, X_L_St#, P_Crit#, Tmp#

k(0) = -0.235694:       k(1) = -0.188838:           k(2) = 0.004:               k(3) = 0.0552466:       k(4) = 0.66918
k(5) = 396.848:         k(6) = 45:                  k(7) = -3.2719E-07: k(8) = 141.699:          k(9) = -0.292631
k(10) = -0.00139991:    k(11) = 1.95965E-06: k(12) = -7.3653E-10: k(13) = 0.904411: k(14) = 0.000769766
k(15) = -1.18658E-06

j(0) = k(0) + k(1) * Exp(-k(2) * T)
j(1) = k(4) + (k(3) - k(4)) / (1 + Exp((T - k(5)) / k(6))) + k(7) * (T + k(8)) ^ 2
j(2) = 0
For ii = 0 To 3
    j(2) = j(2) + k(ii + 9) * T ^ ii
Next ii
j(3) = 0
For ii = 0 To 2
    j(3) = j(3) + k(ii + 13) * T ^ ii
Next ii

If T > 800.7 Then
    P_NaCl = P_Boil(T)
Else
    P_NaCl = P_Subl(T)
End If

x_l = X_L_Sat(T, P)
X_and_P_crit T, 0, P_Crit
'Figure 8A is off by a factor of 10 for salinity, while this code reproduce values made by SoWat
P_Line = (P - P_NaCl) / (P_Crit - P_NaCl)
P_Line = 1 - P_Line
Tmp = Log10(X_L_Sat(T, P_NaCl))
Log_K_Line = 1 + j(0) * P_Line ^ j(1) + j(2) * P_Line + j(3) * P_Line ^ 2 - (1 + j(0) + j(2) + j(3)) * P_Line ^ 3
Log_K_supScr = Log_K_Line * (Log10(P_NaCl / P_Crit) - Tmp) + Tmp

Tmp = Log_K_supScr - Log10(P_NaCl / P)
X_V_Sat = x_l / 10 ^ Tmp


End Function

'Vapor-Liquid-Halite coexistence
'eq 10 from Driesner and Heinrigch 2007
Public Function P_VLH(T) As Double
Dim i As Integer, f(0 To 10) As Double, P_tr_NaCl As Single, T_tr_NaCl As Single
T_tr_NaCl = 800.7
P_tr_NaCl = 0.0005

f(0) = 0.00464:   f(1) = 5E-07:      f(2) = 16.9078:     f(3) = -269.148: f(4) = 7632.04: f(5) = -49563.6
f(6) = 233119#: f(7) = -513556#: f(8) = 549708#: f(9) = -284628#: f(10) = P_tr_NaCl

For i = 0 To 10
    If i <> 10 Then f(10) = f(10) - f(i)
    P_VLH = P_VLH + f(i) * (T / T_tr_NaCl) ^ i
Next i

End Function

'Vapor-liquid field, separated from liquid phase
'Eq 11 from Driesner and Heinrich 2007
Public Function X_VL_Liq(T, P) As Double
Dim d#(1 To 11), XN_Crit#, TmpUnit#, TmpUnit2#, TmpUnit3#, ii%, TH2O_Crit#
Dim h#(1 To 11), G0#, G1#, G2#, P_Crit#

h(1) = 0.00168486:        h(2) = 0.000219379:   h(3) = 438.58:         h(4) = 18.4508
h(5) = -5.6765E-10: h(6) = 6.73704E-06: h(7) = 1.44951E-07: h(8) = 384.904
h(9) = 7.07477:           h(10) = 6.06896E-05: h(11) = 0.00762859

G1 = h(2) + (h(1) - h(2)) / (1 + Exp((T - h(3)) / h(4))) + h(5) * T ^ 2
G2 = h(7) + (h(6) - h(7)) / (1 + Exp((T - h(8)) / h(9))) + h(10) * Exp(-h(11) * T)
X_and_P_crit T, XN_Crit, P_Crit

If T < 800.7 Then
    TmpUnit = P_VLH(T)
    TmpUnit2 = X_L_Sat(T, TmpUnit)
Else
    TmpUnit = P_Boil(T)
    TmpUnit2 = 1
End If
'If P_Crit < P And Abs(P_Crit - P) < 0.0000000001 Then P = P_Crit '0.0000001 Then P = P_Crit
If T < 373.946 Then
    TmpUnit3 = P_H2O_Boiling_Curve(T)
    G0 = (TmpUnit2 + G1 * (TmpUnit - TmpUnit3) + G2 * ((P_Crit - TmpUnit3) ^ 2 - (P_Crit - TmpUnit) ^ 2)) / ((P_Crit - TmpUnit) ^ 0.5 - (P_Crit - TmpUnit3) ^ 0.5)
    X_VL_Liq = G0 * (P_Crit - P) ^ 0.5 - G0 * (P_Crit - TmpUnit3) ^ 0.5 - G1 * (P_Crit - TmpUnit3) - G2 * (P_Crit - TmpUnit3) ^ 2 + G1 * (P_Crit - P) + G2 * (P_Crit - P) ^ 2
Else
    G0 = (TmpUnit2 - XN_Crit - G1 * (P_Crit - TmpUnit) - G2 * (P_Crit - TmpUnit) ^ 2) / (P_Crit - TmpUnit) ^ 0.5
    X_VL_Liq = XN_Crit + G0 * (P_Crit - P) ^ 0.5 + G1 * (P_Crit - P) + G2 * (P_Crit - P) ^ 2
End If
End Function

'Vapor-liquid field, separated from vapor phase
'Eq 13 Eq 14 Eq 15 Eq 16 Eq 17
Public Function X_VL_Vap(T, P) As Double
Dim j#(0 To 3), k#(0 To 15), Log_K_supScr#, Log_K_Line#, P_Line#, P_NaCl#
Dim ii%, X_VL_Lq#, P_Crit#, Tmp#, x#
'On Error Resume Next
k(0) = -0.235694:       k(1) = -0.188838:           k(2) = 0.004:               k(3) = 0.0552466:       k(4) = 0.66918
k(5) = 396.848:         k(6) = 45:                  k(7) = -3.2719E-07: k(8) = 141.699:          k(9) = -0.292631
k(10) = -0.00139991:    k(11) = 1.95965E-06: k(12) = -7.3653E-10: k(13) = 0.904411: k(14) = 0.000769766
k(15) = -1.18658E-06

j(0) = k(0) + k(1) * Exp(-k(2) * T)
j(1) = k(4) + (k(3) - k(4)) / (1 + Exp((T - k(5)) / k(6))) + k(7) * (T + k(8)) ^ 2
j(2) = 0
For ii = 0 To 3
    j(2) = j(2) + k(ii + 9) * T ^ ii
Next ii
j(3) = 0
For ii = 0 To 2
    j(3) = j(3) + k(ii + 13) * T ^ ii
Next ii

If T >= 800.7 Then
    P_NaCl = P_Boil(T)
Else
    P_NaCl = P_Subl(T)
End If

X_VL_Lq = X_VL_Liq(T, P)

X_and_P_crit T, x, P_Crit


P_Line = (P - P_NaCl) / (P_Crit - P_NaCl)
'YK: IF added to avoid an exception to avoid (negative number)^0.5
'CURRENTLY DISABLED
'If Abs(1 - P_Line) < 0.000000001 And P_Line > 1 Then P_Line = 1

Log_K_Line = 1 + j(0) * (1 - P_Line) ^ j(1) + j(2) * (1 - P_Line) + j(3) * (1 - P_Line) ^ 2 - (1 + j(0) + j(2) + j(3)) * (1 - P_Line) ^ 3

Tmp = X_L_Sat(T, P_NaCl)
Log_K_supScr = Log_K_Line * (Log10(P_NaCl / P_Crit) - Log10(Tmp)) + Log10(Tmp)

X_VL_Vap = X_VL_Lq / 10 ^ Log_K_supScr * P_NaCl / P
If P <= P_VLH(T) Then X_VL_Vap = X_L_Sat(T, P) / X_VL_Lq * X_VL_Vap

End Function

' function for low TP region, eq. 17 from Driesner 2007
Public Function V_Extrapol(x_in, T_in, P_in) As Double
Dim o0#, o1#, o2#, T_inv#, mH2O#, mNaCl#, T#, P#, p2#, Rho#, Rho2#, v#, V2#, xNaCl#
Dim o3#, o4#, o5#, Vsat#, Vwat#

If x_in = 0 Then Exit Function

Dim V1000#, dVdP390#

mH2O = 18.015268
mNaCl = 58.4428

T = T_in
xNaCl = x_in
P = P_in

v = X_L_Sat(T, P)
 
If P <= Water_Boiling_Curve(T) And T <= 200 And v - x_in < 0.01 Then
    T = T_Star_V(xNaCl, T, P)
    Vsat = mH2O / Rho_Water_Liq_sat(T) * 1000
    Vwat = mH2O / Rho_Water(T, P) * 1000
    
    If Vsat < Vwat Then

        o2 = 2.0125E-07 + 3.29977E-09 * Exp(-4.31279 * Log10(P)) - 1.17748E-07 * Log10(P) + 7.58009E-08 * (Log10(P)) ^ 2
        v = mH2O / Rho_Water_Liq_sat(T) * 1000
        V2 = mH2O / Rho_Water_Liq_sat(T - 0.005) * 1000
        o1 = (v - V2) / 0.005 - 3 * o2 * T ^ 2
        o0 = v - o1 * T - o2 * T ^ 3
        
        V_Extrapol = o0 + (o2 * T ^ 3) + (o1 * T)
    End If
    
ElseIf P <= 350 And T >= 600 Then
    v = X_VL_Liq(T, P)
    If Round_Down(xNaCl, 5) >= Round_Down(v, 5) Then
'        If P >= New_LV_P(T, P, XNaCl, True) Then
        
            V1000 = (mH2O * (1 - xNaCl) + mNaCl * xNaCl) / Rh_Br_for_V_extr(xNaCl, T, 1000) * 1000
            v = (mH2O * (1 - xNaCl) + mNaCl * xNaCl) / Rh_Br_for_V_extr(xNaCl, T, 390.147) * 1000
            V2 = (mH2O * (1 - xNaCl) + mNaCl * xNaCl) / Rh_Br_for_V_extr(xNaCl, T, 390.137) * 1000
            
            dVdP390 = (v - V2) / (0.01)
            
            o4 = (v - V1000 + dVdP390 * 1609.853) / (LogExp(1390.147 / 2000) - 2390.147 / 1390.147)
            o3 = v - o4 * LogExp(1390.147) - 390.147 * dVdP390 + 390.147 / 1390.147 * o4
            o5 = dVdP390 - o4 / (1390.147)
        
        V_Extrapol = o3 + o4 * LogExp(P + 1000) + o5 * P
'        End If
    End If
Else
    V_Extrapol = 0
End If

End Function
'local function to estimate density of vapor phase neccesery for V_extrapol funtion
Private Function Rh_Br_for_V_extr(xNaCl_frac, T_in_C, P_in_Bar) As Double
Dim T#, P#, T_Star#, V_water#

Dim mH2O#, mNaCl#

mH2O = 18.015268
mNaCl = 58.4428

T = T_in_C
P = P_in_Bar
T_Star = T_Star_V(xNaCl_frac, T, P)

V_water = mH2O / Rho_Water(T_Star, P) * 1000#

Rh_Br_for_V_extr = (mH2O * (1 - xNaCl_frac) + mNaCl * xNaCl_frac) / V_water * 1000#
End Function
Public Function NaCl_Rho_Solid(T_in_C, P_in_Bar)
Dim l#(0 To 5), Rho_Zero#, lParam#, T#, P#, i%
T = T_in_C
P = P_in_Bar
l(0) = 2170.4
l(1) = -0.24599
l(2) = -9.5797E-05
l(3) = 0.005727
l(4) = 0.002715
l(5) = 733.4

For i = 0 To 2
    Rho_Zero = Rho_Zero + l(i) * T ^ i
Next i

lParam = l(3) + l(4) * Exp(T / l(5))

NaCl_Rho_Solid = Rho_Zero + lParam * P
End Function

Public Function NaCl_Rho_Liq(T_in_C, P_in_Bar)
Dim m#(0 To 5), Rho_Zero#, KNaCl#, T#, P#
T = T_in_C
P = P_in_Bar
m(0) = 58443
m(1) = 23.772
m(2) = 0.018639
m(3) = -1.9687E-06
m(4) = -1.5259E-05
m(5) = 5.5058E-08

KNaCl = m(4) + m(5) * T
Rho_Zero = m(0) / (m(1) + m(2) * T + m(3) * T ^ 2)
NaCl_Rho_Liq = Rho_Zero / (1 - 0.1 * LogExp(1 + 10 * P * KNaCl))

End Function

Public Function NaCl_specific_enthalpy#(T_in_C#, P_in_Bar#)
Dim k(), r_rec(), t_trip#, T#, t_t#
T = T_in_C
t_t = T - 800.7

k = Array(226713, 44.6652, -7.41999E-05)
r_rec = Array(1148.81, 0.275774, 8.8103E-05, -0.0017099 - 3.82734E-06 * T / 2 - 8.65455E-09 * T ^ 2 / 3, _
     5.29063E-08 - 9.63084E-11 * T / 2 + 6.50745E-13 * T ^ 2 / 3)
     
NaCl_specific_enthalpy = k(0) + r_rec(0) * t_t + r_rec(1) * t_t ^ 2 + r_rec(2) * t_t ^ 3 + (k(1) + r_rec(3) * T) * P_in_Bar + (k(2) + r_rec(4) * T) * P_in_Bar ^ 2

End Function

Public Function NaCl_isobaric_heat_capacity#(T_in_C#, P_in_Bar#)
Dim r(), t_trip#, T#, t_t#
T = T_in_C
t_t = T - 800.7
r = Array(1148.81, 0.275774, 8.8103E-05, -0.0017099 - 3.82734E-06 * T - 8.65455E-09 * T ^ 2, _
     5.29063E-08 - 9.63084E-11 * T + 6.50745E-13 * T ^ 2)

NaCl_isobaric_heat_capacity = r(0) + 2 * r(1) * t_t + 3 * r(2) * t_t ^ 2 + r(3) * P_in_Bar + r(4) * P_in_Bar ^ 2

End Function
