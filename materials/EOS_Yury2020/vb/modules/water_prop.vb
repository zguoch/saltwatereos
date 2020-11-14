'water constants and equation of state from W. WAGNER AND A. PRUSS
'J. Phys. Chem. Ref. Data, Vol. 31, No. 2, 2002
'
'Huber et al 2009, New International Formulation for the Viscosity of H2O
'Valid up to 1173 K and 1000 MPa or 10 kBar

'WAGNER 2002
'units - MPa
Public Function Water_Pressure_calc#(T_K, Rho)

Dim R_constant As Double, Delta_Rho As Variant, Tau As Variant
R_constant = 0.46151805
Delta_Rho = (Rho / 322): Tau = (647.096 / T_K)

Water_Pressure_calc = ((1 + Delta_Rho * PhiR_Delta(Delta_Rho, Tau)) * Rho * R_constant * T_K) / 1000
End Function

'WAGNER 2002
'Units - Joule per kg
Public Function Water_Enthalpy_calc#(T_in_K, Rho_in_kgm3)
Dim R_constant As Double, Delta_Rho As Variant, T#, Tau As Variant
R_constant = 0.46151805
Delta_Rho = (Rho_in_kgm3 / 322): Tau = (647.096 / T_in_K)

Water_Enthalpy_calc = (1 + Tau * (Phi0_Tau(Delta_Rho, Tau) + PhiR_Tau(Delta_Rho, Tau)) + Delta_Rho * PhiR_Delta(Delta_Rho, Tau)) * R_constant * T_in_K
End Function

'WAGNER 2002
Public Function Water_Isobaric_Heat_capacity_calc#(T, Rho)
Dim R_constant As Double, Delta_Rho As Variant, Tau As Variant, PhirDelta#
R_constant = 0.46151805
Delta_Rho = (Rho / 322): Tau = (647.096 / T)
PhirDelta = PhiR_Delta(Delta_Rho, Tau)

Water_Isobaric_Heat_capacity_calc = Water_Isochoric_Heat_capacity_calc(T, Rho) + R_constant * (1 + Delta_Rho * PhirDelta - Delta_Rho * Tau * PhiR_DeltaTau(Delta_Rho, Tau)) ^ 2 / (1 + 2 * Delta_Rho * PhirDelta + Delta_Rho ^ 2 * PhiR_DeltaDelta(Delta_Rho, Tau))
End Function

'WAGNER 2002
Public Function Water_Isochoric_Heat_capacity_calc#(T, Rho)
Dim R_constant As Double, Delta_Rho As Variant, Tau As Variant
R_constant = 0.46151805
Delta_Rho = (Rho / 322): Tau = (647.096 / T)

Water_Isochoric_Heat_capacity_calc = (-Tau ^ 2 * (Phi0_TauTau(Delta_Rho, Tau) + PhiR_TauTau(Delta_Rho, Tau))) * R_constant
End Function

'WAGNER 2002
Private Function Phi0_TauTau#(Delta_Rho, Tau)
Dim i As Integer
'constants from table 6.1
Dim n0(3 To 8) As Double, gamma0(4 To 8) As Double, Sum1 As Double

'constants from table 6.1
n0(3) = 3.00632:        n0(4) = 0.012436
n0(5) = 0.97315:        n0(6) = 1.2795:       n0(7) = 0.96956: n0(8) = 0.24873

gamma0(4) = 1.28728967: gamma0(5) = 3.53734222: gamma0(6) = 7.74073708
gamma0(7) = 9.24437796: gamma0(8) = 27.5075105

For i = 4 To 8
    Sum1 = Sum1 + n0(i) * gamma0(i) ^ 2 * Exp(-gamma0(i) * Tau) * (1 - Exp(-gamma0(i) * Tau)) ^ -2
Next i
Phi0_TauTau = -n0(3) / Tau ^ 2 - Sum1

End Function
Private Sub LoadConstants(C(), d(), T(), n(), Alpha(), Epsilon(), Beta(), Gamma(), a(), b(), B_Caps(), A_caps(), c_Caps(), D_Caps)

Dim i As Integer
ReDim C(1 To 54), d(1 To 54)
ReDim T(1 To 54), n(1 To 56)
ReDim Alpha(52 To 54), Epsilon(52 To 54)
ReDim Beta(52 To 56), Gamma(52 To 54)
ReDim a(55 To 56), b(55 To 56), B_Caps(55 To 56), A_caps(55 To 56)
ReDim c_Caps(55 To 56), D_Caps(55 To 56)

'table 6.2
For i = 1 To 54
    C(i) = 0
Next i

For i = 8 To 22
C(i) = 1
Next i

For i = 23 To 42
C(i) = 2
Next i

For i = 43 To 46
C(i) = 3
Next i

C(47) = 4

For i = 48 To 51
C(i) = 6
Next i

d(1) = 1:   d(2) = 1:  d(3) = 1:  d(4) = 2:  d(5) = 2:   d(6) = 3:   d(7) = 4
d(8) = 1:   d(9) = 1:  d(10) = 1: d(11) = 2: d(12) = 2:  d(13) = 3:  d(14) = 4
d(15) = 4:  d(16) = 5: d(17) = 7: d(18) = 9: d(19) = 10: d(20) = 11: d(21) = 13
d(22) = 15: d(23) = 1: d(24) = 2: d(25) = 2: d(26) = 2:  d(27) = 3:  d(28) = 4
d(29) = 4:  d(30) = 4: d(31) = 5: d(32) = 6: d(33) = 6:  d(34) = 7:  d(35) = 9
d(36) = 9:  d(37) = 9: d(38) = 9: d(39) = 9: d(40) = 10: d(41) = 10: d(42) = 12
d(43) = 3:  d(44) = 4: d(45) = 4: d(46) = 5: d(47) = 14: d(48) = 3:  d(49) = 6
d(50) = 6:  d(51) = 6: d(52) = 3: d(53) = 3: d(54) = 3

T(1) = -0.5: T(2) = 0.875: T(3) = 1:   T(4) = 0.5: T(5) = 0.75: T(6) = 0.375
T(7) = 1:    T(8) = 4:     T(9) = 6:   T(10) = 12: T(11) = 1:   T(12) = 5
T(13) = 4:   T(14) = 2:    T(15) = 13: T(16) = 9:  T(17) = 3:   T(18) = 4
T(19) = 11:  T(20) = 4:    T(21) = 13: T(22) = 1:  T(23) = 7:   T(24) = 1
T(25) = 9:   T(26) = 10:   T(27) = 10: T(28) = 3:  T(29) = 7:   T(30) = 10
T(31) = 10:  T(32) = 6:    T(33) = 10: T(34) = 10: T(35) = 1:   T(36) = 2
T(37) = 3:   T(38) = 4:    T(39) = 8:  T(40) = 6:  T(41) = 9:   T(42) = 8
T(43) = 16:  T(44) = 22:   T(45) = 23: T(46) = 23: T(47) = 10:  T(48) = 50
T(49) = 44:  T(50) = 46:   T(51) = 50: T(52) = 0:  T(53) = 1:   T(54) = 4

n(1) = 0.012533547935523:      n(2) = 7.8957634722828
n(3) = -8.7803203303561:       n(4) = 0.31802509345418
n(5) = -0.26145533859358:      n(6) = -0.0078199751687981
n(7) = 0.0088089493102134:     n(8) = -0.66856572307965
n(9) = 0.20433810950965:       n(10) = -6.6212605039687E-05
n(11) = -0.19232721156002:     n(12) = -0.25709043003438
n(13) = 0.16074868486251:      n(14) = -0.040092828925807
n(15) = 3.9343422603254E-07:   n(16) = -7.5941377088144E-06
n(17) = 0.00056250979351888:   n(18) = -1.5608652257135E-05
n(19) = 1.1537996422951E-09:   n(20) = 3.6582165144204E-07
n(21) = -1.3251180074668E-12: n(22) = -6.2639586912454E-10
n(23) = -0.10793600908932:     n(24) = 0.017611491008752
n(25) = 0.22132295167546:      n(26) = -0.40247669763528
n(27) = 0.58083399985759:      n(28) = 0.0049969146990806
n(29) = -0.031358700712549:    n(30) = -0.74315929710341
n(31) = 0.4780732991548:       n(32) = 0.020527940895948
n(33) = -0.13636435110343:     n(34) = 0.014180634400617
n(35) = 0.0083326504880713:    n(36) = -0.029052336009585
n(37) = 0.038615085574206:     n(38) = -0.020393486513704
n(39) = -0.0016554050063734:   n(40) = 0.0019955571979541
n(41) = 0.00015870308324157:   n(42) = -1.638856834253E-05
n(43) = 0.043613615723811:     n(44) = 0.034994005463765
n(45) = -0.076788197844621:    n(46) = 0.022446277332006
n(47) = -6.2689710414685E-05:  n(48) = -5.5711118565645E-10
n(49) = -0.19905718354408:     n(50) = 0.31777497330738
n(51) = -0.11841182425981:     n(52) = -31.306260323435
n(53) = 31.546140237781:       n(54) = -2521.3154341695
n(55) = -0.14874640856724:     n(56) = 0.31806110878444

Alpha(52) = 20: Alpha(53) = 20: Alpha(54) = 20

Beta(52) = 150: Beta(53) = 150: Beta(54) = 250: Beta(55) = 0.3: Beta(56) = 0.3

Gamma(52) = 1.21: Gamma(53) = 1.21: Gamma(54) = 1.25

Epsilon(52) = 1: Epsilon(53) = 1: Epsilon(54) = 1

a(55) = 3.5:      a(56) = 3.5:      b(55) = 0.85:      b(56) = 0.95
B_Caps(55) = 0.2: B_Caps(56) = 0.2: c_Caps(55) = 28:   c_Caps(56) = 32
D_Caps(55) = 700: D_Caps(56) = 800: A_caps(55) = 0.32: A_caps(56) = 0.32
'END OF CONSTANT DECLARING.

End Sub

'WAGNER 2002 eq from table 6.5
'Private Function PhiR_TauTau( Delta_Rho,  Tau) As Variant
Private Function PhiR_TauTau#(Delta_Rho, Tau)
Dim i As Integer
Dim C(), d()
Dim T(), n()
Dim Alpha(), Epsilon()
Dim Beta(), Gamma()
Dim a(), b(), B_Caps(), A_caps()
Dim c_Caps(), D_Caps()
'declaring of  components for F0+Fr and experimental PT spaces
Dim Delta#, Theta#, Psi#
'temporary double units
Dim Sum1 As Double, Sum2 As Double, Sum3 As Double, Sum4 As Double, d2DeltaBIDTauTau#, dDeltaBIdTau#, dPsidTau#, d2PsidTauTau#

LoadConstants C(), d(), T(), n(), Alpha(), Epsilon(), Beta(), Gamma(), a(), b(), B_Caps(), A_caps(), c_Caps(), D_Caps()

For i = 1 To 7
    Sum1 = Sum1 + n(i) * T(i) * (T(i) - 1) * Delta_Rho ^ d(i) * Tau ^ (T(i) - 2)
Next i
For i = 8 To 51
    Sum2 = Sum2 + n(i) * T(i) * (T(i) - 1) * Delta_Rho ^ d(i) * Tau ^ (T(i) - 2) * Exp(-Delta_Rho ^ C(i))
Next i
For i = 52 To 54
    Sum3 = Sum3 + n(i) * Delta_Rho ^ d(i) * Tau ^ T(i) * Exp(-Alpha(i) * (Delta_Rho - Epsilon(i)) ^ 2 - Beta(i) * (Tau - Gamma(i)) ^ 2) * ((T(i) / Tau - 2 * Beta(i) * (Tau - Gamma(i))) ^ 2 - T(i) / Tau ^ 2 - 2 * Beta(i))
Next i
For i = 55 To 56
    
    Theta = (1 - Tau) + A_caps(i) * ((Delta_Rho - 1) ^ 2) ^ (1 / (2 * Beta(i)))
    
    Delta = Theta ^ 2 + B_Caps(i) * ((Delta_Rho - 1) ^ 2) ^ a(i)
    
    Psi = Exp(-c_Caps(i) * (Delta_Rho - 1) ^ 2 - D_Caps(i) * (Tau - 1) ^ 2)
    
    d2DeltaBIDTauTau = 2 * b(i) * Delta ^ (b(i) - 1) + 4 * Theta ^ 2 * b(i) * (b(i) - 1) * Delta ^ (b(i) - 2)
    dDeltaBIdTau = -2 * Theta * b(i) * Delta ^ (b(i) - 1)
    dPsidTau = -2 * D_Caps(i) * (Tau - 1) * Psi
    d2PsidTauTau = (2 * D_Caps(i) * (Delta_Rho - 1) ^ 2 - 1) * 2 * D_Caps(i) * Psi
    
    Sum4 = Sum4 + n(i) * Delta_Rho * (d2DeltaBIDTauTau * Psi + 2 * dDeltaBIdTau * dPsidTau + Delta ^ b(i) * d2PsidTauTau)

Next i
PhiR_TauTau = (Sum1 + Sum2 + Sum3 + Sum4)

End Function
'WAGNER 2002 eq from table 6.5
Public Function PhiR_DeltaDelta#(Delta_Rho, Tau)
Dim i As Integer
Dim C(), d()
Dim T(), n()
Dim Alpha(), Epsilon()
Dim Beta(), Gamma()
Dim a(), b(), B_Caps(), A_caps()
Dim c_Caps(), D_Caps()
'declaring of  components for F0+Fr and experimental PT spaces
Dim Delta#, Theta#, Delta2#, Theta2#, Psi#
'temporary double units
Dim Sum1 As Double, Sum2 As Double, Sum3 As Double, Sum4 As Double, TmpUnit#, dPsidTau#, dDeltaddelta#, d2PsidTauDelta#
Dim dDeltaBIdDelta#, dDeltaBIdTau#, dPsidDelta#, d2DeltaBIdRhoTau#, d2Psiddeltadelta#, d2DeltaBIddeltadelta#, d2Deltaddeltadelta#
'table 6.2
LoadConstants C(), d(), T(), n(), Alpha(), Epsilon(), Beta(), Gamma(), a(), b(), B_Caps(), A_caps(), c_Caps(), D_Caps()

For i = 1 To 7
    Sum1 = Sum1 + n(i) * d(i) * (d(i) - 1) * Delta_Rho ^ (d(i) - 2) * Tau ^ T(i)
Next i
For i = 8 To 51
    Sum2 = Sum2 + n(i) * Exp(-Delta_Rho ^ C(i)) * (Delta_Rho ^ (d(i) - 2) * Tau ^ T(i) * ((d(i) - C(i) * Delta_Rho ^ C(i)) * (d(i) - 1 - C(i) * Delta_Rho ^ C(i)) - C(i) ^ 2 * Delta_Rho ^ C(i)))
Next i
For i = 52 To 54
    Sum3 = Sum3 + n(i) * Tau ^ T(i) * Exp(-Alpha(i) * (Delta_Rho - Epsilon(i)) ^ 2 - Beta(i) * (Tau - Gamma(i)) ^ 2) * (-2 * Alpha(i) * Delta_Rho ^ d(i) + 4 * Alpha(i) ^ 2 * Delta_Rho ^ d(i) * (Delta_Rho - Epsilon(i)) ^ 2 - 4 * d(i) * Alpha(i) * Delta_Rho ^ (d(i) - 1) * (Delta_Rho - Epsilon(i)) + d(i) * (d(i) - 1) * Delta_Rho ^ (d(i) - 2))
Next i
For i = 55 To 56
    
    Theta = (1 - Tau) + A_caps(i) * ((Delta_Rho - 1) ^ 2) ^ (1 / (2 * Beta(i)))
    
    Delta = Theta ^ 2 + B_Caps(i) * ((Delta_Rho - 1) ^ 2) ^ a(i)
    
    Psi = Exp(-c_Caps(i) * (Delta_Rho - 1) ^ 2 - D_Caps(i) * (Tau - 1) ^ 2)
    
    dPsidTau = -2 * D_Caps(i) * (Tau - 1) * Psi
    d2PsidTauDelta = 4 * c_Caps(i) * Psi * D_Caps(i) * (Delta_Rho - 1) * (Tau - 1)
    dDeltaddelta = (Delta_Rho - 1) * (A_caps(i) * Theta * 2 / Beta(i) * ((Delta_Rho - 1) ^ 2) ^ (1 / (2 * Beta(i)) - 1) + 2 * B_Caps(i) * a(i) * ((Delta_Rho - 1) ^ 2) ^ (a(i) - 1))
    dDeltaBIdDelta = b(i) * Delta ^ (b(i) - 1) * dDeltaddelta
    dDeltaBIdTau = -2 * Theta * b(i) * Delta ^ (b(i) - 1)
    dPsidDelta = -2 * c_Caps(i) * (Delta_Rho - 1) * Psi
    d2DeltaBIdRhoTau = -A_caps(i) * b(i) * 2 / Beta(i) * Delta ^ (b(i) - 1) * (Delta_Rho - 1) * ((Delta_Rho - 1) ^ 2) ^ (1 / (2 * Beta(i)) - 1) - 2 * Theta * b(i) * (b(i) - 1) * Delta ^ (b(i) - 2) * dDeltaddelta
    d2Psiddeltadelta = (2 * c_Caps(i) * (Delta_Rho - 1) ^ 2 - 1) * 2 * c_Caps(i) * Psi
    d2Deltaddeltadelta = 1 / (Delta_Rho - 1) * dDeltaddelta + (Delta_Rho - 1) ^ 2 * (4 * B_Caps(i) * a(i) * (a(i) - 1) * ((Delta_Rho - 1) ^ 2) ^ (a(i) - 2) + 2 * A_caps(i) ^ 2 * Beta(i) ^ -2 * (((Delta_Rho - 1) ^ 2) ^ (1 / 2 / Beta(i) - 1)) ^ 2 + A_caps(i) * Theta * 4 / Beta(i) * (1 / 2 / Beta(i) - 1) * ((Delta_Rho - 1) ^ 2) ^ (1 / 2 / Beta(i) - 2))
    d2DeltaBIddeltadelta = b(i) * (Delta ^ (b(i) - 1) * d2Deltaddeltadelta + (b(i) - 1) * Delta ^ (b(i) - 2) * dDeltaddelta ^ 2)
    
    Sum4 = Sum4 + n(i) * (Delta ^ b(i) * (2 * dPsidDelta + Delta_Rho * d2Psiddeltadelta) + 2 * dDeltaBIdDelta * (Psi + Delta_Rho * dPsidDelta) + d2DeltaBIddeltadelta * Psi * Delta_Rho)
    

Next i
PhiR_DeltaDelta = (Sum1 + Sum2 + Sum3 + Sum4)

End Function

'WAGNER 2002 eq from table 6.5
Private Function PhiR_DeltaTau#(Delta_Rho, Tau)
Dim i As Integer
Dim C(), d()
Dim T(), n()
Dim Alpha(), Epsilon()
Dim Beta(), Gamma()
Dim a(), b(), B_Caps(), A_caps()
Dim c_Caps(), D_Caps()
'declaring of  components for F0+Fr and experimental PT spaces
Dim Delta#, Theta#, Delta2#, Theta2#, Psi#
'temporary double units
Dim Sum1 As Double, Sum2 As Double, Sum3 As Double, Sum4 As Double, TmpUnit#, dPsidTau#, dDeltaddelta#
Dim d2PsidTauDelta#, dDeltaBIdDelta#, dDeltaBIdTau#, dPsidDelta#, d2DeltaBIdRhoTau#

LoadConstants C(), d(), T(), n(), Alpha(), Epsilon(), Beta(), Gamma(), a(), b(), B_Caps(), A_caps(), c_Caps(), D_Caps()
'END OF CONSTANT DECLARING

For i = 1 To 7
    Sum1 = Sum1 + n(i) * d(i) * T(i) * Delta_Rho ^ (d(i) - 1) * Tau ^ (T(i) - 1)
Next i
For i = 8 To 51
    Sum2 = Sum2 + n(i) * T(i) * Delta_Rho ^ (d(i) - 1) * Tau ^ (T(i) - 1) * (d(i) - C(i) * Delta_Rho ^ C(i)) * Exp(-Delta_Rho ^ C(i))
Next i
For i = 52 To 54
    Sum3 = Sum3 + n(i) * Delta_Rho ^ d(i) * Tau ^ T(i) * Exp(-Alpha(i) * (Delta_Rho - Epsilon(i)) ^ 2 - Beta(i) * (Tau - Gamma(i)) ^ 2) * (d(i) / Delta_Rho - 2 * Alpha(i) * (Delta_Rho - Epsilon(i))) * (T(i) / Tau - 2 * Beta(i) * (Tau - Gamma(i)))
Next i
For i = 55 To 56

    Theta = 0: Delta = 0: Psi = 0: dPsidTau = 0: d2PsidTauDelta = 0: dDeltaddelta = 0: dDeltaBIdDelta = 0: dDeltaBIdTau = 0: dPsidDelta = 0: d2DeltaBIdRhoTau = 0
    
    Theta = (1 - Tau) + A_caps(i) * ((Delta_Rho - 1) ^ 2) ^ (1 / (2 * Beta(i)))
    
    Delta = Theta ^ 2 + B_Caps(i) * ((Delta_Rho - 1) ^ 2) ^ a(i)
    
    Psi = Exp(-c_Caps(i) * (Delta_Rho - 1) ^ 2 - D_Caps(i) * (Tau - 1) ^ 2)
    
    dPsidTau = -2 * D_Caps(i) * (Tau - 1) * Psi
    d2PsidTauDelta = 4 * c_Caps(i) * Psi * D_Caps(i) * (Delta_Rho - 1) * (Tau - 1)
    dDeltaddelta = (Delta_Rho - 1) * (A_caps(i) * Theta * 2 / Beta(i) * ((Delta_Rho - 1) ^ 2) ^ (1 / (2 * Beta(i)) - 1) + 2 * B_Caps(i) * a(i) * ((Delta_Rho - 1) ^ 2) ^ (a(i) - 1))
    dDeltaBIdDelta = b(i) * Delta ^ (b(i) - 1) * dDeltaddelta
    dDeltaBIdTau = -2 * Theta * b(i) * Delta ^ (b(i) - 1)
    dPsidDelta = -2 * c_Caps(i) * (Delta_Rho - 1) * Psi
    d2DeltaBIdRhoTau = -A_caps(i) * b(i) * 2 / Beta(i) * Delta ^ (b(i) - 1) * (Delta_Rho - 1) * ((Delta_Rho - 1) ^ 2) ^ (1 / (2 * Beta(i)) - 1) - 2 * Theta * b(i) * (b(i) - 1) * Delta ^ (b(i) - 2) * dDeltaddelta
    
    Sum4 = Sum4 + n(i) * (Delta ^ b(i) * (dPsidTau + Delta_Rho * d2PsidTauDelta) + Delta_Rho * dDeltaBIdDelta * dPsidTau + dDeltaBIdTau * (Psi + Delta_Rho * dPsidDelta) + d2DeltaBIdRhoTau * Delta_Rho * Psi)
    

Next i
PhiR_DeltaTau = (Sum1 + Sum2 + Sum3 + Sum4)

End Function

'WAGNER 2002 eq from table 6.5
Public Function PhiR_Delta#(Delta_Rho, Tau)
Dim i As Integer

Dim C(), d()
Dim T(), n()
Dim Alpha(), Epsilon()
Dim Beta(), Gamma()
Dim a(), b(), B_Caps(), A_caps()
Dim c_Caps(), D_Caps()
'declaring of  components for F0+Fr and experimental PT spaces
Dim Theta#, Delta#, Psi#, dPsidDelta#, dDeltaddelta#, dDeltaBIdDelta#
'temporary double units
Dim Sum1 As Double, Sum2 As Double, Sum3 As Double, Sum4 As Double

LoadConstants C(), d(), T(), n(), Alpha(), Epsilon(), Beta(), Gamma(), a(), b(), B_Caps(), A_caps(), c_Caps(), D_Caps()

For i = 1 To 7
    Sum1 = Sum1 + n(i) * d(i) * Delta_Rho ^ (d(i) - 1) * (Tau) ^ T(i)
Next i
For i = 8 To 51
    Sum2 = Sum2 + n(i) * Exp(-Delta_Rho ^ C(i)) * (Delta_Rho ^ (d(i) - 1) * Tau ^ T(i) * (d(i) - C(i) * Delta_Rho ^ C(i)))
Next i
For i = 52 To 54
    Sum3 = Sum3 + n(i) * Delta_Rho ^ d(i) * Tau ^ T(i) * Exp(-Alpha(i) * (Delta_Rho - Epsilon(i)) ^ 2 - Beta(i) * (Tau - Gamma(i)) ^ 2) * (d(i) / Delta_Rho - 2 * Alpha(i) * (Delta_Rho - Epsilon(i)))
Next i
For i = 55 To 56
    Theta = (1 - Tau) + A_caps(i) * ((Delta_Rho - 1) ^ 2) ^ (1 / (2 * Beta(i)))
    Delta = Theta ^ 2 + B_Caps(i) * ((Delta_Rho - 1) ^ 2) ^ a(i)
    Psi = Exp(-c_Caps(i) * (Delta_Rho - 1) ^ 2 - D_Caps(i) * (Tau - 1) ^ 2)
    
    dPsidDelta = -2 * c_Caps(i) * (Delta_Rho - 1) * Psi
    dDeltaddelta = (Delta_Rho - 1) * (A_caps(i) * Theta * 2 / Beta(i) * ((Delta_Rho - 1) ^ 2) ^ (1 / (2 * Beta(i)) - 1) + 2 * B_Caps(i) * a(i) * ((Delta_Rho - 1) ^ 2) ^ (a(i) - 1))
    If Delta = 0 Then
        dDeltaBIdDelta = 0
    Else
        dDeltaBIdDelta = dDeltaddelta * b(i) * Delta ^ (b(i) - 1)
    End If

    Sum4 = Sum4 + n(i) * (Delta ^ b(i) * (Psi + Delta_Rho * dPsidDelta) + dDeltaBIdDelta * Delta_Rho * Psi)
Next i
PhiR_Delta = (Sum1 + Sum2 + Sum3 + Sum4)

End Function

Private Function PhiR_Tau#(Delta_Rho, Tau)
Dim i As Integer
Dim C(), d()
Dim T(), n()
Dim Alpha(), Epsilon()
Dim Beta(), Gamma()
Dim a(), b(), B_Caps(), A_caps()
Dim c_Caps(), D_Caps()
'declaring of  components for F0+Fr and experimental PT spaces
Dim Delta#, Theta#, Psi#, dDeltaBIdTau#, dPsidTau#
'temporary double units
Dim Sum1 As Double, Sum2 As Double, Sum3 As Double, Sum4 As Double, Pwr#

LoadConstants C(), d(), T(), n(), Alpha(), Epsilon(), Beta(), Gamma(), a(), b(), B_Caps(), A_caps(), c_Caps(), D_Caps()

For i = 1 To 7
    Sum1 = Sum1 + n(i) * T(i) * Delta_Rho ^ d(i) * Tau ^ (T(i) - 1)
Next i
For i = 8 To 51
    Sum2 = Sum2 + n(i) * T(i) * Delta_Rho ^ d(i) * Tau ^ (T(i) - 1) * Exp(-Delta_Rho ^ C(i))
Next i
For i = 52 To 54
    Sum3 = Sum3 + n(i) * Delta_Rho ^ d(i) * Tau ^ T(i) * Exp(-Alpha(i) * (Delta_Rho - Epsilon(i)) ^ 2 - Beta(i) * (Tau - Gamma(i)) ^ 2) * (T(i) / Tau - 2 * Beta(i) * (Tau - Gamma(i)))
Next i
For i = 55 To 56
    Theta = (1 - Tau) + A_caps(i) * ((Delta_Rho - 1) ^ 2) ^ (1 / (2 * Beta(i)))
    Delta = Theta ^ 2 + B_Caps(i) * ((Delta_Rho - 1) ^ 2) ^ a(i)
    Psi = Exp(-c_Caps(i) * (Delta_Rho - 1) ^ 2 - D_Caps(i) * (Tau - 1) ^ 2)
    
    dDeltaBIdTau = -2 * Theta * b(i) * Delta ^ (b(i) - 1)
    dPsidTau = -2 * D_Caps(i) * (Tau - 1) * Psi
    
    Sum4 = Sum4 + n(i) * Delta_Rho * (dDeltaBIdTau * Psi + Delta ^ b(i) * dPsidTau)
Next i
PhiR_Tau = (Sum1 + Sum2 + Sum3 + Sum4)
End Function
'WAGNER 2002
Private Function Phi0_Tau#(Delta_Rho, Tau)
Dim i%, n0#(2 To 8), gamma0#(4 To 8), sum#

n0(2) = 6.6832105268: n0(3) = 3.00632: n0(4) = 0.012436: n0(5) = 0.97315: n0(6) = 1.2795: n0(7) = 0.96956: n0(8) = 0.24873

gamma0(4) = 1.28728967: gamma0(5) = 3.53734222: gamma0(6) = 7.74073708: gamma0(7) = 9.24437796: gamma0(8) = 27.5075105

For i = 4 To 8
    sum = sum + n0(i) * gamma0(i) * ((1 - Exp(-gamma0(i) * Tau)) ^ -1 - 1)
Next i

Phi0_Tau = n0(2) + n0(3) / Tau + sum
End Function

'WAGNER 2002
Private Function Phi0#(Delta_Rho, Tau)

Dim i As Integer
'constants from table 6.1
Dim n0(1 To 8) As Double, gamma0(4 To 8) As Double, Sum1 As Double

'constants from table 6.1
n0(1) = -8.32044648201: n0(2) = 6.6832105268: n0(3) = 3.00632: n0(4) = 0.012436
n0(5) = 0.97315:        n0(6) = 1.2795:       n0(7) = 0.96956: n0(8) = 0.24873

gamma0(4) = 1.28728967: gamma0(5) = 3.53734222: gamma0(6) = 7.74073708
gamma0(7) = 9.24437796: gamma0(8) = 27.5075105

For i = 4 To 8
    Sum1 = Sum1 + n0(i) * LogExp(1 - Exp(1) ^ (-gamma0(i) * Tau))
Next i
Phi0 = LogExp(Delta_Rho) + n0(1) + n0(2) * Tau + n0(3) * LogExp(Tau) + Sum1

End Function
'WAGNER 2002
Private Function PhiR#(Delta_Rho, Tau)

Dim i As Integer
Dim C(), d()
Dim T(), n()
Dim Alpha(), Epsilon()
Dim Beta(), Gamma()
Dim a(), b(), B_Caps(), A_caps()
Dim c_Caps(), D_Caps()
'declaring of  components for F0+Fr and experimental PT spaces
Dim Delta As Double, Theta As Double, Psi As Double
'temporary double units
Dim Sum1 As Double, Sum2 As Double, Sum3 As Double, Sum4 As Double
LoadConstants C(), d(), T(), n(), Alpha(), Epsilon(), Beta(), Gamma(), a(), b(), B_Caps(), A_caps(), c_Caps(), D_Caps()

For i = 1 To 7
    Sum1 = Sum1 + n(i) * Delta_Rho ^ d(i) * (Tau) ^ T(i)
Next i
For i = 8 To 51
    Sum2 = Sum2 + n(i) * Delta_Rho ^ d(i) * (Tau) ^ T(i) * Exp(1) ^ (-Delta_Rho ^ C(i))
Next i
For i = 52 To 54
    Sum3 = Sum3 + n(i) * Delta_Rho ^ d(i) * (Tau) ^ T(i) * Exp(1) ^ (-Alpha(i) * (Delta_Rho - Epsilon(i)) ^ 2 - Beta(i) * (Tau - Gamma(i)) ^ 2)
Next i
For i = 55 To 56
    Theta = (1 - Tau) + A_caps(i) * ((Delta_Rho - 1) ^ 2) ^ (1 / (2 * Beta(i)))
    Delta = Theta ^ 2 + B_Caps(i) * ((Delta_Rho - 1) ^ 2) ^ a(i)
    Psi = Exp(1) ^ (-c_Caps(i) * (Delta_Rho - 1) ^ 2 - D_Caps(i) * (Tau - 1) ^ 2)
    Sum4 = Sum4 + n(i) * Delta ^ b(i) * Delta_Rho * Psi
Next i
PhiR = (Sum1 + Sum2 + Sum3 + Sum4)
End Function
'Huber et al 2009, New International Formulation for the Viscosity of H2O
'Valid up to 1173 K and 1000 MPa or 10 kBar
' T in K, Rho in kg/m3, output viscosity in micropascal/second (NOT MILIPASCAL/SEC!!!) or 1000 Cp
Public Function Water_Viscosity_calc#(T_Inc, Rho_inc)
Dim i As Integer, j As Integer
Dim T As Single, T_Star As Single, T_ As Double
Dim Rho As Double, Rho_star As Single, Rho_ As Double, Rho2 As Single, Rho_2 As Double
Dim P As Single, P_Star As Single, P_ As Double
Dim Mu As Double, Mu_star As Single, Mu_ As Double
Dim Mu_0 As Double, Mu_1 As Double, Mu_2 As Double
Dim Tmp_Sum1 As Double, Tmp_Sum2 As Double, Tmp_Sum3 As Double
Dim Hi(0 To 3) As Single, Hij(0 To 5, 0 To 6) As Double
Dim Psi_D As Double, l As Double, w As Double
Dim Chi_ As Single, Chi1 As Double, Chi2 As Double
Dim Chi_Mu As Single, qc As Double, qd As Double, Upsilon As Single, Gamma As Single
Dim Xi As Double, Xi0 As Single, gamma0 As Single, T_r As Single
Dim Y As Double

T_Star = 647.096
Rho_star = 322
P_Star = 22.064

Mu_star = 1E-06
Hi(0) = 1.67752:  Hi(1) = 2.20462: Hi(2) = 0.6366564: Hi(3) = -0.241605

Hij(0, 0) = 0.520094:  Hij(1, 0) = 0.0850895:   Hij(2, 0) = -1.08374:  Hij(3, 0) = -0.289555
Hij(0, 1) = 0.222531:  Hij(1, 1) = 0.999115:    Hij(2, 1) = 1.88797:   Hij(3, 1) = 1.26613
Hij(5, 1) = 0.120573:  Hij(0, 2) = -0.281378:   Hij(1, 2) = -0.906851: Hij(2, 2) = -0.772479
Hij(3, 2) = -0.489837: Hij(4, 2) = -0.25704:    Hij(0, 3) = 0.161913:  Hij(1, 3) = 0.257399
Hij(0, 4) = -0.0325372: Hij(3, 4) = 0.0698452:  Hij(4, 5) = 0.00872102: Hij(3, 6) = -0.00435673
Hij(5, 6) = -0.000593264

Chi_Mu = 0.068: qc = 1.9 ^ -1:    qd = 1.1 ^ -1: Upsilon = 0.63: Gamma = 1.239
Xi0 = 0.13:     gamma0 = 0.06: T_r = 1.5
T = T_Inc
T_ = T / T_Star

Rho = Rho_inc
Rho2 = Rho * 0.9995
Rho_ = Rho / Rho_star
Tmp_Sum2 = T_Star * 1.5

Chi1 = ((Rho - Rho2) / (Water_Pressure_calc(T, Rho) - Water_Pressure_calc(T, Rho2))) * P_Star / Rho_star
Chi2 = ((Rho - Rho2) / (Water_Pressure_calc(Tmp_Sum2, Rho) - Water_Pressure_calc(Tmp_Sum2, Rho2))) * P_Star / Rho_star

Chi_ = (Chi1 - Chi2 * Tmp_Sum2 / T) * Rho_

If Chi_ < 0 Then Chi_ = 0

Tmp_Sum2 = 0
For i = 0 To 3
    Tmp_Sum1 = Tmp_Sum1 + Hi(i) / T_ ^ i
Next i

Mu_0 = 100 * T_ ^ 0.5 / Tmp_Sum1
Tmp_Sum1 = 0

For i = 0 To 5
    For j = 0 To 6
        Tmp_Sum1 = Tmp_Sum1 + Hij(i, j) * (Rho_ - 1) ^ j
    Next j
    Tmp_Sum2 = Tmp_Sum2 + ((1 / T_ - 1) ^ i) * Tmp_Sum1
    Tmp_Sum1 = 0
Next i
Mu_1 = Exp(Rho_ * Tmp_Sum2)
Tmp_Sum2 = 0
Tmp_Sum1 = 0

Xi = Xi0 * (Chi_ / gamma0) ^ (Upsilon / Gamma)

If Xi < 0 Then Xi = 0
If Xi >= 0 And Xi <= 0.3817016416 Then
    Y = 0.2 * qc * Xi * (qd * Xi) ^ 5 * (1 - qc * Xi + (qc * Xi) ^ 2 - 765 / 504 * (qd * Xi) ^ 2)
Else
    'Psi_D = WorksheetFunction.Acos((1 + (qd * Xi) ^ 2) ^ -0.5)
    'arccs function have two strings removed from calculation
    Psi_D = Arccs((1 + (qd * Xi) ^ 2) ^ -0.5)
    w = Abs((qc * Xi - 1) / (qc * Xi + 1)) ^ 0.5 * Tan(Psi_D / 2)
    If qc * Xi > 1 Then l = LogExp((1 + w) / (1 - w)) Else l = 2 * Atn(Abs(w))
    Y = 1 / 12 * Sin(3 * Psi_D) - 0.25 / (qc * Xi) * Sin(2 * Psi_D) + 1 / (qc * Xi) ^ 2 * (1 - 1.25 * (qc * Xi) ^ 2) * Sin(Psi_D) - 1 / (qc * Xi) ^ 3 * ((1 - 1.5 * (qc * Xi) ^ 2) * Psi_D - Abs((qc * Xi) ^ 2 - 1) ^ 1.5 * l)
End If
Mu_2 = Exp(Chi_Mu * Y)

Mu = Mu_0 * Mu_1 * Mu_2
Water_Viscosity_calc = Mu

End Function

Public Function Rho_Water_Vap_sat#(Ti)
Dim T#, T_inv#, RhoVapSat#, C#(1 To 6)
'used eq 2.7 to obtain vapor density from Wagner and Pruss, output in bar
T = Ti
If T = 0 Then T = 0.01
T = T + 273.15
T_inv = 1 - T / 647.096
   
C(1) = -2.0315024: C(2) = -2.6830294: C(3) = -5.38626492
C(4) = -17.2991605: C(5) = -44.7586581: C(6) = -63.9201063

Rho_Water_Vap_sat = Exp(C(1) * T_inv ^ (1 / 3) + C(2) * T_inv ^ (2 / 3) + C(3) * T_inv ^ (4 / 3) + C(4) * T_inv ^ 3 + C(5) * T_inv ^ (37 / 6) + C(6) * T_inv ^ (71 / 6)) * 322

End Function
Public Function Rho_Water_Liq_sat#(Ti)
Dim T#, T_inv#, RhoLiqSat#, b#(1 To 6)
'used eq 2.6 to obtain liquid density from Wagner and Pruss, output in bar
T = Ti
'If T > 373.946 Then T = 373.945
If T = 0 Then T = 0.01
T = T + 273.15
T_inv = 1 - T / 647.096

b(1) = 1.99274064:  b(2) = 1.09965342:  b(3) = -0.510839303
b(4) = -1.75493479: b(5) = -45.5170352: b(6) = -674694.45

Rho_Water_Liq_sat = (1 + b(1) * T_inv ^ (1 / 3) + b(2) * T_inv ^ (2 / 3) + b(3) * T_inv ^ (5 / 3) + b(4) * T_inv ^ (16 / 3) + b(5) * T_inv ^ (43 / 3) + b(6) * T_inv ^ (110 / 3)) * 322
    
End Function
Public Function P_H2O_Boiling_Curve#(Ti_C)
'Output is in bar, eq. 2.5
Dim T#, T_inv#, RhoVapSat#, a#(1 To 6)

T = Ti_C
If T = 0 Then T = 0.01
T = T + 273.15
T_inv = 1 - T / 647.096

a(1) = -7.85951783: a(2) = 1.84408259:  a(3) = -11.7866497
a(4) = 22.6807411:  a(5) = -15.9618719: a(6) = 1.80122502

P_H2O_Boiling_Curve = Exp(647.096 / T * (a(1) * T_inv + a(2) * T_inv ^ 1.5 + a(3) * T_inv ^ 3 + a(4) * T_inv ^ 3.5 + a(5) * T_inv ^ 4 + a(6) * T_inv ^ 7.5)) * 220.64

End Function


Public Function Water_Boiling_Curve#(Ti)
Dim T#, T_inv#, RhoVapSat#, C#(1 To 6)
'used eq 2.7 to obtain liquid density from Wagner and Pruss, output in bar
    T = Ti
    If T >= 373.946 Then T = 373.946
    If T = 0 Then T = 0.01
    T = T + 273.15
    T_inv = 1 - T / 647.096
       
    C(1) = -2.0315024: C(2) = -2.6830294: C(3) = -5.38626492
    C(4) = -17.2991605: C(5) = -44.7586581: C(6) = -63.9201063
    
    RhoVapSat = Exp(C(1) * T_inv ^ (1 / 3) + C(2) * T_inv ^ (2 / 3) + C(3) * T_inv ^ (4 / 3) + _
        C(4) * T_inv ^ 3 + C(5) * T_inv ^ (37 / 6) + C(6) * T_inv ^ (71 / 6)) * 322
    
    'SO IF REQUEST RHO BASED ON PRESSURE, VALUES CAN FELT INTO VAPOR STATE
    Water_Boiling_Curve = Water_Pressure_calc(T, RhoVapSat) * 10

End Function
Public Function Water_Melting_Curve#(Ti#, Optional IceIII_andIceV_In_Use As Boolean = False)
'added 2019 input in C, output in bar
'because IceIII and IceV have T range overlapped with IceI
    'the boolean variable implemented, such that if it's on, then results will be for IceIII to iceVII,
    'otherwhise, it's IceI, IceV-IceII
Dim Pn#, Tn#, Theta#
Ti = Ti + 273.15

If IceIII_andIceV_In_Use And Ti < 273.31 Then
    Select Case Ti
    Case 251.165 To 256.164
        'iceIII
        Tn = 251.165
        Pn = 209.9
        Theta = Ti / Tn
        Water_Melting_Curve = Pn * (1 - 0.0295252 * (1 - Theta ^ 60))
        
    Case 256.164 To 273.31
        'iceV
        Tn = 256.164
        Pn = 350.1
        Theta = Ti / Tn
        Water_Melting_Curve = Pn * (1 - 1.18721 * (1 - Theta ^ 8))
        
    End Select
Else
    Select Case Ti
    Case 251.165 To 273.15
        'iceI
        Tn = 273.15
        Pn = 0.000611657
        Theta = Ti / Tn
        Water_Melting_Curve = Pn * (1 - 0.626 * 1000000# * (1 - Theta ^ -3) + 0.197135 * 1000000# * (1 - Theta ^ 21.2))
        
    Case 273.15 To 273.31
        'iceV
        Tn = 256.164
        Pn = 350.1
        Theta = Ti / Tn
        Water_Melting_Curve = Pn * (1 - 1.18721 * (1 - Theta ^ 8))
        
    Case 273.31 To 355
        'iceVI
        Tn = 273.31
        Pn = 632.4
        Theta = Ti / Tn
        Water_Melting_Curve = Pn * (1 - 1.07476 * (1 - Theta ^ 4.6))
    Case 355 To 715
        'iceVII
        Tn = 355
        Pn = 2216
        Theta = Ti / Tn
        Water_Melting_Curve = Pn * Exp(1.73683 * (1 - Theta ^ -1) - 0.0544606 * (1 - Theta ^ 5) + 0.806106 * 1E-07 * (1 - Theta ^ 22))
        
    End Select
End If
'converted to bar:
Water_Melting_Curve# = Water_Melting_Curve# * 10

End Function
Public Function Water_Sublimation_Curve#(Ti)
'added 2019, though it's outside of T limits, so it won't be used
'input in degree C, output in bar
Dim Pn#, Tn#, Theta#
Pn = 0.000611657
Tn = 273.16
Theta = (Ti + 273.15) / Tn

Water_Sublimation_Curve = Pn * Exp(-13.928169 * (1 - Theta ^ -1.5) + 34.7078238 * (1 - Theta ^ -1.25))
'converted to bar:
Water_Sublimation_Curve = Water_Sublimation_Curve * 10

End Function
Public Function Rho_Water#(T_in_C#, P_in_Bar#)

Dim T_K#, P_mPa#, n#, aa#, bb#, Toler#, TH2O_Crit#, RhoH2O_Crit#
Dim Rho1#, Rho2#, RhoAprx#, P_tmp1#, P_tmp2#, P_der#, EndLoop As Boolean

T_K = T_in_C + 273.15
P_mPa = P_in_Bar / 10
TH2O_Crit = 647.096
RhoH2O_Crit = 322

If T_K <= TH2O_Crit Then
    If Round_Down(P_mPa * 10, 3) <= Round_Down(Water_Boiling_Curve(T_in_C), 3) Then
        '0.001021135 stands for density of vapor at T=1000 C and P=0.006 bar
        'Rho1 = 0.001021135
        Rho1 = 1E-06
        Rho2 = Rho_Water_Vap_sat(T_in_C) + 1
    Else
        Rho1 = Rho_Water_Liq_sat(T_in_C) - 1
        Rho2 = 1701
    End If
Else
    'Rho1 = 0.001021135
    Rho1 = 1E-06
    Rho2 = 1701
End If

'Bisscectional aproximation first
Toler = 1
While n <= 1000
    DoEvents
    RhoAprx = (Rho1 + Rho2) / 2
    If Abs(Rho2 - Rho1) / 2 < Toler Then
        RhoAprx = Rho2
        n = 1000
    End If
    n = n + 1
    aa = Water_Pressure_calc(T_K, RhoAprx) - P_mPa
    bb = Water_Pressure_calc(T_K, Rho1) - P_mPa
    If Sgn(aa) = Sgn(bb) Then Rho1 = RhoAprx Else Rho2 = RhoAprx
Wend
If Rho1 < 0 Then Rho1 = 0
If Rho2 < 0 Then Rho2 = 0

'Newthon method to compare desired P_in with pressure,
'and to get accuracy defined as either below 0.01 or 1e-10 bar
If High_Accuracy Then Toler = 1E-05 Else Toler = 0.0001
Rho1 = RhoAprx
Rho2 = Rho1 - Toler

n = 0
While Not EndLoop
    DoEvents
    RhoAprx = Rho2
    P_tmp1 = Water_Pressure_calc(T_K, Rho1) - P_mPa
    P_tmp2 = Water_Pressure_calc(T_K, Rho2) - P_mPa
    P_der = ((P_tmp1 - P_tmp2) / (Rho1 - Rho2))
    Rho1 = Rho1 - P_tmp1 / P_der

    P_tmp1 = P_tmp1 + P_mPa
    Rho2 = Rho1 - Toler
    n = n + 1
    If n >= 10000 Then  'this is a killer for infinite loop, which most likely won't be activated
        Rho1 = 123456789
        EndLoop = True
    End If
    If Abs((1 - P_mPa / P_tmp1)) <= 10 ^ -8 Or Abs(RhoAprx - Rho1) < Toler Then EndLoop = True
Wend

Rho_Water = Rho1
End Function

