
Public Sub Reset_Excel()
' sub to reset excel if it sudenly stopped during operation.
'needs to be executed manually, if SetupExcelForCalc was left on
Application.ScreenUpdating = True
Application.Calculation = xlCalculationAutomatic
Application.EnableEvents = True
Application.DisplayAlerts = True

End Sub

Public Sub NewWheel_GoalSeeker()
'This sub designed to perform built-in Excel Goal Seek analysis on range of selected cells, located in one column.
'Works only for single column which contains function and argument for this function in the neighbouring column
Dim StRow%, SlctCell%
StRow = Selection.Cells.row
Count = 0
SlctCell = Selection.Cells.Column
For Each cell In Selection
    DoEvents
    Application.CutCopyMode = False
    'Selected cells are aimed to reach "specified Goal" (here is 100)
    '"C" should be replaced by the column in which argument needs to be changed. Argument must be number, not function
    Range(CStr(ConvertToLetter(SlctCell) & StRow + Count)).GoalSeek Goal:=100, ChangingCell:=Range(CStr(ConvertToLetter(SlctCell + 1) & StRow + Count))
    Count = Count + 1
Next cell
End Sub


Sub Rho_AND_Visc()
' this is the main sub, that reading data from the table, then
' calculating the values according to user defined list of properties
' and uploading that to the sheet
Dim k&, i&, j%, q2h#, Props(), s#, T#, P#, Tmp#, Tmp2#, xNaCl#

SetupExcelForCalc (True)

i = 2
While Cells(i, 2) <> ""
    i = i + 1
Wend
i = i - 1

Props = Array_Props_Arranging(i)
If Continue_execution = False Then Exit Sub
ProgressInStatusBar True

For k = 2 To i
    DoEvents
    CurTime = Timer
    ' this If is the easiest way for me to allow user's generated datapoints
    ' (not macro generated) be used in calcs
    If Cells(k, 4) <> "Ignored" Then
        q2h = 1
        s = S_Unit_Converter(Cells(k, 1), S_Unit, "WtPer")
        xNaCl = S_Unit_Converter(s * 100, "WtPer", "MolPer")
        If T_Unit <> 0 Then T = Cells(k, 2) - 273.15 Else T = Cells(k, 2)
        If P_Unit <> 1 Then P = Cells(k, 3) * 10 Else P = Cells(k, 3)
        Tmp = Rho_Brine(xNaCl, T, P)
        j = 3
        If Prop_Dens = True Then
            Props(k, j) = Tmp
            j = j + 1
        End If
        
        If Prop_Alpha = True Then
            Props(k, j) = Alpha(xNaCl, T, P, Tmp)
            j = j + 1
        End If
        
        If Prop_Beta = True Then
            Props(k, j) = Beta(xNaCl, T, P, Tmp)
            Tmp2 = Props(k, j)
            j = j + 1
        End If
        
        If Prop_Chi = True Then
            If Tmp2 <> 0 Then
                Props(k, j) = XsiCapital(Tmp2, s, T, P, Tmp)
            Else
                Props(k, j) = XsiCapital(Beta(xNaCl, T, P, Tmp), s, T, P, Tmp)
            End If
            j = j + 1
        End If
        
        If Prop_Cf = True Then
            Props(k, j) = Isob_Heat_cap(xNaCl, T, P)
            j = j + 1
        End If
        
        If Prop_H = True Then
            Tmp2 = T_star_H(xNaCl, q2h, T, P)
            Props(k, j) = Water_Enthalpy_calc(Tmp2, Rho_Water(Tmp2 - 273.15, P))
            j = j + 1
        End If
        
        If Prop_Mu = True Then
            If P <= 5000 And Visc_Inc_data_Fit(T_star_Mu(s, T) - 273.15, P) Then
                Props(k, j) = Viscosity_H2O_NaCl(s, T, P, Tmp)
            Else
                Props(k, j) = 0
            End If
            j = j + 1
        End If
        
        If Prop_MSiO = True Then
            Props(k, j) = Qtz_solubility(xNaCl, T, Tmp)
            j = j + 1
        End If
        
        If Prop_isotM = True Then
            Props(k, j) = dQtzDP_const_T(xNaCl, T, P, Tmp)
            j = j + 1
        End If
        
        If Prop_isobM = True Then
            Props(k, j) = dQtzDT_const_P(xNaCl, T, P, Tmp)
            j = j + 1
        End If
        
        If Prop_Min = True Then
            Props(k, j) = mnslb(Mineral_Name(1), xNaCl, T, Tmp)
            j = j + 1
        End If
    Else
        For i = 3 To UBound(Props, 2)
            Props(k, i) = "Ignored"
            j = 4
        Next i
    End If
    ProgressInStatusBar False, "Step 3 out of 3: Calculation of properties", k - 1, i - 1
    
Next k

MakeHeadersForProps

ActiveWorkbook.Sheets("model").Range("D2:" & ConvertToLetter(j) & k - 1) = Props
ActiveWorkbook.ActiveSheet.Cells(1, 1).Select
SetupExcelForCalc (False)

End Sub

Private Sub MakeHeadersForProps()
Dim j%, Pname$, Tname$, Mname$

j = 4
If P_Unit = 1 Then Pname = "bar" Else Pname = "MPa"
Tname = "K"

With ActiveWorkbook.Sheets("model")
    If Prop_Dens Then .Cells(1, j) = "Density ": j = j + 1
    If Prop_Alpha Then .Cells(1, j) = "alpha 1/K": j = j + 1
    If Prop_Beta Then .Cells(1, j) = "beta 1/" & Pname: j = j + 1
    If Prop_Chi Then .Cells(1, j) = "Chi, dimensionless ": j = j + 1
    If Prop_Cf Then .Cells(1, j) = "Cf, kJ/kg/K": j = j + 1
    If Prop_H Then .Cells(1, j) = "h, kJ/kg": j = j + 1
    If Prop_Mu Then .Cells(1, j) = "viscosity, uPa s": j = j + 1
    If Prop_MSiO Then .Cells(1, j) = "quartz solubility, mol/(kgH2O)": j = j + 1
    If Prop_isotM Then .Cells(1, j) = "qtz sol isothermal coef, mol/(kgH2O)/" & Pname: j = j + 1
    If Prop_isobM Then .Cells(1, j) = "qtz sol isobaric coef, mol/(kgH2O)/" & Tname: j = j + 1
    If Prop_Min Then
        Mname = Mineral_Name(0)
        .Cells(1, j) = Mname & "solubility, mol/(kgH2O)": j = j + 1
    End If
    
End With
End Sub
Private Function Mineral_Name$(short As Boolean)
If short Then
    Select Case Prop_Min_Name
        Case 0
            Mineral_Name = "ap"
        Case 1
            Mineral_Name = "calc"
        Case 2
            Mineral_Name = "cor"
        Case 3
            Mineral_Name = "fl"
        Case 4
            Mineral_Name = "qtz"
        Case 5
            Mineral_Name = "ru"
    End Select
Else
    Select Case Prop_Min_Name
        Case 0
            Mineral_Name = "Apatite (F) "
        Case 1
            Mineral_Name = "Calcite "
        Case 2
            Mineral_Name = "Corrundum "
        Case 3
            Mineral_Name = "Fluorite "
        Case 4
            Mineral_Name = "Quartz "
        Case 5
            Mineral_Name = "Rutile "
    End Select
End If
End Function
Private Function Array_Props_Arranging(TotalRows) As Variant()
Dim i%
On Error GoTo XT
i = 2
If Prop_Dens Then i = i + 1
If Prop_Alpha Then i = i + 1
If Prop_Beta Then i = i + 1
If Prop_Chi Then i = i + 1
If Prop_Cf Then i = i + 1
If Prop_H Then i = i + 1
If Prop_Mu Then i = i + 1
If Prop_MSiO Then i = i + 1
If Prop_isotM Then i = i + 1
If Prop_isobM Then i = i + 1
If Prop_Min Then i = i + 1

ReDim Array_Props_Arranging(2 To TotalRows, 3 To i) As Variant
Continue_execution = True
Exit Function
XT:
Continue_execution = False
End Function
