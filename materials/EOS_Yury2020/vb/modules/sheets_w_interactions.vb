'these variables are used throughout whole macros.
'they store units, list of properties to calculate, and halt/resume variables

'units:
Global S_Unit As String, T_Unit As Single, P_Unit As Single, Dens_Unit As Single, MacInUse As Boolean

'properties:
Global Prop_Dens As Boolean, Prop_Alpha As Boolean, Prop_Beta As Boolean, Prop_Chi As Boolean, _
Prop_Cf As Boolean, Prop_H As Boolean, Prop_Mu As Boolean, Prop_MSiO As Boolean, Prop_isotM As Boolean, _
Prop_isobM As Boolean, Prop_Min As Boolean, Prop_Min_Name As Byte

'supplementary variables - first used to disrupt code execution by request,
'second - to define the accuracy in calculating LV-LH pressure and density.
Global Continue_execution As Boolean, High_Accuracy As Boolean
    
' converts incoming number into alphabet,
' for example
' for 2 result is B (column)
' for 27 = AA
Public Function ConvertToLetter(iCol As Integer) As String

Dim vArr
vArr = Split(Cells(1, iCol).address(True, False), "$")
ConvertToLetter = vArr(0)

End Function

'designed to load all phase boundaries to separate sheet
Public Function CP_And_LV_LH(SSt!, TSt!, SInc!, TInc!, LV_Arr#(), LH_Arr#(), CP_bySalt As Boolean)
Dim i%, Nm$, title$(1 To 5)
If SheetExist("Boundaries") = True Then
    With Sheets("Boundaries").Cells
        .Clear
        .RowHeight = 15
        .ColumnWidth = 8.43
    End With
Else
    Worksheets.Add(After:=Worksheets(1)).Name = "Boundaries"
    With Sheets("Boundaries").Cells
        .Clear
        .RowHeight = 15
        .ColumnWidth = 8.43
    End With
End If

For i = 1 To UBound(LV_Arr, 1)
    ActiveWorkbook.Sheets("Boundaries").Cells(3, i) = (SSt + (i - 1) * SInc) * 100
    ActiveWorkbook.Sheets("Boundaries").Cells(3, i + UBound(LV_Arr, 1) + 1) = (SSt + (i - 1) * SInc) * 100
Next i
Nm = ChrW(&H2190) & "LV pressure" & ChrW(&H2193) & "Temp  " & "LH pressure" & ChrW(&H2192)
ActiveWorkbook.Sheets("Boundaries").Cells(3, UBound(LV_Arr, 1) + 1) = Nm
ActiveWorkbook.Sheets("Boundaries").Cells(i, UBound(LV_Arr, 1) + 1).ColumnWidth = 30
For i = 1 To UBound(LV_Arr, 2)
    ActiveWorkbook.Sheets("Boundaries").Cells(i + 3, UBound(LV_Arr, 1) + 1) = Round((TSt + (i - 1) * TInc) + T_Unit, 3)
Next i

ActiveWorkbook.Sheets("Boundaries").Range("A4:" & ConvertToLetter(UBound(LV_Arr, 1)) & UBound(LV_Arr, 2) + 3) = WorksheetFunction.Transpose(LV_Arr)
ActiveWorkbook.Sheets("Boundaries").Range(ConvertToLetter(UBound(LV_Arr, 1) + 2) & "4:" & ConvertToLetter(2 * UBound(LV_Arr, 1) + 1) & UBound(LH_Arr, 2) + 3) = WorksheetFunction.Transpose(LH_Arr)

CP_Pointer SSt!, TSt!, SInc!, TInc!, UBound(LV_Arr, 1), UBound(LV_Arr, 2), CP_bySalt
ReplaceFives

End Function

Private Sub ReplaceFives()
' used to replace 5000 bar (which refers to absense of single phase at given TX
' so when the calcs are finished, double values replaced by string "N/A" to indicate that
'LV or LH pressure estimation is above limits of the model
If MacInUse Then
    ExecuteExcel4Macro _
        "FORMULA.REPLACE(""5000"",""N/A"",2,1,FALSE,FALSE,,FALSE,FALSE,FALSE,FALSE)"
Else
    ActiveWorkbook.Sheets("Boundaries").Cells.Replace What:="5000", Replacement:="N/A", LookAt:=xlPart, SearchOrder:=xlByRows, MatchCase:=False, SearchFormat:=False, ReplaceFormat:=False
End If
End Sub

Public Static Sub ProgressInStatusBar(Func_Initialisation As Boolean, Optional Label$, Optional Curr_Counter, Optional Max_Counter)
' function to show the progress in excel status bar.
Dim Starting_Time!, Current_Time!, Label_Text$, Progress_Text$, Time_Text$, Time_Left!, Fract_Progress!

If Func_Initialisation = True Then
    Starting_Time = Round(Timer, 0)
    Application.StatusBar = True
Else
    Current_Time = Round(Timer, 0)
    Fract_Progress = Curr_Counter / Max_Counter
    
    Label_Text = Label & ", "
    Progress_Text = "progress: " & Curr_Counter & " of " & Max_Counter & ": " & Round(Fract_Progress * 100, 1) & " %"
    Time_Text = " Time to finish(min) " & Round((Current_Time - Starting_Time) / Fract_Progress * (1 - Fract_Progress) / 60, 0)
    
    Application.StatusBar = Label_Text & Progress_Text & Time_Text
End If

End Sub
Public Function SheetExist(n As String) As Boolean
Dim ws As Worksheet
  SheetExist = False
  For Each ws In Worksheets
    If n = ws.Name Then
      SheetExist = True
      Exit Function
    End If
  Next ws
End Function


Public Sub Matrix_Transpose() '(TP_Space As Boolean)
'this sub should be ready to plot the data in TP or TRho space
' ideally this function converts plain table to matrix
Dim i%, j%, k%, n%, m%, MaxRow%, TotalProps%, TotalTables%, TotalColumns%, tRow%, tCol%
On Error GoTo ExitThisSub
SetupExcelForCalc (True)

GrabDataForMatrix TotalSlts, TotalRows, TotalProps, S_Border, T_array, Prop_Array

If SheetExist("Matrix result") = True Then ActiveWorkbook.Sheets("Matrix result").Delete
Sheets.Add.Name = "Matrix result"

TotalTables = UBound(S_Border, 1) - 1
i = 2
j = 3

While Prop_Array(i, 1) <> Prop_Array(j, 1)
    j = j + 1
Wend
TotalColumns = j - 1
j = i

While ActiveWorkbook.Sheets("model").Cells(S_Border(1), 1) = ActiveWorkbook.Sheets("model").Cells(j, 1)
    j = j + 1
    If Prop_Array(TotalColumns, 1) = Prop_Array(TotalColumns + j, 1) Then MaxRow = MaxRow + 1
Wend

With ActiveWorkbook.Sheets("Matrix result")
    For i = 1 To TotalProps - 3
        For j = 1 To TotalTables
            'headers
            .Cells(1 + (i - 1) * (MaxRow + 3), 1 + (j - 1) * (TotalColumns + 1)) = "Salinity " & ActiveWorkbook.Sheets("model").Cells(S_Border(j), 1)
            .Cells(1 + (i - 1) * (MaxRow + 3), 2 + (j - 1) * (TotalColumns + 1)) = ActiveWorkbook.Sheets("model").Cells(1, i + 3)
            .Cells(1 + (i - 1) * (MaxRow + 3), 3 + (j - 1) * (TotalColumns + 1)) = "T in rows, P in columns"
            'Temperature
            For k = 2 To MaxRow + 1
                .Cells(1 + k + (i - 1) * (MaxRow + 3), 1 + (j - 1) * (TotalColumns + 1)) = T_array(2 + (k - 2) * (TotalColumns - 1))
            Next k
            'Pressure and properties
            For k = 2 To TotalColumns
                m = 0
                tRow = 2 + (i - 1) * (MaxRow + 3)
                tCol = (j - 1) * (TotalColumns + 1) + k
                .Cells(tRow, tCol) = Prop_Array(k, 1)
                    For n = S_Border(j) To S_Border(j + 1) - 1
                        If Prop_Array(n, 1) = Prop_Array(k, 1) Then
                            m = m + 1
                            .Cells(tRow + m, tCol) = Prop_Array(n, i + 1)
                        End If
                    Next n
            Next k
        Next j
    Next i
End With

CollapseMatrixes MaxRow, TotalColumns, TotalTables, TotalProps
ExitThisSub:
SetupExcelForCalc (False)

End Sub
Private Sub CollapseMatrixes(MaxRows%, TotalColumns%, TTablesH%, TTablesW%)
' simplifying matrix navigation by grouping results of matrix calculation through Data-Groups built-in function
Dim i%, j%, k%, l%
j = (TotalColumns - 3) + TTablesH
k = 1
l = 1
With ActiveWorkbook.Sheets("Matrix result")
    For i = 1 To TTablesH
        .Range(CStr(ConvertToLetter(k) & l & ":" & ConvertToLetter(k + TotalColumns - 1) & (l + MaxRows + 1))).Select
        k = k + TotalColumns + 1
        Selection.Columns.Group
        
    Next i
    k = 1
    For i = 1 To TTablesW - 3
        .Range(CStr(ConvertToLetter(1) & k & ":" & ConvertToLetter(6) & (k + MaxRows + 1))).Select
        Selection.Rows.Group
        k = k + MaxRows + 3
    Next i
End With

End Sub

Private Function GrabDataForMatrix(TotalSlts, TotalRows, TotalProps, S_Border As Variant, T_array As Variant, Prop_Array As Variant)
' reads data from the table to upload that to matrix
Dim i%, j%, k%
On Error GoTo ExitThisFunk
'SetupExcelForCalc (True)
TotalRows = 2
TotalProps = 3

With ActiveWorkbook.Sheets("model")
    While .Cells(TotalRows, 1) <> ""
        TotalRows = TotalRows + 1
        If .Cells(TotalRows, 1) <> .Cells(TotalRows - 1, 1) Then TotalSlts = TotalSlts + 1
    Wend
    While .Cells(2, TotalProps) <> ""
        TotalProps = TotalProps + 1
    Wend

    TotalProps = TotalProps - 1
    TotalRows = TotalRows - 1
        
    ReDim S_Border(1 To TotalSlts + 1)
    ReDim T_array(2 To TotalRows)
    ReDim Prop_Array(2 To TotalRows, 1 To TotalProps)
    
    For i = 2 To TotalRows
        T_array(i) = .Cells(i, 2)
        For j = 1 To TotalProps
            Prop_Array(i, j) = .Cells(i, j + 2)
        Next j
    Next i

    k = 1
    j = 2
    S_Border(k) = 2
    S_Border(TotalSlts + 1) = TotalRows + 1
    For i = 2 To TotalRows
        If .Cells(S_Border(k), 1) = .Cells(i, 1) Then
            j = j + 1
        Else
            k = k + 1
            S_Border(k) = i
        End If
    Next i
End With

ExitThisFunk:

End Function

Public Sub SetupExcelForCalc(Start As Boolean)
' This sub speeds up calculations by turning off certain screen/calcualtion updates
If Start Then
    Application.ScreenUpdating = False
    Application.Calculation = xlCalculationManual
    Application.EnableEvents = False
    Application.DisplayAlerts = False
    Application.StatusBar = False
Else
    Application.ScreenUpdating = True
    Application.Calculation = xlCalculationAutomatic
    Application.EnableEvents = True
    Application.DisplayAlerts = True
    Application.StatusBar = True

End If

End Sub

Public Sub Freeze_Code()
'this one and next Sub are designed to stop the execution with user-defined parameters
'and then return to it with manually changed data by clicking "Read from the page" button
    ActiveWorkbook.Sheets("model").Cells(3, 105) = CStr(S_Unit & "@" & T_Unit & "@" & P_Unit & "@" & Dens_Unit _
    & "@" & Prop_Dens & "@" & Prop_Alpha & "@" & Prop_Beta & "@" & Prop_Chi _
    & "@" & Prop_Cf & "@" & Prop_H & "@" & Prop_Mu & "@" & Prop_MSiO & "@" & Prop_isotM & "@" & Prop_isobM & "@" & High_Accuracy _
    & "@" & Prop_Min & "@" & Prop_Min_Name)
    ActiveWorkbook.Sheets("model").Select
    ActiveWorkbook.Sheets("model").Buttons.Add(180, 20, 160, 35).Select
    Selection.OnAction = "Unfreeze_Code"
    Selection.Characters.Text = "Continue calculations"
    SetupExcelForCalc False

    With Selection.Characters(Start:=1, Length:=21).Font
        .Name = "Calibri"
        .FontStyle = "Regular"
        .Size = 11
        .StrikeThrough = False
        .Superscript = False
        .Subscript = False
        .OutlineFont = False
        .Shadow = False
        .Underline = xlUnderlineStyleNone
        .ColorIndex = 1
    End With
    ActiveWorkbook.Sheets("model").Cells(1, 1).Select
    
End Sub

Public Sub Unfreeze_Code()

    Dim vars() As String, Props As String
    Props = ActiveWorkbook.Sheets("model").Cells(3, 105)
    vars() = Split(Props, "@")
    SetupExcelForCalc True
    S_Unit = vars(0)
    T_Unit = Val(vars(1))
    P_Unit = Val(vars(2))
    Dens_Unit = Val(vars(3))

    If vars(4) = "True" Then Prop_Dens = True Else Prop_Dens = False
    If vars(5) = "True" Then Prop_Alpha = True Else Prop_Alpha = False
    If vars(6) = "True" Then Prop_Beta = True Else Prop_Beta = False
    If vars(7) = "True" Then Prop_Chi = True Else Prop_Chi = False
    If vars(8) = "True" Then Prop_Cf = True Else Prop_Cf = False
    If vars(9) = "True" Then Prop_H = True Else Prop_H = False
    If vars(10) = "True" Then Prop_Mu = True Else Prop_Mu = False
    If vars(11) = "True" Then Prop_MSiO = True Else Prop_MSiO = False
    If vars(12) = "True" Then Prop_isotM = True Else Prop_isotM = False
    If vars(13) = "True" Then Prop_isobM = True Else Prop_isobM = False
    If vars(14) = "True" Then High_Accuracy = True Else High_Accuracy = False
    If vars(15) = "True" Then
        Prop_Min = True
        Prop_Min_Name = CByte(Right(Props, 1))
    Else
        Prop_Min = False
    End If

    ActiveWorkbook.Sheets("model").Buttons.Delete
    Rho_AND_Visc
    MainButtonOnSheet
    
    ActiveWorkbook.Sheets("model").Cells(3, 105) = ""
End Sub
Private Function CP_Pointer(SSt!, TSt!, SInc!, TInc!, S_steps%, T_steps%, CalcBySalts As Boolean)
'CP pointer function generates CP_data array in which data stored as
'second index of array stores:
' 0 - salinity in mol %
' 1 - temperature in C
' 2 - pressure in bar
' 3 - density in kg/m3
' 4 - pressure at given density and temperature 1000 C
Dim i%, j%
Dim tmp1#, Tmp2#
Dim CP_data#(), title$()

If CalcBySalts = True Then
    ReDim CP_data(1 To S_steps, 0 To 4)
    For i = 1 To S_steps
        Crit_Array_Filling S_Unit_Converter((SSt + (i - 1) * SInc) * 100, S_Unit, "WtPer"), i, CP_data
    Next i
Else
    ReDim CP_data(1 To T_steps, 0 To 5)
    For i = 1 To T_steps
        X_and_P_crit TSt + (i - 1) * TInc, tmp1, Tmp2
        CP_data(i, 0) = tmp1
        CP_data(i, 1) = TSt + (i - 1) * TInc
        CP_data(i, 2) = Tmp2
        CP_data(i, 3) = Rho_Brine(tmp1, CP_data(i, 1), Tmp2)
        Tmp2 = CP_data(i, 3) * 18.015268 / (18.015268 * (1 - tmp1) + 58.4428 * tmp1)
        CP_data(i, 4) = Water_Pressure_calc(1273.15, Tmp2) * 10
    Next i
End If

tmp1 = 0
j = T_steps + 3

ActiveWorkbook.Sheets("Boundaries").Cells(j + 2, 1) = "NOTE that critical points estimated above 1000 C located outside of the model"
ReDim title$(1 To 5)
title(1) = "CP S ": title(2) = "CP T ":    title(3) = "CP P "
title(4) = "CP density ": title(5) = "P at 1000C and CP density "
For i = 1 To 5
    ActiveWorkbook.Sheets("Boundaries").Cells(j + 3, i) = title(i)
Next i

For i = 1 To UBound(CP_data)
    CP_data(i, 0) = (S_Unit_Converter(CP_data(i, 0) * 100, "MolPer", S_Unit)) * 100
    CP_data(i, 1) = CP_data(i, 1) + T_Unit
    If CP_data(i, 1) < 373 Then
        CP_data(i, 2) = 0
        CP_data(i, 3) = 0
        CP_data(i, 4) = 0
    Else
        CP_data(i, 2) = CP_data(i, 2) / P_Unit
        CP_data(i, 3) = CP_data(i, 3) * Dens_Unit
        CP_data(i, 4) = CP_data(i, 4) / P_Unit
    End If

Next i

ActiveWorkbook.Sheets("Boundaries").Range("A" & j + 4 & ":E" & j + UBound(CP_data, 1) + 3) = CP_data

End Function

Public Sub Crit_Array_Filling(S_Tmp, SaltRun, CP_array)
Dim n As Integer, mH2O As Single, mNaCl As Single, xNaCl As Double, T_min As Single, T_Max As Single

Dim i As Integer, j As Integer
Dim Sum1 As Double, PH2O_Crit As Double, TH2O_Crit As Double

Dim C(1 To 14) As Double, CA(1 To 11) As Single, d(1 To 11) As Double
Dim T As Double, P#, P_Crit#
Dim x_crit As Double, X_Crit_der As Double, t2 As Double, X_Crit2 As Double

mH2O = 18.015268
mNaCl = 58.4428
xNaCl = S_Tmp / mNaCl / (S_Tmp / mNaCl + (1 - S_Tmp) / mH2O)
TH2O_Crit = 373.973 - 0.027

T = 0
'estimation of T
For j = 1 To 80
    If xNaCl = 0 Then
        T = TH2O_Crit
        Exit For
    End If
    T = j * 25
    x_crit = 0: P_Crit = 0
    X_and_P_crit T, x_crit, P_Crit

    If x_crit >= xNaCl Then Exit For
Next j

j = 0
'Newtons algorythm to obtain X_Crit
While Abs(x_crit - xNaCl) > 1E-06
    x_crit = 0: X_Crit2 = 0: X_Crit_der = 0
    X_and_P_crit T, x_crit, P_Crit
    x_crit = x_crit - xNaCl

    t2 = T - 0.001
    X_and_P_crit t2, X_Crit2, P_Crit
    X_Crit2 = X_Crit2 - xNaCl

    X_Crit_der = (x_crit - X_Crit2) / (T - t2)
    T = T - x_crit / X_Crit_der
    x_crit = x_crit + xNaCl

    If xNaCl = 0 Then x_crit = 0
    j = j + 1
    If j > 1000 Then Stop
Wend

If T = 0 Then T = 0.01
X_and_P_crit T, x_crit, P_Crit
    'TEMPERATURE AND PRESSURE CORRECTED FOR CRITICAL POINT OF PURE WATER,
    'ACCORDING TO ITS POSITION IN IAPWS-95


'Savings data in array
CP_array(SaltRun, 0) = xNaCl
CP_array(SaltRun, 1) = T
CP_array(SaltRun, 2) = P_Crit

CP_array(SaltRun, 3) = Rho_Brine(xNaCl, T, P_Crit)
' 1000 C and pressure at this temperature for drawing isochores
'CP_array(SaltRun, 4) = 1000

'!!!!!!!!!This part has not been properly tested!!!!!!!!
Sum1 = CP_array(SaltRun, 3) * mH2O / (mH2O * (1 - xNaCl) + mNaCl * xNaCl)
CP_array(SaltRun, 4) = Water_Pressure_calc(1273.15, Sum1) * 10


End Sub


Public Sub MainButtonOnSheet()
'bringning back buttons that are wiped out from sheet during cleaning
Dim btn As Button
Dim btn2 As Button
Dim btn3 As Button
Dim btn4 As Button
Dim btn5 As Button

ActiveWorkbook.Sheets("model").Select
Set btn = ActiveWorkbook.Sheets("model").Buttons.Add(750, 15, 85, 30)
btn.OnAction = "StartButtonClick"
Set btn2 = ActiveWorkbook.Sheets("model").Buttons.Add(750, 50, 105, 23.25)
btn2.OnAction = "PageReaderRun"
Set btn3 = ActiveWorkbook.Sheets("model").Buttons.Add(750, 80, 105, 23.25)
btn3.OnAction = "LoadMicrothermometry"
Set btn4 = ActiveWorkbook.Sheets("model").Buttons.Add(750, 110, 105, 23.25)
btn4.OnAction = "SupSubs_List_of_User_Functions_for_Mac"
Set btn5 = ActiveWorkbook.Sheets("model").Buttons.Add(750, 140, 105, 23.25)
btn5.OnAction = "Model_Page_setup"



If MacInUse = False Then
    btn.Characters.Text = "Run model"
    With btn.Characters(Start:=1, Length:=9).Font
        .Name = "Calibri"
        .FontStyle = "Regular"
        .Size = 11
        .StrikeThrough = False
        .Superscript = False
        .Subscript = False
        .OutlineFont = False
        .Shadow = False
        .Underline = xlUnderlineStyleNone
        .ColorIndex = 1
    End With
    
    btn2.Characters.Text = "Read from the sheet"
    With btn2.Characters(Start:=1, Length:=18).Font
        .Name = "Calibri"
        .FontStyle = "Regular"
        .Size = 11
        .StrikeThrough = False
        .Superscript = False
        .Subscript = False
        .OutlineFont = False
        .Shadow = False
        .Underline = xlUnderlineStyleNone
        .ColorIndex = 1
    End With

    btn3.Characters.Text = "Microthermometry"
    With btn3.Characters(Start:=1, Length:=18).Font
        .Name = "Calibri"
        .FontStyle = "Regular"
        .Size = 11
        .StrikeThrough = False
        .Superscript = False
        .Subscript = False
        .OutlineFont = False
        .Shadow = False
        .Underline = xlUnderlineStyleNone
        .ColorIndex = 1
    End With
    
    btn4.Characters.Text = "Functions"
    With btn3.Characters(Start:=1, Length:=18).Font
        .Name = "Calibri"
        .FontStyle = "Regular"
        .Size = 11
        .StrikeThrough = False
        .Superscript = False
        .Subscript = False
        .OutlineFont = False
        .Shadow = False
        .Underline = xlUnderlineStyleNone
        .ColorIndex = 1
    End With

    btn5.Characters.Text = "Clear workbook"
    With btn3.Characters(Start:=1, Length:=18).Font
        .Name = "Calibri"
        .FontStyle = "Regular"
        .Size = 11
        .StrikeThrough = False
        .Superscript = False
        .Subscript = False
        .OutlineFont = False
        .Shadow = False
        .Underline = xlUnderlineStyleNone
        .ColorIndex = 1
    End With
Else
    btn.Select
    Selection.Characters.Text = "Run model"
    With Selection.Characters(Start:=1, Length:=9).Font
        .Name = "Helvetica"
        .FontStyle = "Regular"
        .Size = 11
        .StrikeThrough = False
        .Superscript = False
        .Subscript = False
        .OutlineFont = False
        .Shadow = False
        .Underline = xlUnderlineStyleNone
        .ColorIndex = xlAutomatic
    End With

    btn2.Select
    Selection.Characters.Text = "Read from the sheet"
    With Selection.Characters(Start:=1, Length:=18).Font
        .Name = "Helvetica"
        .FontStyle = "Regular"
        .Size = 11
        .StrikeThrough = False
        .Superscript = False
        .Subscript = False
        .OutlineFont = False
        .Shadow = False
        .Underline = xlUnderlineStyleNone
        .ColorIndex = xlAutomatic
    End With
    
    btn3.Select
    btn3.Characters.Text = "Microthermometry"
    With Selection.Characters(Start:=1, Length:=18).Font
        .Name = "Helvetica"
        .FontStyle = "Regular"
        .Size = 11
        .StrikeThrough = False
        .Superscript = False
        .Subscript = False
        .OutlineFont = False
        .Shadow = False
        .Underline = xlUnderlineStyleNone
        .ColorIndex = xlAutomatic
    End With

    btn4.Select
    btn4.Characters.Text = "Functions"
    With Selection.Characters(Start:=1, Length:=18).Font
        .Name = "Helvetica"
        .FontStyle = "Regular"
        .Size = 11
        .StrikeThrough = False
        .Superscript = False
        .Subscript = False
        .OutlineFont = False
        .Shadow = False
        .Underline = xlUnderlineStyleNone
        .ColorIndex = xlAutomatic
    End With
    
    btn5.Select
    btn5.Characters.Text = "Clear workbook"
    With Selection.Characters(Start:=1, Length:=18).Font
        .Name = "Helvetica"
        .FontStyle = "Regular"
        .Size = 11
        .StrikeThrough = False
        .Superscript = False
        .Subscript = False
        .OutlineFont = False
        .Shadow = False
        .Underline = xlUnderlineStyleNone
        .ColorIndex = xlAutomatic
    End With
    
End If

Range("A1").Select
End Sub

'Sub to show user the initial parameters that have been selected
Sub IncomingData(s_s, s_e, s_i, t_s, t_e, t_i, p_s, p_e, p_i)
Dim Tmp As String
With ActiveWorkbook.Sheets("model")

    .Cells(3, 19) = "Salinity, " & CStr(S_Unit)
    If T_Unit = 0 Then Tmp = "C" Else Tmp = "K"
    .Cells(4, 19) = "Temperature, " & Tmp
    If P_Unit = 10 Then Tmp = "MPa" Else Tmp = "bar"
    .Cells(5, 19) = "Pressure, " & Tmp
    .Cells(2, 21) = "Start"
    .Cells(2, 22) = "End"
    .Cells(2, 23) = "Increment"
    
    .Cells(3, 21) = s_s * 100
    .Cells(3, 22) = s_e * 100
    .Cells(3, 23) = s_i * 100
    
    .Cells(4, 21) = t_s
    .Cells(4, 22) = t_e
    .Cells(4, 23) = t_i
    
    .Cells(5, 21) = p_s
    .Cells(5, 22) = p_e
    .Cells(5, 23) = p_i
    
    Range(Cells(1, 19), Cells(1, 23)).Select
    With Selection
        .HorizontalAlignment = xlCenter
        .VerticalAlignment = xlBottom
        .WrapText = False
        .Orientation = 0
        .AddIndent = False
        .IndentLevel = 0
        .ShrinkToFit = False
        .ReadingOrder = xlContext
        .MergeCells = False
    End With
    Selection.Merge
    ActiveCell.FormulaR1C1 = "Incoming data"

End With

End Sub

Sub Model_Page_setup()
'rearranging columns, formats etc
Dim Ttxt As String, Ptxt As String, Dtxt As String
Reset_Excel
SetupExcelForCalc (True)

For Each Worksheet In ActiveWorkbook.Worksheets
    If Worksheet.Name <> "model" Then
        Worksheet.Activate
        ActiveWorkbook.ActiveSheet.Delete
    End If
Next


With Sheets("model").Cells
    .Clear
    .RowHeight = 15
    .ColumnWidth = 8.43
End With

If T_Unit = 273.15 Then Ttxt = "K" Else Ttxt = "C"
If P_Unit = 10 Then Ptxt = "MPa" Else Ptxt = "bar"
If Dens_Unit = 1 Then Dtxt = "kg*m-3" Else Dtxt = "g*cm-3"

Cells(1, 1) = "x, " & S_Unit
Cells(1, 2) = "T, " & Ttxt
Cells(1, 3) = "P, " & Ptxt

Rows("2:2").Select
With ActiveWindow
    .SplitColumn = 0
    .SplitRow = 1
End With
ActiveWindow.FreezePanes = True
Columns("A:F").Select
Columns("A:C").ColumnWidth = 7
Columns("D:N").ColumnWidth = 10#
Range("F1").Activate
Selection.AutoFilter
Columns("D:D").Select
Selection.NumberFormat = "0.0"
Columns("E:N").Select
Selection.NumberFormat = "0.00E+0"
'Columns("C:C").Select
'Selection.NumberFormat = "0.0"
ActiveWorkbook.Sheets("model").Buttons.Delete
MainButtonOnSheet

SetupExcelForCalc (False)

End Sub

Private Sub StartButtonClick()
H2O_NaCl_model.Show
End Sub

Private Sub PageReaderRun()
Dim i%, s#, T#, P#
Continue_execution = False
MsgBox "Data outside of the single phase field are excluded from calculation", vbOKOnly, "Notification"
Unit_config_From_Page.Show
If Continue_execution = False Then Exit Sub
i = 2

With ActiveWorkbook.Sheets("model")
    While .Cells(i, 2) <> ""
        s = .Cells(i, 1)
        T = .Cells(i, 2)
        P = .Cells(i, 3)
        If .Cells(i, 4) = "Ignored" Then .Cells(i, 4) = ""
        If s < 0 Or s > 100 Or (T < 0 And T_Unit = 0) Or (T < 273.15 And T_Unit <> 0) Or P < 0.01 Or P > 5000 Then
            .Cells(i, 4) = "Ignored"
        ElseIf New_LV_P(T, P, 0, S_Unit_Converter(CDbl(s), S_Unit, "MolPer"), True) / P_Unit > P Then
            .Cells(i, 4) = "Ignored"
        End If
        
        i = i + 1
    Wend
End With

props_IDs.Show
Rho_AND_Visc


End Sub
