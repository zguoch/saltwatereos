

'Four following subs are hiding part of the form interface for constant PTx
Private Sub None_const_Click()
    Form_Text_Reset
End Sub

Private Sub P_const_Click()
    Visibility_switcher
End Sub

Private Sub T_Const_Click()
    Visibility_switcher
End Sub

Private Sub S_Const_Click()
    Visibility_switcher
End Sub

Private Sub ShowPhases_Click()
If ShowPhases.Value = True Then
    SortBySalinity.Locked = False
    SortBySalinity.Value = True
    SortBySalinity.ForeColor = &H80000012
Else
    SortBySalinity.Locked = True
    SortBySalinity.Value = False
    SortBySalinity.ForeColor = &H8000000A
End If
End Sub

Private Sub Form_Text_Reset()

S_Label_From.Caption = "From"
T_Label_From.Caption = "From"
P_Label_From.Caption = "From"

S_Label_To.Visible = True
T_Label_To.Visible = True
P_Label_To.Visible = True

S_End.Visible = True
T_End.Visible = True
P_End.Visible = True

Label4.Visible = True
Label6.Visible = True
Label7.Visible = True

Label11.Visible = True
Label13.Visible = True

Salt_incr.Visible = True
T_incr.Visible = True
P_incr.Visible = True
End Sub

Private Sub Visibility_switcher()

If None_const.Value = True Then
    Form_Text_Reset
End If

If S_Const.Value = True Then
    Form_Text_Reset
    S_Label_From.Caption = "   Is"
    S_End.Value = S_start.Value
    S_Label_To.Visible = False
    S_End.Visible = False
    Label4.Visible = False
    Salt_incr.Visible = False
End If

If T_Const.Value = True Then
    Form_Text_Reset
    T_Label_From.Caption = "   Is"
    T_End.Value = T_Start.Value
    T_Label_To.Visible = False
    T_End.Visible = False
    Label6.Visible = False
    Label11.Visible = False
    T_incr.Visible = False
End If

If P_const.Value = True Then
    Form_Text_Reset
    P_Label_From.Caption = "   Is"
    P_End.Value = P_start.Value
    P_Label_To.Visible = False
    P_End.Visible = False
    Label7.Visible = False
    Label13.Visible = False
    P_incr.Visible = False
End If

End Sub
Private Sub RecordPropsToListen()

Prop_Dens = Prop_Dens_Show.Value
Prop_Alpha = Prop_Alpha_Show.Value
Prop_Beta = Prop_Beta_Show.Value
Prop_Chi = Prop_Chi_Show.Value
Prop_Cf = Prop_Cf_Show.Value
Prop_H = Prop_Enth_Show.Value
Prop_Mu = Prop_Visc_Show.Value
Prop_MSiO = Prop_Qtz_Show.Value
Prop_isotM = Prop_Qtz_Isot_Show.Value
Prop_isobM = Prop_Qtz_Isob_Show.Value
Prop_Min = B_M_solub.Value
Select Case MinName.Value
    Case "Ap"
        Prop_Min_Name = 0
    Case "Calc"
        Prop_Min_Name = 1
    Case "Cor"
        Prop_Min_Name = 2
    Case "Fl"
        Prop_Min_Name = 3
    Case "Qtz"
        Prop_Min_Name = 4
    Case "Ru"
        Prop_Min_Name = 5
End Select
End Sub


Private Sub Run_Model_Click()

'''Workbook reset to initial state
Visibility_switcher
RecordPropsToListen
If SheetExist("model") = True Then
    If MsgBox("Reset sheet?", vbOKCancel) = vbOK Then
        Model_Page_setup
    Else
        Exit Sub
    End If
Else
    Worksheets.Add(After:=Worksheets(1)).Name = "model"
    Model_Page_setup
End If
ActiveWorkbook.Sheets("model").Select
SetupExcelForCalc True

Dim S_St!, S_En!, S_Inc!, T_St!, T_En!, T_Inc!, P_St!, P_En!, P_Inc!, tmp1!, Tmp2!, Tmp3!, Tmp4!, Row_Counter&, Total_Counter&
Dim S_index%, T_index%, P_index%, i%, j%, k%, l%
Dim LV_boundary_Pres#(), LH_boundary_Pres#(), DataPointsToWorkWith() As Boolean, Data_Array!()

Dim LH_Pressurizer#

S_St = Val(S_start.Value) / 100: S_En = Val(S_End.Value) / 100: S_Inc = Val(Salt_incr.Value) / 100
T_St = Val(T_Start.Value): T_En = Val(T_End.Value): T_Inc = Val(T_incr.Value)
P_St = Val(P_start.Value): P_En = Val(P_End.Value): P_Inc = Val(P_incr.Value)

'first estimate of total datapoints to analyze
tmp1 = S_St
While Round(tmp1, 5) <= Round(S_En, 5)
    tmp1 = tmp1 + S_Inc
    Total_Counter = Total_Counter + 1
Wend

Tmp2 = Total_Counter
S_index = Total_Counter
Total_Counter = 0
tmp1 = T_St

While Round(tmp1, 5) <= Round(T_En, 5)
    tmp1 = tmp1 + T_Inc
    Total_Counter = Total_Counter + 1
Wend

T_index = Total_Counter
Tmp2 = Total_Counter * Tmp2
Total_Counter = 0
tmp1 = P_St

While Round(tmp1, 5) <= Round(P_En, 5)
    tmp1 = tmp1 + P_Inc
    Total_Counter = Total_Counter + 1
Wend

P_index = Total_Counter
Total_Counter = Total_Counter * Tmp2
Tmp3 = Total_Counter / Tmp2

tmp1 = 0
Tmp2 = 0
Unload H2O_NaCl_model

IncomingData S_St, S_En, S_Inc, T_St, T_En, T_Inc, P_St, P_En, P_Inc
ProgressInStatusBar True

ReDim LV_boundary_Pres#(1 To S_index, 1 To T_index)
ReDim LH_boundary_Pres#(1 To S_index, 1 To T_index)
ReDim DataPointsToWorkWith(1 To S_index, 1 To T_index, 1 To P_index) As Boolean

'here are the exact estimates of data points' amount
'i is counter for X
    'j for T
    'k for P
For i = 1 To S_index
    tmp1 = S_St + (i - 1) * S_Inc
    For j = 1 To T_index
        DoEvents
        Tmp2 = T_St + (j - 1) * T_Inc
        For k = 1 To P_index
            DoEvents
            Tmp3 = (k - 1) * P_Inc + P_St
            If tmp1 = 0 Then
                DataPointsToWorkWith(i, j, k) = True
                Tmp4 = Tmp4 + 1
            Else
                If Single_Phase_Checker(Tmp2 - T_Unit, Tmp3 / P_Unit, WtToMol(tmp1 * 100) / 100, False) Then
                    DataPointsToWorkWith(i, j, k) = True
                    Tmp4 = Tmp4 + 1
                End If
            End If
            Row_Counter = Row_Counter + 1
        Next k
    ProgressInStatusBar False, "Step 1 from 3: Phase boundaries estimation", Row_Counter, Total_Counter
    Next j
Next i
ReDim Data_Array!(0 To Tmp4, 1 To 3)
tmp1 = 0

Row_Counter = 2
For i = 1 To S_index
    For j = 1 To T_index
        For k = 1 To P_index
            If DataPointsToWorkWith(i, j, k) = True Then
                Data_Array(tmp1, 1) = Round((S_St + (i - 1) * S_Inc) * 100, 3)
                Data_Array(tmp1, 2) = Round((T_St + (j - 1) * T_Inc), 3)
                Data_Array(tmp1, 3) = Round((P_St + (k - 1) * P_Inc), 3)
                tmp1 = tmp1 + 1
                ProgressInStatusBar False, "Step 2 from 3: Recording points to the table", tmp1, Total_Counter
            End If
        Next k
    Next j
Next i
ProgressInStatusBar False, "Step 2 from 3: finishing. Can take a few minutes", 1, 1
Range("A2:C" & UBound(Data_Array, 1) + 1) = Data_Array

If ShowPhases.Value = True Then CP_And_LV_LH S_St, T_St, S_Inc, T_Inc, LV_boundary_Pres, LH_boundary_Pres, SortBySalinity.Value
If PauseWorkBefore.Value = True Then
    Freeze_Code
    Exit Sub
End If

ActiveWorkbook.Sheets("model").Activate
Rho_AND_Visc
If Matrix_Report_page = True Then Matrix_Transpose

SetupExcelForCalc False
End Sub
Private Sub P_End_AfterUpdate()
    If Not IsNumeric(P_End) Or PTx_inc_checker(Val(T_Start), Val(T_End)) Then
        P_End.Value = ""
        MsgBox "Please correct final pressure", vbOKOnly
        P_End.SetFocus
    End If
End Sub
Private Sub P_Start_AfterUpdate()
    If P_const Then P_End.Value = P_start.Value
    If Not IsNumeric(P_start) Or PTx_inc_checker(Val(P_start), Val(P_End)) Then
        P_start.Value = ""
        MsgBox "Please correct initial pressure", vbOKOnly
        P_start.SetFocus
    End If
End Sub
Private Sub S_End_AfterUpdate()
    If Not IsNumeric(S_End) Or PTx_inc_checker(Val(S_start), Val(S_End)) Then
        S_End.Value = ""
        MsgBox "Please correct final salinity", vbOKOnly
        S_End.SetFocus
    End If
End Sub
Private Sub S_Start_AfterUpdate()
    If S_Const Then S_End.Value = S_start.Value
    If Not IsNumeric(S_start) Or PTx_inc_checker(Val(S_start), Val(S_End)) Then
        S_start.Value = ""
        MsgBox "Please correct initial salinity", vbOKOnly
        S_start.SetFocus
    End If
End Sub
Private Sub T_End_AfterUpdate()
    If Not IsNumeric(T_End) Or PTx_inc_checker(Val(T_Start), Val(T_End)) Then
        T_End.Value = ""
        MsgBox "Please correct final temperature", vbOKOnly
        T_End.SetFocus
    End If
End Sub
Private Sub T_Start_AfterUpdate()
    If T_Const Then T_End.Value = T_Start.Value
    If Not IsNumeric(T_Start) Or PTx_inc_checker(Val(T_Start), Val(T_End)) Then
        T_Start.Value = ""
        MsgBox "Please correct initial temperature", vbOKOnly
        T_Start.SetFocus
    End If
End Sub
Private Function PTx_inc_checker(Init#, Final#) As Boolean

If Final < Init Or Init < 0 Then
    PTx_inc_checker = True
End If

End Function
Private Sub Units_Click()
Unit_config.Show
End Sub

Private Sub UserForm_Initialize()

S_Unit = "WtPer"
P_Unit = 1
T_Unit = 0
Dens_Unit = 1

With Salt_incr
    .Value = 5
    .AddItem 0.1
    .AddItem 0.2
    .AddItem 0.5
    .AddItem 1
    .AddItem 2
    .AddItem 3
    .AddItem 5
    .AddItem 10
    .AddItem 20
    .AddItem 25
End With
With T_incr
    .Value = 100
    .AddItem 0.2
    .AddItem 0.5
    .AddItem 1
    .AddItem 2
    .AddItem 3
    .AddItem 5
    .AddItem 10
    .AddItem 25
    .AddItem 50
    .AddItem 100
    .AddItem 250
End With
With P_incr
    .Value = 100
    .AddItem 0.2
    .AddItem 0.5
    .AddItem 1
    .AddItem 2
    .AddItem 3
    .AddItem 5
    .AddItem 10
    .AddItem 25
    .AddItem 50
    .AddItem 100
    .AddItem 250
End With
With MinName
    .Value = "Qtz"
    .AddItem "Ap"
    .AddItem "Calc"
    .AddItem "Cor"
    .AddItem "Fl"
    .AddItem "Qtz"
    .AddItem "Ru"
End With

End Sub

