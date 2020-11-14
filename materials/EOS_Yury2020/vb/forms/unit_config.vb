Private Sub Cancel_Units_Click()
    Unit_config.Hide
End Sub

Public Sub Ok_Click()

If Temp.Value = " K" Then T_Unit = 273.15 Else T_Unit = 0
If Pres.Value = "MPa" Then P_Unit = 10 Else P_Unit = 1
If Dens.Value = "kg*m-3" Then Dens_Unit = 1 Else Dens_Unit = 0.001
If Prec_check Then High_Accuracy = True Else High_Accuracy = False

If NaCl.Value = "Wt. % NaCl" Then
    S_Unit = "WtPer"
ElseIf NaCl.Value = "mol. % NaCl" Then
    S_Unit = "MolPer"
ElseIf NaCl.Value = "vol. % NaCl" Then
    S_Unit = "VolPer"
Else
    S_Unit = "Molal"
End If

With H2O_NaCl_model
    .Salt_Frame.Caption = "Salinity, in " + NaCl
    .Temp_Frame.Caption = "Temperature, in " + Temp
    .Label11.Caption = Temp
    .Pres_Frame.Caption = "Pressure, in " + Pres
    .Label13.Caption = Pres
    
    If T_Unit <> 0 Then
        .T_Start.Value = Round(.T_Start.Value + T_Unit, 2)
        .T_End.Value = Round(.T_End.Value + T_Unit, 2)
    End If
    
    If P_Unit <> 1 Then
        .P_start.Value = .P_start.Value / P_Unit
        .P_End.Value = .P_End.Value / P_Unit
        
        .P_incr.Clear
        .P_incr.AddItem 0.2 / P_Unit
        .P_incr.AddItem 0.5 / P_Unit
        .P_incr.AddItem 1 / P_Unit
        .P_incr.AddItem 2 / P_Unit
        .P_incr.AddItem 3 / P_Unit
        .P_incr.AddItem 5 / P_Unit
        .P_incr.AddItem 10 / P_Unit
        .P_incr.AddItem 25 / P_Unit
        .P_incr.AddItem 50 / P_Unit
        .P_incr.AddItem 100 / P_Unit
        'If .P_incr.Value = "" Then .P_incr.Value = 10
        .P_incr.Value = .P_incr.Value / P_Unit
    End If
End With

Unit_config.Hide

End Sub

Private Sub UserForm_Initialize()

With NaCl
    .AddItem "Wt. % NaCl"
    .AddItem "mol. % NaCl"
    .AddItem "vol. % NaCl"
    .AddItem "molality"
End With

With Temp
    .AddItem " C"
    .AddItem " K"
End With

With Pres
    .AddItem "Bar"
    .AddItem "MPa"
End With

With Dens
    .AddItem "kg*m-3"
    .AddItem "g*cm-3"
End With

End Sub
