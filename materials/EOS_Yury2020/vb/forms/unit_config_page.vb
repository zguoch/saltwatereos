
Private Sub Cancel_Units_Click()
    Unit_config_From_Page.Hide
    Continue_execution = False
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

Unit_config_From_Page.Hide
Continue_execution = True

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
