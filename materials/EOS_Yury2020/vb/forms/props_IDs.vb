Private Sub CommandButton1_Click()

Columns("E:M").Clear

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


Unload props_IDs

End Sub

Private Sub Min2Frame_Click()

End Sub

Private Sub UserForm_Initialize()

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
