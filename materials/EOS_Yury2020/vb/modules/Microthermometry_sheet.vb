'This module creates an example sheet with microthermometrical data
'to demonstrate proper usage of functions.
'Example includes incorrectly typed data to demonstrate result of such mistakes
Sub LoadMicrothermometry()
Dim SheetName As String
SetupExcelForCalc True

SheetName = "microthermometry"
If SheetExist(SheetName) = True Then ActiveWorkbook.Sheets(SheetName).Delete
Sheets.Add.Name = SheetName

ActiveWorkbook.Sheets(SheetName).Activate
Range("A1:P48") = SheetContents()
SheetFormat

SetupExcelForCalc False
End Sub

Sub SheetFormat()
Dim StNm As String

Range("D:D,I:I,M:M").Select
Selection.ColumnWidth = 2
With Selection.Interior
    .Pattern = xlSolid
    .PatternColorIndex = xlAutomatic
    .ThemeColor = xlThemeColorDark1
    .TintAndShade = -0.499984740745262
    .PatternTintAndShade = 0
End With
Range("E1:H1").Merge
Range("J1:L1").Merge
Range("N1:P1").Merge
Range("A1:C1").Merge
Columns("J:L").ColumnWidth = 9.43
Columns("E:E").ColumnWidth = 35
Columns("E:E").WrapText = True
Columns("C:C").ColumnWidth = 15
Columns("L:L").NumberFormat = "0.0"
Range("K2:L48").Select
ActiveSheet.Shapes.AddChart2(240, xlXYScatterSmoothNoMarkers, 750, 100, 500, 300).Select

StNm = ActiveSheet.Name
With ActiveChart
    .SetSourceData Source:=Range((StNm) & "!$K$2:$L$48")
    Application.CutCopyMode = False
    .SeriesCollection.NewSeries
    .FullSeriesCollection(2).Name = "=""Isochore"""
    .FullSeriesCollection(2).XValues = "=" & (StNm) & "!$O$3:$O$4"
    .FullSeriesCollection(2).Values = "=" & (StNm) & "!$P$3:$P$4"
    .Axes(xlCategory, xlPrimary).HasTitle = True
    .Axes(xlCategory, xlPrimary).AxisTitle.Characters.Text = "Temperature (°C)"
    .Axes(xlValue, xlPrimary).HasTitle = True
    .Axes(xlValue, xlPrimary).AxisTitle.Characters.Text = "Pressure (bar)"
    .ChartTitle.Text = _
        "Isochore (red line) for FI with 22.4 salinity." & Chr(13) & "Blue line separates single phase (above) " & Chr(13) & "from L+V two phases field (below)"
    Selection.Format.TextFrame2.TextRange.Characters.Text = _
        "Isochore (red line) for FI with 22.4 salinity." & Chr(13) & "Blue line separates single phase (above) " & Chr(13) & "from L+V two phases field (below)"

    
End With
Rows("3:3").Select
ActiveWindow.FreezePanes = True
Range("A1").Select

End Sub

Private Function SheetContents() As Variant
Dim Val() As Variant
Dim i%, j%, refToTm$, refToTh$, RefToX$

ReDim Val(48, 16)
Val(0, 0) = "Measured data":   Val(0, 4) = "Interpreted data":    Val(0, 9) = "Phase boundary for first row": Val(0, 13) = "Isochore for first row"
Val(1, 0) = "T melt":          Val(1, 2) = "phase":               Val(1, 1) = "T hom": Val(1, 4) = "Bulk salinity": Val(1, 5) = "P hom": Val(1, 6) = "Density": Val(1, 7) = "dP/dT": Val(1, 9) = "S": Val(1, 10) = "T": Val(1, 11) = "P": Val(1, 13) = "x": Val(1, 14) = "T": Val(1, 15) = "P"

Val(2, 0) = "-20":     Val(2, 2) = "ice":     Val(2, 13) = "=E3": Val(2, 14) = "=B3": Val(2, 15) = "=F3"
Val(3, 0) = "-20":     Val(3, 2) = "ICE":     Val(3, 13) = "=N3": Val(3, 14) = "355.66": Val(3, 15) = "=H3*(O4-O3)+P3"
Val(4, 0) = "-20":     Val(4, 2) = "":
Val(5, 0) = "-20":     Val(5, 2) = "i":
Val(6, 0) = "-25":     Val(6, 2) = "ice":
Val(7, 0) = "-10":     Val(7, 2) = "hydrohalite":
Val(8, 0) = "-10":     Val(8, 2) = "HYDROHALITE":
Val(9, 0) = "-10":     Val(9, 2) = "Hydrohalite":
Val(10, 0) = "-10":    Val(10, 2) = "hh":
Val(11, 0) = "-25":    Val(11, 2) = "hydrohalite":
Val(12, 0) = "25":     Val(12, 2) = "halite":
Val(13, 0) = "200":    Val(13, 2) = "halite":
Val(14, 0) = "200":    Val(14, 2) = "HALITE":
Val(15, 0) = "200":    Val(15, 2) = "Halite":
Val(16, 0) = "253":    Val(16, 2) = "h":
Val(17, 0) = "460":    Val(17, 2) = "h":
Val(18, 0) = "950":    Val(18, 2) = "h":
Val(19, 0) = "-25":    Val(19, 2) = "h":


For i = 2 To 19
    Val(i, 1) = "250"
    refToTm = "A" & i + 1
    refToTh = ",B" & i + 1
    RefToX = "E" & i + 1
    
    If i = 4 Then
        Val(i, 4) = "=FlInc_Salinity(" & refToTm & refToTh & ")"
    Else
        Val(i, 4) = "=FlInc_Salinity(" & refToTm & refToTh & ",C" & i + 1 & ")"
    End If
    Val(i, 5) = "=Phases_Single_ph_pressure_at_TX(E" & i + 1 & ",max(" & refToTm & refToTh & "),true)"
    Val(i, 6) = "=Brine_Density(" & RefToX & ",max(" & refToTm & refToTh & "),F" & i + 1 & ",true)"
    Val(i, 7) = "=Brine_dP_dT_for_Isochore(" & RefToX & ",max(" & refToTm & refToTh & "),F" & i + 1 & ",G" & i + 1 & ")"
    If i = 2 Then
        Val(i, 9) = "=E3"
        Val(i, 10) = "100"
    Else
        Val(i, 9) = "=J" & i
        Val(i, 10) = "=K" & i & "+20"
    End If
    Val(i, 11) = "=Phases_Single_ph_pressure_at_TX(J" & i + 1 & ",K" & i + 1 & ",true)"
Next i

For i = 19 To 47
    Val(i, 9) = "=J" & i
    Val(i, 10) = "=K" & i & "+20"
    Val(i, 11) = "=Phases_Single_ph_pressure_at_TX(J" & i + 1 & ",K" & i + 1 & ",true)"
Next i

For i = 0 To 48
    For j = 0 To 16
        If IsEmpty(Val(i, j)) Then
            Val(i, j) = ""
        Else
            If Right(Val(i, j), 1) = ")" Then
                Val(i, j) = "=iferror(" & Right(Val(i, j), Len(Val(i, j)) - 1) & _
                     "," & Chr(34) & "n/a" & Chr(34) & ")"
            End If
        End If
    Next j
Next i
SheetContents = Val

End Function

