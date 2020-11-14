'Excel VBA math functions are missing few important math functions
'and even though they can be substituted with built-in Excel functions,
'separate functions (below) perform calculation significantly faster

Public Function LogExp#(x)
Dim Tmp As Double
Tmp = Exp(1)
LogExp = Log(x) / Log(Tmp#)
End Function

Public Function Log10#(x)
Log10 = Log(x) / Log(10#)
End Function

Public Function Arccos#(x)
If Round(x, 8) = 1# Then Arccos = 0#: Exit Function
If Round(x, 8) = -1# Then Arccos = Pi: Exit Function
Arccos = Atn(-x / Sqr(-x * x + 1)) + 2 * Atn(1)
End Function

Public Function Round_Down#(Nmbr, DecPl%)
Round_Down = Sgn(Nmbr) * Fix(Abs(Nmbr) * 10 ^ (DecPl)) / 10 ^ DecPl
End Function

Public Function Round_Up#(Nmbr, DecPl%)
Round_Up = -Sgn(Nmbr) * Int(-Abs(Nmbr) * 10 ^ (DecPl)) / 10 ^ DecPl
End Function

Public Function Arccs#(x)
If Round(x, 12) = 1# Then Arccs = 0#: Exit Function
If Round(x, 12) = -1# Then Arccs = Pi: Exit Function
Arccs = Atn(-x / Sqr(-x * x + 1)) + 2 * Atn(1)
End Function

Public Function Divisible%(i%)
Divisible = Round_Down((i Mod 30), 0)
End Function
