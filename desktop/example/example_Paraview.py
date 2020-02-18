
from paraview.simple import *

# create a new 'Legacy VTK Reader'
xHvtk = LegacyVTKReader(FileNames=['XH.vtk'])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

xHvtkDisplay = Show(xHvtk, renderView1)

xHvtkDisplay.Representation = 'Surface'

renderView1.AxesGrid.Visibility = 1
 
xHvtkDisplay.Scale = [1.0, 2.0, 1.0]
renderView1.AxesGrid.DataScale = [1.0, 2, 1.0]

renderView1.AxesGrid.XTitle = 'Enthalpy (kJ/kg)'
renderView1.AxesGrid.YTitle = 'Salinity'
renderView1.AxesGrid.ZTitle = 'Pressure (bar)'

renderView1.AxesGrid.XTitleFontSize = 16
renderView1.AxesGrid.XTitleBold = 1
renderView1.AxesGrid.YTitleFontSize = 16
renderView1.AxesGrid.YTitleBold = 1
# renderView1.AxesGrid.ZTitleFontSize = 16
# renderView1.AxesGrid.ZYTitleBold = 1

renderView1.ResetCamera()