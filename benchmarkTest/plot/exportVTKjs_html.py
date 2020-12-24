#!/Applications/ParaView5.9.0-RC1.app/Contents/bin/pvpython
# -*-coding:utf-8-*-
# trace generated using paraview version 5.9.0-RC1
import sys
import os
#### import the simple module from the paraview
from paraview.simple import *

def export2VTKJS_HTML(path_vtkjs):
    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    # load state
    LoadState('model_geometry_geology_exportVTKJS.pvsm', DataDirectory='.')

    # find view
    renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')

    # set active view
    SetActiveView(renderView1)

    # export view
    ExportView(path_vtkjs, view=renderView1, ParaViewGlanceHTML='/Applications/ParaView5.9.0-RC1.app/Contents/Resources/web/glance/ParaViewGlance.html')

def main(argv):
    argc=len(argv)
    if(not argc==2):
        print('Example: ./exportVTKjs_html.py my.vtkjs')
        exit(0)
    path_vtkjs=argv[1]
    os.system('rm %s'%(path_vtkjs.replace('.vtkjs','.*')))
    export2VTKJS_HTML(path_vtkjs)
    os.system('rm %s'%(path_vtkjs))
if __name__ == '__main__':
    sys.exit(main(sys.argv))