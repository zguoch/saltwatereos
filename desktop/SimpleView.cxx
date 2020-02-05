/*
 * Copyright 2007 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */


#include "ui_SimpleView.h"
#include "SimpleView.h"

#include <vtkDataObjectToTable.h>
#include <vtkElevationFilter.h>
#include "vtkGenericOpenGLRenderWindow.h"
#include <vtkNew.h>
#include <vtkPolyDataMapper.h>
#include <vtkQtTableView.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include "vtkSmartPointer.h"
#include <vtkVectorText.h>

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkLookupTable.h>
#include <vtkNamedColors.h>
#include <vtkPolyData.h>
#include <vtkPolyhedron.h>
#include <vtkProperty.h>
#include <vtkRenderWindowInteractor.h>

#include "vtkConeSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderWindow.h"
#include "vtkCamera.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include <vtkRenderWindowInteractor.h>

namespace
{
vtkSmartPointer<vtkPolyhedron> MakeDodecahedron();
}
vtkSmartPointer<vtkActor> test()
{
  vtkSmartPointer<vtkNamedColors> colors =
    vtkSmartPointer<vtkNamedColors>::New();

  vtkSmartPointer<vtkPolyhedron> dodecahedron = MakeDodecahedron();

  // Visualize
  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData(dodecahedron->GetPolyData());

  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetColor(
    colors->GetColor3d("PapayaWhip").GetData());
  return actor;
  // vtkSmartPointer<vtkRenderer> renderer =
  //   vtkSmartPointer<vtkRenderer>::New();
  // vtkSmartPointer<vtkRenderWindow> renderWindow =
  //   vtkSmartPointer<vtkRenderWindow>::New();
  // renderWindow->SetWindowName("Dodecahedron");
  // renderWindow->AddRenderer(renderer);
  // vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
  //   vtkSmartPointer<vtkRenderWindowInteractor>::New();
  // renderWindowInteractor->SetRenderWindow(renderWindow);

  // renderer->AddActor(actor);
  // renderer->SetBackground(colors->GetColor3d("CadetBlue").GetData());
  // renderer->GetActiveCamera()->Azimuth(30);
  // renderer->GetActiveCamera()->Elevation(30);

  // renderer->ResetCamera();

  // renderWindow->Render();
  // renderWindowInteractor->Start();
}
// Constructor
SimpleView::SimpleView()
{
  this->ui = new Ui_SimpleView;
  this->ui->setupUi(this);

  // Qt Table View
  this->TableView = vtkSmartPointer<vtkQtTableView>::New();

  // Place the table view in the designer form
  // this->ui->layout()->addWidget(this->TableView->GetWidget());

  // Geometry
  vtkNew<vtkVectorText> text;
  text->SetText("Salt Water EOS");
  vtkNew<vtkElevationFilter> elevation;
  elevation->SetInputConnection(text->GetOutputPort());
  elevation->SetLowPoint(0,0,0);
  elevation->SetHighPoint(10,0,0);

  // Mapper
  vtkNew<vtkPolyDataMapper> mapper;
  mapper->SetInputConnection(elevation->GetOutputPort());

  // Actor in scene
  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper);

  // VTK Renderer
  vtkNew<vtkRenderer> ren;


vtkConeSource *cone = vtkConeSource::New();
  cone->SetHeight( 3.0 );
  cone->SetRadius( 1.0 );
  cone->SetResolution( 10 );

  vtkPolyDataMapper *coneMapper = vtkPolyDataMapper::New();
  coneMapper->SetInputConnection( cone->GetOutputPort() );

  vtkActor *coneActor = vtkActor::New();
  coneActor->SetMapper( coneMapper );


  // Add Actor to renderer
  ren->AddActor(actor);
  // ren->AddActor(test());
  ren->AddActor(coneActor);

  ren->SetBackground( 0.1, 0.2, 0.4 );


  // VTK/Qt wedded
  vtkNew<vtkGenericOpenGLRenderWindow> renderWindow;
  this->ui->qvtkWidget->SetRenderWindow(renderWindow);
  this->ui->qvtkWidget->GetRenderWindow()->AddRenderer(ren);

  // Just a bit of Qt interest: Culling off the
  // point data and handing it to a vtkQtTableView
  // vtkNew<vtkDataObjectToTable> toTable;
  // toTable->SetInputConnection(elevation->GetOutputPort());
  // toTable->SetFieldType(vtkDataObjectToTable::POINT_DATA);

  // Here we take the end of the VTK pipeline and give it to a Qt View
  // this->TableView->SetRepresentationFromInputConnection(toTable->GetOutputPort());

  // Set up action signals and slots
  connect(this->ui->actionOpenFile, SIGNAL(triggered()), this, SLOT(slotOpenFile()));
  connect(this->ui->actionExit, SIGNAL(triggered()), this, SLOT(slotExit()));

};

SimpleView::~SimpleView()
{
  // The smart pointers should clean up for up

}

// Action to be taken upon file open
void SimpleView::slotOpenFile()
{

}

void SimpleView::slotExit() {
  qApp->exit();
}

namespace
{

  vtkSmartPointer<vtkPolyhedron>MakeDodecahedron()
  {
    vtkSmartPointer<vtkPolyhedron> aDodecahedron =
      vtkSmartPointer<vtkPolyhedron>::New();

    for (int i = 0; i < 20; ++i)
    {
      aDodecahedron->GetPointIds()->InsertNextId(i);
    }

    aDodecahedron->GetPoints()->InsertNextPoint(1.21412,    0,          1.58931);
    aDodecahedron->GetPoints()->InsertNextPoint(0.375185,   1.1547,     1.58931);
    aDodecahedron->GetPoints()->InsertNextPoint(-0.982247,  0.713644,   1.58931);
    aDodecahedron->GetPoints()->InsertNextPoint(-0.982247,  -0.713644,  1.58931);
    aDodecahedron->GetPoints()->InsertNextPoint(0.375185,   -1.1547,    1.58931);
    aDodecahedron->GetPoints()->InsertNextPoint(1.96449,    0,          0.375185);
    aDodecahedron->GetPoints()->InsertNextPoint(0.607062,   1.86835,    0.375185);
    aDodecahedron->GetPoints()->InsertNextPoint(-1.58931,   1.1547,     0.375185);
    aDodecahedron->GetPoints()->InsertNextPoint(-1.58931,   -1.1547,    0.375185);
    aDodecahedron->GetPoints()->InsertNextPoint(0.607062,   -1.86835,   0.375185);
    aDodecahedron->GetPoints()->InsertNextPoint(1.58931,    1.1547,     -0.375185);
    aDodecahedron->GetPoints()->InsertNextPoint(-0.607062,  1.86835,    -0.375185);
    aDodecahedron->GetPoints()->InsertNextPoint(-1.96449,   0,          -0.375185);
    aDodecahedron->GetPoints()->InsertNextPoint(-0.607062,  -1.86835,   -0.375185);
    aDodecahedron->GetPoints()->InsertNextPoint(1.58931,    -1.1547,    -0.375185);
    aDodecahedron->GetPoints()->InsertNextPoint(0.982247,   0.713644,   -1.58931);
    aDodecahedron->GetPoints()->InsertNextPoint(-0.375185,  1.1547,     -1.58931);
    aDodecahedron->GetPoints()->InsertNextPoint(-1.21412,   0,          -1.58931);
    aDodecahedron->GetPoints()->InsertNextPoint(-0.375185,  -1.1547,    -1.58931);
    aDodecahedron->GetPoints()->InsertNextPoint(0.982247,   -0.713644,  -1.58931);

    vtkIdType faces[73] =
      {12,                   // number of faces
      5, 0, 1, 2, 3, 4,     // number of ids on face, ids
      5, 0, 5, 10, 6, 1,
      5, 1, 6, 11, 7, 2,
      5, 2, 7, 12, 8, 3,
      5, 3, 8, 13, 9, 4,
      5, 4, 9, 14, 5, 0,
      5, 15, 10, 5, 14, 19,
      5, 16, 11, 6, 10, 15,
      5, 17, 12, 7, 11, 16,
      5, 18, 13, 8, 12, 17,
      5, 19, 14, 9, 13, 18,
      5, 19, 18, 17, 16, 15};

    aDodecahedron->SetFaces(faces);
    aDodecahedron->Initialize();

    return aDodecahedron;
  }
}