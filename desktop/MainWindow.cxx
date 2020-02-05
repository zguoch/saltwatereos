/*
 * Copyright 2007 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */

#include "MainWindow.h"
#include "ui_MainWindow.h"

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

// Constructor
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , m_calculationMode(CALCULATION_SINGLE_POINT)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

  // Qt Table View
  this->TableView = vtkSmartPointer<vtkQtTableView>::New();

  // Place the table view in the designer form
    m_IndependentVar1_old=ui->doubleSpinBox->value()*1e5;
    m_IndependentVar2_old=ui->doubleSpinBox_2->value();
    m_IndependentVar3_old=ui->doubleSpinBox_3->value();
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

  // Add Actor to renderer
  ren->AddActor(actor);
  // ren->AddActor(test());

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

MainWindow::~MainWindow()
{
  // The smart pointers should clean up for up

}

// Action to be taken upon file open
void MainWindow::slotOpenFile()
{

}

void MainWindow::slotExit() {
  qApp->exit();
}


void MainWindow::on_pushButton_2_clicked()
{
     double p=ui->doubleSpinBox->value()*1e5;
     double T=ui->doubleSpinBox_2->value();
     double X=ui->doubleSpinBox_3->value();
     QString color_value_1="Black", color_value_2="Black", color_value_3="Black";
     color_value_1=(m_IndependentVar1_old==p ? color_value_1="Black": color_value_1="Green");
     color_value_2=(m_IndependentVar2_old==T ? color_value_2="Black": color_value_2="Green");
     color_value_3=(m_IndependentVar3_old==X ? color_value_3="Black": color_value_3="Green");

    SWEOS::cH2ONaCl eos(p, T, X);
    eos.Calculate();

    QString result_str;
    // T, P, X
    result_str="<font color=Purple>Pressure</font> = <font color="+color_value_1+">"+QString::number(p)
            +"</font> Pa, <font color=Purple>Perssure</font> = <font color="+color_value_2+">"+QString::number(T)
            +"</font> deg. C, <font color=Purple>Salinity</font>  = <font color="+color_value_3+">"+QString::number(X*100)
            +"</font> wt. % NaCl<br>";
    result_str+="======================================================================<br>";
    // Region
    result_str+="<font color=Blue>Phase Region</font>: "+QString::fromStdString(eos.m_phaseRegion_name[eos.m_prop.Region])+"<br>";
    // Xl, Xv
    result_str+="<font color=Blue>Xl: </font> &nbsp; "+QString::number(eos.m_prop.X_l)+
              ",&nbsp;&nbsp;&nbsp; <font color=Blue>Xv: </font> &nbsp; "+QString::number(eos.m_prop.X_v)+"<br>";
    // Bulk rho, bulk H, bulk mu
    result_str+="<font color=Blue>Bulk density: </font> &nbsp; "+QString::number(eos.m_prop.Rho)+
              ",&nbsp;&nbsp;&nbsp; <font color=Blue>Bulk enthalpy: </font> &nbsp; "+QString::number(eos.m_prop.H)+
              ",&nbsp;&nbsp;&nbsp; <font color=Blue>Bulk viscosity: </font> &nbsp; "+QString::number(eos.m_prop.H)+"<br>";
    // rhol, rhov, rhoh
    result_str+="<font color=Blue>Liquid density: </font> &nbsp; "+QString::number(eos.m_prop.Rho_l)+
              ",&nbsp;&nbsp;&nbsp; <font color=Blue>Vapour density: </font> &nbsp; "+QString::number(eos.m_prop.Rho_v)+
              ",&nbsp;&nbsp;&nbsp; <font color=Blue>Halite density: </font> &nbsp; "+QString::number(eos.m_prop.Rho_h)+"<br>";

    // Hl, Hv, Hh
    result_str+="<font color=Blue>Liquid enthalpy: </font> &nbsp; "+QString::number(eos.m_prop.H_l)+
              ",&nbsp;&nbsp;&nbsp; <font color=Blue>Vapour enthalpy: </font> &nbsp; "+QString::number(eos.m_prop.H_v)+
              ",&nbsp;&nbsp;&nbsp; <font color=Blue>Halite enthalpy: </font> &nbsp; "+QString::number(eos.m_prop.H_h)+"<br>";
    // Mul, Muv
    result_str+="<font color=Blue>Liquid viscosity: </font> &nbsp; "+QString::number(eos.m_prop.Mu_l)+
              ",&nbsp;&nbsp;&nbsp; <font color=Blue>Vapour viscosity: </font> &nbsp; "+QString::number(eos.m_prop.Mu_v)+"<br>";
    // Time stamp
    QDateTime current_date_time =QDateTime::currentDateTime();
    QString current_date =current_date_time.toString("yyyy.MM.dd hh:mm:ss ddd");
    result_str+="<br> ----------------------------------- <font color=Grey>"+current_date+" </font> ----------------------------------- <br><br>";


    ui->textEdit->setText(ui->textEdit->toPlainText());
    ui->textEdit->moveCursor(QTextCursor::Start);
    ui->textEdit->insertHtml(result_str);


    //update old value
    m_IndependentVar1_old=p;
    m_IndependentVar2_old=T;
    m_IndependentVar3_old=X;
}

void MainWindow::on_radioButton_pressed()
{

}

void MainWindow::updateCalculationModelSelection(bool isSinglePoint)
{
    ui->comboBox->setEnabled(isSinglePoint);
    ui->doubleSpinBox->setEnabled(isSinglePoint);
    ui->doubleSpinBox_2->setEnabled(isSinglePoint);
    ui->doubleSpinBox_3->setEnabled(isSinglePoint);
}
void MainWindow::on_radioButton_clicked()
{
    if(ui->radioButton->isChecked())
    {
        m_calculationMode=CALCULATION_SINGLE_POINT;
        ui->pushButton_2->setText("Calculatet");
        updateCalculationModelSelection(true);
    }
}

void MainWindow::on_radioButton_2_clicked()
{
    if(ui->radioButton_2->isChecked())
    {
        m_calculationMode=CALCULATION_MULTI_POINTS;
        ui->pushButton_2->setText("Open File");
        updateCalculationModelSelection(false);
    }
}
