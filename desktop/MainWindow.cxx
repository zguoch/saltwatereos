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

#include <vtkRendererCollection.h>


// Constructor
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , m_calculationMode(CALCULATION_SINGLE_POINT)
    , m_dimension(1)
    ,m_calculationMode_123Dim(oneDim_Temperature)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

//  // Qt Table View
//  this->TableView = vtkSmartPointer<vtkQtTableView>::New();

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

// ------2D line


  // Set up the view
//  vtkSmartPointer<vtkContextView> view =
//    vtkSmartPointer<vtkContextView>::New();
  m_vtkChartView=vtkSmartPointer<vtkContextView>::New();
  m_vtkChartView->GetRenderer()->SetBackground(1.0, 1.0, 1.0);
  // Add multiple line plots, setting the colors etc
  vtkSmartPointer<vtkChartXY> chart = 
    vtkSmartPointer<vtkChartXY>::New();
  m_vtkChartView->GetScene()->AddItem(chart);
//  vtkPlot *line = chart->AddPlot(vtkChart::LINE);
//  line->SetInputData(table, 0, 1);
//  line->SetColor(0, 255, 0, 255);
//  line->SetWidth(1.0);
//  line = chart->AddPlot(vtkChart::LINE);
//  line->SetInputData(table, 0, 2);
//  line->SetColor(255, 0, 0, 255);
//  line->SetWidth(5.0);
  // VTK/Qt wedded
  vtkNew<vtkGenericOpenGLRenderWindow> renderWindow;
  this->ui->qvtkWidget->SetRenderWindow(renderWindow);
//this->ui->qvtkWidget->GetRenderWindow()->AddRenderer(ren);

  // 2. key for 2D charts
   m_vtkChartView->SetRenderWindow(this->ui->qvtkWidget->GetRenderWindow()); //must be here




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

void MainWindow::on_tabWidget_currentChanged(int index)
{
    switch (index) {
    case 0:
        ui->textEdit->setEnabled(true);
//        ui->qvtkWidget->setEnabled(false);
        break;
    case 1:
        ui->textEdit->setEnabled(false);
//        ui->qvtkWidget->setEnabled(true);
        break;

    }

}

void MainWindow::on_pushButton_clicked()
{
    switch (m_dimension) {
    case 1:
    {
        QString varName="Temperature (C)";
        QString propName="Density (kg/m3)";
        double p=ui->doubleSpinBox_4->value()*1e5;
        double X=ui->doubleSpinBox_6->value();
        double Tmin=ui->doubleSpinBox_10->value();
        double Tmax=ui->doubleSpinBox_11->value();
        double dT=ui->doubleSpinBox_12->value();
        int numPoints=(int)((Tmax-Tmin)/dT+1);


        // Create a table with some points in it
          vtkSmartPointer<vtkTable> table =
            vtkSmartPointer<vtkTable>::New();

          vtkSmartPointer<vtkFloatArray> arrX =
            vtkSmartPointer<vtkFloatArray>::New();
          arrX->SetName("Temperature (C)");
          table->AddColumn(arrX);

          vtkSmartPointer<vtkFloatArray> arrC =
            vtkSmartPointer<vtkFloatArray>::New();
          arrC->SetName("Density (kg/m3)");
          table->AddColumn(arrC);

//          vtkSmartPointer<vtkFloatArray> arrS =
//            vtkSmartPointer<vtkFloatArray>::New();
//          arrS->SetName("Sine");
//          table->AddColumn(arrS);

          // Fill in the table with some example values
//          int numPoints = 69;
//          float inc = 7.5 / (numPoints-1);
          table->SetNumberOfRows(numPoints);
//          for (int i = 0; i < numPoints; ++i)
//          {
//            table->SetValue(i, 0, i * inc);
//            table->SetValue(i, 1, cos(i * inc));
//            table->SetValue(i, 2, sin(i * inc));
//          }
//          vtkSmartPointer<vtkChartXY> chart;
//          m_vtkChartLine->SetInputData(table);
//          m_vtkChartView->GetRenderer()->SetBackground(1.0, 0, 1.0);
//           m_vtkChartView->GetRenderer()->Render();
          for (int i=0;i<numPoints;++i) {
              double T=Tmin+i*dT;
              SWEOS::cH2ONaCl eos(p, T, X);
              eos.Calculate();
              table->SetValue(i, 0, T);
              table->SetValue(i, 1, eos.m_prop.Rho);
          }
          int num_oldItems=m_vtkChartView->GetScene()->GetNumberOfItems();
          for (int i=0;i<num_oldItems;i++) {
              m_vtkChartView->GetScene()->RemoveItem(m_vtkChartView->GetScene()->GetItem(0));
          }

           vtkSmartPointer<vtkChartXY> chart =
             vtkSmartPointer<vtkChartXY>::New();
           m_vtkChartView->GetScene()->AddItem(chart);

//           (vtkChartXY*) chart0=m_vtkChartView->GetScene()->GetItem(0);
           vtkPlot *line = chart->AddPlot(vtkChart::LINE);
           line->SetInputData(table, 0, 1);
           line->SetColor(0, 255, 0, 255);
           line->SetWidth(1.0);
//        ui->textEdit->append(QString::number(m_vtkChartView->GetScene()->GetNumberOfItems()));
          m_vtkChartView->GetRenderWindow()->Render();
    }
        ui->textEdit->append("1D is comming soon");

        break;
    case 2:
        ui->textEdit->append("2D is comming soon");
        break;
    case 3:
        ui->textEdit->append("3D is comming soon");
        break;

    }
}

void MainWindow::on_radioButton_3_clicked()
{
    if(ui->radioButton_3->isChecked())
    {
        m_dimension=1;
    }
}

void MainWindow::on_radioButton_4_clicked()
{
    if(ui->radioButton_4->isChecked())
    {
        m_dimension=2;
    }
}

void MainWindow::on_radioButton_5_clicked()
{
    if(ui->radioButton_5->isChecked())
    {
        m_dimension=3;
    }
}
