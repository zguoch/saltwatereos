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

    m_IndependentVar1_old=ui->doubleSpinBox->value()*1e5;
    m_IndependentVar2_old=ui->doubleSpinBox_2->value();
    m_IndependentVar3_old=ui->doubleSpinBox_3->value();

    //default set 1D UI
    update1dUI("Temperature");


    m_vtkTable=vtkSmartPointer<vtkTable>::New();
    m_structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
    // ------------renderwindow------------
    m_renderWindow=vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    this->ui->qvtkWidget2->SetRenderWindow(m_renderWindow);
    initRenderWindow();

    // ------2D chart---------
    // Set up the view
  m_vtkChartView=vtkSmartPointer<vtkContextView>::New();
  m_vtkChartView->GetRenderer()->SetBackground(1.0, 1.0, 1.0);
  vtkSmartPointer<vtkChartXY> chart = 
    vtkSmartPointer<vtkChartXY>::New();
  m_vtkChartView->GetScene()->AddItem(chart);
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
void MainWindow::initRenderWindow()
{
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
    this->ui->qvtkWidget2->GetRenderWindow()->AddRenderer(ren);
    this->ui->qvtkWidget2->GetRenderWindow()->Render();
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
              ",&nbsp;&nbsp;&nbsp; <font color=Blue>Bulk viscosity: </font> &nbsp; "+QString::number(eos.m_prop.Mu)+"<br>";
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
        QString varName;
        int index_var=0;
        vector<double>arrP,arrT,arrX;
        switch (ui->comboBox_2->currentIndex()) {
        case 0:   //temperature
        {
            varName="Temperature (C)";
            index_var=0;
            double dT=ui->doubleSpinBox_12->value();
            double Tmin=ui->doubleSpinBox_10->value();
            double Tmax=ui->doubleSpinBox_11->value();
            double p0=ui->doubleSpinBox_4->value()*1e5;
            double X0=ui->doubleSpinBox_6->value();
            for (double T=Tmin;T<Tmax;T=T+dT) {
                arrT.push_back(T);
                arrP.push_back(p0);
                arrX.push_back(X0);
            }

        }
            break;
        case 1:   //pressure
        {
            varName="Pressure (bar)";
            index_var=1;
            double dP=ui->doubleSpinBox_12->value()*1e5;
            double Pmin=ui->doubleSpinBox_10->value()*1e5;
            double Pmax=ui->doubleSpinBox_11->value()*1e5;
            double T0=ui->doubleSpinBox_4->value();
            double X0=ui->doubleSpinBox_6->value();
            for (double P=Pmin;P<Pmax;P=P+dP) {
                arrT.push_back(T0);
                arrP.push_back(P);
                arrX.push_back(X0);
            }
        }
            break;
        case 2:  //salinity
        {
            varName="Salinity";
            index_var=2;
            double dX=ui->doubleSpinBox_12->value();
            double Xmin=ui->doubleSpinBox_10->value();
            double Xmax=ui->doubleSpinBox_11->value();
            double P0=ui->doubleSpinBox_4->value()*1e5;
            double T0=ui->doubleSpinBox_6->value();
            for (double X=Xmin;X<Xmax;X=X+dX) {
                arrT.push_back(T0);
                arrP.push_back(P0);
                arrX.push_back(X);
            }
        }
            break;
        }
        // Create a table with some points in it
    //          vtkSmartPointer<vtkTable> table =
    //            table=vtkSmartPointer<vtkTable>::New();
        vtkSmartPointer<vtkFloatArray> varT =
          vtkSmartPointer<vtkFloatArray>::New();
        varT->SetName("Temperature(C)");
        m_vtkTable->AddColumn(varT);

        vtkSmartPointer<vtkFloatArray> varP =
          vtkSmartPointer<vtkFloatArray>::New();
        varP->SetName("Pressure (bar)");
        m_vtkTable->AddColumn(varP);

        vtkSmartPointer<vtkFloatArray> varX =
          vtkSmartPointer<vtkFloatArray>::New();
        varX->SetName("Salinity");
        m_vtkTable->AddColumn(varX);

        vtkSmartPointer<vtkFloatArray> prop_density =
          vtkSmartPointer<vtkFloatArray>::New();
        prop_density->SetName("Density (kg/m3)");
        m_vtkTable->AddColumn(prop_density);

        vtkSmartPointer<vtkFloatArray> prop_enthalpy =
          vtkSmartPointer<vtkFloatArray>::New();
        prop_enthalpy->SetName("Enthalpy (J/kg)");
        m_vtkTable->AddColumn(prop_enthalpy);

        vtkSmartPointer<vtkFloatArray> prop_viscosity =
          vtkSmartPointer<vtkFloatArray>::New();
        prop_viscosity->SetName("Viscosity (Pa s)");
        m_vtkTable->AddColumn(prop_viscosity);

        m_vtkTable->SetNumberOfRows(arrT.size());
        for (size_t i=0;i<arrT.size();++i) {
            SWEOS::cH2ONaCl eos(arrP[i],arrT[i],arrX[i]);
            eos.Calculate();
            m_vtkTable->SetValue(i, 0, arrT[i]);
            m_vtkTable->SetValue(i, 1, arrP[i]);
            m_vtkTable->SetValue(i, 2, arrX[i]);
            m_vtkTable->SetValue(i, 3, eos.m_prop.Rho);
            m_vtkTable->SetValue(i, 4, eos.m_prop.H);
            m_vtkTable->SetValue(i, 5, eos.m_prop.Mu_l);
        }

        int index_prop=ui->comboBox_3->currentIndex()+3;

          int num_oldItems=m_vtkChartView->GetScene()->GetNumberOfItems();
          for (int i=0;i<num_oldItems;i++) {
              m_vtkChartView->GetScene()->RemoveItem(m_vtkChartView->GetScene()->GetItem(0));
          }

           vtkSmartPointer<vtkChartXY> chart =
             vtkSmartPointer<vtkChartXY>::New();
           m_vtkChartView->GetScene()->AddItem(chart);

           vtkPlot *line = chart->AddPlot(vtkChart::LINE);
           line->SetInputData(m_vtkTable, index_var, index_prop);
           line->SetColor(0, 255, 0, 255);
           line->SetWidth(3.0);

            //            chart->GetAxis(1)->SetRange(0, 11);
            //            chart->GetAxis(1)->SetBehavior(vtkAxis::FIXED);
            chart->GetAxis(1)->SetTitle(m_vtkTable->GetColumn(index_var)->GetName());

            chart->SetShowLegend(true);

            chart->GetAxis(0)->SetTitle(m_vtkTable->GetColumn(index_prop)->GetName());

          m_vtkChartView->GetRenderWindow()->Render();

          //print results into textedit box
    //          ui->textEdit->clear();
    }

        break;
    case 2:
    {
        double X0 = ui->doubleSpinBox_6->value();
        double Tmin=ui->doubleSpinBox_10->value();
        double Tmax=ui->doubleSpinBox_11->value();
        double dT=ui->doubleSpinBox_12->value();
        double Pmin=ui->doubleSpinBox_15->value();
        double Pmax=ui->doubleSpinBox_14->value();
        double dP=ui->doubleSpinBox_13->value();
        vector<double> vectorT, vectorP;
        for (double T=Tmin; T<Tmax;T=T+dT)
        {
          vectorT.push_back(T);
        }
        for (double P=Pmin; P<Pmax; P=P+dP)
        {
          vectorP.push_back(P);
        }
        
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

          // unsigned int gridSize = 8;
          // unsigned int counter = 0;
          // Create a 5x5 grid of points
          for(size_t j = 0; j < vectorT.size(); j++)
          {
            for(size_t i = 0; i < vectorP.size(); i++)
            {
              SWEOS::cH2ONaCl eos(vectorP[i],vectorT[j],X0);
              eos.Calculate();
              points->InsertNextPoint(i, j, eos.m_prop.Rho);
              // if(i == 3 && j == 3) // Make one point higher than the rest
              // {
              //   points->InsertNextPoint(i, j, 2);
              //   std::cout << "The different point is number " << counter << std::endl;
              // }
              // else
              // {
              //   points->InsertNextPoint(i, j, 0); // Make most of the points the same height
              // }
//              counter++;
            }
          }

          // Specify the dimensions of the grid
          m_structuredGrid->SetDimensions(vectorP.size(),vectorT.size(),1);

          m_structuredGrid->SetPoints(points);

          // m_structuredGrid->BlankPoint(27);
          // m_structuredGrid->Modified();

          // Create a mapper and actor
          vtkSmartPointer<vtkDataSetMapper> gridMapper = vtkSmartPointer<vtkDataSetMapper>::New();
          gridMapper->SetInputData(m_structuredGrid);

          vtkSmartPointer<vtkActor> gridActor = vtkSmartPointer<vtkActor>::New();
          gridActor->SetMapper(gridMapper);
          gridActor->GetProperty()->EdgeVisibilityOn();
          gridActor->GetProperty()->SetEdgeColor(0,0,1);

          // Create a renderer, render window, and interactor
          vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();

          // Add the actor to the scene
          renderer->AddActor(gridActor);
          renderer->SetBackground(.3, .6, .3); // Background color green

        // before adding new renderer, remove all the old renderer, always keep only renderer in m_renderwindow
        m_renderWindow->GetRenderers()->RemoveAllItems();
        this->ui->qvtkWidget2->GetRenderWindow()->AddRenderer(renderer);
        // Render and interact
        this->ui->qvtkWidget2->GetRenderWindow()->Render();
        // ui->textEdit->append(QString::number(m_renderWindow->GetRenderers()->GetNumberOfItems()));
    }
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
        ui->vtkWindowTab->setCurrentIndex(0);
    }
}

void MainWindow::on_radioButton_4_clicked()
{
    if(ui->radioButton_4->isChecked())
    {
        m_dimension=2;
        ui->vtkWindowTab->setCurrentIndex(1);
    }
}

void MainWindow::on_radioButton_5_clicked()
{
    if(ui->radioButton_5->isChecked())
    {
        m_dimension=3;
        ui->vtkWindowTab->setCurrentIndex(1);
    }
}

void MainWindow::update1dUI(QString arg)
{
    if(arg=="Temperature")
    {
        ui->label_5->setText("Pressure (bar)");
        ui->doubleSpinBox_4->setDecimals(2);
        ui->doubleSpinBox_4->setRange(SWEOS::PMIN/1e5, SWEOS::PMAX/1e5);
        ui->doubleSpinBox_4->setSingleStep(1);
        ui->doubleSpinBox_4->setValue(316);
        ui->label_7->setText("Salinity");
        ui->doubleSpinBox_6->setDecimals(4);
        ui->doubleSpinBox_6->setRange(SWEOS::XMIN, SWEOS::XMAX);
        ui->doubleSpinBox_6->setSingleStep(0.001);
        ui->doubleSpinBox_6->setValue(0.032);

        ui->label_14->setText("dT");
        ui->doubleSpinBox_10->setDecimals(2);
        ui->doubleSpinBox_10->setRange(SWEOS::TMIN, SWEOS::TMAX);
        ui->doubleSpinBox_10->setSingleStep(1);
        ui->doubleSpinBox_10->setValue(SWEOS::TMIN);

        ui->doubleSpinBox_11->setDecimals(2);
        ui->doubleSpinBox_11->setRange(SWEOS::TMIN, SWEOS::TMAX);
        ui->doubleSpinBox_11->setSingleStep(1);
        ui->doubleSpinBox_11->setValue(SWEOS::TMAX);

        ui->doubleSpinBox_12->setDecimals(2);
        ui->doubleSpinBox_12->setRange(0.1,100);
        ui->doubleSpinBox_12->setSingleStep(1);
        ui->doubleSpinBox_12->setValue(1);

    }else if(arg=="Pressure")
    {
        ui->label_5->setText("Temperature (C)");
        ui->doubleSpinBox_4->setDecimals(2);
        ui->doubleSpinBox_4->setRange(SWEOS::TMIN, SWEOS::TMAX);
        ui->doubleSpinBox_4->setSingleStep(1);
        ui->doubleSpinBox_4->setValue(100);
        ui->label_7->setText("Salinity");
        ui->doubleSpinBox_6->setDecimals(4);
        ui->doubleSpinBox_6->setRange(SWEOS::XMIN, SWEOS::XMAX);
        ui->doubleSpinBox_6->setSingleStep(0.001);
        ui->doubleSpinBox_6->setValue(0.032);

        ui->label_14->setText("dP (bar)");
        ui->doubleSpinBox_10->setDecimals(2);
        ui->doubleSpinBox_10->setRange(SWEOS::PMIN/1E5, SWEOS::PMAX/1E5);
        ui->doubleSpinBox_10->setSingleStep(1);
        ui->doubleSpinBox_10->setValue(SWEOS::PMIN/1E5);

        ui->doubleSpinBox_11->setDecimals(2);
        ui->doubleSpinBox_11->setRange(SWEOS::PMIN/1E5, SWEOS::PMAX/1E5);
        ui->doubleSpinBox_11->setSingleStep(1);
        ui->doubleSpinBox_11->setValue(SWEOS::PMAX/1E5);

        ui->doubleSpinBox_12->setDecimals(2);
        ui->doubleSpinBox_12->setRange(0.1,100);
        ui->doubleSpinBox_12->setSingleStep(1);
        ui->doubleSpinBox_12->setValue(1);

    }else if(arg=="Salinity")
    {
        ui->label_5->setText("Pressure (bar)");
        ui->doubleSpinBox_4->setDecimals(2);
        ui->doubleSpinBox_4->setRange(SWEOS::PMIN/1e5, SWEOS::PMAX/1e5);
        ui->doubleSpinBox_4->setSingleStep(1);
        ui->doubleSpinBox_4->setValue(316);
        ui->label_7->setText("Temperature (C)");
        ui->doubleSpinBox_6->setDecimals(2);
        ui->doubleSpinBox_6->setRange(SWEOS::TMIN, SWEOS::TMAX);
        ui->doubleSpinBox_6->setSingleStep(1);
        ui->doubleSpinBox_6->setValue(100);

        ui->label_14->setText("dX");
        ui->doubleSpinBox_10->setDecimals(4);
        ui->doubleSpinBox_10->setRange(SWEOS::XMIN, SWEOS::XMAX);
        ui->doubleSpinBox_10->setSingleStep(0.001);
        ui->doubleSpinBox_10->setValue(0.0001);

        ui->doubleSpinBox_11->setDecimals(4);
        ui->doubleSpinBox_11->setRange(SWEOS::XMIN, SWEOS::XMAX);
        ui->doubleSpinBox_11->setSingleStep(0.001);
        ui->doubleSpinBox_11->setValue(SWEOS::XMAX);

        ui->doubleSpinBox_12->setDecimals(4);
        ui->doubleSpinBox_12->setRange(0.0001,0.9);
        ui->doubleSpinBox_12->setSingleStep(0.001);
        ui->doubleSpinBox_12->setValue(0.1);
    }else
    {
        std::cout<<"error: update1dUI, no such item: "<<arg.toStdString()<<std::endl;
    }
}
void MainWindow::on_comboBox_2_activated(const QString &arg1)
{
    switch (m_dimension) {
        case 1:
            update1dUI(arg1);
        break;

    }
}
