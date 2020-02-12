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
    m_geometry_Groupbox_variables=ui->groupBox_Variables->geometry();
    //update1dUI("Temperature");
    updateUI(m_dimension);


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


    // init plot
  on_pushButton_clicked();
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
//    switch (m_dimension) {
//    case 1:
//    {
//        QString varName;
//        int index_var=0;
//        vector<double>arrP,arrT,arrX;
//        switch (ui->comboBox_2->currentIndex()) {
//        case 0:   //temperature
//        {
//            varName="Temperature (C)";
//            index_var=0;
//            double dT=ui->doubleSpinBox_12->value();
//            double Tmin=ui->doubleSpinBox_10->value();
//            double Tmax=ui->doubleSpinBox_11->value();
//            double p0=ui->doubleSpinBox_4->value()*1e5;
//            double X0=ui->doubleSpinBox_6->value();
//            for (double T=Tmin;T<Tmax;T=T+dT) {
//                arrT.push_back(T);
//                arrP.push_back(p0);
//                arrX.push_back(X0);
//            }

//        }
//            break;
//        case 1:   //pressure
//        {
//            varName="Pressure (bar)";
//            index_var=1;
//            double dP=ui->doubleSpinBox_12->value()*1e5;
//            double Pmin=ui->doubleSpinBox_10->value()*1e5;
//            double Pmax=ui->doubleSpinBox_11->value()*1e5;
//            double T0=ui->doubleSpinBox_4->value();
//            double X0=ui->doubleSpinBox_6->value();
//            for (double P=Pmin;P<Pmax;P=P+dP) {
//                arrT.push_back(T0);
//                arrP.push_back(P);
//                arrX.push_back(X0);
//            }
//        }
//            break;
//        case 2:  //salinity
//        {
//            varName="Salinity";
//            index_var=2;
//            double dX=ui->doubleSpinBox_12->value();
//            double Xmin=ui->doubleSpinBox_10->value();
//            double Xmax=ui->doubleSpinBox_11->value();
//            double P0=ui->doubleSpinBox_4->value()*1e5;
//            double T0=ui->doubleSpinBox_6->value();
//            for (double X=Xmin;X<Xmax;X=X+dX) {
//                arrT.push_back(T0);
//                arrP.push_back(P0);
//                arrX.push_back(X);
//            }
//        }
//            break;
//        }
//        // Create a table with some points in it
//    //          vtkSmartPointer<vtkTable> table =
//    //            table=vtkSmartPointer<vtkTable>::New();
//        vtkSmartPointer<vtkFloatArray> varT =
//          vtkSmartPointer<vtkFloatArray>::New();
//        varT->SetName("Temperature(C)");
//        m_vtkTable->AddColumn(varT);

//        vtkSmartPointer<vtkFloatArray> varP =
//          vtkSmartPointer<vtkFloatArray>::New();
//        varP->SetName("Pressure (bar)");
//        m_vtkTable->AddColumn(varP);

//        vtkSmartPointer<vtkFloatArray> varX =
//          vtkSmartPointer<vtkFloatArray>::New();
//        varX->SetName("Salinity");
//        m_vtkTable->AddColumn(varX);

//        vtkSmartPointer<vtkFloatArray> prop_density =
//          vtkSmartPointer<vtkFloatArray>::New();
//        prop_density->SetName("Density (kg/m3)");
//        m_vtkTable->AddColumn(prop_density);

//        vtkSmartPointer<vtkFloatArray> prop_enthalpy =
//          vtkSmartPointer<vtkFloatArray>::New();
//        prop_enthalpy->SetName("Enthalpy (J/kg)");
//        m_vtkTable->AddColumn(prop_enthalpy);

//        vtkSmartPointer<vtkFloatArray> prop_viscosity =
//          vtkSmartPointer<vtkFloatArray>::New();
//        prop_viscosity->SetName("Viscosity (Pa s)");
//        m_vtkTable->AddColumn(prop_viscosity);

//        m_vtkTable->SetNumberOfRows(arrT.size());
//        for (size_t i=0;i<arrT.size();++i) {
//            SWEOS::cH2ONaCl eos(arrP[i],arrT[i],arrX[i]);
//            eos.Calculate();
//            m_vtkTable->SetValue(i, 0, arrT[i]);
//            m_vtkTable->SetValue(i, 1, arrP[i]);
//            m_vtkTable->SetValue(i, 2, arrX[i]);
//            m_vtkTable->SetValue(i, 3, eos.m_prop.Rho);
//            m_vtkTable->SetValue(i, 4, eos.m_prop.H);
//            m_vtkTable->SetValue(i, 5, eos.m_prop.Mu_l);
//        }

//        int index_prop=ui->comboBox_3->currentIndex()+3;

//          int num_oldItems=m_vtkChartView->GetScene()->GetNumberOfItems();
//          for (int i=0;i<num_oldItems;i++) {
//              m_vtkChartView->GetScene()->RemoveItem(m_vtkChartView->GetScene()->GetItem(0));
//          }

//           vtkSmartPointer<vtkChartXY> chart =
//             vtkSmartPointer<vtkChartXY>::New();
//           m_vtkChartView->GetScene()->AddItem(chart);

//           vtkNew<vtkNamedColors> colors;
//           vtkPlot *line = chart->AddPlot(vtkChart::LINE);
//           line->SetInputData(m_vtkTable, index_var, index_prop);
//           vtkColor3d color3d = colors->GetColor3d("banana");
//           line->SetColor(color3d.GetRed(), color3d.GetGreen(), color3d.GetBlue());
//           line->SetWidth(3.0);

//           int fontsize=20;

//            //            chart->GetAxis(1)->SetRange(0, 11);
//            //            chart->GetAxis(1)->SetBehavior(vtkAxis::FIXED);
//            chart->GetAxis(1)->SetTitle(m_vtkTable->GetColumn(index_var)->GetName());
//            // chart->GetAxis(1)->GetLabelProperties()->SetColor(1,0,0);
//            chart->GetAxis(1)->GetLabelProperties()->SetFontSize(fontsize);
//            chart->GetAxis(1)->GetLabelProperties()->SetFontFamilyToArial();
////            chart->GetAxis(1)->SetMinimum(273);
////            chart->GetAxis(1)->SetMaximum(573);
////            chart->GetAxis(1)->SetNumberOfTicks(5);
//            chart->GetAxis(1)->GetTitleProperties()->SetFontSize(fontsize);
//            chart->GetAxis(1)->GetTitleProperties()->SetBold(false);
//            chart->GetAxis(1)->GetTitleProperties()->SetFontFamilyToArial();
////            chart->GetAxis(1)->SetTickLabelAlgorithm(vtkAxis::TICK_SIMPLE);
////            chart->GetAxis(1)->RecalculateTickSpacing();
////            chart->GetAxis(1)->SetBehavior(2);


////            vtkAxisActor2D *axesX = vtkAxisActor2D::New();
////            axesX->SetTitle("Test");
////            axesX->SetTickLength(2);
////            axesX->SetRange(0, 10);
////            axesX->SetPoint1(0, 0);
////            axesX->SetPoint2(10, 0);
////            m_vtkChartView->AddViewProp(axesX);

//            chart->GetAxis(0)->GetLabelProperties()->SetFontSize(fontsize);
//            chart->GetAxis(0)->GetLabelProperties()->SetFontFamilyToArial();
////            chart->GetAxis(0)->SetLabelFormat("%.0f");
////            chart->GetAxis(0)->SetTickLabelAlgorithm(vtkAxis::TICK_SIMPLE);
////            chart->GetAxis(0)->SetMinimum(273);
////            chart->GetAxis(0)->SetMaximum(573);
////            chart->GetAxis(0)->SetNumberOfTicks(5);
//            chart->GetAxis(0)->GetTitleProperties()->SetFontSize(fontsize);
//            chart->GetAxis(0)->GetTitleProperties()->SetBold(false);
//            chart->GetAxis(0)->GetTitleProperties()->SetFontFamilyToArial();

//            chart->GetTooltip()->GetTextProperties()->SetFontSize(fontsize);


//            chart->SetShowLegend(true);
//            chart->GetLegend()->GetLabelProperties()->SetFontSize(fontsize);


//            chart->GetAxis(0)->SetTitle(m_vtkTable->GetColumn(index_prop)->GetName());

//          m_vtkChartView->GetRenderWindow()->Render();

//          //print results into textedit box
//    //          ui->textEdit->clear();
//    }

//        break;
//    case 2:
//    {
//        double X0 = ui->doubleSpinBox_6->value();
//        double Tmin=ui->doubleSpinBox_10->value();
//        double Tmax=ui->doubleSpinBox_11->value();
//        double dT=ui->doubleSpinBox_12->value();
//        double Pmin=ui->doubleSpinBox_15->value();
//        double Pmax=ui->doubleSpinBox_14->value();
//        double dP=ui->doubleSpinBox_13->value();
//        std::string xlabel="Temperature (K)";
//        std::string ylabel="Pressure (bar)";
//        std::string zlabel="Salinity";
//        std::string propleName="Density (kg/m3)";
//        vector<double> vectorT, vectorP;
//        for (double T=Tmin; T<Tmax;T=T+dT)
//        {
//          vectorT.push_back(T);
//        }
//        for (double P=Pmin; P<Pmax; P=P+dP)
//        {
//          vectorP.push_back(P);
//        }

//        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
//        vtkSmartPointer<vtkDoubleArray> doublevalue = vtkSmartPointer<vtkDoubleArray>::New();
//        doublevalue->SetNumberOfComponents(1);
//        doublevalue->SetName("Density");
//        for(size_t j = 0; j < vectorP.size(); j++)
//        {
//        for(size_t i = 0; i < vectorT.size(); i++)
//        {
//        SWEOS::cH2ONaCl eos(vectorP[j]*1e5,vectorT[i],X0);
//        eos.Calculate();
//        points->InsertNextPoint(vectorT[i],vectorP[j],X0);
//        doublevalue->InsertNextValue(eos.m_prop.Rho_l);
//        }
//        }
//        // Specify the dimensions of the grid
//        m_structuredGrid->SetDimensions(vectorT.size(),vectorP.size(),1);
//        m_structuredGrid->SetPoints(points);
//        m_structuredGrid->GetPointData()->SetScalars(doublevalue);


//        vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();
//        writer->SetFileName("/Users/zguo/Downloads/rho.vtk");
//        writer->SetInputData(m_structuredGrid);
//        writer->Write();


//        // Create a mapper and actor
//        vtkSmartPointer<vtkDataSetMapper> gridMapper = vtkSmartPointer<vtkDataSetMapper>::New();
//        gridMapper->SetInputData(m_structuredGrid);
//        gridMapper->SetScalarRange(m_structuredGrid->GetScalarRange());

//        // lookup table
//        vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
//        lut->SetNumberOfColors(256);
//        lut->SetHueRange(0.0,0.667);
//        //        lut->SetTableRange(100,1000);
//        //        lut->Build();
//        gridMapper->SetLookupTable(lut);
//        vtkSmartPointer<vtkActor> gridActor = vtkSmartPointer<vtkActor>::New();
//        gridActor->SetMapper(gridMapper);


//        // Create a renderer, render window, and interactor
//        vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
//        // Add the actor to the scene
//        renderer->AddActor(gridActor);
//        gridActor->SetScale((SWEOS::PMAX-SWEOS::PMIN)/1e5/(SWEOS::TMAX-SWEOS::TMIN),1,(SWEOS::PMAX-SWEOS::PMIN)/1e5/(SWEOS::XMAX-SWEOS::XMIN));

//        //axis actor
//        vtkSmartPointer<vtkCubeAxesActor> axis=vtkSmartPointer<vtkCubeAxesActor>::New();
//        axis->SetCamera(renderer->GetActiveCamera());
//        axis->SetBounds(gridActor->GetBounds());

//        //color scale
//        vtkSmartPointer<vtkScalarBarActor> scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
//          scalarBar->SetLookupTable(gridMapper->GetLookupTable());
//          scalarBar->SetLookupTable( lut );
//          scalarBar->SetTitle(propleName.c_str());
//          scalarBar->SetDrawFrame(true);
//          scalarBar->SetWidth(scalarBar->GetWidth()/2);
////          scalarBar->GetLabelTextProperty()->SetFontSize(scalarBar->GetLabelTextProperty()->GetFontSize()*6);
////          scalarBar->SetTextureGridWidth(20);
////          scalarBar->SetNumberOfLabels(4);
//          renderer->AddActor2D(scalarBar);

//        renderer->AddActor(axis);
//        InitCubeAxes(axis,vtkBoundingBox(gridActor->GetBounds()),vtkBoundingBox(m_structuredGrid->GetBounds()),xlabel,ylabel,zlabel);
//        renderer->SetBackground(0,0,0); // Background color green
//        SetCamera(renderer,vtkBoundingBox(gridActor->GetBounds()));
//        axis->SetUse2DMode(1); //set font size
//        // before adding new renderer, remove all the old renderer, always keep only renderer in m_renderwindow
//        m_renderWindow->GetRenderers()->RemoveAllItems();
//        this->ui->qvtkWidget2->GetRenderWindow()->AddRenderer(renderer);
//        renderer->Render();
//        // Render and interact
//        vtkSmartPointer<vtkInteractorStyleImage> style = vtkSmartPointer<vtkInteractorStyleImage>::New();
////        vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
//        vtkSmartPointer<vtkRenderWindowInteractor> iren=vtkSmartPointer<vtkRenderWindowInteractor>::New();
//        this->ui->qvtkWidget2->GetRenderWindow()->GetInteractor()->SetInteractorStyle(style);
//        this->ui->qvtkWidget2->GetRenderWindow()->Render();

//    }
//        ui->textEdit->append("2D is comming soon");
//        break;
//    case 3:
//        ui->textEdit->append("3D is comming soon");
//        break;

//    }
}
int MainWindow::InitCubeAxes(vtkCubeAxesActor* axes, vtkBoundingBox boundingbox, vtkBoundingBox rangebox, std::string xlabel, std::string ylabel, std::string zlabel,int fontsize)
{
    axes->SetFlyModeToClosestTriad();
//    axes->XAxisMinorTickVisibilityOff();
//    axes->YAxisMinorTickVisibilityOff();
//    axes->ZAxisMinorTickVisibilityOff();
    double bounds[6];
    boundingbox.GetBounds(bounds);
    axes->SetBounds(bounds);
    double ranges[6];
    rangebox.GetBounds(ranges);
    axes->SetXAxisRange(ranges[0], ranges[1]);//坐标轴上显示的坐标值
    axes->SetYAxisRange(ranges[2], ranges[3]);
    axes->SetZAxisRange(ranges[4], ranges[5]);
    //m_CubeAxes->SetXLabelFormat("%6.4f");
    axes->SetXTitle(xlabel.c_str());
    axes->SetYTitle(ylabel.c_str());
    axes->SetZTitle(zlabel.c_str());
    //font color
    for (int i = 0; i < 3; i++)
    {
        axes->GetTitleTextProperty(i)->SetFontSize(fontsize);
        axes->GetTitleTextProperty(i)->SetFontFamilyToTimes();
        axes->GetLabelTextProperty(i)->SetFontFamilyToTimes();
        axes->GetLabelTextProperty(i)->SetFontSize(fontsize);
    }
    return 0;
}
int MainWindow::SetCamera(vtkSmartPointer<vtkRenderer> renderer, vtkBoundingBox boundingbox, int type)
{
    double center[3];
    boundingbox.GetCenter(center);
    double bounds[6];
    boundingbox.GetBounds(bounds);
    double xlength = bounds[1] - bounds[0];
    double ylength = bounds[3] - bounds[2];
//	double zlength = bounds[5] - bounds[4];
//	double xyBiZhi = ylength / xlength;
//	double viewup_x = 0.3;
    switch (type)
    {
    case ID_CAMERA_FRONT:
        center[1] = bounds[2];
        renderer->GetActiveCamera()->SetFocalPoint(center);//焦点
        renderer->GetActiveCamera()->SetPosition(center[0], bounds[2] - ylength, center[2]);//相机位置
        renderer->GetActiveCamera()->SetViewUp(0.0, 0, 1);//相机“上”方向
        break;
    case ID_CAMERA_BACK:
        center[1] = bounds[3];
        renderer->GetActiveCamera()->SetFocalPoint(center);//焦点
        renderer->GetActiveCamera()->SetPosition(center[0], bounds[3] + ylength, center[2]);//相机位置
        renderer->GetActiveCamera()->SetViewUp(0.0, 0, 1);//相机“上”方向
        break;
    case ID_CAMERA_LEFT:
        center[0] = bounds[0];
        renderer->GetActiveCamera()->SetFocalPoint(center);//焦点
        renderer->GetActiveCamera()->SetPosition(center[0] - xlength, center[1], center[2]);//相机位置
        renderer->GetActiveCamera()->SetViewUp(0.0, 0, 1);//相机“上”方向
        break;
    case ID_CAMERA_RIGHT:
        center[0] = bounds[1];
        renderer->GetActiveCamera()->SetFocalPoint(center);//焦点
        renderer->GetActiveCamera()->SetPosition(center[0] + xlength, center[1], center[2]);//相机位置
        renderer->GetActiveCamera()->SetViewUp(0.0, 0, 1);//相机“上”方向
        break;
    case ID_CAMERA_UP:
        center[2] = bounds[5];
        renderer->GetActiveCamera()->SetFocalPoint(center);//焦点
        renderer->GetActiveCamera()->SetPosition(center[0], center[1], center[2] + 2 * (xlength > ylength ? xlength : ylength));//相机位置
        renderer->GetActiveCamera()->SetViewUp(0.0, 1, 0);//相机“上”方向
        break;
    case ID_CAMERA_DOWN:
        center[2] = bounds[4];
        renderer->GetActiveCamera()->SetFocalPoint(center);//焦点
        renderer->GetActiveCamera()->SetPosition(center[0], center[1], center[2] - 2 * (xlength > ylength ? xlength : ylength));//相机位置
        renderer->GetActiveCamera()->SetViewUp(0.0, 1, 0);//相机“上”方向
        break;
    default:
        center[2] = bounds[5];
        renderer->GetActiveCamera()->SetFocalPoint(center);//焦点
        renderer->GetActiveCamera()->SetPosition(center[0], bounds[2] - 1.5*ylength, bounds[5] + ylength / 2.0);//相机位置
        renderer->GetActiveCamera()->SetViewUp(0.0, 0, 1);//相机“上”方向
        break;
    }
    return 0;
}

void MainWindow::on_radioButton_3_clicked()
{
    if(ui->radioButton_3->isChecked())
    {
        m_dimension=1;
        ui->vtkWindowTab->setCurrentIndex(0);
        updateUI(m_dimension);
    }
}

void MainWindow::on_radioButton_4_clicked()
{
    if(ui->radioButton_4->isChecked())
    {
        m_dimension=2;
        ui->vtkWindowTab->setCurrentIndex(1);
        updateUI(m_dimension);
    }
}

void MainWindow::on_radioButton_5_clicked()
{
    if(ui->radioButton_5->isChecked())
    {
        m_dimension=3;
        ui->vtkWindowTab->setCurrentIndex(1);
        updateUI(m_dimension);
    }
}
void MainWindow::updateUI(int dim)
{
    switch (dim) {
    case 1:
        updateUILayout(false,false, true, true,true,1/3.0);
        update1dUI(ui->comboBox_selectVariable->currentText());
        break;
    case 2:
        updateUILayout(true,false, true, false,true,2/3.0);
        update2dUI(ui->comboBox_selectVariable->currentText());
        break;
    case 3:
        updateUILayout(true,true, true, false,false,1);
        update3dUI(ui->comboBox_selectVariable->currentText());
        break;
    }
}
void MainWindow::updateUILayout(bool show_secondVariable, bool show_thirdVariable, bool show_firstFixedVar, bool show_secondFixedVar,bool show_groupbox_fixedVars, double shinkWidth_groupbox_Vars)
{
    if(show_thirdVariable)
    {
        ui->comboBox_selectVariable->clear();
        ui->comboBox_selectVariable->addItem("PTX");
        ui->comboBox_selectVariable->addItem("PHX");
    }else if(show_secondVariable)
    {
        ui->comboBox_selectVariable->clear();
        ui->comboBox_selectVariable->addItem("PT");
        ui->comboBox_selectVariable->addItem("PX");
        ui->comboBox_selectVariable->addItem("TX");
    }else
    {
        ui->comboBox_selectVariable->clear();
        ui->comboBox_selectVariable->addItem("Temperature");
        ui->comboBox_selectVariable->addItem("Pressure");
        ui->comboBox_selectVariable->addItem("Salinity");
    }
    //second variable
    ui->label_max_secondVar->setVisible(show_secondVariable);
    ui->label_delta_secondVar->setVisible(show_secondVariable);
    ui->label_min_secondVar->setVisible(show_secondVariable);
    ui->doubleSpinBox_max_secondVar->setVisible(show_secondVariable);
    ui->doubleSpinBox_delta_secondVar->setVisible(show_secondVariable);
    ui->doubleSpinBox_min_secondVar->setVisible(show_secondVariable);

    //third variable
    ui->label_max_thirdVar->setVisible(show_thirdVariable);
    ui->label_delta_thirdVar->setVisible(show_thirdVariable);
    ui->label_min_thirdVar->setVisible(show_thirdVariable);
    ui->doubleSpinBox_max_thirdVar->setVisible(show_thirdVariable);
    ui->doubleSpinBox_delta_thirdVar->setVisible(show_thirdVariable);
    ui->doubleSpinBox_min_thirdVar->setVisible(show_thirdVariable);

    // groupbox of variables
    QRect geometry=ui->groupBox_Variables->geometry();
    ui->groupBox_Variables->setGeometry(geometry.x(),geometry.y(),m_geometry_Groupbox_variables.width()*shinkWidth_groupbox_Vars,geometry.height());

    // first fixed variable
    ui->label_fixed_firsVar->setVisible(show_firstFixedVar);
    ui->doubleSpinBox_fixed_firstVar->setVisible(show_firstFixedVar);

    //second fixed variable
    ui->label_fixed_secondVar->setVisible(show_secondFixedVar);
    ui->doubleSpinBox_fixed_secondVar->setVisible(show_secondFixedVar);

    // groupbox of fixed variables
    QRect geometry2=ui->groupBox_fixed_Var->geometry();
    ui->groupBox_fixed_Var->setGeometry(ui->groupBox_Variables->geometry().x()+ui->groupBox_Variables->geometry().width()+10,geometry2.y(),geometry2.width(),geometry2.height());
    ui->groupBox_fixed_Var->setVisible(show_groupbox_fixedVars);
    if(show_groupbox_fixedVars)
    {
        // button of calculation
        QRect geometry3=ui->pushButton->geometry();
        ui->pushButton->setGeometry(ui->groupBox_fixed_Var->geometry().x()+ui->groupBox_fixed_Var->geometry().width()+10, geometry3.y(),geometry3.width(),geometry3.height());
    }else
    {
        // button of calculation
        QRect geometry3=ui->pushButton->geometry();
        ui->pushButton->setGeometry(ui->groupBox_Variables->geometry().x()+ui->groupBox_Variables->geometry().width()+10, geometry3.y(),geometry3.width(),geometry3.height());
    }

}

void MainWindow::update1dUI(QString arg)
{
    if(arg=="Temperature")
    {
        //fixed vars
        UpdateUI_fixedP(ui->label_fixed_firsVar,ui->doubleSpinBox_fixed_firstVar);
        UpdateUI_fixedX(ui->label_fixed_secondVar,ui->doubleSpinBox_fixed_secondVar);
        //independent variable
        UpdateUI_T(ui->label_delta_firstVar, ui->doubleSpinBox_delta_firstVar, ui->doubleSpinBox_max_firstVar, ui->doubleSpinBox_min_firstVar);

    }else if(arg=="Pressure")
    {
        //fixed vars
        UpdateUI_fixedT(ui->label_fixed_firsVar,ui->doubleSpinBox_fixed_firstVar);
        UpdateUI_fixedX(ui->label_fixed_secondVar,ui->doubleSpinBox_fixed_secondVar);
        //independent variable
        UpdateUI_P(ui->label_delta_firstVar, ui->doubleSpinBox_delta_firstVar, ui->doubleSpinBox_max_firstVar, ui->doubleSpinBox_min_firstVar);

    }else if(arg=="Salinity")
    {
        //fixed vars
        UpdateUI_fixedP(ui->label_fixed_firsVar,ui->doubleSpinBox_fixed_firstVar);
        UpdateUI_fixedT(ui->label_fixed_secondVar,ui->doubleSpinBox_fixed_secondVar);
        //independent variable
        UpdateUI_X(ui->label_delta_firstVar, ui->doubleSpinBox_delta_firstVar, ui->doubleSpinBox_max_firstVar, ui->doubleSpinBox_min_firstVar);
    }else
    {
        std::cout<<"error: update1dUI, no such item: "<<arg.toStdString()<<std::endl;
    }
}
void MainWindow::UpdateUI_fixedP(QLabel* label, QDoubleSpinBox* box, double defaultValue)
{
    label->setText("Pressure (bar)");
    box->setDecimals(2);
    box->setRange(SWEOS::PMIN/1e5, SWEOS::PMAX/1e5);
    box->setSingleStep(1);
    box->setValue(defaultValue);
}
void MainWindow::UpdateUI_fixedX(QLabel* label, QDoubleSpinBox* box, double defaultValue)
{
    label->setText("Salinity");
    box->setDecimals(4);
    box->setRange(SWEOS::XMIN, SWEOS::XMAX);
    box->setSingleStep(0.001);
    box->setValue(defaultValue);
}
void MainWindow::UpdateUI_fixedT(QLabel* label, QDoubleSpinBox* box, double defaultValue)
{
    label->setText("Temperature (K)");
    box->setDecimals(2);
    box->setRange(SWEOS::TMIN, SWEOS::TMAX);
    box->setSingleStep(1);
    box->setValue(defaultValue);
}
void MainWindow::UpdateUI_X(QLabel* label, QDoubleSpinBox* deltaBox, QDoubleSpinBox* maxBox, QDoubleSpinBox* minBox)
{
    label->setText("dX:");
    deltaBox->setDecimals(4);
    deltaBox->setRange(SWEOS::XMIN, SWEOS::XMAX);
    deltaBox->setSingleStep(0.001);
    deltaBox->setValue(0.0001);

    maxBox->setDecimals(4);
    maxBox->setRange(SWEOS::XMIN, SWEOS::XMAX);
    maxBox->setSingleStep(0.001);
    maxBox->setValue(SWEOS::XMAX);

    minBox->setDecimals(4);
    minBox->setRange(0.0001,0.9);
    minBox->setSingleStep(0.001);
    minBox->setValue(0.01);
}
void MainWindow::UpdateUI_T(QLabel* label, QDoubleSpinBox* deltaBox, QDoubleSpinBox* maxBox, QDoubleSpinBox* minBox)
{
    label->setText("dT(K):");
    deltaBox->setDecimals(2);
    deltaBox->setRange(0.01, 100);
    deltaBox->setSingleStep(1);
    deltaBox->setValue(1);

    maxBox->setDecimals(2);
    maxBox->setRange(SWEOS::TMIN, SWEOS::TMAX);
    maxBox->setSingleStep(1);
    maxBox->setValue(SWEOS::TMAX);

    minBox->setDecimals(2);
    minBox->setRange(SWEOS::TMIN,SWEOS::TMAX);
    minBox->setSingleStep(1);
    minBox->setValue(SWEOS::TMIN);
}
void MainWindow::UpdateUI_P(QLabel* label, QDoubleSpinBox* deltaBox, QDoubleSpinBox* maxBox, QDoubleSpinBox* minBox)
{
    label->setText("dP (bar):");
    deltaBox->setDecimals(2);
    deltaBox->setRange(SWEOS::PMIN/1E5, SWEOS::PMAX/1E5);
    deltaBox->setSingleStep(1);
    deltaBox->setValue(SWEOS::PMIN/1E5);

    maxBox->setDecimals(2);
    maxBox->setRange(SWEOS::PMIN/1E5, SWEOS::PMAX/1E5);
    maxBox->setSingleStep(1);
    maxBox->setValue(SWEOS::PMAX/1E5);

    minBox->setDecimals(2);
    minBox->setRange(0.1,100);
    minBox->setSingleStep(1);
    minBox->setValue(1);
}
void MainWindow::update2dUI(QString arg)
{
    if(arg=="PT")
    {
        //fixed vars
        UpdateUI_fixedX(ui->label_fixed_secondVar,ui->doubleSpinBox_fixed_secondVar);
        //independent variable
        UpdateUI_P(ui->label_delta_firstVar, ui->doubleSpinBox_delta_firstVar, ui->doubleSpinBox_max_firstVar, ui->doubleSpinBox_min_firstVar);
        UpdateUI_T(ui->label_delta_secondVar, ui->doubleSpinBox_delta_secondVar, ui->doubleSpinBox_max_secondVar, ui->doubleSpinBox_min_secondVar);

    }else if(arg=="PX")
    {
        //fixed vars
        UpdateUI_fixedT(ui->label_fixed_secondVar,ui->doubleSpinBox_fixed_secondVar);
        //independent variable
        UpdateUI_P(ui->label_delta_firstVar, ui->doubleSpinBox_delta_firstVar, ui->doubleSpinBox_max_firstVar, ui->doubleSpinBox_min_firstVar);
        UpdateUI_X(ui->label_delta_secondVar, ui->doubleSpinBox_delta_secondVar, ui->doubleSpinBox_max_secondVar, ui->doubleSpinBox_min_secondVar);

    }else if(arg=="TX")
    {
        //fixed vars
        UpdateUI_fixedP(ui->label_fixed_secondVar,ui->doubleSpinBox_fixed_secondVar);
        //independent variable
        UpdateUI_T(ui->label_delta_firstVar, ui->doubleSpinBox_delta_firstVar, ui->doubleSpinBox_max_firstVar, ui->doubleSpinBox_min_firstVar);
        UpdateUI_X(ui->label_delta_secondVar, ui->doubleSpinBox_delta_secondVar, ui->doubleSpinBox_max_secondVar, ui->doubleSpinBox_min_secondVar);
    }else
    {
        std::cout<<"error: update2dUI, no such item: "<<arg.toStdString()<<std::endl;
    }
}
void MainWindow::update3dUI(QString arg)
{
    if(arg=="Temperature")
    {
        ui->label_fixed_firsVar->setText("Pressure (bar)");
        ui->doubleSpinBox_fixed_firstVar->setDecimals(2);
        ui->doubleSpinBox_fixed_firstVar->setRange(SWEOS::PMIN/1e5, SWEOS::PMAX/1e5);
        ui->doubleSpinBox_fixed_firstVar->setSingleStep(1);
        ui->doubleSpinBox_fixed_firstVar->setValue(316);
        ui->label_fixed_secondVar->setText("Salinity");
        ui->doubleSpinBox_fixed_secondVar->setDecimals(4);
        ui->doubleSpinBox_fixed_secondVar->setRange(SWEOS::XMIN, SWEOS::XMAX);
        ui->doubleSpinBox_fixed_secondVar->setSingleStep(0.001);
        ui->doubleSpinBox_fixed_secondVar->setValue(0.032);

        ui->label_delta_firstVar->setText("dT(K):");
        ui->doubleSpinBox_delta_firstVar->setDecimals(2);
        ui->doubleSpinBox_delta_firstVar->setRange(SWEOS::TMIN, SWEOS::TMAX);
        ui->doubleSpinBox_delta_firstVar->setSingleStep(1);
        ui->doubleSpinBox_delta_firstVar->setValue(SWEOS::TMIN);

        ui->doubleSpinBox_max_firstVar->setDecimals(2);
        ui->doubleSpinBox_max_firstVar->setRange(SWEOS::TMIN, SWEOS::TMAX);
        ui->doubleSpinBox_max_firstVar->setSingleStep(1);
        ui->doubleSpinBox_max_firstVar->setValue(SWEOS::TMAX);

        ui->doubleSpinBox_min_firstVar->setDecimals(2);
        ui->doubleSpinBox_min_firstVar->setRange(0.1,100);
        ui->doubleSpinBox_min_firstVar->setSingleStep(1);
        ui->doubleSpinBox_min_firstVar->setValue(SWEOS::TMIN);

    }else if(arg=="Pressure")
    {
        ui->label_fixed_firsVar->setText("Temperature (K)");
        ui->doubleSpinBox_fixed_firstVar->setDecimals(2);
        ui->doubleSpinBox_fixed_firstVar->setRange(SWEOS::TMIN, SWEOS::TMAX);
        ui->doubleSpinBox_fixed_firstVar->setSingleStep(1);
        ui->doubleSpinBox_fixed_firstVar->setValue(373);
        ui->label_fixed_secondVar->setText("Salinity");
        ui->doubleSpinBox_fixed_secondVar->setDecimals(4);
        ui->doubleSpinBox_fixed_secondVar->setRange(SWEOS::XMIN, SWEOS::XMAX);
        ui->doubleSpinBox_fixed_secondVar->setSingleStep(0.001);
        ui->doubleSpinBox_fixed_secondVar->setValue(0.032);

        ui->label_delta_firstVar->setText("dP (bar):");
        ui->doubleSpinBox_delta_firstVar->setDecimals(2);
        ui->doubleSpinBox_delta_firstVar->setRange(SWEOS::PMIN/1E5, SWEOS::PMAX/1E5);
        ui->doubleSpinBox_delta_firstVar->setSingleStep(1);
        ui->doubleSpinBox_delta_firstVar->setValue(SWEOS::PMIN/1E5);

        ui->doubleSpinBox_max_firstVar->setDecimals(2);
        ui->doubleSpinBox_max_firstVar->setRange(SWEOS::PMIN/1E5, SWEOS::PMAX/1E5);
        ui->doubleSpinBox_max_firstVar->setSingleStep(1);
        ui->doubleSpinBox_max_firstVar->setValue(SWEOS::PMAX/1E5);

        ui->doubleSpinBox_min_firstVar->setDecimals(2);
        ui->doubleSpinBox_min_firstVar->setRange(0.1,100);
        ui->doubleSpinBox_min_firstVar->setSingleStep(1);
        ui->doubleSpinBox_min_firstVar->setValue(1);

    }else if(arg=="Salinity")
    {
        ui->label_fixed_firsVar->setText("Pressure (bar)");
        ui->doubleSpinBox_fixed_firstVar->setDecimals(2);
        ui->doubleSpinBox_fixed_firstVar->setRange(SWEOS::PMIN/1e5, SWEOS::PMAX/1e5);
        ui->doubleSpinBox_fixed_firstVar->setSingleStep(1);
        ui->doubleSpinBox_fixed_firstVar->setValue(316);
        ui->label_fixed_secondVar->setText("Temperature (K)");
        ui->doubleSpinBox_fixed_secondVar->setDecimals(2);
        ui->doubleSpinBox_fixed_secondVar->setRange(SWEOS::TMIN, SWEOS::TMAX);
        ui->doubleSpinBox_fixed_secondVar->setSingleStep(1);
        ui->doubleSpinBox_fixed_secondVar->setValue(373);

        ui->label_delta_firstVar->setText("dX:");
        ui->doubleSpinBox_delta_firstVar->setDecimals(4);
        ui->doubleSpinBox_delta_firstVar->setRange(SWEOS::XMIN, SWEOS::XMAX);
        ui->doubleSpinBox_delta_firstVar->setSingleStep(0.001);
        ui->doubleSpinBox_delta_firstVar->setValue(0.0001);

        ui->doubleSpinBox_max_firstVar->setDecimals(4);
        ui->doubleSpinBox_max_firstVar->setRange(SWEOS::XMIN, SWEOS::XMAX);
        ui->doubleSpinBox_max_firstVar->setSingleStep(0.001);
        ui->doubleSpinBox_max_firstVar->setValue(SWEOS::XMAX);

        ui->doubleSpinBox_min_firstVar->setDecimals(4);
        ui->doubleSpinBox_min_firstVar->setRange(0.0001,0.9);
        ui->doubleSpinBox_min_firstVar->setSingleStep(0.001);
        ui->doubleSpinBox_min_firstVar->setValue(0.01);
    }else
    {
        std::cout<<"error: update1dUI, no such item: "<<arg.toStdString()<<std::endl;
    }
}
void MainWindow::on_comboBox_selectVariable_activated(const QString &arg1)
{
    switch (m_dimension) {
    case 1:
        update1dUI(arg1);
    break;
    case 2:
        update2dUI(arg1);
    break;
    case 3:
        update3dUI(arg1);
    break;
    }
}
