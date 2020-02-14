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
    ,m_vtkFontSize(25)
    ,m_index_var(0)
    ,m_xlabel("xlabel")
    ,m_ylabel("ylabel")
    ,m_zlabel("zlabel")
    ,m_showScatter_1Dchart(false)
    ,m_alphaPhaseRegion(0.4)
    ,m_vtkLineWidth(5)
    ,m_showPhaseRegion_1Dchart(true)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    m_vtkColorSeries_PhaseRegion=vtkSmartPointer<vtkColorSeries>::New();
    m_vtkColorSeries_PhaseRegion->SetColorScheme(vtkColorSeries::BREWER_QUALITATIVE_SET3);

    m_vtkColorSeries_Lines_1Dchart=vtkSmartPointer<vtkColorSeries>::New();
    m_vtkColorSeries_Lines_1Dchart->SetColorScheme(0); //SPECTRUM: 7 colors

    //progressbar
    m_progressBar = new QProgressBar();
    m_progressBar->setRange(0, 100);
    m_progressBar->setValue(0);
    m_progressBar->setTextVisible(true);
    m_progressBar->setFormat("");
    ui->statusbar->addPermanentWidget(m_progressBar, 2);
    //  // Qt Table View
    ui->tabWidget->setCurrentIndex(0);
    ui->vtkWindowTab->setCurrentIndex(0);

    m_IndependentVar1_old=ui->doubleSpinBox->value()*1e5;
    m_IndependentVar2_old=ui->doubleSpinBox_2->value();
    m_IndependentVar3_old=ui->doubleSpinBox_3->value();

    // used to set actor to a nice width/height ratio
    m_actorScale_T=(SWEOS::PMAX-SWEOS::PMIN)/1e5/(SWEOS::TMAX-SWEOS::TMIN);
    m_actorScale_P=1;
    m_actorScale_X=(SWEOS::PMAX-SWEOS::PMIN)/1e5/(SWEOS::XMAX-SWEOS::XMIN);
    m_actorScale[0]=m_actorScale_T;
    m_actorScale[1]=m_actorScale_P;
    m_actorScale[2]=m_actorScale_X;

    //default set 1D UI
    m_geometry_Groupbox_variables=ui->groupBox_Variables->geometry();
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
//  on_pushButton_clicked();
};

MainWindow::~MainWindow()
{
  // The smart pointers should clean up for up
  if(m_progressBar) delete m_progressBar;
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
    switch (m_calculationMode) {
        case CALCULATION_SINGLE_POINT:
        {
            SinglePointCalculation(ui->comboBox->currentIndex());
        }
        break;
        case CALCULATION_MULTI_POINTS:
        {
            QString fileName;
            fileName = QFileDialog::getOpenFileName(this, tr("Open T(C) P(bar) X File: three columns separated by spaces"), "", tr("Text File (*.txt)"));

            if (!fileName.isNull())
            {
                ifstream fin(fileName.toStdString());
                if(!fin)
                {
                    std::cout<<"error: open file failed, "<<fileName.toStdString()<<endl;
                    return;
                }
                vector<double>arrT,arrP,arrX;
                double T,P,X;
                while (!fin.eof()) {
                    fin>>T>>P>>X;
                    arrT.push_back(T);
                    arrP.push_back(P*1e5);
                    arrX.push_back(X);
                }
                fin.close();

                // calculate
                CalculateProps_PTX(arrT,arrP,arrX,m_vtkTable);

                // display in textedit
                ui->textEdit->clear();
                // head
                QString str_line="";
                for (int j=0;j<m_vtkTable->GetNumberOfColumns();j++) {
                    QString str;
                    str.sprintf("%30s, ",m_vtkTable->GetColumn(j)->GetName());
                    str_line+=str;
                }
                ui->textEdit->append(str_line);
                for (int i=0;i<m_vtkTable->GetNumberOfRows();i++) {
                    QString str_line="";
                    for (int j=0;j<m_vtkTable->GetNumberOfColumns();j++) {
                        QString str;
                        str.sprintf("%30.2f, ",m_vtkTable->GetValue(i,j).ToDouble());
//                        str_line+=QString::number(m_vtkTable->GetValue(i,j).ToDouble(), 'f', 2)+", ";
                        str_line+=str;
                    }
                    ui->textEdit->append(str_line);
                }
            }
            else
            {

            }
        }
        break;
    }
}
void MainWindow::SinglePointCalculation(int index_varsSelection)
{
    double p=ui->doubleSpinBox->value()*1e5;
    double X=ui->doubleSpinBox_3->value();
    double T_H=ui->doubleSpinBox_2->value();
    QString color_value_1="Black", color_value_2="Black", color_value_3="Black";
    color_value_1=(m_IndependentVar1_old==p ? color_value_1="Black": color_value_1="Green");
    color_value_2=(m_IndependentVar2_old==T_H ? color_value_2="Black": color_value_2="Green");
    color_value_3=(m_IndependentVar3_old==X ? color_value_3="Black": color_value_3="Green");
    QString name_T_H, name_unit_T_H, name_T_H_display;
    double value_T_H_display;
    SWEOS::cH2ONaCl eos;
    switch (index_varsSelection) {
    case 0: //PTX
    {
        eos.prop_pTX(p, T_H+SWEOS::Kelvin, X);
        name_T_H="Temperature";
        name_unit_T_H="deg. C";
        name_T_H_display="Bulk enthalpy";
        value_T_H_display=eos.m_prop.H;
    }
        break;
    case 1: //PHX
    {
        eos.prop_pHX(p, T_H*1000, X); //Enthalpy unit in UI is kJ/kg
        name_T_H="Enthalpy";
        name_unit_T_H="kJ/kg";
        name_T_H_display="Temperature";
        value_T_H_display=eos.m_prop.T;
    }
        break;
    }
   QString result_str;
   // T, P, X
   result_str="<font color=Purple>Pressure</font> = <font color="+color_value_1+">"+QString::number(p)
           +"</font> Pa, <font color=Purple>"+name_T_H+"</font> = <font color="+color_value_2+">"+QString::number(T_H)
           +"</font> "+name_unit_T_H+", <font color=Purple>Salinity</font>  = <font color="+color_value_3+">"+QString::number(X*100)
           +"</font> wt. % NaCl<br>";
   result_str+="======================================================================<br>";
   // Region
   result_str+="<font color=Blue>Phase Region</font>: "+QString::fromStdString(eos.m_phaseRegion_name[eos.m_prop.Region])+"<br>";
   // Xl, Xv
   result_str+="<font color=Blue>Xl: </font> &nbsp; "+QString::number(eos.m_prop.X_l)+
             ",&nbsp;&nbsp;&nbsp; <font color=Blue>Xv: </font> &nbsp; "+QString::number(eos.m_prop.X_v)+"<br>";
   // Bulk rho, bulk H, bulk mu
   result_str+="<font color=Blue>Bulk density: </font> &nbsp; "+QString::number(eos.m_prop.Rho)+
             ",&nbsp;&nbsp;&nbsp; <font color=Blue>"+name_T_H_display+": </font> &nbsp; "+QString::number(value_T_H_display)+"<br>";
//             ",&nbsp;&nbsp;&nbsp; <font color=Blue>Bulk viscosity: </font> &nbsp; "+QString::number(eos.m_prop.Mu)
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
   m_IndependentVar2_old=T_H;
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
        ui->vtkWindowTab->setEnabled(false);
        break;
    case 1:
        ui->textEdit->setEnabled(false);
        ui->vtkWindowTab->setEnabled(true);
        break;

    }

}

void MainWindow::Calculate_Diagram1D()
{
    QString varName;
    //            int index_var=0;
    vector<double>arrP,arrT,arrX;
    switch (ui->comboBox_selectVariable->currentIndex()) {
    case 0:   //temperature
    {
        varName="Temperature (C)";
        m_index_var=0;
        double dT=ui->doubleSpinBox_delta_firstVar->value();
        double Tmin=ui->doubleSpinBox_min_firstVar->value();
        double Tmax=ui->doubleSpinBox_max_firstVar->value();
        double p0=ui->doubleSpinBox_fixed_firstVar->value()*1e5;
        double X0=ui->doubleSpinBox_fixed_secondVar->value();
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
        m_index_var=1;
        double dP=ui->doubleSpinBox_delta_firstVar->value()*1e5;
        double Pmin=ui->doubleSpinBox_min_firstVar->value()*1e5;
        double Pmax=ui->doubleSpinBox_max_firstVar->value()*1e5;
        double T0=ui->doubleSpinBox_fixed_firstVar->value();
        double X0=ui->doubleSpinBox_fixed_secondVar->value();
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
        m_index_var=2;
        double dX=ui->doubleSpinBox_delta_firstVar->value();
        double Xmin=ui->doubleSpinBox_min_firstVar->value();
        double Xmax=ui->doubleSpinBox_max_firstVar->value();
        double P0=ui->doubleSpinBox_fixed_firstVar->value()*1e5;
        double T0=ui->doubleSpinBox_fixed_secondVar->value();
        for (double X=Xmin;X<Xmax;X=X+dX) {
            arrT.push_back(T0);
            arrP.push_back(P0);
            arrX.push_back(X);
        }
    }
        break;
    }

    //calculate
    CalculateProps_PTX(arrT, arrP,arrX, m_vtkTable);
    //display
    ShowProps_1D();

}

void MainWindow::CalculateProps_PTX(std::vector<double> arrT,std::vector<double> arrP, std::vector<double> arrX, vtkSmartPointer<vtkTable> table )
{
    // add columns to table
    //ui->textEdit->append(QString::number(m_vtkTable->GetNumberOfColumns()));
//    vtkSmartPointer<vtkTable> table=vtkSmartPointer<vtkTable>::New();

    vtkSmartPointer<vtkFloatArray> varT = vtkSmartPointer<vtkFloatArray>::New();
    varT->SetName("Temperature(C)");
    table->AddColumn(varT);

    vtkSmartPointer<vtkFloatArray> varP = vtkSmartPointer<vtkFloatArray>::New();
    varP->SetName("Pressure (bar)");
    table->AddColumn(varP);

    vtkSmartPointer<vtkFloatArray> varX = vtkSmartPointer<vtkFloatArray>::New();
    varX->SetName("Salinity");
    table->AddColumn(varX);

    //phase region
    vtkSmartPointer<vtkIntArray> prop_region = vtkSmartPointer<vtkIntArray>::New();
    prop_region->SetName("Phase Region");
    table->AddColumn(prop_region);
    //density
    vtkSmartPointer<vtkFloatArray> prop_density = vtkSmartPointer<vtkFloatArray>::New();
    prop_density->SetName("Bulk density (kg/m3)");
    table->AddColumn(prop_density);
    vtkSmartPointer<vtkFloatArray> prop_rho_l = vtkSmartPointer<vtkFloatArray>::New();
    prop_rho_l->SetName("Liquid density (kg/m3)");
    table->AddColumn(prop_rho_l);
    vtkSmartPointer<vtkFloatArray> prop_rho_v = vtkSmartPointer<vtkFloatArray>::New();
    prop_rho_v->SetName("Vapour density (kg/m3)");
    table->AddColumn(prop_rho_v);
    vtkSmartPointer<vtkFloatArray> prop_rho_h = vtkSmartPointer<vtkFloatArray>::New();
    prop_rho_h->SetName("Halite density (kg/m3)");
    table->AddColumn(prop_rho_h);
    //Enthalpy
    vtkSmartPointer<vtkFloatArray> prop_enthalpy = vtkSmartPointer<vtkFloatArray>::New();
    prop_enthalpy->SetName("Bulk enthalpy (J/kg)");
    table->AddColumn(prop_enthalpy);
    vtkSmartPointer<vtkFloatArray> prop_enthalpy_l = vtkSmartPointer<vtkFloatArray>::New();
    prop_enthalpy_l->SetName("Liquid enthalpy (J/kg)");
    table->AddColumn(prop_enthalpy_l);
    vtkSmartPointer<vtkFloatArray> prop_enthalpy_v = vtkSmartPointer<vtkFloatArray>::New();
    prop_enthalpy_v->SetName("Vapour enthalpy (J/kg)");
    table->AddColumn(prop_enthalpy_v);
    vtkSmartPointer<vtkFloatArray> prop_enthalpy_h = vtkSmartPointer<vtkFloatArray>::New();
    prop_enthalpy_h->SetName("Halite enthalpy (J/kg)");
    table->AddColumn(prop_enthalpy_h);

    //Saturation
    vtkSmartPointer<vtkFloatArray> prop_S_l = vtkSmartPointer<vtkFloatArray>::New();
    prop_S_l->SetName("Liquid saturation");
    table->AddColumn(prop_S_l);
    vtkSmartPointer<vtkFloatArray> prop_S_v = vtkSmartPointer<vtkFloatArray>::New();
    prop_S_v->SetName("Vapour saturation");
    table->AddColumn(prop_S_v);
    vtkSmartPointer<vtkFloatArray> prop_S_h = vtkSmartPointer<vtkFloatArray>::New();
    prop_S_h->SetName("Halite saturation");
    table->AddColumn(prop_S_h);

    //viscosity
    vtkSmartPointer<vtkFloatArray> prop_Mu_l = vtkSmartPointer<vtkFloatArray>::New();
    prop_Mu_l->SetName("Liquid viscosity (Pa s)");
    table->AddColumn(prop_Mu_l);
    vtkSmartPointer<vtkFloatArray> prop_Mu_v = vtkSmartPointer<vtkFloatArray>::New();
    prop_Mu_v->SetName("Vapour viscosity (Pa s)");
    table->AddColumn(prop_Mu_v);

    //salinity
    vtkSmartPointer<vtkFloatArray> prop_X_l = vtkSmartPointer<vtkFloatArray>::New();
    prop_X_l->SetName("Liquid salinity");
    table->AddColumn(prop_X_l);
    vtkSmartPointer<vtkFloatArray> prop_X_v = vtkSmartPointer<vtkFloatArray>::New();
    prop_X_v->SetName("Vapour salinity");
    table->AddColumn(prop_X_v);

    //calculate and set value to vtktable
    table->SetNumberOfRows(arrT.size());
    for (size_t i=0;i<arrT.size();++i) {
        SWEOS::cH2ONaCl eos;
        eos.prop_pTX(arrP[i],arrT[i]+SWEOS::Kelvin,arrX[i]);
        table->SetValue(i, 0, arrT[i]);
        table->SetValue(i, 1, arrP[i]/1e5);
        table->SetValue(i, 2, arrX[i]);
        table->SetValue(i, 3, eos.m_prop.Region);

        table->SetValue(i, 4, eos.m_prop.Rho);
        table->SetValue(i, 5, eos.m_prop.Rho_l);
        table->SetValue(i, 6, eos.m_prop.Rho_v);
        table->SetValue(i, 7, eos.m_prop.Rho_h);

        table->SetValue(i, 8, eos.m_prop.H);
        table->SetValue(i, 9, eos.m_prop.H_l);
        table->SetValue(i, 10, eos.m_prop.H_v);
        table->SetValue(i, 11, eos.m_prop.H_h);

        table->SetValue(i, 12, eos.m_prop.S_l);
        table->SetValue(i, 13, eos.m_prop.S_v);
        table->SetValue(i, 14, eos.m_prop.S_h);

        table->SetValue(i, 15, eos.m_prop.Mu_l);
        table->SetValue(i, 16, eos.m_prop.Mu_v);

        table->SetValue(i, 17, eos.m_prop.X_l);
        table->SetValue(i, 18, eos.m_prop.X_v);
    }
}

void MainWindow::ShowProps_1D()
{
    int index_var=m_index_var;
    int index_prop_combox=ui->comboBox_selectProps->currentIndex();
    if(m_vtkTable->GetNumberOfRows()==0)return;
    switch (index_prop_combox) {
        case 0: //phase region
        {
            std::vector<int> index_props={3};
            std::vector<bool> showcomponents={true};
            update1dChart(index_var, "Phase region", index_props, showcomponents, m_vtkColorSeries_Lines_1Dchart);
        }
            break;
        case 1: //Density, rho_l, rho_v, rho_h
        {
            std::vector<int> index_props={4,5,6,7};
            std::vector<bool> showcomponents={ui->checkBox_2->isChecked(), ui->checkBox_3->isChecked(), ui->checkBox_4->isChecked(), ui->checkBox_5->isChecked()};
            update1dChart(index_var, "Density (kg/m3)", index_props,showcomponents, m_vtkColorSeries_Lines_1Dchart);
        }
            break;
        case 2: //Enthalpy
        {
            std::vector<int> index_props={8, 9, 10, 11};
            std::vector<bool> showcomponents={ui->checkBox_2->isChecked(), ui->checkBox_3->isChecked(), ui->checkBox_4->isChecked(), ui->checkBox_5->isChecked()};
            update1dChart(index_var, "Specific enthalpy (J/kg)", index_props,showcomponents, m_vtkColorSeries_Lines_1Dchart);
        }
            break;
        case 3://saturation
        {
            std::vector<int> index_props={12,13,14};
            std::vector<bool> showcomponents={ui->checkBox_2->isChecked(), ui->checkBox_3->isChecked(), ui->checkBox_4->isChecked(), ui->checkBox_5->isChecked()};
            update1dChart(index_var, "Saturation", index_props,showcomponents, m_vtkColorSeries_Lines_1Dchart);
        }
            break;
        case 4://viscosity
        {
            std::vector<int> index_props={15, 16};
            std::vector<bool> showcomponents={ui->checkBox_2->isChecked(), ui->checkBox_3->isChecked(), ui->checkBox_4->isChecked(), ui->checkBox_5->isChecked()};
            update1dChart(index_var, "Dynamic viscosity (Pa s)", index_props,showcomponents, m_vtkColorSeries_Lines_1Dchart);
        }
            break;
        case 5://salinity
        {
            std::vector<int> index_props={17, 18};
            std::vector<bool> showcomponents={ui->checkBox_2->isChecked(), ui->checkBox_3->isChecked(), ui->checkBox_4->isChecked(), ui->checkBox_5->isChecked()};
            update1dChart(index_var, "Salinity", index_props,showcomponents, m_vtkColorSeries_Lines_1Dchart);
        }
            break;
    }
}

void MainWindow::GetMaxMin_vtkTableColumn(const vtkSmartPointer<vtkTable> table, int index_col, double& min, double& max)
{
    min=1e30;
    max=-1e30;
    double value_tab;
    for (int i=0;i<table->GetNumberOfRows();i++) {
        value_tab=table->GetValue(i,index_col).ToDouble();
        min=(value_tab<min ? value_tab: min);
        max=(value_tab>max ? value_tab: max);
    }
}

void MainWindow::update1dChart(int index_var, std::string name_prop, std::vector<int> index_props, std::vector<bool> components, vtkSmartPointer<vtkColorSeries> colorseries)
{
    //clean old charts
    int num_oldItems=m_vtkChartView->GetScene()->GetNumberOfItems();
    for (int i=0;i<num_oldItems;i++) {
        m_vtkChartView->GetScene()->RemoveItem(m_vtkChartView->GetScene()->GetItem(0));
    }
    // add new chart
    vtkSmartPointer<vtkChartXY> chart = vtkSmartPointer<vtkChartXY>::New();
    m_vtkChartView->GetScene()->AddItem(chart);


    vtkNew<vtkNamedColors> colors;
    vtkColor3d color3d;
    double min_allProps=1e20, max_allProps=-1e20;

    for (size_t i=0;i<index_props.size();i++) {
        if(!components[i]) continue;
        double min,max;
        GetMaxMin_vtkTableColumn(m_vtkTable, index_props[i], min,max);
        min_allProps=(min< min_allProps ? min: min_allProps);
        max_allProps=(max>max_allProps ? max: max_allProps);
    }
    //extract phase region and plot
    if(m_showPhaseRegion_1Dchart)
    {
        const int index_col_table_phaseRegion=3;
        vector<int> indexInTable_phaseRegion;
        int phaseRegion_start=m_vtkTable->GetValue(0,index_col_table_phaseRegion).ToInt();
        indexInTable_phaseRegion.push_back(0);
        for (int i=1;i<m_vtkTable->GetNumberOfRows();i++) {
            if(phaseRegion_start!=m_vtkTable->GetValue(i,index_col_table_phaseRegion).ToInt())
            {
                indexInTable_phaseRegion.push_back(i);
                phaseRegion_start=m_vtkTable->GetValue(i,index_col_table_phaseRegion).ToInt();
            }
        }
        if(indexInTable_phaseRegion[indexInTable_phaseRegion.size()-1]!=(m_vtkTable->GetNumberOfRows()-1))
            indexInTable_phaseRegion.push_back(m_vtkTable->GetNumberOfRows()-1);
        // plot region area
        min_allProps=min_allProps-(max_allProps-min_allProps)*0.02;
        max_allProps=max_allProps+(max_allProps-min_allProps)*0.02;
        SWEOS::cH2ONaCl eos;
        bool PhaseRegion_present[8]={false, false, false, false,false, false, false, false};
        for (size_t i=1;i<indexInTable_phaseRegion.size();i++)
        {
            int phaseRegion=m_vtkTable->GetValue(indexInTable_phaseRegion[i-1],index_col_table_phaseRegion).ToInt();
            vtkNew<vtkTable> table_region;
              vtkNew<vtkFloatArray> arrX;
              arrX->SetName("X Axis");
              table_region->AddColumn(arrX);
              vtkNew<vtkFloatArray> arrBottom;
              arrBottom->SetName(eos.m_phaseRegion_name[phaseRegion].c_str());
              table_region->AddColumn(arrBottom);
              vtkNew<vtkFloatArray> arrTop;
              arrTop->SetName("Top");
              table_region->AddColumn(arrTop);
              int numPoints = indexInTable_phaseRegion[i]-indexInTable_phaseRegion[i-1] +1;
              table_region->SetNumberOfRows(numPoints);
              for (int j = indexInTable_phaseRegion[i-1]; j <= indexInTable_phaseRegion[i]; j++)
              {
                int ind=j-indexInTable_phaseRegion[i-1];
                table_region->SetValue(ind, 0, m_vtkTable->GetValue(j,index_var).ToDouble());
                table_region->SetValue(ind, 1, min_allProps);
                table_region->SetValue(ind, 2, max_allProps);
              }
              // Add multiple line plots, setting the colors etc
              vtkPlotArea* area = dynamic_cast<vtkPlotArea*>(chart->AddPlot(vtkChart::AREA));
              area->SetInputData(table_region);
              area->SetInputArray(0, "X Axis");
              area->SetInputArray(1, arrBottom->GetName());
              area->SetInputArray(2, "Top");
              area->GetBrush()->SetColorF(m_vtkColorSeries_PhaseRegion->GetColor(phaseRegion).GetRed()/255.0,
                                          m_vtkColorSeries_PhaseRegion->GetColor(phaseRegion).GetGreen()/255.0,
                                          m_vtkColorSeries_PhaseRegion->GetColor(phaseRegion).GetBlue()/255.0,
                                          m_alphaPhaseRegion);
              if(PhaseRegion_present[phaseRegion])area->LegendVisibilityOff();
              if(!PhaseRegion_present[phaseRegion]) PhaseRegion_present[phaseRegion]=true; //used to mark duplicated phase region legend
        }
    }

    // plot properties curve
    for (size_t i=0;i<index_props.size();i++) {
        if(!components[i]) continue;
        vtkPlot *line = chart->AddPlot(vtkChart::LINE);
        line->SetInputData(m_vtkTable, index_var, index_props[i]);
        line->SetColor(color3d.GetRed(), color3d.GetGreen(), color3d.GetBlue());
        line->SetColor(colorseries->GetColor(i).GetRed()/255.0,
                       colorseries->GetColor(i).GetGreen()/255.0,
                       colorseries->GetColor(i).GetBlue()/255.0);
        line->SetWidth(m_vtkLineWidth);
        // scatter
        if(m_showScatter_1Dchart)
        {
            vtkPlot *points = chart->AddPlot(vtkChart::POINTS);
            points->SetInputData(m_vtkTable, index_var, index_props[i]);
            points->SetColor(0, 0, 0);
            points->SetWidth(m_vtkLineWidth*2);
            points->LegendVisibilityOff();
            dynamic_cast<vtkPlotPoints*>(points)->SetMarkerStyle(vtkPlotPoints::PLUS);
        }
    }


    chart->GetAxis(1)->SetTitle(m_vtkTable->GetColumn(index_var)->GetName());
    chart->GetAxis(1)->GetLabelProperties()->SetFontSize(m_vtkFontSize);
    chart->GetAxis(1)->GetLabelProperties()->SetFontFamilyToArial();
    chart->GetAxis(1)->GetTitleProperties()->SetFontSize(m_vtkFontSize);
    chart->GetAxis(1)->GetTitleProperties()->SetBold(false);
    chart->GetAxis(1)->GetTitleProperties()->SetFontFamilyToArial();

    chart->GetAxis(0)->SetTitle(name_prop);
    chart->GetAxis(0)->GetLabelProperties()->SetFontSize(m_vtkFontSize);
    chart->GetAxis(0)->GetLabelProperties()->SetFontFamilyToArial();
    chart->GetAxis(0)->GetTitleProperties()->SetFontSize(m_vtkFontSize);
    chart->GetAxis(0)->GetTitleProperties()->SetBold(false);
    chart->GetAxis(0)->GetTitleProperties()->SetFontFamilyToArial();

    chart->GetTooltip()->GetTextProperties()->SetFontSize(m_vtkFontSize);

    chart->SetShowLegend(true);
    chart->GetLegend()->GetLabelProperties()->SetFontSize(m_vtkFontSize);
    chart->GetLegend()->SetSymbolWidth(60);

    m_vtkChartView->GetRenderWindow()->Render();
}

void MainWindow::Calculate_Diagram2D()
{
//    vector<double>arrP,arrT,arrX;
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    //phase region
    vtkSmartPointer<vtkDoubleArray> arrPhaseRegion = vtkSmartPointer<vtkDoubleArray>::New();
    arrPhaseRegion->SetNumberOfComponents(1);
    arrPhaseRegion->SetName("PhaseRegion");
    //density
    vtkSmartPointer<vtkDoubleArray> arrDensity = vtkSmartPointer<vtkDoubleArray>::New();
    arrDensity->SetNumberOfComponents(1);
    arrDensity->SetName("Bulk density (kg/m3)");
    vtkSmartPointer<vtkDoubleArray> arrDensity_l = vtkSmartPointer<vtkDoubleArray>::New();
    arrDensity_l->SetNumberOfComponents(1);
    arrDensity_l->SetName("Liquid density (kg/m3)");
    vtkSmartPointer<vtkDoubleArray> arrDensity_v = vtkSmartPointer<vtkDoubleArray>::New();
    arrDensity_v->SetNumberOfComponents(1);
    arrDensity_v->SetName("Vapour density(kg/m3)");
    vtkSmartPointer<vtkDoubleArray> arrDensity_h = vtkSmartPointer<vtkDoubleArray>::New();
    arrDensity_h->SetNumberOfComponents(1);
    arrDensity_h->SetName("Halite density (kg/m3)");
    //enthalpy
    vtkSmartPointer<vtkDoubleArray> arrH = vtkSmartPointer<vtkDoubleArray>::New();
    arrH->SetNumberOfComponents(1);
    arrH->SetName("Bulk enthalpy (J/kg)");
    vtkSmartPointer<vtkDoubleArray> arrH_l = vtkSmartPointer<vtkDoubleArray>::New();
    arrH_l->SetNumberOfComponents(1);
    arrH_l->SetName("Liquid enthalpy (J/kg)");
    vtkSmartPointer<vtkDoubleArray> arrH_v = vtkSmartPointer<vtkDoubleArray>::New();
    arrH_v->SetNumberOfComponents(1);
    arrH_v->SetName("Vapour enthalpy (J/kg)");
    vtkSmartPointer<vtkDoubleArray> arrH_h = vtkSmartPointer<vtkDoubleArray>::New();
    arrH_h->SetNumberOfComponents(1);
    arrH_h->SetName("Halite enthalpy (J/kg)");
    //saturation
    vtkSmartPointer<vtkDoubleArray> arrS_l = vtkSmartPointer<vtkDoubleArray>::New();
    arrS_l->SetNumberOfComponents(1);
    arrS_l->SetName("Liquid saturation");
    vtkSmartPointer<vtkDoubleArray> arrS_v = vtkSmartPointer<vtkDoubleArray>::New();
    arrS_v->SetNumberOfComponents(1);
    arrS_v->SetName("Vapour saturation");
    vtkSmartPointer<vtkDoubleArray> arrS_h = vtkSmartPointer<vtkDoubleArray>::New();
    arrS_h->SetNumberOfComponents(1);
    arrS_h->SetName("Halite saturation");
    //viscosity
    vtkSmartPointer<vtkDoubleArray> arrMu_l = vtkSmartPointer<vtkDoubleArray>::New();
    arrMu_l->SetNumberOfComponents(1);
    arrMu_l->SetName("Liquid viscosity (Pa s)");
    vtkSmartPointer<vtkDoubleArray> arrMu_v = vtkSmartPointer<vtkDoubleArray>::New();
    arrMu_v->SetNumberOfComponents(1);
    arrMu_v->SetName("Vapour viscosity (Pa s)");
    //salinity
    vtkSmartPointer<vtkDoubleArray> arrX_l = vtkSmartPointer<vtkDoubleArray>::New();
    arrX_l->SetNumberOfComponents(1);
    arrX_l->SetName("Liquid salinity");
    vtkSmartPointer<vtkDoubleArray> arrX_v = vtkSmartPointer<vtkDoubleArray>::New();
    arrX_v->SetNumberOfComponents(1);
    arrX_v->SetName("Vapour salinity");

    switch (ui->comboBox_selectVariable->currentIndex()) {
        case 0:   //PT
        {
            m_xlabel="Temperature (C)";
            m_ylabel="Pressure (bar)";
            m_actorScale[0]=m_actorScale_T;
            m_actorScale[1]=m_actorScale_P;
            m_actorScale[2]=m_actorScale_X;
            double dP=ui->doubleSpinBox_delta_firstVar->value();
            double Pmin=ui->doubleSpinBox_min_firstVar->value();
            double Pmax=ui->doubleSpinBox_max_firstVar->value();

            double dT=ui->doubleSpinBox_delta_secondVar->value();
            double Tmin=ui->doubleSpinBox_min_secondVar->value();
            double Tmax=ui->doubleSpinBox_max_secondVar->value();

            double X0=ui->doubleSpinBox_fixed_firstVar->value();

            vector<double> vectorT, vectorP;
            for (double T=Tmin; T<Tmax;T=T+dT)
            {
              vectorT.push_back(T);
            }
            for (double P=Pmin; P<Pmax; P=P+dP)
            {
              vectorP.push_back(P);
            }
            //calculate
            for(size_t j = 0; j < vectorP.size(); j++)
            {
                for(size_t i = 0; i < vectorT.size(); i++)
                {
                    SWEOS::cH2ONaCl eos;
                    eos.prop_pTX(vectorP[j]*1e5,vectorT[i]+SWEOS::Kelvin,X0);
                    points->InsertNextPoint(vectorT[i],vectorP[j],X0);
                    arrPhaseRegion->InsertNextValue(eos.m_prop.Region);
                    arrDensity->InsertNextValue(eos.m_prop.Rho);
                    arrDensity_l->InsertNextValue(eos.m_prop.Rho_l);
                    arrDensity_v->InsertNextValue(eos.m_prop.Rho_v);
                    arrDensity_h->InsertNextValue(eos.m_prop.Rho_h);
                    arrH->InsertNextValue(eos.m_prop.H);
                    arrH_l->InsertNextValue(eos.m_prop.H_l);
                    arrH_v->InsertNextValue(eos.m_prop.H_v);
                    arrH_h->InsertNextValue(eos.m_prop.H_h);
                    arrS_l->InsertNextValue(eos.m_prop.S_l);
                    arrS_v->InsertNextValue(eos.m_prop.S_v < 1e-30 ? 0 : eos.m_prop.S_v);
                    arrS_h->InsertNextValue(eos.m_prop.S_h < 1e-30 ? 0: eos.m_prop.S_h);
                    arrMu_l->InsertNextValue(eos.m_prop.Mu_l);
                    arrMu_v->InsertNextValue(eos.m_prop.Mu_v);
                    arrX_l->InsertNextValue(eos.m_prop.X_l);
                    arrX_v->InsertNextValue(eos.m_prop.X_v);
//                    QApplication::processEvents();
                }
                m_progressBar->setValue(j/(vectorP.size()-1)*100);
                m_progressBar->setFormat("Calculating");
            }
            // Specify the dimensions of the grid
            m_structuredGrid->SetDimensions(vectorT.size(),vectorP.size(),1);
        }
        break;
        case 1:   //PX
        {
            m_xlabel="Salinity";
            m_ylabel="Pressure (bar)";
            m_actorScale[0]=m_actorScale_X;
            m_actorScale[1]=m_actorScale_P;
            m_actorScale[2]=m_actorScale_T;
            double dP=ui->doubleSpinBox_delta_firstVar->value();
            double Pmin=ui->doubleSpinBox_min_firstVar->value();
            double Pmax=ui->doubleSpinBox_max_firstVar->value();

            double dX=ui->doubleSpinBox_delta_secondVar->value();
            double Xmin=ui->doubleSpinBox_min_secondVar->value();
            double Xmax=ui->doubleSpinBox_max_secondVar->value();

            double T0=ui->doubleSpinBox_fixed_firstVar->value();

            vector<double> vectorX, vectorP;
            for (double X=Xmin; X<Xmax;X=X+dX)
            {
              vectorX.push_back(X);
            }
            for (double P=Pmin; P<Pmax; P=P+dP)
            {
              vectorP.push_back(P);
            }
            //calculate
            for(size_t j = 0; j < vectorP.size(); j++)
            {
                for(size_t i = 0; i < vectorX.size(); i++)
                {
                    SWEOS::cH2ONaCl eos;
                    eos.prop_pTX(vectorP[j]*1e5,T0+SWEOS::Kelvin,vectorX[i]);
                    points->InsertNextPoint(vectorX[i],vectorP[j],T0);
                    arrPhaseRegion->InsertNextValue(eos.m_prop.Region);
                    arrDensity->InsertNextValue(eos.m_prop.Rho);
                    arrDensity_l->InsertNextValue(eos.m_prop.Rho_l);
                    arrDensity_v->InsertNextValue(eos.m_prop.Rho_v);
                    arrDensity_h->InsertNextValue(eos.m_prop.Rho_h);
                    arrH->InsertNextValue(eos.m_prop.H);
                    arrH_l->InsertNextValue(eos.m_prop.H_l);
                    arrH_v->InsertNextValue(eos.m_prop.H_v);
                    arrH_h->InsertNextValue(eos.m_prop.H_h);
                    arrS_l->InsertNextValue(eos.m_prop.S_l);
                    arrS_v->InsertNextValue(eos.m_prop.S_v < 1e-30 ? 0 : eos.m_prop.S_v);
                    arrS_h->InsertNextValue(eos.m_prop.S_h < 1e-30 ? 0: eos.m_prop.S_h);
                    arrMu_l->InsertNextValue(eos.m_prop.Mu_l);
                    arrMu_v->InsertNextValue(eos.m_prop.Mu_v);
                    arrX_l->InsertNextValue(eos.m_prop.X_l);
                    arrX_v->InsertNextValue(eos.m_prop.X_v);
                }
            }
            // Specify the dimensions of the grid
            m_structuredGrid->SetDimensions(vectorX.size(),vectorP.size(),1);
        }
        break;
        case 2:  //TX
        {
            m_xlabel="Temperature (C)";
            m_ylabel="Salinity";
            m_actorScale[0]=m_actorScale_T;
            m_actorScale[1]=m_actorScale_X;
            m_actorScale[2]=m_actorScale_P;
            double dT=ui->doubleSpinBox_delta_firstVar->value();
            double Tmin=ui->doubleSpinBox_min_firstVar->value();
            double Tmax=ui->doubleSpinBox_max_firstVar->value();

            double dX=ui->doubleSpinBox_delta_secondVar->value();
            double Xmin=ui->doubleSpinBox_min_secondVar->value();
            double Xmax=ui->doubleSpinBox_max_secondVar->value();

            double P0=ui->doubleSpinBox_fixed_firstVar->value();

            vector<double> vectorX, vectorT;
            for (double X=Xmin; X<Xmax;X=X+dX)
            {
              vectorX.push_back(X);
            }
            for (double T=Tmin; T<Tmax; T=T+dT)
            {
              vectorT.push_back(T);
            }
            //calculate
            for(size_t j = 0; j < vectorX.size(); j++)
            {
                for(size_t i = 0; i < vectorT.size(); i++)
                {
                    SWEOS::cH2ONaCl eos;
                    eos.prop_pTX(P0*1e5,vectorT[i]+SWEOS::Kelvin,vectorX[j]);
                    points->InsertNextPoint(vectorT[i],vectorX[j],P0);
                    arrPhaseRegion->InsertNextValue(eos.m_prop.Region);
                    arrDensity->InsertNextValue(eos.m_prop.Rho);
                    arrDensity_l->InsertNextValue(eos.m_prop.Rho_l);
                    arrDensity_v->InsertNextValue(eos.m_prop.Rho_v);
                    arrDensity_h->InsertNextValue(eos.m_prop.Rho_h);
                    arrH->InsertNextValue(eos.m_prop.H);
                    arrH_l->InsertNextValue(eos.m_prop.H_l);
                    arrH_v->InsertNextValue(eos.m_prop.H_v);
                    arrH_h->InsertNextValue(eos.m_prop.H_h);
                    arrS_l->InsertNextValue(eos.m_prop.S_l);
                    arrS_v->InsertNextValue(eos.m_prop.S_v < 1e-30 ? 0 : eos.m_prop.S_v);
                    arrS_h->InsertNextValue(eos.m_prop.S_h < 1e-30 ? 0: eos.m_prop.S_h);
                    arrMu_l->InsertNextValue(eos.m_prop.Mu_l);
                    arrMu_v->InsertNextValue(eos.m_prop.Mu_v);
                    arrX_l->InsertNextValue(eos.m_prop.X_l);
                    arrX_v->InsertNextValue(eos.m_prop.X_v);
                }
            }
            // Specify the dimensions of the grid
            m_structuredGrid->SetDimensions(vectorT.size(),vectorX.size(),1);
        }
        break;
    }
    m_structuredGrid->SetPoints(points);
    //        m_structuredGrid->GetPointData()->SetScalars(doublevalue);
    m_structuredGrid->GetPointData()->AddArray(arrPhaseRegion);
    m_structuredGrid->GetPointData()->AddArray(arrDensity);
    m_structuredGrid->GetPointData()->AddArray(arrDensity_l);
    m_structuredGrid->GetPointData()->AddArray(arrDensity_v);
    m_structuredGrid->GetPointData()->AddArray(arrDensity_h);
    m_structuredGrid->GetPointData()->AddArray(arrH);
    m_structuredGrid->GetPointData()->AddArray(arrH_l);
    m_structuredGrid->GetPointData()->AddArray(arrH_v);
    m_structuredGrid->GetPointData()->AddArray(arrH_h);
    m_structuredGrid->GetPointData()->AddArray(arrS_l);
    m_structuredGrid->GetPointData()->AddArray(arrS_v);
    m_structuredGrid->GetPointData()->AddArray(arrS_h);
    m_structuredGrid->GetPointData()->AddArray(arrMu_l);
    m_structuredGrid->GetPointData()->AddArray(arrMu_v);
    m_structuredGrid->GetPointData()->AddArray(arrX_l);
    m_structuredGrid->GetPointData()->AddArray(arrX_v);

    ShowProps_2D(ui->comboBox_selectProps->currentIndex(), m_xlabel, m_ylabel, m_zlabel,m_actorScale);

}

void MainWindow::ShowProps_2D(int index_prop, std::string xlabel,std::string ylabel, std::string zlabel , double scale_actor[3])
{
    if(m_structuredGrid->GetPointData()->GetNumberOfArrays()==0)return;

    vtkSmartPointer<vtkStructuredGrid> grid=vtkSmartPointer<vtkStructuredGrid>::New();
    grid->SetDimensions(m_structuredGrid->GetDimensions());
    grid->SetPoints(m_structuredGrid->GetPoints());
    grid->GetPointData()->SetScalars(m_structuredGrid->GetPointData()->GetArray(index_prop));

    // Create a mapper and actor
    vtkSmartPointer<vtkDataSetMapper> gridMapper = vtkSmartPointer<vtkDataSetMapper>::New();
    gridMapper->SetInputData(grid);
    gridMapper->SetScalarRange(grid->GetScalarRange());

    // lookup table
    vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
    lut->SetNumberOfColors(256);
    lut->SetHueRange(0.0,0.667);
    gridMapper->SetLookupTable(lut);
    vtkSmartPointer<vtkActor> gridActor = vtkSmartPointer<vtkActor>::New();
    gridActor->SetMapper(gridMapper);

    // Create a renderer, render window, and interactor
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->AddActor(gridActor);
    //    gridActor->SetScale((SWEOS::PMAX-SWEOS::PMIN)/1e5/(SWEOS::TMAX-SWEOS::TMIN),1,(SWEOS::PMAX-SWEOS::PMIN)/1e5/(SWEOS::XMAX-SWEOS::XMIN));
    gridActor->SetScale(scale_actor);

    //axis actor
    vtkSmartPointer<vtkCubeAxesActor> axis=vtkSmartPointer<vtkCubeAxesActor>::New();
    axis->SetCamera(renderer->GetActiveCamera());
    axis->SetBounds(gridActor->GetBounds());
    renderer->AddActor(axis);
    axis->SetUse2DMode(1); //set font size

    //color scale
    vtkSmartPointer<vtkScalarBarActor> scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
    scalarBar->SetLookupTable(gridMapper->GetLookupTable());
    scalarBar->SetLookupTable( lut );
    scalarBar->SetTitle(grid->GetPointData()->GetArray(0)->GetName());
    scalarBar->SetDrawFrame(true);
    scalarBar->SetWidth(scalarBar->GetWidth()/2);
    renderer->AddActor2D(scalarBar);

    InitCubeAxes(axis,vtkBoundingBox(gridActor->GetBounds()),vtkBoundingBox(m_structuredGrid->GetBounds()),xlabel,ylabel,zlabel);
    renderer->SetBackground(0,0,0); //
    SetCamera(renderer,vtkBoundingBox(gridActor->GetBounds()));

    // before adding new renderer, remove all the old renderer, always keep only renderer in m_renderwindow
    m_renderWindow->GetRenderers()->RemoveAllItems();
    this->ui->qvtkWidget2->GetRenderWindow()->AddRenderer(renderer);
    renderer->Render();
    // Render and interact
    vtkSmartPointer<vtkInteractorStyleImage> style = vtkSmartPointer<vtkInteractorStyleImage>::New();
    // vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
    vtkSmartPointer<vtkRenderWindowInteractor> iren=vtkSmartPointer<vtkRenderWindowInteractor>::New();
    this->ui->qvtkWidget2->GetRenderWindow()->GetInteractor()->SetInteractorStyle(style);
    this->ui->qvtkWidget2->GetRenderWindow()->Render();
}

void MainWindow::on_pushButton_clicked()
{
    switch (m_dimension) {
        case 1:
        {
            Calculate_Diagram1D();
        }
        break;
        case 2:
        {
            Calculate_Diagram2D();
        }
        break;
    case 3:
        ui->textEdit->append("3D is comming soon");
        break;

    }
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
        update1dUI(ui->comboBox_selectVariable->currentText(), ui->comboBox_selectProps->currentIndex());
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
        if(m_dimension==1)
        {
            // chart options
            ui->groupBox_chartOptions->setVisible(true);
            QRect geometry=ui->groupBox_chartOptions->geometry();
            ui->groupBox_chartOptions->setGeometry(ui->groupBox_fixed_Var->geometry().x()+ui->groupBox_fixed_Var->geometry().width()+10, geometry.y(),geometry.width(),geometry.height());
            // button of calculation
            QRect geometry3=ui->pushButton->geometry();
            ui->pushButton->setGeometry(ui->groupBox_chartOptions->geometry().x()+ui->groupBox_chartOptions->geometry().width()+10, geometry3.y(),geometry3.width(),geometry3.height());
        }else
        {
            // chart options
            ui->groupBox_chartOptions->setVisible(false);
            QRect geometry3=ui->pushButton->geometry();
            ui->pushButton->setGeometry(ui->groupBox_fixed_Var->geometry().x()+ui->groupBox_fixed_Var->geometry().width()+10, geometry3.y(),geometry3.width(),geometry3.height());
        }
    }else
    {
        if(m_dimension==1)
        {
            //chart options
            ui->groupBox_chartOptions->setVisible(true);
            QRect geometry=ui->groupBox_chartOptions->geometry();
            ui->groupBox_chartOptions->setGeometry(ui->groupBox_Variables->geometry().x()+ui->groupBox_Variables->geometry().width()+10, geometry.y(),geometry.width(),geometry.height());
            // button of calculation
            QRect geometry3=ui->pushButton->geometry();
            ui->pushButton->setGeometry(ui->groupBox_chartOptions->geometry().x()+ui->groupBox_chartOptions->geometry().width()+10, geometry3.y(),geometry3.width(),geometry3.height());
        }else
        {
            // chart options
            ui->groupBox_chartOptions->setVisible(false);
            // button of calculation
            QRect geometry3=ui->pushButton->geometry();
            ui->pushButton->setGeometry(ui->groupBox_Variables->geometry().x()+ui->groupBox_Variables->geometry().width()+10, geometry3.y(),geometry3.width(),geometry3.height());
        }
    }

}
void MainWindow::update1dUI_chartOptions(int index_propSelection)
{
    //update UI of chart options for 1D chart
    ui->checkBox->setVisible(true);
    ui->checkBox->setText("Point marker");
    ui->checkBox_2->setVisible(false);
    ui->checkBox_3->setVisible(false);
    ui->checkBox_4->setVisible(false);
    ui->checkBox_5->setVisible(false);
    switch (index_propSelection) {
        case 0: //phase region
        {
            ui->checkBox->setText("Point marker");
        }
        break;
        case 1: //density
        {
            ui->checkBox_2->setVisible(true);
            ui->checkBox_2->setChecked(true);
            ui->checkBox_2->setText("Bulk density");
            ui->checkBox_3->setVisible(true);
            ui->checkBox_3->setChecked(true);
            ui->checkBox_3->setText("Liquid density");
            ui->checkBox_4->setVisible(true);
            ui->checkBox_4->setChecked(true);
            ui->checkBox_4->setText("Vapour density");
            ui->checkBox_5->setVisible(true);
            ui->checkBox_5->setChecked(true);
            ui->checkBox_5->setText("Halite density");
        }
        break;
        case 2: //enthalpy
        {
            ui->checkBox_2->setVisible(true);
            ui->checkBox_2->setChecked(true);
            ui->checkBox_2->setText("Bulk enthalpy");
            ui->checkBox_3->setVisible(true);
            ui->checkBox_3->setChecked(true);
            ui->checkBox_3->setText("Liquid enthalpy");
            ui->checkBox_4->setVisible(true);
            ui->checkBox_4->setChecked(true);
            ui->checkBox_4->setText("Vapour enthalpy");
            ui->checkBox_5->setVisible(true);
            ui->checkBox_5->setChecked(true);
            ui->checkBox_5->setText("Halite enthalpy");
        }
        break;
        case 3: //saturation
        {
            ui->checkBox_2->setVisible(true);
            ui->checkBox_2->setChecked(true);
            ui->checkBox_2->setText("Liquid saturation");
            ui->checkBox_3->setVisible(true);
            ui->checkBox_3->setChecked(true);
            ui->checkBox_3->setText("Vapour saturation");
            ui->checkBox_4->setVisible(true);
            ui->checkBox_4->setChecked(true);
            ui->checkBox_4->setText("Halite saturation");
        }
        break;
        case 4: //viscosity
        {
            ui->checkBox_2->setVisible(true);
            ui->checkBox_2->setChecked(true);
            ui->checkBox_2->setText("Liquid viscosity");
            ui->checkBox_3->setVisible(true);
            ui->checkBox_3->setChecked(true);
            ui->checkBox_3->setText("Vapour viscosity");
        }
        break;
        case 5: //viscosity
        {
            ui->checkBox_2->setVisible(true);
            ui->checkBox_2->setChecked(true);
            ui->checkBox_2->setText("Liquid salinity");
            ui->checkBox_3->setVisible(true);
            ui->checkBox_3->setChecked(true);
            ui->checkBox_3->setText("Vapour salinity");
        }
        break;
    }
}
void MainWindow::update1dUI(QString arg, int index_propSelection)
{
    //update properties selection
    ui->comboBox_selectProps->clear();
    ui->comboBox_selectProps->addItem("Phase Region");
    ui->comboBox_selectProps->addItem("Density");
    ui->comboBox_selectProps->addItem("Enthalpy");
    ui->comboBox_selectProps->addItem("Saturation");
    ui->comboBox_selectProps->addItem("Viscosity");
    ui->comboBox_selectProps->addItem("Salinity");
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

    // chart options
    update1dUI_chartOptions(index_propSelection);
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
    label->setText("Temperature (C)");
    box->setDecimals(2);
    box->setRange(SWEOS::TMIN-SWEOS::Kelvin, SWEOS::TMAX-SWEOS::Kelvin);
    box->setSingleStep(1);
    box->setValue(defaultValue);
}
void MainWindow::UpdateUI_fixedH(QLabel* label, QDoubleSpinBox* box, double defaultValue)
{
    label->setText("Enthalpy (kJ/kg)");
    box->setDecimals(4);
    box->setRange(1, 5000);
    box->setSingleStep(1);
    box->setValue(defaultValue);
}
void MainWindow::UpdateUI_X(QLabel* label, QDoubleSpinBox* deltaBox, QDoubleSpinBox* maxBox, QDoubleSpinBox* minBox)
{
    label->setText("dX:");
    deltaBox->setDecimals(4);
    deltaBox->setRange(0.0001, 1);
    deltaBox->setSingleStep(0.001);
    deltaBox->setValue(0.01);

    maxBox->setDecimals(4);
    maxBox->setRange(SWEOS::XMIN, SWEOS::XMAX);
    maxBox->setSingleStep(0.001);
    maxBox->setValue(SWEOS::XMAX);

    minBox->setDecimals(4);
    minBox->setRange(SWEOS::XMIN, SWEOS::XMAX);
    minBox->setSingleStep(0.001);
    minBox->setValue(0.0001);
}
void MainWindow::UpdateUI_T(QLabel* label, QDoubleSpinBox* deltaBox, QDoubleSpinBox* maxBox, QDoubleSpinBox* minBox)
{
    label->setText("dT(C):");
    deltaBox->setDecimals(2);
    deltaBox->setRange(0.01, 100);
    deltaBox->setSingleStep(1);
    deltaBox->setValue(1);

    maxBox->setDecimals(2);
    maxBox->setRange(SWEOS::TMIN-SWEOS::Kelvin, SWEOS::TMAX-SWEOS::Kelvin);
    maxBox->setSingleStep(1);
    maxBox->setValue(673);

    minBox->setDecimals(2);
    minBox->setRange(SWEOS::TMIN-SWEOS::Kelvin,SWEOS::TMAX-SWEOS::Kelvin);
    minBox->setSingleStep(1);
    minBox->setValue(SWEOS::TMIN-SWEOS::Kelvin);
}
void MainWindow::UpdateUI_P(QLabel* label, QDoubleSpinBox* deltaBox, QDoubleSpinBox* maxBox, QDoubleSpinBox* minBox)
{
    label->setText("dP (bar):");
    deltaBox->setDecimals(2);
    deltaBox->setRange(0.1, 100);
    deltaBox->setSingleStep(1);
    deltaBox->setValue(1);

    maxBox->setDecimals(2);
    maxBox->setRange(SWEOS::PMIN/1E5, SWEOS::PMAX/1E5);
    maxBox->setSingleStep(1);
    maxBox->setValue(400);

    minBox->setDecimals(2);
    minBox->setRange(SWEOS::PMIN/1E5, SWEOS::PMAX/1E5);
    minBox->setSingleStep(1);
    minBox->setValue(SWEOS::PMIN/1E5);
}
void MainWindow::update2dUI(QString arg)
{
    //update properties selection
    ui->comboBox_selectProps->clear();
    ui->comboBox_selectProps->addItem("Phase Region");
    ui->comboBox_selectProps->addItem("Bulk density");
    ui->comboBox_selectProps->addItem("Liquid density");
    ui->comboBox_selectProps->addItem("Vapour density");
    ui->comboBox_selectProps->addItem("Halite density");
    ui->comboBox_selectProps->addItem("Bulk enthalpy");
    ui->comboBox_selectProps->addItem("Liquid enthalpy");
    ui->comboBox_selectProps->addItem("Vapour enthalpy");
    ui->comboBox_selectProps->addItem("Halite enthalpy");
    ui->comboBox_selectProps->addItem("Liquid saturation");
    ui->comboBox_selectProps->addItem("Vapour saturation");
    ui->comboBox_selectProps->addItem("Halite saturation");
    ui->comboBox_selectProps->addItem("Liquid viscosity");
    ui->comboBox_selectProps->addItem("Vapour viscosity");
    ui->comboBox_selectProps->addItem("Liquid salinity");
    ui->comboBox_selectProps->addItem("Vapour salinity");
    if(arg=="PT")
    {
        //fixed vars
        UpdateUI_fixedX(ui->label_fixed_firsVar,ui->doubleSpinBox_fixed_firstVar);
        //independent variable
        UpdateUI_P(ui->label_delta_firstVar, ui->doubleSpinBox_delta_firstVar, ui->doubleSpinBox_max_firstVar, ui->doubleSpinBox_min_firstVar);
        UpdateUI_T(ui->label_delta_secondVar, ui->doubleSpinBox_delta_secondVar, ui->doubleSpinBox_max_secondVar, ui->doubleSpinBox_min_secondVar);

    }else if(arg=="PX")
    {
        //fixed vars
        UpdateUI_fixedT(ui->label_fixed_firsVar,ui->doubleSpinBox_fixed_firstVar);
        //independent variable
        UpdateUI_P(ui->label_delta_firstVar, ui->doubleSpinBox_delta_firstVar, ui->doubleSpinBox_max_firstVar, ui->doubleSpinBox_min_firstVar);
        UpdateUI_X(ui->label_delta_secondVar, ui->doubleSpinBox_delta_secondVar, ui->doubleSpinBox_max_secondVar, ui->doubleSpinBox_min_secondVar);

    }else if(arg=="TX")
    {
        //fixed vars
        UpdateUI_fixedP(ui->label_fixed_firsVar,ui->doubleSpinBox_fixed_firstVar);
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
    //update properties selection
    ui->comboBox_selectProps->clear();
    ui->comboBox_selectProps->addItem("Phase Region");
    ui->comboBox_selectProps->addItem("Bulk density");
    ui->comboBox_selectProps->addItem("Liquid density");
    ui->comboBox_selectProps->addItem("Vapour density");
    ui->comboBox_selectProps->addItem("Halite density");
    ui->comboBox_selectProps->addItem("Bulk enthalpy");
    ui->comboBox_selectProps->addItem("Liquid enthalpy");
    ui->comboBox_selectProps->addItem("Vapour enthalpy");
    ui->comboBox_selectProps->addItem("Halite enthalpy");
    ui->comboBox_selectProps->addItem("Liquid saturation");
    ui->comboBox_selectProps->addItem("Vapour saturation");
    ui->comboBox_selectProps->addItem("Halite saturation");
    ui->comboBox_selectProps->addItem("Liquid viscosity");
    ui->comboBox_selectProps->addItem("Vapour viscosity");
    ui->comboBox_selectProps->addItem("Liquid salinity");
    ui->comboBox_selectProps->addItem("Vapour salinity");
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

        ui->label_delta_firstVar->setText("dT(C):");
        ui->doubleSpinBox_delta_firstVar->setDecimals(2);
        ui->doubleSpinBox_delta_firstVar->setRange(SWEOS::TMIN-SWEOS::Kelvin, SWEOS::TMAX-SWEOS::Kelvin);
        ui->doubleSpinBox_delta_firstVar->setSingleStep(1);
        ui->doubleSpinBox_delta_firstVar->setValue(SWEOS::TMIN-SWEOS::Kelvin);

        ui->doubleSpinBox_max_firstVar->setDecimals(2);
        ui->doubleSpinBox_max_firstVar->setRange(SWEOS::TMIN-SWEOS::Kelvin, SWEOS::TMAX-SWEOS::Kelvin);
        ui->doubleSpinBox_max_firstVar->setSingleStep(1);
        ui->doubleSpinBox_max_firstVar->setValue(SWEOS::TMAX-SWEOS::Kelvin);

        ui->doubleSpinBox_min_firstVar->setDecimals(2);
        ui->doubleSpinBox_min_firstVar->setRange(0.1,100);
        ui->doubleSpinBox_min_firstVar->setSingleStep(1);
        ui->doubleSpinBox_min_firstVar->setValue(SWEOS::TMIN-SWEOS::Kelvin);

    }else if(arg=="Pressure")
    {
        ui->label_fixed_firsVar->setText("Temperature (C)");
        ui->doubleSpinBox_fixed_firstVar->setDecimals(2);
        ui->doubleSpinBox_fixed_firstVar->setRange(SWEOS::TMIN-SWEOS::Kelvin, SWEOS::TMAX-SWEOS::Kelvin);
        ui->doubleSpinBox_fixed_firstVar->setSingleStep(1);
        ui->doubleSpinBox_fixed_firstVar->setValue(100);
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
        ui->label_fixed_secondVar->setText("Temperature (C)");
        ui->doubleSpinBox_fixed_secondVar->setDecimals(2);
        ui->doubleSpinBox_fixed_secondVar->setRange(SWEOS::TMIN-SWEOS::Kelvin, SWEOS::TMAX-SWEOS::Kelvin);
        ui->doubleSpinBox_fixed_secondVar->setSingleStep(1);
        ui->doubleSpinBox_fixed_secondVar->setValue(100);

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
        update1dUI(arg1, ui->comboBox_selectProps->currentIndex());
    break;
    case 2:
        update2dUI(arg1);
    break;
    case 3:
        update3dUI(arg1);
    break;
    }
}

void MainWindow::on_actionSave_triggered()
{
    int index_tab=ui->tabWidget->currentIndex();
    switch (index_tab) {
        case 0:
        {
            if(m_calculationMode==CALCULATION_MULTI_POINTS)
            {
                std::string filter_ext, title_dlg;
                QString fileName;
                filter_ext="CSV File (*.csv);;Delimited Text File (*.txt)";
                title_dlg="Save Multiple Points Calculation Results to File: Scatter";
                fileName = QFileDialog::getSaveFileName(this, tr(title_dlg.c_str()), "", tr(filter_ext.c_str()));
                if (!fileName.isNull())
                {
                    vtkSmartPointer<vtkDelimitedTextWriter> writer =  vtkSmartPointer<vtkDelimitedTextWriter>::New();
                    writer->SetFileName(fileName.toStdString().c_str());
                    writer->SetInputData(m_vtkTable);
                    writer->Write();
                }
            }

        }
        break;
        case 1:
        {
            std::string filter_ext, title_dlg;
            QString fileName;
            switch (m_dimension) {
                case 1:
                {
                    filter_ext="CSV File (*.csv);;Delimited Text File (*.txt)";
                    title_dlg="Save Diagram Calculation Results to File: 1D";
                    fileName = QFileDialog::getSaveFileName(this, tr(title_dlg.c_str()), "", tr(filter_ext.c_str()));
                    if (!fileName.isNull())
                    {
                        vtkSmartPointer<vtkDelimitedTextWriter> writer =  vtkSmartPointer<vtkDelimitedTextWriter>::New();
                        writer->SetFileName(fileName.toStdString().c_str());
                        writer->SetInputData(m_vtkTable);
                        writer->Write();
                    }
                    else
                    {

                    }
                }
                break;
            case 2:
            {
                filter_ext="VTK File (*.vtk)";
                title_dlg="Save Diagram Calculation Results to File: 2D";
                fileName = QFileDialog::getSaveFileName(this, tr(title_dlg.c_str()), "", tr(filter_ext.c_str()));
                if (!fileName.isNull())
                {
                    vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();
                    writer->SetFileName(fileName.toStdString().c_str());
                    writer->SetInputData(m_structuredGrid);
                    writer->Write();
                }
                else
                {

                }
            }
                break;
            }
        }
        break;
    }

}

void MainWindow::on_comboBox_selectProps_activated(const QString &)
{
    switch (m_dimension) {
    case 1:
        update1dUI_chartOptions(ui->comboBox_selectProps->currentIndex());
        ShowProps_1D();
        break;
    case 2:
        ShowProps_2D(ui->comboBox_selectProps->currentIndex(),m_xlabel,m_ylabel,m_zlabel,m_actorScale);
        break;
    }
}

void MainWindow::on_checkBox_stateChanged(int )
{
    m_showScatter_1Dchart=ui->checkBox->isChecked();
    ShowProps_1D();
}

void MainWindow::on_checkBox_2_stateChanged(int )
{
    ShowProps_1D();
}

void MainWindow::on_checkBox_3_stateChanged(int )
{
    ShowProps_1D();
}

void MainWindow::on_checkBox_4_stateChanged(int )
{
    ShowProps_1D();
}

void MainWindow::on_checkBox_5_stateChanged(int )
{
    ShowProps_1D();
}
void MainWindow::on_checkBox_6_stateChanged(int )
{
    m_showPhaseRegion_1Dchart=ui->checkBox_6->isChecked();
    ShowProps_1D();
}

void MainWindow::on_comboBox_activated(const QString &)
{
    updateScatterCalculationUI(ui->comboBox->currentIndex());
}

void MainWindow::updateScatterCalculationUI(int index_varsSelection)
{
    switch (index_varsSelection) {
    case 0:
    {
        UpdateUI_fixedT(ui->label_3,ui->doubleSpinBox_2);
    }
        break;
    case 1:
    {
        UpdateUI_fixedH(ui->label_3, ui->doubleSpinBox_2);
    }
        break;
    }
}
