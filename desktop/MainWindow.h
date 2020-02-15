/*=========================================================================

  Program:   Visualization Toolkit
  Module:    MainWindow.h
  Language:  C++

  Copyright 2009 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/
#ifndef MainWindow_H
#define MainWindow_H

#include "vtkSmartPointer.h"    // Required for smart pointer internal ivars.
#include <QMainWindow>
#include "H2ONaCl.H"
#include <QDateTime>
#include <QTime>
#include <QMessageBox>
#include <QFileDialog>
#include <QProgressBar>
#include <QThread>
#include <QtConcurrent>
#include <QFuture>

#include <vector>
#include <fstream>
using namespace std;
#include "MainWindow.h"
#include "ui_MainWindow.h"

#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>
#include <vtkChartXY.h>
#include <vtkTable.h>
#include <vtkDoubleArray.h>
#include <vtkPlot.h>
#include <vtkFloatArray.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>
#include <vtkPen.h>
#include <vtkAxis.h>
#include <vtkTextProperty.h>
#include <vtkAxisActor2D.h>
#include <vtkViewport.h>
#include <vtkAxisActor.h>
#include <vtkXYPlotActor.h>
#include <vtkFieldData.h>
#include <vtkNamedColors.h>
#include <vtkStructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkTooltipItem.h>
#include <vtkChartLegend.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkWarpScalar.h>
#include <vtkPolyDataNormals.h>
#include <vtkClipPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkStructuredGridWriter.h>
#include <vtkStructuredData.h>
#include <vtkPointData.h>
#include <vtkImageMapToColors.h>
#include <vtkCubeAxesActor.h>
#include <vtkCubeAxesActor2D.h>
#include <vtkAxisActor.h>
#include <vtkBoundingBox.h>
#include <vtkAxesActor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkInteractorStyleImage.h>
#include <vtkScalarBarActor.h>
#include <vtkColorSeries.h>
#include <vtkDelimitedTextWriter.h>
#include <vtkPlotPoints.h>
#include <vtkPlotArea.h>
#include <vtkBrush.h>
#include <vtkColorSeries.h>
#include <vtkLegendBoxActor.h>
#include <vtkCubeSource.h>

#define CALCULATION_SINGLE_POINT 1
#define CALCULATION_MULTI_POINTS 2

//definition for 1,2,3Dimension
#define oneDim_Temperature 1
#define oneDim_Pressure 2
#define oneDim_Salinity 3

//define using PTX or PHX
#define USING_PTX 0
#define USING_PHX 1

//define range of enthalpy kJ/kg
#define MIN_ENTHALPY 1
#define MAX_ENTHALPY 5000

//define props index in vtkTable
#define INDEX_T_VTKTABLE 0
#define INDEX_P_VTKTABLE 1
#define INDEX_X_VTKTABLE 2
#define INDEX_REGION_VTKTABLE 3
#define INDEX_RHO_VTKTABLE 4
#define INDEX_RHO_L_VTKTABLE 5
#define INDEX_RHO_V_VTKTABLE 6
#define INDEX_RHO_H_VTKTABLE 7
#define INDEX_H_VTKTABLE 8
#define INDEX_H_L_VTKTABLE 9
#define INDEX_H_V_VTKTABLE 10
#define INDEX_H_H_VTKTABLE 11
#define INDEX_S_L_VTKTABLE 12
#define INDEX_S_V_VTKTABLE 13
#define INDEX_S_H_VTKTABLE 14
#define INDEX_MU_L_VTKTABLE 15
#define INDEX_MU_V_VTKTABLE 16
#define INDEX_X_L_VTKTABLE 17
#define INDEX_X_V_VTKTABLE 18

//camera view
#define ID_CAMERA_GENERAL	1
#define ID_CAMERA_FRONT	2
#define ID_CAMERA_BACK	3
#define ID_CAMERA_LEFT	4
#define ID_CAMERA_RIGHT	5
#define ID_CAMERA_UP	6
#define ID_CAMERA_DOWN	7

// Forward Qt class declarations
class Ui_MainWindow;

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

// Forward VTK class declarations
class vtkQtTableView;


class MainWindow : public QMainWindow
{
  Q_OBJECT

public:

  // Constructor/Destructor
  MainWindow(QWidget *parent = nullptr);
  ~MainWindow() override;

public slots:

  virtual void slotOpenFile();
  virtual void slotExit();

protected:
  int m_calculationMode;
  double m_IndependentVar1_old, m_IndependentVar2_old, m_IndependentVar3_old;
  int m_dimension, m_calculationMode_123Dim;
//    vtk variable
  vtkSmartPointer<vtkContextView> m_vtkChartView;
  vtkSmartPointer<vtkGenericOpenGLRenderWindow> m_renderWindow;
    vtkSmartPointer<vtkTable> m_vtkTable;
  vtkSmartPointer<vtkStructuredGrid> m_structuredGrid;
  int m_vtkFontSize;
  int m_index_var;
  double m_actorScale_T, m_actorScale_P, m_actorScale_X;
  std::string m_xlabel, m_ylabel, m_zlabel;
  double m_actorScale[3];
  QProgressBar* m_progressBar;
  bool m_showScatter_1Dchart;
  vtkColor3d m_VTKbackgroundColor;
  vtkSmartPointer<vtkColorSeries> m_vtkColorSeries_Lines_1Dchart, m_vtkColorSeries_PhaseRegion;
  double m_alphaPhaseRegion;
  int m_vtkLineWidth;
  bool m_showPhaseRegion_1Dchart;
  bool m_vtkCameraInitialized;
  bool m_resetChartRange;
protected slots:

private slots:
  void updateCalculationModelSelection(bool isSinglePoint);
  int SetCamera(vtkSmartPointer<vtkRenderer> renderer, vtkBoundingBox boundingbox, int type=ID_CAMERA_UP);
  int InitCubeAxes(vtkCubeAxesActor* axes, vtkBoundingBox boundingbox, vtkBoundingBox rangebox, std::string xlabel, std::string ylabel, std::string zlabel,int fontsize=30);
  void Calculate_Diagram1D();
  void Calculate_Diagram2D();
  void on_pushButton_2_clicked();

  void on_radioButton_pressed();

  void on_radioButton_clicked();

  void on_radioButton_2_clicked();

  void on_tabWidget_currentChanged(int index);

  void on_pushButton_clicked();

  void on_radioButton_3_clicked();

  void on_radioButton_4_clicked();

  void on_radioButton_5_clicked();

  void on_comboBox_selectVariable_activated(const QString &arg1);
    void UpdateUI_fixedT(QLabel* label, QDoubleSpinBox* box, double defaultValue=100);
    void UpdateUI_fixedH(QLabel* label, QDoubleSpinBox* box, double defaultValue=2000);
    void UpdateUI_fixedP(QLabel* label, QDoubleSpinBox* box, double defaultValue=316);
    void UpdateUI_fixedX(QLabel* label, QDoubleSpinBox* box, double defaultValue=0.032);
    void updateScatterCalculationUI(int index_varsSelection); //PTX, PHX
    void update1dUI(QString arg, int index_propSelection);
    void update1dUI_chartOptions(int index_propSelection);
    void update2dUI(QString arg);
    void update3dUI(QString arg);
    void initRenderWindow();
    void updateUI(int dim); //update UI of phase diagram
    void updateUILayout(bool show_secondVariable, bool show_thirdVariable, bool show_firstFixedVar, bool show_secondFixedVar,bool show_groupbox_fixedVars, double shinkWidth_groupbox_Vars);
    void on_actionSave_triggered();
    void update1dChart(int index_var, std::string name_prop, std::vector<int> index_props, std::vector<bool> components, vtkSmartPointer<vtkColorSeries> colorseries);

    void ShowProps_1D();
    void ShowProps_2D(int index_prop, std::string xlabel,std::string ylabel, std::string zlabel , double scale_actor[3]);
    void SinglePointCalculation(int index_varsSelection);
    void on_comboBox_selectProps_activated(const QString &arg1);
    void CalculateProps_PTX_PHX(int PTX_PHX, std::vector<double> arrT_H,std::vector<double> arrP, std::vector<double> arrX, vtkSmartPointer<vtkTable> table);

    void on_checkBox_stateChanged(int arg1);

    void on_checkBox_2_stateChanged(int arg1);

    void on_checkBox_3_stateChanged(int arg1);

    void on_checkBox_4_stateChanged(int arg1);

    void on_checkBox_5_stateChanged(int arg1);
    void GetMaxMin_vtkTableColumn(const vtkSmartPointer<vtkTable> table, int index_col, double& min, double& max);

    void on_checkBox_6_stateChanged(int arg1);

    void on_comboBox_activated(const QString &arg1);

    void on_doubleSpinBox_fixed_firstVar_valueChanged(double arg1);

    void on_doubleSpinBox_fixed_secondVar_valueChanged(double arg1);
    void updateEnthalpyRange_UI_1D();

    void on_vtkCameraControl_reset_triggered();

private:
    QRect m_geometry_Groupbox_variables;
    void UpdateUI_P(QLabel* label, QDoubleSpinBox* deltaBox, QDoubleSpinBox* maxBox, QDoubleSpinBox* minBox);
    void UpdateUI_T(QLabel* label, QDoubleSpinBox* deltaBox, QDoubleSpinBox* maxBox, QDoubleSpinBox* minBox);
    void UpdateUI_H(QLabel* label, QDoubleSpinBox* deltaBox, QDoubleSpinBox* maxBox, QDoubleSpinBox* minBox, double pMinMax[2], double XMinMax[2]);
    void UpdateUI_X(QLabel* label, QDoubleSpinBox* deltaBox, QDoubleSpinBox* maxBox, QDoubleSpinBox* minBox);

  vtkSmartPointer<vtkQtTableView>         TableView;

  // Designer form
  Ui::MainWindow *ui;
};

#endif // MainWindow_H



