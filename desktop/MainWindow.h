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
#include <vector>
using namespace std;

#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>
#include <vtkChartXY.h>
#include <vtkTable.h>
#include <vtkPlot.h>
#include <vtkFloatArray.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>
#include <vtkPen.h>
#include <vtkAxis.h>

#define CALCULATION_SINGLE_POINT 1
#define CALCULATION_MULTI_POINTS 2

//definition for 1,2,3Dimension
#define oneDim_Temperature 1
#define oneDim_Pressure 2
#define oneDim_Salinity 3

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
    vtkSmartPointer<vtkTable> m_vtkTable;
protected slots:

private slots:
  void updateCalculationModelSelection(bool isSinglePoint);
  void on_pushButton_2_clicked();

  void on_radioButton_pressed();

  void on_radioButton_clicked();

  void on_radioButton_2_clicked();

  void on_tabWidget_currentChanged(int index);

  void on_pushButton_clicked();



  void on_radioButton_3_clicked();

  void on_radioButton_4_clicked();

  void on_radioButton_5_clicked();

  void on_comboBox_2_activated(const QString &arg1);
    void update1dUI(QString arg);
private:

  vtkSmartPointer<vtkQtTableView>         TableView;

  // Designer form
  Ui::MainWindow *ui;
};

#endif // MainWindow_H
