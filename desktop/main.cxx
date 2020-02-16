/*
 * Copyright 2007 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
// QT includes
#include <QApplication>
#include <QSurfaceFormat>

#include "QVTKOpenGLWidget.h"
#include "MainWindow.h"
#include "SWEOSbash.h"

extern int qInitResources_icons();

int main( int argc, char** argv )
{
 if(argc==2 || argc==1)
 {
    // needed to ensure appropriate OpenGL context is created for VTK rendering.
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLWidget::defaultFormat());   // must be here
      // QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
    // QT Stuff
    QApplication app( argc, argv );

    QApplication::setStyle("fusion");

    qInitResources_icons();

    MainWindow myMainWindow;
  //  myMainWindow.showMaximized();
    myMainWindow.show();

    return app.exec();
 }else
 {
   SWEOSbash::bash_run(argc, argv);
 }
}
