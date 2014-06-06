#pragma once

#include <QtGui>
#include <QtWidgets>

#include "ui_mainwin.h"
#include "ogrewidget.hpp"

class MainWind : public QMainWindow {

    Q_OBJECT

public:
    MainWind(QWidget *parent = 0);

    void initGui();

private:
    OgreWidget * _w;
    Ui::MainWind _ui;
};