#include "pprojectwidget.h"

PProjectWidget::PProjectWidget(QWidget * parent) 
    : QWidget(parent) {
    
    QVBoxLayout * vbox = new QVBoxLayout;
    vbox->setMargin(1);
    setLayout(vbox);
    _tabWidget = new QTabWidget;
    vbox->addWidget(_tabWidget);
    _tabWidget->setTabShape(QTabWidget::TabShape::Rounded);
    _tabWidget->setTabPosition(QTabWidget::TabPosition::West);
    _tabWidget->addTab(new QWidget, tr("Original"));
    _tabWidget->addTab(new QWidget, tr("Reconstruction"));

}

PProjectWidget::~PProjectWidget(){

}