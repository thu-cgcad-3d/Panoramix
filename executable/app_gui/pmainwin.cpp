
#include "pprojectwidget.h"
#include "pmainwin.h"

#include "../src/Panoramix.Core"
#include "../src/Panoramix.Vis"

PMainWin::PMainWin(QWidget * parent)
    : QMainWindow(parent) {

    // init gui
    // mdi
    _mdiArea = new QMdiArea(this);
    QLinearGradient grad(0, 0, 0, 1600);
    grad.setColorAt(0, qRgb(200, 200, 255));
    grad.setColorAt(0.2, qRgb(240, 240, 240));
    grad.setColorAt(1.0, Qt::black);
    _mdiArea->setBackground(grad);
    setCentralWidget(_mdiArea);
    QMdiSubWindow * swin = _mdiArea->addSubWindow(new PProjectWidget, Qt::Widget);
    swin->resize(400, 400);
    // menu bar
    
    // status bar
    

    // style sheet    
    QFile cssFile(":/qss/gui.css");
    if (cssFile.open(QFile::ReadOnly)){
        setStyleSheet(tr(cssFile.readAll()));
    }

    setWindowTitle(tr("GUI for Panoramix"));
}

PMainWin::~PMainWin () {

}