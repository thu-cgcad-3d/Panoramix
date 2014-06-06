
#include "mainwin.hpp"


MainWind::MainWind(QWidget * parent) : QMainWindow(parent) {
    _ui.setupUi(this);
    _w = new OgreWidget;
    setCentralWidget(_w);
    initGui();
}


void MainWind::initGui() {

    // init actions


}

