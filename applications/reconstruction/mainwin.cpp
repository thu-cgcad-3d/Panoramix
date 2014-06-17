
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

void MainWind::on_actionInsert_Panorama_triggered() {
    QString fn = QFileDialog::getOpenFileName(this, tr("Choose a panoramic image file"), tr(PROJECT_DATA_DIR_STR),
        tr("Image file (*.png;*.jpg);;All Files (*.*)"));
    if (fn.isEmpty())
        return;
    _w->setupPanorama(fn);
}

void MainWind::on_actionInsert_View_triggered() {
    qDebug() << "new view inserted!";
}

