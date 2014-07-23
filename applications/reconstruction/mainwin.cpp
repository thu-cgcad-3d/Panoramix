
#include "mainwin.hpp"


MainWind::MainWind(QWidget * parent) : QMainWindow(parent) {
    _ui.setupUi(this);
    _w = new OgreWidget;
    setCentralWidget(_w);
    _thread = new WorkThread(this);
    initGui();
}


void MainWind::initGui() {

    // init actions
    _progressBar = new QProgressBar(this);
    statusBar()->addPermanentWidget(_progressBar);
    _progressBar->hide();
    _progressBar->setRange(0, 0);
    _progressBar->setFixedHeight(15);
    _progressBar->setFixedWidth(300);

    connect(_thread, SIGNAL(started()), _progressBar, SLOT(show()));
    connect(_thread, SIGNAL(finished()), _progressBar, SLOT(hide()));
}

void MainWind::on_actionInsert_Panorama_triggered() {
    QString fn = QFileDialog::getOpenFileName(this, tr("Choose a panoramic image file"), tr(PROJECT_DATA_DIR_STR),
        tr("Image file (*.png;*.jpg);;All Files (*.*)"));
    if (fn.isEmpty())
        return;
    _thread->start();
}

void MainWind::on_actionInsert_View_triggered() {
    qDebug() << "new view inserted!";
}

