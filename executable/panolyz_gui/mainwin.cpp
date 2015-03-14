#include "widgets.hpp"
#include "workthread.hpp"
#include "project.hpp"
#include "mainwin.hpp"

MainWin::MainWin(QWidget *parent) : QMainWindow(parent) {
   
    //// menus
    auto menuFile = menuBar()->addMenu(tr("&File"));
    auto menuView = menuBar()->addMenu(tr("&View"));
    auto menuSettings = menuBar()->addMenu(tr("&Settings"));
    auto menuHelp = menuBar()->addMenu(tr("Help"));

    //// actions
    // file
    auto actionNewProjPanorama = menuFile->addAction(tr("Create A &New Panolyz Project from Panorama ..."));
    actionNewProjPanorama->setShortcut(QKeySequence::New);
    connect(actionNewProjPanorama, &QAction::triggered, [this](){
        QStringList filenames = QFileDialog::getOpenFileNames(this, tr("Select a Panoramic Image File"), tr("."),
            tr("Image Files (*.png;*.jpg;*.jpeg);;All Files (*.*)"));
        for (auto & filename : filenames)
            installProject(new Project(filename, true, this));
    });
    auto actionNewProjNormal = menuFile->addAction(tr("Create A &New Panolyz Project from Normal Photo ..."));
    connect(actionNewProjNormal, &QAction::triggered, [this](){
        QStringList filenames = QFileDialog::getOpenFileNames(this, tr("Select a Normal Image File"), tr("."),
            tr("Image Files (*.png;*.jpg;*.jpeg);;All Files (*.*)"));
        for (auto & filename : filenames)
            installProject(new Project(filename, false, this));
    });
    auto actionOpenProj = menuFile->addAction(tr("&Open An Existing Panolyz Project ..."));
    actionOpenProj->setShortcut(QKeySequence::Open);
    connect(actionOpenProj, &QAction::triggered, [this](){
        QStringList filenames = QFileDialog::getOpenFileNames(this, tr("Select Panolyz Project Files"), tr("."),
            tr("Panolyz Project File (*.panoproj)"));
        for (auto & filename : filenames)
            installProject(new Project(filename, this));
    });

    menuFile->addSeparator();
    
    auto actionSaveProj = menuFile->addAction(tr("&Save Current Panolyz Project"));
    actionSaveProj->setShortcut(QKeySequence::Save);
    auto actionSaveAsProj = menuFile->addAction(tr("&Save Current Panolyz Project As ..."));
    actionSaveAsProj->setShortcut(QKeySequence::SaveAs);

    menuFile->addSeparator();

    menuFile->addAction(tr("&Quit"), this, SLOT(close()), QKeySequence::Quit);

    // view


    // settings



    // help
    auto actionHelpManual = menuHelp->addAction(tr("&Check Manual"));
    auto actionAbout = menuHelp->addAction(tr("&About"));
    actionAbout->setShortcut(QKeySequence::HelpContents);
    QObject::connect(actionAbout, &QAction::triggered, [this](){
        QMessageBox::about(this, tr("About This"),
            tr("This is an experimental project.\n Author: Yang Hao \n Contact: yangh2007@gmail.com"));
    });

    
    _mdiArea = new QMdiArea;
    setCentralWidget(_mdiArea);
    _mdiArea->setBackground(Qt::lightGray);
    _workThread = new WorkThread(this);


    ////auto actionLoadPanorama = 

    //connect(ui.actionLoad_Meshes, SIGNAL(triggered()), this, SLOT(loadMesh()));
    //connect(ui.actionExport_Mesh, SIGNAL(triggered()), this, SLOT(exportMesh()));
    //connect(ui.actionNew, SIGNAL(triggered()), this, SLOT(newFile()));
    //connect(ui.actionOpen, SIGNAL(triggered()), this, SLOT(openFile()));
    //connect(ui.actionSave, SIGNAL(triggered()), this, SLOT(saveFile()));
    //connect(ui.actionSave_As, SIGNAL(triggered()), this, SLOT(saveAsFile()));
    //connect(ui.actionIntersection, SIGNAL(triggered()), this, SLOT(intersectAll()));
    //connect(ui.actionUnion, SIGNAL(triggered()), this, SLOT(uniteAll()));
    //connect(ui.actionSubtraction, SIGNAL(triggered()), this, SLOT(subtractAll()));
    //connect(ui.actionSubtractionReverse, SIGNAL(triggered()), this, SLOT(subtractReverse()));
    //connect(ui.actionNormalize, SIGNAL(triggered()), this, SLOT(normalizeGraphs()));
    //connect(ui.actionClear, SIGNAL(triggered()), this, SLOT(removeSelected()));
    //connect(ui.actionReset_View, SIGNAL(triggered()), this, SLOT(resetView()));
    //connect(ui.actionSelect_All, SIGNAL(triggered()), this, SLOT(selectAll()));
    //connect(ui.actionDeselect_All, SIGNAL(triggered()), this, SLOT(deselectAll()));
    //connect(ui.actionReverse_Orientation, SIGNAL(triggered()), this, SLOT(reverseOrientation()));
    //connect(ui.actionAbout, SIGNAL(triggered()), this, SLOT(about()));
    //connect(ui.actionTime_Recording_Enabled, SIGNAL(toggled(bool)), this, SLOT(setTimeRecordingEnabled(bool)));
    //connect(ui.actionFocus_At_Selection, SIGNAL(triggered()), this, SLOT(focusAt()));
    //connect(ui.actionToggle_Appearance, SIGNAL(triggered()), this, SLOT(toggleApearance()));

    //t_ = new WorkThread(this);

    //p_ = new QProgressBar(this);
    //p_->setFixedHeight(15);
    //p_->setFixedWidth(300);
    //statusBar()->addPermanentWidget(p_);
    //p_->hide();
    //p_->setRange(0, 0);

    //connect(t_, SIGNAL(started()), p_, SLOT(show()));
    //connect(t_, SIGNAL(finished()), p_, SLOT(hide()));
    //connect(t_, SIGNAL(finished()), this, SLOT(postWork()));

    //if (ui.toolBarBottom->actions().empty())
    //    ui.toolBarBottom->hide();

    //newFile();

    //time_recording_enabled_ = false;
    setAcceptDrops(true);
}

void MainWin::installProject(Project * proj) {
    _projects.append(proj);
    // add widgets and actions
    QList<QMdiSubWindow*> ws;
    for (QWidget * w : proj->widgets()){
        ws << _mdiArea->addSubWindow(w);
        ws.last()->hide();
    }
    _subwins << ws;
    QList<QAction*> as;
    for (QAction * a : proj->actions()){
        as << a;
    }
    _actions << as;
    _activeProject = _projects.size() - 1;
}