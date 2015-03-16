#include "widgets.hpp"
#include "workthread.hpp"
#include "project.hpp"
#include "mainwin.hpp"

MainWin::MainWin(QWidget *parent) : QMainWindow(parent) {
   
    //// menus
    auto menuFile = menuBar()->addMenu(tr("&File"));
    auto menuView = menuBar()->addMenu(tr("&View"));
    auto menuSettings = menuBar()->addMenu(tr("&Settings"));
    auto menuTools = menuBar()->addMenu(tr("&Tools"));
    auto menuHelp = menuBar()->addMenu(tr("Help"));

    //// actions
    // file
    auto actionNewProjPanorama = menuFile->addAction(tr("Create A &New Panolyz Project from Panorama ..."));
    actionNewProjPanorama->setShortcut(QKeySequence::New);
    connect(actionNewProjPanorama, &QAction::triggered, [this](){
        QStringList filenames = QFileDialog::getOpenFileNames(this, tr("Select a Panoramic Image File"), 
            tr(PROJECT_TEST_DATA_DIR_STR),
            tr("Image Files (*.png;*.jpg;*.jpeg);;All Files (*.*)"));
        for (auto & filename : filenames)
            installProject(new Project(filename, true, this));
    });
    auto actionNewProjNormal = menuFile->addAction(tr("Create A &New Panolyz Project from Normal Photo ..."));
    connect(actionNewProjNormal, &QAction::triggered, [this](){
        QStringList filenames = QFileDialog::getOpenFileNames(this, tr("Select a Normal Image File"),
            tr(PROJECT_TEST_DATA_DIR_STR),
            tr("Image Files (*.png;*.jpg;*.jpeg);;All Files (*.*)"));
        for (auto & filename : filenames)
            installProject(new Project(filename, false, this));
    });
    auto actionOpenProj = menuFile->addAction(tr("&Open An Existing Panolyz Project ..."));
    actionOpenProj->setShortcut(QKeySequence::Open);
    connect(actionOpenProj, &QAction::triggered, [this](){
        QStringList filenames = QFileDialog::getOpenFileNames(this, tr("Select Panolyz Project Files"), 
            tr(PROJECT_TEST_DATA_DIR_STR),
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


    // tools


    // help
    auto actionHelpManual = menuHelp->addAction(tr("&Check Manual"));
    auto actionAbout = menuHelp->addAction(tr("&About"));
    actionAbout->setShortcut(QKeySequence::HelpContents);
    QObject::connect(actionAbout, &QAction::triggered, [this](){
        QMessageBox::about(this, tr("About This"),
            tr("This is an experimental project.\n Author: Yang Hao \n Contact: yangh2007@gmail.com"));
    });

    
    _tabWidget = new QTabWidget;
    _tabWidget->setTabsClosable(true);
    _tabWidget->setTabShape(QTabWidget::TabShape::Rounded);
    QObject::connect(_tabWidget, &QTabWidget::tabCloseRequested, [this](int index){
        _subwins.removeAt(index);
        _actions.removeAt(index);
        _projects[index]->deleteLater();
        _projects.removeAt(index);
        _tabWidget->removeTab(index);
    });
    QObject::connect(_tabWidget, &QTabWidget::currentChanged, [this](int index){
        switchToProject(index); 
    });
    _tabWidget->setTabPosition(QTabWidget::TabPosition::West);
    setCentralWidget(_tabWidget);
    _workThread = new WorkThread(this);

    _progressBar = new QProgressBar(this);
    _progressBar->setFixedHeight(15);
    _progressBar->setFixedWidth(300);
    statusBar()->addPermanentWidget(_progressBar);
    _progressBar->hide();
    _progressBar->setRange(0, 0);

    connect(_workThread, SIGNAL(started()), _progressBar, SLOT(show()));
    connect(_workThread, SIGNAL(finished()), _progressBar, SLOT(hide()));
    //connect(_workThread, SIGNAL(finished()), this, SLOT(postWork()));

    setAcceptDrops(true);
}

void MainWin::installProject(Project * proj) {
    _projects.append(proj);
    QMdiArea * mdiArea = new QMdiArea;
    int tabId = _tabWidget->addTab(mdiArea, proj->projectFileInfo().fileName());
    // add widgets and actions
    QList<QMdiSubWindow*> ws;
    for (QWidget * w : proj->widgets()){
        ws << mdiArea->addSubWindow(w);
        ws.last()->hide();
    }
    _subwins << ws;
    QList<QAction*> as;
    for (QAction * a : proj->actions()){
        as << a;
    }
    _actions << as;
    Q_ASSERT(tabId == _projects.size() - 1);
    Q_ASSERT(tabId == _subwins.size() - 1);
    Q_ASSERT(tabId == _actions.size() - 1);
}

void MainWin::switchToProject(int index){
    qDebug() << "switch to " << index;
}