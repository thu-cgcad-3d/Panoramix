#include "project.hpp"
#include "mainwin.hpp"


class SimpleThread : public QThread {
public:
    explicit SimpleThread(std::function<void(void)> && fun, QObject * parent = nullptr) 
        : QThread(parent), _fun(std::move(fun)) {}
protected:
    virtual void run() override { _fun(); }
private:
    std::function<void(void)> _fun;
};

class ThreadPool : public QObject {
public:
    explicit ThreadPool(int maxThreads, QObject * parent = nullptr) 
        : QObject(parent), _maxThreads(maxThreads) {}
    
    QPair<QProgressBar *, SimpleThread *> attach(std::function<void(void)> && fun, QWidget * barParent = nullptr){
        SimpleThread * t = new SimpleThread(std::move(fun), this);
        QProgressBar * bar = new QProgressBar(barParent);
        bar->hide();
        connect(t, SIGNAL(started()), bar, SLOT(show()));
        connect(t, SIGNAL(finished()), bar, SLOT(hide()));
        connect(t, SIGNAL(finished()), t, SLOT(deleteLater()));
        connect(t, SIGNAL(finished()), bar, SLOT(deleteLater()));
        return qMakePair(bar, t);
    }

private:
    int _maxThreads;
};



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
            installProject(Project::createProjectFromImage(filename, true, this));
    });
    auto actionNewProjNormal = menuFile->addAction(tr("Create A &New Panolyz Project from Normal Photo ..."));
    connect(actionNewProjNormal, &QAction::triggered, [this](){
        QStringList filenames = QFileDialog::getOpenFileNames(this, tr("Select a Normal Image File"),
            tr(PROJECT_TEST_DATA_DIR_STR),
            tr("Image Files (*.png;*.jpg;*.jpeg);;All Files (*.*)"));
        for (auto & filename : filenames)
            installProject(Project::createProjectFromImage(filename, false, this));
    });
    auto actionOpenProj = menuFile->addAction(tr("&Open An Existing Panolyz Project ..."));
    actionOpenProj->setShortcut(QKeySequence::Open);
    connect(actionOpenProj, &QAction::triggered, [this](){
        QStringList filenames = QFileDialog::getOpenFileNames(this, tr("Select Panolyz Project Files"), 
            tr(PROJECT_TEST_DATA_DIR_STR),
            tr("Panolyz Project File (*.panoproj)"));
        for (auto & filename : filenames)
            installProject(Project::loadProjectFromDisk(filename, this));
    });

    menuFile->addSeparator();
    
    auto actionSaveProj = menuFile->addAction(tr("&Save Current Panolyz Project"));
    actionSaveProj->setShortcut(QKeySequence::Save);
    auto actionSaveAsProj = menuFile->addAction(tr("&Save Current Panolyz Project As ..."));
    actionSaveAsProj->setShortcut(QKeySequence::SaveAs);

    menuFile->addSeparator();

    menuFile->addAction(tr("&Quit"), this, SLOT(close()), QKeySequence::Quit);

    // view
    auto actionViewCascade = menuView->addAction(tr("Cascade Views"));
    connect(actionViewCascade, &QAction::triggered, [this](){
        if (_tabWidget->currentIndex() < 0)
            return;
        qobject_cast<QMdiArea*>(_tabWidget->currentWidget())->cascadeSubWindows();
    });
    auto actionViewTile = menuView->addAction(tr("Tile Views"));
    connect(actionViewTile, &QAction::triggered, [this](){
        if (_tabWidget->currentIndex() < 0)
            return;
        qobject_cast<QMdiArea*>(_tabWidget->currentWidget())->tileSubWindows();
    });


    // settings


    // tools
    auto actionRerunAll = menuTools->addAction(tr("Rerun All StepsDAG"));
    QObject::connect(actionRerunAll, &QAction::triggered, [this](){
        int index = _tabWidget->currentIndex();
        if (index < 0)
            return;
        updateProject(index, true);
    });
    auto actionRerunThoseNeedUpdate = menuTools->addAction(tr("Rerun Those Need Update"));
    QObject::connect(actionRerunThoseNeedUpdate, &QAction::triggered, [this](){
        int index = _tabWidget->currentIndex();
        if (index < 0)
            return;
        updateProject(index, false);
    });


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


    _threadPool = new ThreadPool(4, this);

    setAcceptDrops(true);
}

void MainWin::installProject(Project * proj) {
    _projects.append(proj);
    QMdiArea * mdiArea = new QMdiArea;
    int tabId = _tabWidget->addTab(mdiArea, proj->projectFileInfo().fileName());
    // add widgets and actions
    QList<QMdiSubWindow*> ws;
    for (int i = 0; i < proj->widgets().size(); i++){
        auto w = proj->widgets().at(i);
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
    //updateProject(tabId, true);
}

void MainWin::switchToProject(int index){
    qDebug() << "switch to " << index;

}

void MainWin::updateProject(int index, bool forceSourceStepUpdate) {
    if (index < 0)
        return;
    auto tb = _threadPool->attach([this, index, forceSourceStepUpdate](){
        _projects[index]->update(forceSourceStepUpdate);
    });
    connect(_tabWidget->widget(index), SIGNAL(destroyed()), tb.second, SLOT(terminate()));
    auto b = tb.first;
    b->setFixedHeight(15);
    b->setFixedWidth(300);
    statusBar()->addPermanentWidget(b);
    b->hide();
    b->setRange(0, 0);
    tb.second->start();
}