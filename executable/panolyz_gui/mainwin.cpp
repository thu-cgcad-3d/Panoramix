#include "../../src/gui/singleton.hpp"

#include "project.hpp"
#include "mainwin.hpp"

using namespace panoramix;

class WorkThread : public QThread {
public:
    explicit WorkThread(std::function<void(void)> && fun, QObject * parent = nullptr) 
        : QThread(parent), _fun(std::move(fun)) {}
protected:
    virtual void run() override { _fun(); }
private:
    std::function<void(void)> _fun;
};

class WorkProgressBar : public QProgressBar {
    Q_OBJECT
public:
    explicit WorkProgressBar(QWidget * parent = nullptr) : QProgressBar(parent) {}
    Q_SLOT void setText(const QString & t) { _text = t; }
    virtual QString text() const override { return fontMetrics().elidedText(_text, Qt::TextElideMode::ElideRight, width() / 3); }
private:
    QString _text;
};

class ThreadPool : public QObject {
public:
    explicit ThreadPool(int maxThreads, QObject * parent = nullptr) 
        : QObject(parent), _maxThreads(maxThreads) {}
    
    QPair<WorkProgressBar *, WorkThread *> attach(std::function<void(void)> && fun, QWidget * barParent = nullptr){
        WorkThread * t = new WorkThread(std::move(fun), this);
        WorkProgressBar * bar = new WorkProgressBar(barParent);
        bar->hide();
        connect(t, SIGNAL(started()), bar, SLOT(show()), Qt::QueuedConnection);
        connect(t, SIGNAL(finished()), bar, SLOT(hide()), Qt::QueuedConnection);
        connect(t, SIGNAL(finished()), t, SLOT(deleteLater()), Qt::QueuedConnection);
        connect(t, SIGNAL(finished()), bar, SLOT(deleteLater()), Qt::QueuedConnection);
        return qMakePair(bar, t);
    }

private:
    int _maxThreads;
};

class SubWindow : public QMdiSubWindow {
    Q_OBJECT
public:
    explicit SubWindow(QWidget * parent = nullptr) : QMdiSubWindow(parent, Qt::SubWindow | Qt::Tool) {}
    Q_SIGNAL void visibilityChanged(bool);
protected:
    virtual void showEvent(QShowEvent * event) override {
        emit visibilityChanged(true);
        QMdiSubWindow::showEvent(event);
    }
    virtual void hideEvent(QHideEvent * event) override{
        emit visibilityChanged(false);
        QMdiSubWindow::hideEvent(event);
    }
};

class MdiArea : public QMdiArea {
    Q_OBJECT
public:
    explicit MdiArea(QWidget * parent = nullptr) : QMdiArea(parent) {}
    Q_SIGNAL void closed();
protected:
    virtual void closeEvent(QCloseEvent * event) override {
        emit closed();
        QMdiArea::closeEvent(event);
    }
};

class ConfigurationModel : public QAbstractItemModel {
public:
    ConfigurationModel(Project * p = nullptr, QObject * parent = nullptr) : QAbstractItemModel(parent), _project(p) {}
    
    virtual QVariant data(const QModelIndex & index, int role) const override {
        if (!index.isValid() || !_project)
            return QVariant();
        auto names = _project->configurations().keys();
        if (role == Qt::DisplayRole || role == Qt::EditRole){
            int col = index.column();
            int row = index.row();
            return col == 0 ? names.at(row) : _project->conf(names.at(row));
        }
        else if (role == Qt::DecorationRole && index.column() == 0){
            return gui::Singleton::DefaultIcon();
        }
        else if (role == Qt::BackgroundRole){
            //return QBrush(Qt::darkGray);
        }
        else if (role == Qt::ForegroundRole){
            //return QBrush(Qt::white);
        }
        return QVariant();
    }
    
    bool setData(const QModelIndex &index, const QVariant &value, int role) {
        if (!_project)
            return false;
        if (role == Qt::EditRole && index.isValid() && index.column() == 1){
            auto names = _project->configurations().keys();
            _project->setConf(names.at(index.row()), value);
            emit dataChanged(index, index);
            return true;
        }
        return false;
    }

    Qt::ItemFlags flags(const QModelIndex &index) const{
        if (!index.isValid() || !_project)
            return 0;
        auto names = _project->configurations().keys();
        if (index.column() == 0)
            return (Qt::ItemIsEnabled | Qt::ItemIsSelectable);
        return (Qt::ItemIsEnabled | Qt::ItemIsSelectable | Qt::ItemIsEditable);
    }

    QVariant headerData(int section, Qt::Orientation orientation,
        int role = Qt::DisplayRole) const{
        if (orientation == Qt::Horizontal && role == Qt::DisplayRole)
            return section == 0 ? tr("Name") : tr("Value");
        return QVariant();
    }

    QModelIndex index(int row, int column,
        const QModelIndex &parent = QModelIndex()) const{
        if (!hasIndex(row, column, parent))
            return QModelIndex();
        if (parent.isValid())
            return QModelIndex();
        return createIndex(row, column);
    }

    QModelIndex parent(const QModelIndex &index) const{
        return QModelIndex();
    }

    int rowCount(const QModelIndex &parent = QModelIndex()) const {
        if (!parent.isValid() && _project)
            return _project->configurations().size();
        return 0;
    }

    int columnCount(const QModelIndex &parent = QModelIndex()) const {
        return 2;
    }

private:
    Project * _project;
};



#include "mainwin.moc"



static void CustomItemEditorFactory() {
    struct DoubleEditorCreator : QItemEditorCreatorBase {
        virtual QWidget *createWidget(QWidget *parent) const {
            QDoubleSpinBox * sb = new QDoubleSpinBox(parent);
            sb->setRange(0.0, 100.0);
            sb->setDecimals(10);
            sb->setSingleStep(0.005);
            return sb;
        }
        virtual QByteArray valuePropertyName() const {
            return "value";
        }
    };

    QItemEditorFactory * f = new QItemEditorFactory(*QItemEditorFactory::defaultFactory());
    f->registerEditor(QVariant::Double, new DoubleEditorCreator);
    QItemEditorFactory::setDefaultFactory(f);
}



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
            tr("F:\\DataSets\\YorkUrbanDB\\data"),
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
    auto actionViewCascade = new QAction(tr("Cascade Subwindows"), this);
    menuView->addAction(actionViewCascade);
    connect(actionViewCascade, &QAction::triggered, [this](){
        if (_tabWidget->currentIndex() < 0)
            return;
        qobject_cast<QMdiArea*>(_tabWidget->currentWidget())->cascadeSubWindows();
    });
    auto actionViewTile = new QAction(tr("Tile Subwindows"), this);
    menuView->addAction(actionViewTile);
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
        _subwinShowActions.removeAt(index);
        _actions.removeAt(index);
        _confModels.removeAt(index);
        _tabWidget->removeTab(index);
        _projects[index]->deleteLater();
        _projects.removeAt(index);
        qDebug() << "removed: " << index;
    });

    _confDock = new QDockWidget(tr("Configurations"), this);
    QTreeView * vw = new QTreeView;
    vw->setAlternatingRowColors(true);
    vw->setAutoScroll(true);
    vw->setIndentation(0);
    _confDock->setWidget(vw);
    addDockWidget(Qt::DockWidgetArea::RightDockWidgetArea, _confDock);
    _confDock->hide();
    // item editor factory
    CustomItemEditorFactory();
    QObject::connect(_tabWidget, &QTabWidget::currentChanged, [this, menuView, actionViewCascade, actionViewTile, vw](int index){
        switchToProject(index); 
        // set menu actions
        menuView->clear();
        menuView->addAction(actionViewCascade);
        menuView->addAction(actionViewTile);
        if (index != -1){
            menuView->addSeparator();
            menuView->addActions(_subwinShowActions[index]);
            vw->setModel(_confModels.at(index));
        }
        else{
            _confDock->hide();
        }
    });
    _tabWidget->setTabPosition(QTabWidget::TabPosition::West);
    setCentralWidget(_tabWidget);


    _threadPool = new ThreadPool(4, this);

    setAcceptDrops(true);
}



void MainWin::installProject(Project * proj) {
    _projects.append(proj);
    MdiArea * mdiArea = new MdiArea;
    QLinearGradient gradient(0, 0, 800.0, 1000.0);
    gradient.setColorAt(0, QColor(80, 80, 80));
    gradient.setColorAt(1, Qt::black);
    mdiArea->setBackground(Qt::white);

    // add widgets and actions
    QList<QMdiSubWindow*> ws;
    QList<QAction*> wShowActions;
    for (int i = 0; i < proj->widgets().size(); i++){
        auto w = proj->widgets().at(i);
        SubWindow * subwin = new SubWindow;
        subwin->setWidget(w);
        mdiArea->addSubWindow(subwin);
        //auto subwin = mdiArea->addSubWindow(w);
        subwin->hide();
        subwin->setStyleSheet(
            "QMdiSubWindow { border: 5px solid darkGray; background: gray }"
            "QMdiSubWindow:title{ background: white; }"
        );
        QAction * a = new QAction(subwin->windowTitle(), subwin);
        connect(a, SIGNAL(triggered(bool)), subwin, SLOT(setVisible(bool)));
        connect(subwin, SIGNAL(visibilityChanged(bool)), a, SLOT(setChecked(bool)));
        a->setCheckable(true);
        a->setChecked(false);
        ws << subwin;
        wShowActions << a;
    }
    _subwins << ws;
    _subwinShowActions << wShowActions;

    QList<QAction*> as;
    for (QAction * a : proj->actions()){
        as << a;
    }
    _actions << as;

    _confModels << (new ConfigurationModel(proj, proj));

    int tabId = _tabWidget->addTab(mdiArea, proj->projectFileInfo().fileName());
    Q_ASSERT(tabId == _projects.size() - 1);
    Q_ASSERT(tabId == _subwins.size() - 1);
    Q_ASSERT(tabId == _actions.size() - 1);
    _tabWidget->setCurrentIndex(tabId);
    _confDock->show();

    updateProject(tabId, true);
}

void MainWin::switchToProject(int index){
    qDebug() << "switch to " << index;

}

void MainWin::updateProject(int index, bool forceSourceStepUpdate) {
    if (index < 0)
        return;
    QPair<WorkProgressBar*, WorkThread*> tb = _threadPool->attach([this, index, forceSourceStepUpdate](){
        _projects[index]->update(forceSourceStepUpdate);
    });
    connect(_projects[index], SIGNAL(messageUpdated(QString)), tb.first, SLOT(setText(QString)));
    connect(_tabWidget->widget(index), SIGNAL(closed()), tb.second, SLOT(terminate()));
    auto b = tb.first;
    b->setFixedHeight(15);
    b->setFixedWidth(300);
    statusBar()->addPermanentWidget(b);
    b->hide();
    b->setRange(0, 0);
    tb.second->start();
}