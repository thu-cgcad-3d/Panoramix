#ifndef MAINWIN_HPP
#define MAINWIN_HPP

#include <functional>
#include <QtWidgets>

class Project;
class WorkThread;

class ThreadPool : public QObject {
    Q_OBJECT

    void start(std::function<void(void)> && fun);

};

class MainWin : public QMainWindow {
	Q_OBJECT
	
public:
    MainWin(QWidget *parent = 0);
    ~MainWin(){}

    void installProject(Project * proj);

public slots:
    void switchToProject(int index);
    void updateProject(int index);

private:
    QList<Project *> _projects;
    QList<QList<QMdiSubWindow*>> _subwins;
    QList<QList<QAction*>> _actions;
    QTabWidget * _tabWidget;
    WorkThread * _workThread;
    QProgressBar * _progressBar;
};
 
#endif