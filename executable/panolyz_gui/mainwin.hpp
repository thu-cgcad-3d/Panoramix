#ifndef MAINWIN_HPP
#define MAINWIN_HPP

#include <functional>
#include <QtWidgets>

class Project;
class WorkThread;
class ThreadPool;
class MainWin : public QMainWindow {
	Q_OBJECT
	
public:
    MainWin(QWidget *parent = 0);
    ~MainWin(){}

    void installProject(Project * proj);

public slots:
    void switchToProject(int index);
    void updateProject(int index, bool forceSourceStepUpdate = false);

private:
    QList<Project *> _projects;
    QList<QList<QMdiSubWindow*>> _subwins;
    QList<QList<QAction*>> _actions;
    QList<QList<QAction*>> _subwinShowActions;
    QTabWidget * _tabWidget;
    ThreadPool * _threadPool;
};
 
#endif