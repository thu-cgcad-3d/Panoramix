#ifndef MAINWIN_HPP
#define MAINWIN_HPP
 
#include <QtWidgets>

class Project;
class WorkThread;
class MainWin : public QMainWindow {
	Q_OBJECT
	
public:
    MainWin(QWidget *parent = 0);
    ~MainWin(){}

    void installProject(Project * proj);

private:
    QMdiArea * _mdiArea;
    QList<Project *> _projects;
    QList<QList<QMdiSubWindow*>> _subwins;
    QList<QList<QAction*>> _actions;
    int _activeProject;
    WorkThread * _workThread;
};
 
#endif