#ifndef PROJECT_HPP
#define PROJECT_HPP
 
#include <memory>
#include <QtCore>
#include <QtGui>

class Steps;
class Project : public QObject {    
    Q_OBJECT

public:
    explicit Project(QObject * parent = 0);
    explicit Project(const QString & projFile, QObject * parent = 0);
    virtual ~Project();

    void saveToDisk(const QString & filename) const;
    void loadFromDisk(const QString & filename);

    const QFileInfo & projectFileInfo() const { return _projectFileInfo; }
    const QList<QWidget*> & widgets() const { return _widgets; }
    const QList<QAction*> & actions() const { return _actions; }

public slots:
    void update(bool forceSourceStepUpdate = false);

public:
    static Project * createProjectFromImage(const QString & image, bool isPano, QObject * parent = 0);
    static Project * loadProjectFromDisk(const QString & filename, QObject * parent = 0);

protected:
    QFileInfo _projectFileInfo;
    QList<QWidget*> _widgets;
    QList<QAction*> _actions;
    std::unique_ptr<Steps> _steps;
};


 
#endif