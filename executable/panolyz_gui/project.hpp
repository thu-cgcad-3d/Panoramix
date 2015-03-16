#ifndef PROJECT_HPP
#define PROJECT_HPP
 
#include <memory>
#include <QtCore>
#include <QtGui>

namespace panoramix {namespace vis {
    class Steps;
}}


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

    Q_SLOT void update();

    static Project * createProject(const QString & image, bool isPano, QObject * parent = 0);

protected:
    QFileInfo _projectFileInfo;
    QList<QWidget*> _widgets;
    QList<QAction*> _actions;
    std::unique_ptr<panoramix::vis::Steps> _steps;
};


 
#endif