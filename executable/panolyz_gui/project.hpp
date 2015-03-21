#ifndef PROJECT_HPP
#define PROJECT_HPP
 
#include <memory>
#include <QtCore>
#include <QtGui>

class StepsDAG;
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

    const QVariantMap & configurations() const { return _configurations; }
    QVariant conf(const QString & name) const { return _configurations.value(name); }
    void setConf(const QString & name, const QVariant & v) { _configurations.insert(name, v); }

signals:
    void messageUpdated(QString);

public slots:
    void update(bool forceSourceStepUpdate = false);

public:
    static Project * createProjectFromImage(const QString & image, bool isPano, QObject * parent = 0);
    static Project * loadProjectFromDisk(const QString & filename, QObject * parent = 0);

protected:
    QFileInfo _projectFileInfo;
    QList<QWidget*> _widgets;
    QList<QAction*> _actions;
    StepsDAG * _steps;
    QVariantMap _configurations;
    QReadWriteLock _confLock;
};


 
#endif