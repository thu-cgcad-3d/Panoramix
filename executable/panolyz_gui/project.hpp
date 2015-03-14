#ifndef PROJECT_HPP
#define PROJECT_HPP
 
#include <memory>
#include <QtCore>
#include <QtGui>

namespace panoramix {namespace vis {
    class ProjectCore;
}}


class Project : public QObject {
    Q_OBJECT
public:
    explicit Project(QObject * parent = 0);
    explicit Project(const QString & image, bool isPano, QObject * parent = 0);
    explicit Project(const QString & projFile, QObject * parent = 0);
    virtual ~Project();

    void initialize(const QString & image, bool isPano);

    void saveToDisk(const QString & filename) const;
    void loadFromDisk(const QString & filename);

    const QFileInfo & projectFileInfo() const { return _projectFileInfo; }
    const QList<QWidget*> & widgets() const { return _widgets; }
    const QList<QAction*> & actions() const { return _actions; }

    Q_SLOT void update();

private:
    QFileInfo _projectFileInfo;
    QFileInfo _imageFileInfo;
    bool _isPanorama;
    QList<QWidget*> _widgets;
    QList<QAction*> _actions;
    std::unique_ptr<panoramix::vis::ProjectCore> _core;
};

 
#endif