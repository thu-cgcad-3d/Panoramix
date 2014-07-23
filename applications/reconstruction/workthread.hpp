#pragma once

#include <functional>
#include <QtCore>



class WorkThread : public QThread {

    Q_OBJECT

public:
    WorkThread(QObject *parent = nullptr);
    ~WorkThread();

    const QList<QVector3D> & camDirections() const { return _camDirections; }
    QList<QVector3D> & camDirections() { return _camDirections; }

protected:
    void run() override;

private:
    void reconstructPanorama(const QImage & pan);

private:
    QList<QVector3D> _camDirections;
};