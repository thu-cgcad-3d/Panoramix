#ifndef WORKTHREAD_HPP
#define WORKTHREAD_HPP

#include <functional>
#include <QtCore>


class WorkThread : public QThread
{
    Q_OBJECT

public:
    WorkThread(QObject *parent);
    ~WorkThread();

protected:
    void run() override;

private:
};

 
#endif