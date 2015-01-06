#pragma once

#include <QtWidgets>

class PProjectWidget : public QWidget {

    Q_OBJECT

public:
    explicit PProjectWidget(QWidget * parent = 0);
    ~PProjectWidget();

private:
    QTabWidget * _tabWidget;

};