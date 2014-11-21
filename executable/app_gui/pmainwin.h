#pragma once

#include <QtGui>
#include <QtWidgets>

class PMainWin : public QMainWindow {

	Q_OBJECT

public:
    explicit PMainWin(QWidget * parent = 0);
	~ PMainWin();

private:
    QMdiArea * _mdiArea;

};