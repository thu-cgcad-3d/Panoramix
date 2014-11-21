
#include <QtGui>
#include "pmainwin.h"


int main(int argc, char ** argv) {
    Q_INIT_RESOURCE(gui);
	QApplication app(argc, argv);
	PMainWin win;
	win.resize(800, 600);
	win.show();
	return app.exec();
}