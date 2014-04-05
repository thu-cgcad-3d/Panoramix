#include "qt_resources.hpp"

inline void _InitResources() {
	Q_INIT_RESOURCE(shaders);
	Q_INIT_RESOURCE(textures);
}

namespace panoramix {
    namespace vis {

        QApplication* InitGui(int argc, char ** argv) {
            if (qApp)
                return qApp;
            QApplication* app = new QApplication(argc, argv);
        	_InitResources();
            app->setQuitOnLastWindowClosed(true);
        	return app;
        }

        static char * appName = "Qt Gui";
        QApplication* InitGui(){
            return InitGui(1, &appName);
        }

        void ContinueGui() {            
            if (!qApp){
                qDebug() << "call InitGui first!";
                return;
            }
            qApp->setQuitOnLastWindowClosed(true);
            qApp->exec();
        }

    }
}