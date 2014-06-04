#include "qt_glue.hpp"
#include "qt_opengl_object.hpp"
#include "singleton.hpp"

inline void _InitResources() {
    Q_INIT_RESOURCE(shaders);
    Q_INIT_RESOURCE(textures);
}

namespace panoramix {
    namespace vis {

        static char * appName = "Qt Gui";

        QApplication* Singleton::InitGui(int argc, char ** argv){
            if (qApp)
                return qApp;
            QApplication* app = new QApplication(argc, argv);
            _InitResources();
            app->setQuitOnLastWindowClosed(true);
            return app;
        }

        QApplication* Singleton::InitGui(){
            return InitGui(1, &appName);
        }

        void Singleton::ContinueGui(){
            if (!qApp){
                qDebug() << "call InitGui first!";
                return;
            }
            qApp->setQuitOnLastWindowClosed(true);
            qApp->exec();
        }


    }
}