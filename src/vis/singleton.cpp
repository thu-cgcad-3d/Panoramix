#include "qt_glue.hpp"
#include "renderable_object.hpp"
#include "singleton.hpp"

namespace panoramix {
    namespace vis {

        static char * appName = "Qt Gui";

        QGuiApplication* Singleton::InitGui(int argc, char ** argv) {
            if (qApp)
                return qApp;
            QGuiApplication* app = new QGuiApplication(argc, argv);
            app->setQuitOnLastWindowClosed(true);
            return app;
        }

        QGuiApplication* Singleton::InitGui() {
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