#include <QtGui>
#include <QtOpenGL>
#include <QApplication>

#include "qt_glue.hpp"
#include "renderable_object.hpp"
#include "singleton.hpp"

namespace panoramix {
    namespace vis {

        static QIcon defaultIcon;
        const QIcon & Singleton::DefaultIcon(){
            return defaultIcon;
        }

        static char * appName = "Qt Gui";

        QApplication* Singleton::InitGui(int argc, char ** argv) {
            if (qApp)
                return qApp;
            Q_INIT_RESOURCE(vis);
            QApplication* app = new QApplication(argc, argv);
            defaultIcon = QIcon(":/icons/icon.png");
            Q_ASSERT(!defaultIcon.isNull());
            QApplication::setWindowIcon(defaultIcon);
            app->setQuitOnLastWindowClosed(true);
            QGLFormat glf = QGLFormat::defaultFormat();
            qDebug("OpenGL version: %d.%d", glf.majorVersion(), glf.minorVersion());
            glf.setSampleBuffers(true);
            glf.setSamples(16);
            QGLFormat::setDefaultFormat(glf);
            return app;
        }

        QApplication* Singleton::InitGui() {
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