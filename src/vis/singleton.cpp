#include <QtGui>
#include <QtOpenGL>
#include <QApplication>

#include "qt_glue.hpp"
#include "renderable_object.hpp"
#include "singleton.hpp"

namespace panoramix {
    namespace vis {

        static Singleton::Configuration defaultConfiguration;
        const Singleton::Configuration & Singleton::DefaultConfiguration(){
            return defaultConfiguration;
        }

        static char * appName = "Qt Gui";

        QApplication* Singleton::InitGui(int argc, char ** argv) {
            if (qApp)
                return qApp;
            Q_INIT_RESOURCE(vis);
            QApplication* app = new QApplication(argc, argv);
            
            defaultConfiguration.icon = QIcon(":/icons/icon.png");
            Q_ASSERT(!defaultConfiguration.icon.isNull());
            QFile file(":/css/vis_win.css");
            bool opened = file.open(QFile::ReadOnly);
            Q_ASSERT(opened);
            defaultConfiguration.css = QTextStream(&file).readAll();
            QApplication::setWindowIcon(defaultConfiguration.icon);
            app->setStyleSheet(defaultConfiguration.css);

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