#include <QtGui>
#include <QtOpenGL>
#include <QApplication>

#include "qt_glue.hpp"
#include "singleton.hpp"

namespace panoramix {
    namespace gui {

        static QIcon defaultIcon;
        const QIcon & Singleton::DefaultIcon(){
            return defaultIcon;
        }
        static QString defaultCSS;
        const QString & Singleton::DefaultCSS(){
            return defaultCSS;
        }

        static char * appName = "Qt Gui";

        QApplication* Singleton::InitGui(int argc, char ** argv) {
            if (qApp)
                return qApp;
            Q_INIT_RESOURCE(gui);
            QApplication* app = new QApplication(argc, argv);
            
            defaultIcon = QIcon(":/icons/icon.png");
            Q_ASSERT(!defaultIcon.isNull());
            QFile file(":/css/gui_win.css");
            bool opened = file.open(QFile::ReadOnly);
            Q_ASSERT(opened);
            defaultCSS = QTextStream(&file).readAll();
            QApplication::setWindowIcon(defaultIcon);
            app->setStyleSheet(defaultCSS);

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

        int Singleton::ContinueGui(){
            if (!qApp){
                qDebug() << "call InitGui first!";
                return 0;
            }
            qApp->setQuitOnLastWindowClosed(true);            
            return qApp->exec();
        }


    }
}