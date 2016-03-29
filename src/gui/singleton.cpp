#include <QApplication>
#include <QtGui>
#include <QtOpenGL>
//#include <QtPlugin>

#include "qttools.hpp"
#include "singleton.hpp"

namespace pano {
namespace gui {

static QIcon defaultIcon;
const QIcon &Singleton::DefaultIcon() { return defaultIcon; }
static QString defaultCSS;
const QString &Singleton::DefaultCSS() { return defaultCSS; }

static int _argc = 1;
static char **_argv;
static char **_envp;

void Singleton::SetCmdArgs(int argc, char **argv, char **envp) {
  _argc = argc;
  _argv = argv;
  _envp = envp;
}

QApplication *Singleton::InitGui(int argc, char **argv) {
  if (qApp)
    return qApp;

  //Q_IMPORT_PLUGIN(QWindowsIntegrationPlugin);
  Q_INIT_RESOURCE(gui);
  _argc = argc;
  _argv = argv;
  QApplication *app = new QApplication(argc, argv);

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

QApplication *Singleton::InitGui() { return InitGui(_argc, _argv); }

int Singleton::ContinueGui() {
  if (!qApp) {
    qDebug() << "call InitGui first!";
    return 0;
  }
  qApp->setQuitOnLastWindowClosed(true);
  return qApp->exec();
}
}
}