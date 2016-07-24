#include "pch.hpp"

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

  auto p = qgetenv("QT_QPA_PLATFORM_PLUGIN_PATH");
  std::cout << "QT_QPA_PLATFORM_PLUGIN_PATH=" << p.toStdString() << std::endl;

  //Q_IMPORT_PLUGIN(QWindowsIntegrationPlugin);
  Q_INIT_RESOURCE(gui);
  //Q_IMPORT_PLUGIN(QWindowsIntegrationPlugin)
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
	QSurfaceFormat sf = QSurfaceFormat::defaultFormat();
	sf.setSamples(16);
	 qDebug("OpenGL version: %d.%d", sf.majorVersion(), sf.minorVersion());
	QSurfaceFormat::setDefaultFormat(sf);

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