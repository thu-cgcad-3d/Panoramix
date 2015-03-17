#include "../../src/core/version.hpp"
#include "../../src/vis/singleton.hpp"
#include "mainwin.hpp"

int main(int argc, char ** argv) {
    panoramix::vis::Singleton::InitGui(argc, argv);
    QApplication::setApplicationName(QObject::tr("PANOLYZ"));
    QApplication::setApplicationVersion(QString::fromStdString(panoramix::core::GetVersion().toString()));
    QApplication::setQuitOnLastWindowClosed(true);
    MainWin mwin;
    mwin.resize(900, 800);
    mwin.show();
    return panoramix::vis::Singleton::ContinueGui();
}