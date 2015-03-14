#include "../../src/core/version.hpp"
#include "../../src/vis/singleton.hpp"
#include "mainwin.hpp"
#include "widgets.hpp"
#include "workthread.hpp"

int main(int argc, char ** argv) {
    panoramix::vis::Singleton::InitGui(argc, argv);
    QApplication::setApplicationName(QObject::tr("PANOLYZ"));
    QApplication::setApplicationVersion(QString::fromStdString(panoramix::core::GetVersion().toString()));
    MainWin mwin;
    mwin.resize(900, 800);
    mwin.show();
    return panoramix::vis::Singleton::ContinueGui();
}