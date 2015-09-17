#include <QApplication>


#include "../../src/gui/singleton.hpp"

#include "mainwin_panoindoor.hpp"

namespace category = panolyz::PanoramaIndoor;

int main(int argc, char ** argv) {
    pano::gui::Singleton::InitGui(argc, argv);
    QApplication::setApplicationName(QObject::tr("PANOLYZ_LABEL"));
    QApplication::setQuitOnLastWindowClosed(true);

    QString dir = "H:\\DataSet";
    auto filenames = QFileDialog::getOpenFileNames(nullptr, QObject::tr("Select an image file"),
        dir,
        QObject::tr("Image Files (*.png;*.jpg;*.jpeg);;All Files (*.*)"));

    for (auto & fname : filenames) {
        category::MainWin mwin;
        mwin.selectFile(fname);
        mwin.resize(900, 800);
        mwin.show();
    }

    return pano::gui::Singleton::ContinueGui();
}