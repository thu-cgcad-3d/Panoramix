#include <QApplication>
#include <QFileDialog>
#include <QtCore>

#include "../../src/gui/singleton.hpp"
#include "../../src/experimental/pi_graph_annotation.hpp"


int main(int argc, char ** argv) {
    
    pano::gui::Singleton::InitGui(argc, argv);

    QApplication::setApplicationName(QObject::tr("PANOLYZ_LABEL"));
    QApplication::setQuitOnLastWindowClosed(true);

    QString dir = "H:\\DataSet\\pi";
    auto filenames = QFileDialog::getOpenFileNames(nullptr, QObject::tr("Select an image file"),
        dir,
        QObject::tr("Image Files (*.png;*.jpg;*.jpeg);;All Files (*.*)"));

    for (auto & fname : filenames) {
        auto anno = pano::experimental::LoadOrInitializeNewAnnotation(fname.toStdString());
        pano::experimental::EditAnnotation(anno);
        pano::experimental::SaveAnnotation(fname.toStdString(), anno);
    }
}