#include <QApplication>
#include <QFileDialog>
#include <QtCore>

#include "../../src/gui/singleton.hpp"

#include "../../src/experimental/pi_graph_annotation.hpp"
#include "../../src/experimental/pi_graph_solve.hpp"
#include "../../src/experimental/pi_graph_vis.hpp"


int main(int argc, char ** argv) {
    
    pano::gui::Singleton::InitGui(argc, argv);

    QApplication::setApplicationName(QObject::tr("PANOLYZ_LABEL"));
    QApplication::setQuitOnLastWindowClosed(true);

    QString dir = "H:\\DataSet\\pi";
    auto filenames = QFileDialog::getOpenFileNames(nullptr, QObject::tr("Select an image file"),
        dir,
        QObject::tr("Image Files (*.png;*.jpg;*.jpeg);;All Files (*.*)"));

    for (auto & fname : filenames) {
        auto anno = pano::experimental::LoadOrInitializeNewLayoutAnnotation(fname.toStdString());
        pano::experimental::EditLayoutAnnotation(anno);
        pano::experimental::ReconstructLayoutAnnotation(anno);
        pano::experimental::VisualizeLayoutAnnotation(anno);
        pano::experimental::SaveLayoutAnnotation(fname.toStdString(), anno);
    }

}