#include <QApplication>
#include <QFileDialog>
#include <QtCore>

#include "../../src/gui/singleton.hpp"
#include "../../src/gui/utility.hpp"

#include "../../src/experimental/pi_graph_annotation.hpp"
#include "../../src/experimental/pi_graph_solve.hpp"
#include "../../src/experimental/pi_graph_vis.hpp"


int main(int argc, char ** argv) {
    
    pano::gui::Singleton::InitGui(argc, argv);

    QApplication::setApplicationName(QObject::tr("PANOLYZ_LABEL"));
    QApplication::setQuitOnLastWindowClosed(true);

    QString dir = "H:\\DataSet\\PanoContext";
    auto filenames = QFileDialog::getOpenFileNames(nullptr, QObject::tr("Select an image file"),
        dir,
        QObject::tr("Image Files (*.jpg;*.jpeg);;All Files (*.*)"));

    pano::misc::Matlab matlab;

    for (auto & fname : filenames) {
        auto anno = pano::experimental::LoadOrInitializeNewLayoutAnnotation(fname.toStdString());
        while (true) {
            pano::experimental::EditLayoutAnnotation(fname.toStdString(), anno);
            pano::experimental::ReconstructLayoutAnnotation3(anno, matlab);
            pano::experimental::VisualizeLayoutAnnotation(anno);
            int selected = pano::gui::SelectFrom({ "Accept", "Edit Again", "Abandon" }, 
                "Your decision?", 
                "Accept the edit, or edit it again, or just abandon the edit this time?", 0, 2);
            if (selected == 0) {
                pano::experimental::SaveLayoutAnnotation(fname.toStdString(), anno);
                break;
            } else if (selected == 2) {
                break;
            }
        }
    }

}