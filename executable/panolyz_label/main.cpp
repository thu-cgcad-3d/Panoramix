#include <QApplication>
#include <QFileDialog>
#include <QtCore>

#include "../../src/gui/singleton.hpp"
#include "../../src/gui/utility.hpp"

#include "../../src/experimental/pi_graph_annotation.hpp"
#include "../../src/experimental/pi_graph_solve.hpp"
#include "../../src/experimental/pi_graph_vis.hpp"

using namespace pano;
using namespace pano::core;
using namespace pano::experimental;

int main(int argc, char **argv) {

  gui::Singleton::InitGui();
  misc::Matlab matlab;

  std::vector<std::string> impaths;
  //gui::PickImages("H:\\DataSet\\pi\\dataset\\selected\\", &impaths);
  
  impaths.push_back("");
  gui::PickAnImage("F:\\CVPR2016", &(impaths[0]));

  for (auto &impath : impaths) {
    auto anno = pano::experimental::LoadOrInitializeNewLayoutAnnotation(impath);
    
    while (true) {
      pano::experimental::EditLayoutAnnotation(impath, anno);
      pano::experimental::ReconstructLayoutAnnotation(anno, matlab);
      pano::experimental::VisualizeLayoutAnnotation(anno, 0.08);
      int selected = pano::gui::SelectFrom(
          {"Accept", "Edit Again", "Abandon"}, "Your decision?",
          "Accept the edit, or edit it again, or just abandon the edit this "
          "time?",
          0, 2);
      if (selected == 0) {
        pano::experimental::SaveLayoutAnnotation(impath, anno);
        break;
      } else if (selected == 2) {
        break;
      }
    }
  }

  return 0;
}