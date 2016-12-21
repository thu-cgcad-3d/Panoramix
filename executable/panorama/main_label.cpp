#include "panorama_reconstruction.hpp"

int main_label(int argc, char **argv,
               const std::map<std::string, std::string> &options) {

  std::string dataDir = PANORAMIX_TEST_DATA_DIR_STR;
  if (options.count("data")) {
    dataDir = options.at("data");
  }

  gui::UI::InitGui(argc, argv);
  misc::Matlab matlab;

  std::string impath;
  gui::FileDialog::PickAnImage(dataDir, &impath);

  if (impath.empty()) {
    return 0;
  }

  auto anno = LoadOrInitializeNewLayoutAnnotation(impath);

  while (true) {
    EditLayoutAnnotation(impath, anno);
    ReconstructLayoutAnnotation(anno, matlab);
    VisualizeLayoutAnnotation(anno, 0.08);
    int selected = pano::gui::SelectFrom(
        {"Accept", "Edit Again", "Abandon"}, "Your decision?",
        "Accept the edit, or edit it again, or just abandon the edit this "
        "time?",
        0, 2);
    if (selected == 0) {
      SaveLayoutAnnotation(impath, anno);
      break;
    } else if (selected == 2) {
      break;
    }
  }

  return 0;
}