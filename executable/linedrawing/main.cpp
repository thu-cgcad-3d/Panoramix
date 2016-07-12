#include "cache.hpp"
#include "singleton.hpp"

#include "line_drawing_reconstruction.hpp"
#include "line_drawing_tools.hpp"

using namespace pano;
using namespace pano::core;
using namespace pano::experimental;

int main(int argc, char **argv, char **env) {
  gui::Singleton::SetCmdArgs(argc, argv, env);
  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath(PANORAMIX_CACHE_DATA_DIR_STR "\\LineDrawing\\");

  enum TaskName { EditImageAnnotation, TrainEnergyWeights, RunReconstruction };

  TaskName task = TrainEnergyWeights;

  if (task == EditImageAnnotation) { // edit annotation on image
    std::string im_file_name;
    Image3ub im = gui::FileDialog::PickAnImage(PANORAMIX_TEST_DATA_DIR_STR,
                                               &im_file_name);
    ResizeToMakeHeightUnder(im, 600);
    std::string folder = misc::FolderOfFile(im_file_name);
    std::string name = misc::NameOfFile(im_file_name);
    std::string anno_file = folder + "\\" + name + ".ldanno.cereal";
    while (true) {
      LineDrawingWidget w(im);
      LoadFromDisk(anno_file, w.projection);
      w.projection.image = im;
      w.show();

      gui::Singleton::ContinueGui();
      int choice = gui::SelectFrom({"save", "abandon"}, "make your choice",
                                   "your choice?", 0, 1);
      if (choice == 0) {
        SaveToDisk(anno_file, w.projection);
      }
    }
  } else if (task == TrainEnergyWeights) {

    std::vector<std::string> model_names = {
        "hex",       "triangle", "twotriangles", "stool", "plane",
        "towers",    "tower",    "towerx",       "car",   "gate",
        "smallgate", "castle",   "desk"};

    for (auto &model_name : model_names) {
      std::string folder = PANORAMIX_TEST_DATA_DIR_STR "\\linedrawing\\objs\\" +
                           model_name + "\\";
      auto line_drawing =
          LineDrawing3FromObjFile(folder + model_name + "_w_intf.obj");
      // randomly generate cameras

      // todo
    }

  } else if (task == RunReconstruction) {
    LineDrawingReconstructionSteps steps;
    steps.preprocess.rerun = false;
    steps.find_vps.rerun = steps.preprocess.rerun || false;
    steps.estimate_orientations.rerun =
        steps.find_vps.rerun || false;
    steps.reconstruction.rerun =
        steps.estimate_orientations.rerun || false;

    std::vector<LineDrawingInput> inputs = {
        LineDrawingInput::FromObjFile("hex", "cam1"),          //
        LineDrawingInput::FromObjFile("triangle", "cam1"),     //
        LineDrawingInput::FromObjFile("twotriangles", "cam1"), //
        LineDrawingInput::FromObjFile("stool", "cam1"),        //
        LineDrawingInput::FromObjFile("plane", "cam1"),        //
        LineDrawingInput::FromObjFile("towers", "cam1"),       // skewed a bit
        LineDrawingInput::FromObjFile("tower", "cam1"),        //
        LineDrawingInput::FromObjFile("towerx", "cam1"),       //
        LineDrawingInput::FromObjFile("car", "cam1"),          //
        LineDrawingInput::FromObjFile("gate", "cam1"),         // skewed a bit
        LineDrawingInput::FromObjFile("smallgate", "cam1"),    //
        LineDrawingInput::FromObjFile("castle", "cam1"),       //
        LineDrawingInput::FromObjFile("desk", "cam1")          //
    };

    for (auto &input : inputs) {
      RunLineDrawingReconstruction(input, steps, GUILevelNecessary);
    }
  }

  return 0;
}