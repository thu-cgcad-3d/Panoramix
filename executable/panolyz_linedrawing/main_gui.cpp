#include "line_drawing_widget.hpp"
#include "cache.hpp"
#include "cameras.hpp"
#include "canvas.hpp"
#include "line_drawing.hpp"
#include "scene.hpp"
#include "singleton.hpp"
#include "optimization.hpp"
#include "parallel.hpp"
#include "clock.hpp"
#include "tools.hpp"

using namespace pano;
using namespace pano::core;
using namespace pano::experimental;

int main(int argc, char **argv, char **env) {
  gui::Singleton::SetCmdArgs(argc, argv, env);
  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath(PANORAMIX_CACHE_DATA_DIR_STR "\\LineDrawing\\");

  std::string im_file_name;
  Image3ub im = gui::PickAnImage(PANORAMIX_TEST_DATA_DIR_STR, &im_file_name);
  ResizeToMakeHeightUnder(im, 600);
  std::string folder = misc::FolderOfFile(im_file_name);
  std::string name = misc::NameOfFile(im_file_name);
  std::string anno_file = folder + "\\" + name + ".ldanno.cereal";

  while (true) {
    LineDrawingWidget w(im);
    LoadFromDisk(anno_file, w.anno);
    w.show();

    gui::Singleton::ContinueGui();

    int choice = gui::SelectFrom({"save", "abandon"}, "make your choice",
                                 "your choice?", 0, 1);
    if (choice == 0) {
      SaveToDisk(anno_file, w.anno);
    }
  }

  return 0;
}