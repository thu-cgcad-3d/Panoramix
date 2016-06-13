#include "cache.hpp"
#include "cameras.hpp"
#include "canvas.hpp"
#include "clock.hpp"
#include "line_drawing.hpp"
#include "line_drawing_widget.hpp"
#include "optimization.hpp"
#include "parallel.hpp"
#include "scene.hpp"
#include "singleton.hpp"
#include "tools.hpp"

using namespace pano;
using namespace pano::core;
using namespace pano::experimental;

void RoutineEditImageAnno() {
  std::string im_file_name;
  Image3ub im = gui::PickAnImage(PANORAMIX_TEST_DATA_DIR_STR, &im_file_name);
  ResizeToMakeHeightUnder(im, 600);
  std::string folder = misc::FolderOfFile(im_file_name);
  std::string name = misc::NameOfFile(im_file_name);
  std::string anno_file = folder + "\\" + name + ".ldanno.cereal";

  while (true) {
    LineDrawingWidget w(im);
    LoadFromDisk(anno_file, w.anno);
    w.anno.image = im;
    w.show();

    gui::Singleton::ContinueGui();
    int choice = gui::SelectFrom({"save", "abandon"}, "make your choice",
                                 "your choice?", 0, 1);
    if (choice == 0) {
      SaveToDisk(anno_file, w.anno);
    }
  }
}