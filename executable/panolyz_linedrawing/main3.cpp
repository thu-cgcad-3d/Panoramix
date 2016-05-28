#include "factor_graph.hpp"
#include "parallel.hpp"
#include "cache.hpp"
#include "clock.hpp"

#include "canvas.hpp"
#include "qttools.hpp"
#include "scene.hpp"
#include "singleton.hpp"
#include "gui_util.hpp"

#include "line_drawing.hpp"
#include "mesh_advanced_util.hpp"

using namespace pano;
using namespace pano::core;
using namespace pano::experimental;

int main(int argc, char **argv, char **env) {
  gui::Singleton::SetCmdArgs(argc, argv, env);
  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath("D:\\Panoramix\\LineDrawing\\");
  misc::Matlab matlab;

  // std::string name = "gate";
  // std::string name = "towerx";
  std::string name = "tower";
  // std::string name = "hex";
  // std::string name = "triangle";
  // std::string name = "twotriangles";
  // std::string name = "bridge";
  std::string cam_name = "cam1";

  auto lineDrawing = LoadLineDrawingFromObjFile("F:\\LineDrawings\\manifold\\" +
                                                name + "\\" + name + ".obj");
  


  return 0;
}