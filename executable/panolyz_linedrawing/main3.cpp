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

int main3(int argc, char **argv, char **env) {
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
  PerspectiveCamera cam;
  LoadFromDisk("F:\\LineDrawings\\manifold\\" + name + "\\" + name + ".obj." +
                   cam_name + ".cereal",
               cam);

  {
    gui::SceneBuilder sb;
    sb.add(lineDrawing);
    sb.show(true, true, gui::RenderOptions()
                            .renderMode(gui::Lines)
                            .fixUpDirectionInCameraMove(false));
  }

  // decompose line drawings
  {
	  if (lineDrawing.topo.maybeManifold()) {

	  } else {
		  NOT_IMPLEMENTED_YET();
	  }
  }


  // estimate pp focal candidates using lines

  // vp labeling

  // reconstruction loop
  // 1. optimize with vps
  // 2. optimize without vps


  return 0;
}