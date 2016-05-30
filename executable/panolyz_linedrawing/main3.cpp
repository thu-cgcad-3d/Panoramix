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

struct AlgorithmInput {
  LineDrawingTopo lineDrawingTopo;
  std::vector<Point2> corners2d;
  std::vector<Point3> cornersGT;
  PerspectiveCamera cameraGT;
};

// ParseInput
AlgorithmInput ParseInput(const std::string &modelName,
                          const std::string &camName) {
  auto lineDrawingGT = LoadLineDrawingFromObjFile(
      "F:\\LineDrawings\\manifold\\" + modelName + "\\" + modelName + ".obj");
  assert(lineDrawingGT.ncorners() > 0);
  PerspectiveCamera cam;
  std::string camPath = "F:\\LineDrawings\\manifold\\" + modelName + "\\" +
                        modelName + ".obj." + camName + ".cereal";
  bool succ = LoadFromDisk(camPath, cam);
  if (!succ) {
    gui::SceneBuilder sb;
    sb.add(lineDrawingGT);
    cam = sb.show(true, true, gui::RenderOptions()
                                  .renderMode(gui::Lines)
                                  .fixUpDirectionInCameraMove(false))
              .camera;
    SaveToDisk(camPath, cam);
  }
  
  std::vector<Point2> corners2d(lineDrawingGT.corners.size());
  for (int i = 0; i < corners2d.size(); i++) {
    corners2d[i] = cam.toScreen(lineDrawingGT.corners[i]);
  }
  return AlgorithmInput{std::move(lineDrawingGT.topo), std::move(corners2d),
                        std::move(lineDrawingGT.corners), std::move(cam)};
}

int main3(int argc, char **argv, char **env) {
  gui::Singleton::SetCmdArgs(argc, argv, env);
  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath("D:\\Panoramix\\LineDrawing\\");
  misc::Matlab matlab;

  AlgorithmInput input = ParseInput("tower", "cam1");

  // group line drawing entities
  {
	  
  }


  // estimate pp focal candidates using lines

  // vp labeling

  // reconstruction loop
  // 1. optimize with vps
  // 2. optimize without vps


  return 0;
}