#include "../../src/core/parallel.hpp"
#include "../../src/misc/cache.hpp"
#include "../../src/misc/clock.hpp"

#include "../../src/gui/canvas.hpp"
#include "../../src/gui/scene.hpp"
#include "../../src/gui/singleton.hpp"
#include "../../src/gui/utility.hpp"

#include "../../src/experimental/pi_graph_annotation.hpp"
#include "../../src/experimental/pi_graph_cg.hpp"
#include "../../src/experimental/pi_graph_control.hpp"
#include "../../src/experimental/pi_graph_occlusion.hpp"
#include "../../src/experimental/pi_graph_optimize.hpp"
#include "../../src/experimental/pi_graph_vis.hpp"

using namespace pano;
using namespace pano::core;
using namespace pano::experimental;

template <class PointT> struct LineDrawing {
  std::vector<PointT> vertPositions;
  std::vector<std::pair<int, int>> line2verts;
  std::vector<std::vector<int>> face2verts;
  DenseMatd vertAdjMat;
};

LineDrawing<Point3> LoadLineDrawing(const std::string &filename,
                                    const std::string &gtfilename) {
  LineDrawing<Point3> drawing;

  std::ifstream ifs(filename);
  if (!ifs.is_open()) {
    return drawing;
  }

  int lineNum = 0;
  ifs >> lineNum;
  std::vector<Line2> lineMatrix(lineNum);
  for (int i = 0; i < lineNum; i++) {
    ifs >> lineMatrix[i].first[0] >> lineMatrix[i].first[1] >>
        lineMatrix[i].second[0] >> lineMatrix[i].second[1];
  }

  int vertexNum = 0;
  ifs >> vertexNum;
  drawing.vertAdjMat = DenseMat<bool>::zeros(vertexNum, vertexNum);
  for (int i = 0; i < vertexNum; i++) {
    for (int j = 0; j < vertexNum; j++) {
      ifs >> drawing.vertAdjMat(i, j);
    }
  }

  int faceNum = 0;
  ifs >> faceNum;
  drawing.face2verts.resize(faceNum);
  std::string line;
  std::getline(ifs, line);
  for (int i = 0; i < faceNum; i++) {
    std::getline(ifs, line);
    std::istringstream ss(line);
    int vnum = 0;
    ss >> vnum;
    if (vnum == 0) {
      continue;
    }
    int id = 0;
    for (int k = 0; k < vnum; k++) {
      ss >> id;
      drawing.face2verts[i].push_back(id);
    }
  }

  int dummy = 0;
  ifs >> dummy >> dummy;
  drawing.line2verts.resize(lineNum);
  for (int i = 0; i < lineNum; i++) {
    ifs >> drawing.line2verts[i].first >> drawing.line2verts[i].second;
  }

  ifs.close();
  ifs.open(gtfilename);
  if (!ifs.is_open()) {
    return drawing;
  }
  drawing.vertPositions.resize(vertexNum);
  for (int k = 0; k < 3; k++) {
    for (int i = 0; i < vertexNum; i++) {
      ifs >> drawing.vertPositions[i][k];
    }
  }

  return drawing;
}

int main(int argc, char **argv) {
  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath("D:\\Cache\\Panoramix\\LineDrawing\\");
  misc::Matlab matlab;

  std::string folder = "F:\\DataSets\\linedrawing_"
                       "dataset\\DatabaseFile\\bigHouse\\";
  auto drawing =
      LoadLineDrawing(folder + "bigHouse_cy.3dr", folder + "bigHouse.gt");

  gui::SceneBuilder sb;
  sb.installingOptions().defaultShaderSource =
      gui::OpenGLShaderSourceDescriptor::XLines;
  sb.installingOptions().discretizeOptions.color = gui::Black;
  sb.installingOptions().lineWidth = 3.0;
  double scale =
      BoundingBoxOfContainer(drawing.vertPositions).outerSphere().radius;
  for (auto &l : drawing.line2verts) {
    Line3 line3(drawing.vertPositions[l.first],
                drawing.vertPositions[l.second]);
    sb.add(line3 / scale);
  }
  sb.show(
      true, true,
      gui::RenderOptions().backgroundColor(gui::White).renderMode(gui::Lines));
}