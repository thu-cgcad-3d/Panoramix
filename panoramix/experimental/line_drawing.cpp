#include "eigen.hpp"
#include <Eigen/SVD>

#include "line_drawing.hpp"
#include "mesh.hpp"

namespace pano {
namespace experimental {

LineDrawingTopo::LineDrawingTopo(const std::vector<std::pair<int, int>> &e2cs,
                                 const std::vector<std::vector<int>> &f2cs)
    : face2corners(f2cs) {

  size_t nfaces = f2cs.size();
  size_t ncorners = 0;

  for (int edge = 0; edge < e2cs.size(); edge++) {
    auto cpair = MakeOrderedPair(e2cs[edge]);
    ncorners =
        std::max({ncorners, (size_t)cpair.first + 1, (size_t)cpair.second + 1});
    if (!Contains(corners2edge, cpair)) {
      corners2edge[cpair] = corners2edge.size();
    }
  }

  for (int face = 0; face < nfaces; face++) {
    auto &cs = face2corners[face];
    assert(cs.size() >= 3);
    for (int i = 0; i < cs.size(); i++) {
      int c1 = cs[i];
      int c2 = cs[(i + 1) % cs.size()];
      ncorners = std::max(ncorners, (size_t)c1 + 1);
      auto cpair = MakeOrderedPair(c1, c2);
      if (!Contains(corners2edge, cpair)) {
        corners2edge[cpair] = corners2edge.size();
      }
    }
  }

  size_t nedges = corners2edge.size();
  edge2corners.resize(nedges);
  corner2edges.resize(ncorners);
  for (auto &p : corners2edge) {
    edge2corners[p.second] = p.first;
    corner2edges[p.first.first].push_back(p.second);
    corner2edges[p.first.second].push_back(p.second);
  }

  face2edges.resize(nfaces);
  edge2faces.resize(nedges);
  corner2faces.resize(ncorners);
  for (int face = 0; face < nfaces; face++) {
    auto &cs = face2corners[face];
    assert(cs.size() >= 3);
    for (int i = 0; i < cs.size(); i++) {
      int c1 = cs[i];
      int c2 = cs[(i + 1) % cs.size()];
      corner2faces[c1].push_back(face);
      auto cpair = MakeOrderedPair(c1, c2);
      int edge = corners2edge.at(cpair);
      face2edges[face].push_back(edge);
      edge2faces[edge].push_back(face);
    }
  }
}

bool LineDrawingTopo::maybeManifold() const {
  for (int c = 0; c < ncorners(); c++) {
    if (corner2edges[c].size() <= 1) {
      return false;
    }
    if (corner2faces[c].size() <= 1) {
      return false;
    }
    if (corner2edges[c].size() != corner2faces[c].size()) {
      return false;
    }
  }
  for (int e = 0; e < nedges(); e++) {
    if (edge2corners[e].first == edge2corners[e].second) {
      return false;
    }
    if (edge2faces[e].size() != 2) {
      return false;
    }
  }
  for (int f = 0; f < nfaces(); f++) {
    if (face2corners[f].size() <= 2) {
      return false;
    }
    if (face2edges[f].size() <= 2) {
      return false;
    }
    if (face2corners[f].size() != face2edges[f].size()) {
      return false;
    }
  }
  return true;
}

LineDrawing<Point3> LoadLineDrawingFromObjFile(const std::string &fname) {
  std::ifstream ifs(fname);
  if (ifs.is_open()) {
    std::string line;

    std::vector<Point3> corners;
    std::vector<std::vector<int>> face2corners;

    while (std::getline(ifs, line)) {
      if (line.empty()) {
        continue;
      }
      std::istringstream ss(line);
      std::string token;
      ss >> token;
      if (token == "v") {
        Point3 pos;
        ss >> pos[0] >> pos[1] >> pos[2];
        corners.push_back(pos);
      } else if (token == "f") {
        std::vector<int> corners;
        while (ss >> token) {
          if (token.empty()) {
            continue;
          }
          int vid = -1;
          size_t p = token.find_first_of('/');
          if (p == std::string::npos) {
            vid = std::stoi(token);
          } else {
            vid = std::stoi(token.substr(0, p));
          }
          assert(vid != -1);
          corners.push_back(vid - 1);
        }
        if (!corners.empty()) {
          face2corners.push_back(std::move(corners));
        }
      }
    }

    return LineDrawing<Point3>({}, face2corners, std::move(corners));
  }
  return LineDrawing<Point3>();
}
// LineDrawing<Point2> LoadLineDrawing(const std::string &filename) {
//  std::ifstream ifs(filename);
//  if (!ifs.is_open()) {
//    return LineDrawing<Point2>();
//  }
//
//  int lineNum = 0;
//  ifs >> lineNum;
//  std::vector<Line2> lines(lineNum);
//  for (int i = 0; i < lineNum; i++) {
//    ifs >> lines[i].first[0] >> lines[i].first[1] >> lines[i].second[0] >>
//        lines[i].second[1];
//  }
//
//  int vertexNum = 0;
//  ifs >> vertexNum;
//  for (int i = 0; i < vertexNum; i++) {
//    for (int j = 0; j < vertexNum; j++) {
//      int dummy = 0;
//      ifs >> dummy;
//    }
//  }
//
//  int faceNum = 0;
//  ifs >> faceNum;
//  drawing.face2corners.resize(faceNum);
//  std::string line;
//  std::getline(ifs, line);
//  for (int i = 0; i < faceNum; i++) {
//    std::getline(ifs, line);
//    std::istringstream ss(line);
//    int vnum = 0;
//    ss >> vnum;
//    if (vnum == 0) {
//      continue;
//    }
//    int id = 0;
//    for (int k = 0; k < vnum; k++) {
//      ss >> id;
//      drawing.face2corners[i].push_back(id);
//    }
//  }
//
//  int dummy = 0;
//  ifs >> dummy >> dummy;
//  drawing.edge2corners.resize(lineNum);
//  for (int i = 0; i < lineNum; i++) {
//    ifs >> drawing.edge2corners[i].first >> drawing.edge2corners[i].second;
//  }
//
//  ifs.close();
//
//  drawing.corners.resize(vertexNum);
//  for (int i = 0; i < lineNum; i++) {
//    drawing.corners[drawing.edge2corners[i].first] = lines[i].first;
//    drawing.corners[drawing.edge2corners[i].second] = lines[i].second;
//  }
//
//  return drawing;
//}
//
// LineDrawing<Point3> LoadLineDrawing(const std::string &filename,
//                                    const std::string &gtfilename) {
//  LineDrawing<Point3> drawing;
//
//  std::ifstream ifs(filename);
//  if (!ifs.is_open()) {
//    return drawing;
//  }
//
//  int lineNum = 0;
//  ifs >> lineNum;
//  std::vector<Line2> lines(lineNum);
//  for (int i = 0; i < lineNum; i++) {
//    ifs >> lines[i].first[0] >> lines[i].first[1] >> lines[i].second[0] >>
//        lines[i].second[1];
//  }
//
//  int vertexNum = 0;
//  ifs >> vertexNum;
//  for (int i = 0; i < vertexNum; i++) {
//    for (int j = 0; j < vertexNum; j++) {
//      int dummy = 0;
//      ifs >> dummy;
//    }
//  }
//
//  int faceNum = 0;
//  ifs >> faceNum;
//  drawing.face2corners.resize(faceNum);
//  std::string line;
//  std::getline(ifs, line);
//  for (int i = 0; i < faceNum; i++) {
//    std::getline(ifs, line);
//    std::istringstream ss(line);
//    int vnum = 0;
//    ss >> vnum;
//    if (vnum == 0) {
//      continue;
//    }
//    int id = 0;
//    for (int k = 0; k < vnum; k++) {
//      ss >> id;
//      drawing.face2corners[i].push_back(id);
//    }
//  }
//
//  int dummy = 0;
//  ifs >> dummy >> dummy;
//  drawing.edge2corners.resize(lineNum);
//  for (int i = 0; i < lineNum; i++) {
//    ifs >> drawing.edge2corners[i].first >> drawing.edge2corners[i].second;
//  }
//
//  ifs.close();
//  ifs.open(gtfilename);
//  if (!ifs.is_open()) {
//    return drawing;
//  }
//  drawing.corners.resize(vertexNum);
//  for (int k = 0; k < 3; k++) {
//    for (int i = 0; i < vertexNum; i++) {
//      ifs >> drawing.corners[i][k];
//    }
//  }
//
//  return drawing;
//}
}
}