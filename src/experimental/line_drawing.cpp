#include "line_drawing.hpp"

namespace pano {
namespace experimental {

LineDrawing<Point2> LoadLineDrawing(const std::string &filename) {
  LineDrawing<Point2> drawing;
  std::ifstream ifs(filename);
  if (!ifs.is_open()) {
    return drawing;
  }

  int lineNum = 0;
  ifs >> lineNum;
  std::vector<Line2> lines(lineNum);
  for (int i = 0; i < lineNum; i++) {
    ifs >> lines[i].first[0] >> lines[i].first[1] >> lines[i].second[0] >>
        lines[i].second[1];
  }

  int vertexNum = 0;
  ifs >> vertexNum;
  for (int i = 0; i < vertexNum; i++) {
    for (int j = 0; j < vertexNum; j++) {
      int dummy = 0;
      ifs >> dummy;
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

  drawing.vertPositions.resize(vertexNum);
  for (int i = 0; i < lineNum; i++) {
    drawing.vertPositions[drawing.line2verts[i].first] = lines[i].first;
    drawing.vertPositions[drawing.line2verts[i].second] = lines[i].second;
  }

  return drawing;
}

LineDrawing<Point3> LoadLineDrawing(const std::string &filename,
                                    const std::string &gtfilename) {
  LineDrawing<Point3> drawing;

  std::ifstream ifs(filename);
  if (!ifs.is_open()) {
    return drawing;
  }

  int lineNum = 0;
  ifs >> lineNum;
  std::vector<Line2> lines(lineNum);
  for (int i = 0; i < lineNum; i++) {
    ifs >> lines[i].first[0] >> lines[i].first[1] >> lines[i].second[0] >>
        lines[i].second[1];
  }

  int vertexNum = 0;
  ifs >> vertexNum;
  for (int i = 0; i < vertexNum; i++) {
    for (int j = 0; j < vertexNum; j++) {
      int dummy = 0;
      ifs >> dummy;
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

}
}