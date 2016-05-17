#include "mesh_util.hpp"

namespace pano {
namespace core {
Mesh<Point3> LoadFromObjFile(const std::string &fname) {
  Mesh<Point3> mesh;
  std::ifstream ifs(fname);
  if (ifs.is_open()) {
    std::string line;
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
        mesh.addVertex(pos);
      } else if (token == "f") {
        std::vector<VertHandle> vhs;
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
          vhs.push_back(VertHandle(vid - 1));
        }
        if (!vhs.empty()) {
          auto fh = mesh.addFace(vhs, true);
          assert(fh.valid());
        }
      }
    }
  }
  return mesh;
}
}
}
