#pragma once

#include "../core/basic_types.hpp"
#include "../core/mesh.hpp"

namespace pano {
namespace experimental {

using namespace pano::core;

// LineDrawing
enum LineType { SolidLine, DashedLine };
template <class PointT> struct LineDrawing {
  std::vector<PointT> corners;
  std::vector<std::pair<int, int>> line2corners;
  std::vector<LineType> line2type;
  std::vector<std::vector<int>> face2corners;
  template <class ArchiverT> void serialize(ArchiverT &ar) {
    ar(corners, line2corners, line2type, face2corners);
  }
};

// LoadLineDrawing
LineDrawing<Point2> LoadLineDrawing(const std::string &filename);
LineDrawing<Point3> LoadLineDrawing(const std::string &filename,
                                    const std::string &gtfilename);

// ToMesh
template <class T, class PointConvertT = std::identity<T>>
Mesh<std::decay_t<typename FunctionTraits<PointConvertT>::ResultType>>
ToMesh(const LineDrawing<T> &ld, PointConvertT cvtFun = PointConvertT()) {
  Mesh<std::decay_t<typename FunctionTraits<PointConvertT>::ResultType>> mesh;
  using VHandle = VertHandle;
  mesh.internalVertices().reserve(ld.corners.size());
  mesh.internalHalfEdges().reserve(ld.line2corners.size() * 2);
  mesh.internalFaces().reserve(ld.face2corners.size());
  for (int i = 0; i < ld.corners.size(); i++) {
    mesh.addVertex(cvtFun(ld.corners[i]));
  }
  for (int i = 0; i < ld.line2corners.size(); i++) {
    mesh.addEdge(VHandle(ld.line2corners[i].first),
                 VHandle(ld.line2corners[i].second));
  }
  if (ld.face2corners.empty()) {
    return mesh;
  }
  std::vector<int> faceInsertOrder = {0};
  std::set<int> notInsertedFace;
  for (int i = 1; i < ld.face2corners.size(); i++) {
    notInsertedFace.insert(i);
  }
  while (!notInsertedFace.empty()) {
  LabelInsertFace:
    for (int candFace : notInsertedFace) {
      for (int insertedFace : faceInsertOrder) {
        auto &candCorners = ld.face2corners[candFace];
        auto &insertedCorners = ld.face2corners[insertedFace];
        for (int i1 = 0; i1 < candCorners.size(); i1++) {
          int from1 = candCorners[i1];
          int to1 = candCorners[(i1 + 1) % candCorners.size()];
          for (int i2 = 0; i2 < insertedCorners.size(); i2++) {
            int from2 = insertedCorners[i2];
            int to2 = insertedCorners[(i2 + 1) % insertedCorners.size()];
            if (from1 == from2 && to1 == to2 || from1 == to2 && to1 == from2) {
              // has common edge!
              notInsertedFace.erase(candFace);
              faceInsertOrder.push_back(candFace);
              goto LabelInsertFace;
            }
          }
        }
      }
    }
  }

  for (int i : faceInsertOrder) {
    std::vector<VHandle> vhs(ld.face2corners[i].size());
    for (int j = 0; j < vhs.size(); j++) {
      vhs[j] = VHandle(ld.face2corners[i][j]);
    }
    mesh.addFace(vhs, true);
  }

  AssertEdgesAreStiched(mesh);
  return mesh;
}

// ToLineDrawing
template <class T, class H, class F, class VertConvertT = std::identity<T>>
LineDrawing<std::decay_t<typename FunctionTraits<VertConvertT>::ResultType>>
ToLineDrawing(const Mesh<T, H, F> &mesh, VertConvertT cvtFun = VertConvertT()) {
  LineDrawing<std::decay_t<typename FunctionTraits<VertConvertT>::ResultType>>
      ld;
  using VHandle = VertHandle;
  std::vector<int> vh2corner(mesh.internalVertices().size(), -1);
  for (auto &v : mesh.vertices()) {
    ld.corners.push_back(cvtFun(v.data));
    vh2corner[v.topo.hd.id] = ld.corners.size() - 1;
  }
  std::vector<bool> hh2used(mesh.internalHalfEdges().size(), false);
  for (auto &h : mesh.halfedges()) {
    if (hh2used[h.topo.hd.id] || hh2used[mesh.topo(h.topo.hd).opposite.id]) {
      continue;
    }
    ld.line2corners.emplace_back(vh2corner[mesh.topo(h.topo.hd).from().id],
                                 vh2corner[mesh.topo(h.topo.hd).to().id]);
    hh2used[h.topo.hd.id] = true;
  }
  ld.line2type.resize(ld.line2corners.size(), SolidLine);
  for (auto &f : mesh.faces()) {
    std::vector<int> corners(f.topo.halfedges.size());
    for (int j = 0; j < corners.size(); j++) {
      corners[j] = vh2corner[mesh.topo(f.topo.halfedges[j]).to().id];
    }
    ld.face2corners.push_back(std::move(corners));
  }
  return ld;
}

// Transform
template <class T, class FunT>
auto Transform(const LineDrawing<T> &in, const FunT &fun) {
  LineDrawing<typename FunctionTraits<FunT>::ResultType> out;
  out.corners.resize(in.corners.size());
  std::transform(in.corners.begin(), in.corners.end(), out.corners.begin(),
                 fun);
  out.line2corners = in.line2corners;
  out.line2type = in.line2type;
  out.face2corners = in.face2corners;
  return out;
}

template <class T, class FunT>
auto Transform(LineDrawing<T> &&in, const FunT &fun) {
  LineDrawing<typename FunctionTraits<FunT>::ResultType> out;
  out.corners.resize(in.corners.size());
  std::transform(in.corners.begin(), in.corners.end(), out.corners.begin(),
                 fun);
  out.line2corners = std::move(in.line2corners);
  out.line2type = std::move(in.line2type);
  out.face2corners = std::move(in.face2corners);
  return out;
}

// SubMesh
struct SubMesh {
  std::unordered_set<VertHandle> vhs;
  std::unordered_set<HalfHandle> hhs;
  std::unordered_set<FaceHandle> fhs;
  int drfub;
};

// ExtractSubMeshes
template <class VertDataT, class HalfDataT, class FaceDataT>
std::vector<SubMesh>
ExtractSubMeshes(const Mesh<VertDataT, HalfDataT, FaceDataT> &mesh,
                 int drfTryNum = 10) {
  using MeshT = Mesh<VertDataT, HalfDataT, FaceDataT>;
  std::unordered_map<int, std::unordered_set<VertHandle>> cc2vhs;
  HandledTable<VertHandle, int> vh2ccid(mesh.internalVertices().size(), -1);
  int nccs = ConnectedComponents(
      mesh, [&cc2vhs, &vh2ccid](const MeshT &, VertHandle vh, int ccid) {
        cc2vhs[ccid].insert(vh);
        vh2ccid[vh] = ccid;
      });

  std::vector<SubMesh> subMeshes(nccs);
  for (int i = 0; i < nccs; i++) {
    subMeshes[i].vhs = std::move(cc2vhs[i]);
  }

  HandledTable<HalfHandle, int> hh2ccid(mesh.internalHalfEdges().size(), -1);
  for (auto &h : mesh.halfedges()) {
    assert(vh2ccid[h.topo.from()] == vh2ccid[h.topo.to()]);
    int ccid = vh2ccid[h.topo.from()];
    subMeshes[ccid].hhs.insert(h.topo.hd);
    hh2ccid[h.topo.hd] = ccid;
  }

  for (auto &f : mesh.faces()) {
    int ccid = hh2ccid[f.topo.halfedges.front()];
    assert(std::all_of(
        f.topo.halfedges.begin(), f.topo.halfedges.end(),
        [ccid, &hh2ccid](HalfHandle hh) { return hh2ccid[hh] == ccid; }));
    subMeshes[ccid].fhs.insert(f.topo.hd);
  }

  // compute drfub
  std::default_random_engine rng;
  for (int i = 0; i < nccs; i++) {
    std::vector<FaceHandle> fhs(subMeshes[i].fhs.begin(),
                                subMeshes[i].fhs.end());
    subMeshes[i].drfub = 0;
    for (int k = 0; k < drfTryNum; k++) {
      std::shuffle(fhs.begin(), fhs.end(), rng);
      int ub = FindUpperBoundOfDRF(mesh, fhs.begin(), fhs.end(),
                                   [](auto hhsBegin, auto hhsEnd) -> bool {
                                     if (hhsBegin == hhsEnd) {
                                       return true;
                                     }
                                     // todo
                                   });
      if (subMeshes[i].drfub < ub) {
        subMeshes[i].drfub = ub;
      }
    }
  }

  return subMeshes;
}

// Reconstruct (TPAMI 2008, Plane Based Optimization ... )
// - PlaneObjectiveFunT: (((FaceHandle)[int]->double) fh2planeeq, double*
// scores) -> void
// - VertDataToDirectionFunT: (VertDataT) -> Vec3
template <class VertDataT, class HalfDataT, class FaceDataT,
          class PlaneObjectiveFunT,
          class VertDataToDirectionFunT = std::identity<VertDataT>>
std::unordered_map<FaceHandle, Plane3>
Reconstruct(const Mesh<VertDataT, HalfDataT, FaceDataT> &mesh,
            const SubMesh &sub, PlaneObjectiveFunT planes2scores, int scoresNum,
            VertDataToDirectionFunT vert2dir = VertDataToDirectionFunT(),
            int drfSlack = 2) {

  using namespace Eigen;

  int nverts = sub.vhs.size();
  int nhalves = sub.hhs.size();
  int nedges = nhalves / 2;
  int nfaces = sub.fhs.size();

  int nequations = nedges * 2 - nverts;
  int nvars = nfaces * 3;

  // order faces
  std::vector<FaceHandle> fhs(sub.fhs.begin(), sub.fhs.end());
  std::unordered_map<FaceHandle, int> fh2varStart;
  int varid = 0;
  for (FaceHandle fh : fhs) {
    fh2varStart[fh] = varid;
    varid += 3;
  }

  // construct projection matrix P
  MatrixXd P = MatrixXd::Zero(nequations, nvars);
  int eid = 0;
  for (VertHandle vh : sub.vhs) {
    auto &hhs = mesh.topo(vh).halfedges;
    for (auto it = std::next(hhs.begin()); it != hhs.end()++ it) {
      auto hh = *it;
      assert(Contains(sub.hhs, hh));
      HalfHandle oppohh = mesh.topo(hh).opposite;
      FaceHandle fh1 = mesh.topo(hh).face;
      FaceHandle fh2 = mesh.topo(oppohh).face;
      int varStartId1 = fh2varStart.at(fh1);
      int varStartId2 = fh2varStart.at(fh2);
      assert(varStartId1 != -1 && varStartId2 != -1);

      auto p = normalize(vert2dir(mesh.data(vh)));
      for (int i = 0; i < 3; i++) {
        P(eid, varStartId1 + i) = p[i];
        P(eid, varStartId2 + i) = -p[i];
      }
      eid++;
    }
  }

  // perform SVD on P
  // P = USV^t
  JacobiSVD<MatrixXd> svd(P, ComputeFullV);
  auto deltas = svd.singularValues(); // decreasing order defaultly
  auto V = svd.matrixV();
  assert(V.cols() == nvars);
  assert(sub.drfub + drfSlack <= nvars);
  // see Algorithm 2 of the paper
  MatrixXd H = V.rightCols(sub.drfub + drfSlack);
  assert(H.cols() == sub.drfub + drfSlack && H.rows() == nvars);

  auto functor = misc::MakeGenericNumericDiffFunctor<double>(
      [&planes2scores, &H, &fh2varStart](const VectorXd &x, VectorXd &v) {
        VectorXd planeEqVec = H * x;
        planes2scores(
            [&planeEqVec, &fh2varStart](FaceHandle fh) -> const double * {
              return planeEqVec.data() + fh2varStart.at(fh);
            },
            v.data());
      },
      H.cols(), scoresNum);
  LevenbergMarquardt<decltype(functor), double> lm(functor);

  // initialize vars
  VectorXd X = VectorXd::Ones(H.cols());

  // solve
  lm.minimize(X);

  // convert to planes
  std::unordered_map<FaceHandle, Plane3> planes;
  VectorXd planeEqVec = H * X;
  for (FaceHandle fh : fhs) {
    int varStart = fh2varStart[fh];
    planes[fh] =
        Plane3FromEquation(planeEqVec[varStart], planeEqVec[varStart + 1],
                           planeEqVec[varStart + 2]);
  }
  return planes;
}
}

namespace core {
template <class PointT>
inline auto
BoundingBox(const pano::experimental::LineDrawing<PointT> &drawing) {
  return BoundingBoxOfContainer(drawing.corners);
}
}
}
