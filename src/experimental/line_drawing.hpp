#pragma once

#include "../core/basic_types.hpp"
#include "../core/mesh.hpp"
#include "../core/mesh_util.hpp"

#include "../misc/eigen.hpp"

#include "optimization.hpp"

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
  bool contains(VertHandle h) const {return Contains(vhs, h); }
  bool contains(HalfHandle h) const {return Contains(hhs, h); }
  bool contains(FaceHandle h) const {return Contains(fhs, h); }

  std::unordered_set<VertHandle> fundamental_vhs;
  std::vector<FaceHandle> ordered_fhs;
  template <class VT, class HT, class FT, class VertColinearFunT>
  void ComputeVertexFaceDependencies(const Mesh<VT, HT, FT> &mesh,
                                     VertColinearFunT vhs_colinear_fun) {
    ordered_fhs = std::vector<FaceHandle>(fhs.begin(), fhs.end());
    fundamental_vhs = SortFacesWithPlanarityDependency(
        mesh, ordered_fhs.begin(), ordered_fhs.end(), vhs_colinear_fun);
  }

  template <class ArchiverT> void serialize(ArchiverT &ar) {
    ar(vhs, hhs, fhs, drfub);
  }
};

// ExtractSubMeshes
// - HalfEdgeColinearFunT: (HalfEdgeIterT hhsBegin, HalfEdgeIterT hhsEnd) ->
// bool
template <class VertDataT, class HalfDataT, class FaceDataT,
          class HalfEdgeColinearFunT>
std::vector<SubMesh>
ExtractSubMeshes(const Mesh<VertDataT, HalfDataT, FaceDataT> &mesh,
                 HalfEdgeColinearFunT colinearFun, int drfTryNum = 10) {
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
      int ub = FindUpperBoundOfDRF(mesh, fhs.begin(), fhs.end(), colinearFun);
      if (subMeshes[i].drfub < ub) {
        subMeshes[i].drfub = ub;
      }
    }
  }

  return subMeshes;
}


// Reconstruct (TPAMI 2008, Plane Based Optimization ... )
// - PlaneObjectiveFunT: ( ((FaceHandle)[int]->double) fh2planeeq) -> double
// - VertDataToDirectionFunT: (VertDataT) -> Vec3
template <class VertDataT, class HalfDataT, class FaceDataT,
          class PlaneObjectiveFunT,
          class VertDataToDirectionFunT = std::identity<VertDataT>>
std::unordered_map<FaceHandle, Plane3>
Reconstruct(const Mesh<VertDataT, HalfDataT, FaceDataT> &mesh,
            const SubMesh &sub,
            PlaneObjectiveFunT planes2energy,
            VertDataToDirectionFunT vert2dir = VertDataToDirectionFunT(),
            int drfSlack = 2, int maxIters = 10000) {

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
    for (auto it = std::next(hhs.begin()); it != hhs.end(); ++it) {
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

  // initialize vars
  VectorXd X = VectorXd::Ones(H.cols());

  static const bool use_simulated_annealing = false;
  //if (use_simulated_annealing) {
    std::default_random_engine rng;
    SimulatedAnnealing(
        X,
        [&planes2energy, &H, &fh2varStart](const VectorXd &curX) -> double {
          VectorXd planeEqVec = H * curX;
          planeEqVec /= planeEqVec.cwiseAbs().maxCoeff();
          return planes2energy(
              [&planeEqVec, &fh2varStart](FaceHandle fh) -> const double * {
                return planeEqVec.data() + fh2varStart.at(fh);
              });
        },
        [](int iter) {
          return std::max(1.0 / log(log(iter + 2)), 1e-10);
        }, // temperature
        [maxIters](VectorXd curX, int iter, auto &&forEachNeighborFun) {
          if (iter >= maxIters) {
            return;
          }
          double step = std::max(1.0 / log(iter + 2), 1e-10);
          for (int i = 0; i < curX.size(); i++) {
            double curXHere = curX[i];
            curX[i] = curXHere + step;
            forEachNeighborFun(curX);
            curX[i] = curXHere - step;
            forEachNeighborFun(curX);
            curX[i] = curXHere;
          }
        },
        rng);
  //} else {
    auto functor = misc::MakeGenericNumericDiffFunctor<double>(
        [&planes2energy, &H, &fh2varStart](const VectorXd &curX, VectorXd &e) {
          VectorXd planeEqVec = H * curX;
          planeEqVec /= planeEqVec.cwiseAbs().maxCoeff();
          e[0] = planes2energy(
              [&planeEqVec, &fh2varStart](FaceHandle fh) -> const double * {
                return planeEqVec.data() + fh2varStart.at(fh);
              });
        },
        X.size(), 1);
    LevenbergMarquardt<decltype(functor)> lm(std::move(functor));
    lm.minimize(X);
  //}

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

// ConstructPerspectiveTransformMatrices
template <class VertT, class FaceT, class Face2RelatedVertsFunT,
          class Vert2RowVector3dFunT>
void ConstructPerspectiveTransformMatrices(
    const std::vector<VertT> &fundamental_verts,
    const std::vector<FaceT> &ordered_faces,
    Face2RelatedVertsFunT face2related_verts,
    Vert2RowVector3dFunT vert2direction,
    std::map<VertHandle, Eigen::RowVectorXd>
        &v2matrix, // inversed_depth of vh = matrix * X
    std::map<FaceHandle, Eigen::Matrix3Xd>
        &f2matrix // plane equation of fh = matrix * X
    ) {
  v2matrix.clear();
  f2matrix.clear();

  using namespace Eigen;
  size_t nvars = fundamental_verts.size();

  std::unordered_map<VertT, int> v2position;
  for (int pos = 0; pos < nvars; pos++) {
    v2position[fundamental_verts[pos]] = pos;
  }

  // build matrix for each vertex and face in current sub using the dependency
  // information

  // insert the matrices of fundamental verts
  for (const VertT &v : fundamental_verts) {
    auto &mat = v2matrix[v];
    mat = RowVectorXd::Zero(1, nvars);
    mat(0, v2position[v]) = 1;
  }

  for (const FaceT &f : ordered_faces) {
    // pick three determined adjacent verts
    auto related_verts = face2related_verts(f);

    Matrix3d best_matrix = Matrix3d::Zero();
    double best_matrix_abs_det = 0.0;
    VertT best_dependencies[3];

    bool found = false;
    for (int i1 = 0; !found && i1 < related_verts.size(); i1++) {
      const VertT &v1 = related_verts[i1];
      if (!Contains(v2matrix, v1)) {
        continue;
      }
      RowVector3d dir1 = vert2direction(v1).normalized();
      for (int i2 = i1 + 1; !found && i2 < related_verts.size(); i2++) {
        const VertT &v2 = related_verts[i2];
        if (!Contains(v2matrix, v2)) {
          continue;
        }
        RowVector3d dir2 = vert2direction(v2).normalized();
        for (int i3 = i2 + 1; !found && i3 < related_verts.size(); i3++) {
          const VertT &v3 = related_verts[i3];
          if (!Contains(v2matrix, v3)) {
            continue;
          }
          RowVector3d dir3 = vert2direction(v3).normalized();

          Matrix3d matrix = Matrix3d::Zero();
          matrix << dir1, dir2, dir3;

          double abs_det = abs(matrix.determinant());
          if (abs_det > best_matrix_abs_det) {
            best_matrix = matrix;
            best_matrix_abs_det = abs_det;
            best_dependencies[0] = v1;
            best_dependencies[1] = v2;
            best_dependencies[2] = v3;
            if (best_matrix_abs_det > 0.1) {
              found = true;
            }
          }
        } // for (int i3 = i2 + 1; i3 < hhs.size(); i3++)
      }   // for (int i2 = i1 + 1; i2 < hhs.size(); i2++)
    }     // for (int i1 = 0; i1 < hhs.size(); i1++)

    // matrix A = [x1 y1 z1]
    //			[x2 y2 z2]
    //			[x3 y3 z3]
    // [plane eq] = A^-1 * [id1 id2 id3]^T

    assert(best_matrix_abs_det > 0);

    // collect the matrices of three vhs
    Matrix3Xd vmatrices = Matrix3Xd::Zero(3, nvars);
    vmatrices << v2matrix.at(best_dependencies[0]),
        v2matrix.at(best_dependencies[1]), v2matrix.at(best_dependencies[2]);

    for (int row = 0; row < 3; row++) {
      assert(vmatrices.row(row) == v2matrix.at(best_dependencies[row]));
    }

    // 3 x nvars
    f2matrix[f] = best_matrix.inverse() * vmatrices;
    // best_matrix.fullPivLu().solve(vh_matrices);

    // update all the related vh matrices
    for (const VertT &v : related_verts) {
      if (Contains(v2matrix, v)) {
        continue;
      }
      RowVector3d dir = vert2direction(v).normalized();
      // 1 x nvars
      v2matrix[v] = dir * f2matrix[f];
    }
  } //  for (FaceHandle fh_proxy : sub.ordered_fhs)
}

// EstimatePerspectiveDepths
// Vert2InitialDepthFunT: (VertT vert) -> double
// Face2RelatedVertsFunT: (FaceT face) -> std::vector<VertT>
// Vert2RowVector3dFunT: (VertT vert) -> Eigen::RowVector3d
// EnergyFunT: (((FaceT face) -> Eigen::Vector3d) face2equation, ((VertT vert) -> double) vert2depth) -> Eigen::VectorXd
template <class VertT, class FaceT, class Vert2InitialDepthFunT,
          class Face2RelatedVertsFunT, class Vert2RowVector3dFunT,
          class EnergyFunT>
Eigen::VectorXd EstimatePerspectiveDepths(
    const std::vector<VertT> &fundamental_verts,
    const std::vector<FaceT> &ordered_faces,
    Vert2InitialDepthFunT vert2initial_depth,
    Face2RelatedVertsFunT face2related_verts,
    Vert2RowVector3dFunT vert2direction, EnergyFunT energy_fun,
    std::map<VertT, Eigen::RowVectorXd> *v2matrix_ptr = nullptr,
    std::map<FaceT, Eigen::Matrix3Xd> *f2matrix_ptr = nullptr) {

  using namespace Eigen;
  size_t nvars = fundamental_verts.size();

  VectorXd inversed_depths = VectorXd::Zero(fundamental_verts.size());
  for (int pos = 0; pos < nvars; pos++) {
    inversed_depths[pos] = 1.0 / vert2initial_depth(fundamental_verts[pos]);
  }

  std::unordered_map<VertT, int> v2position;
  for (int pos = 0; pos < nvars; pos++) {
    v2position[fundamental_verts[pos]] = pos;
  }

  // build matrix for each vertex and face in current sub using the dependency
  // information
  std::map<VertT, RowVectorXd> v2matrix; // inversed_depth of vh = matrix * X
  std::map<FaceT, Matrix3Xd> f2matrix;   // plane equation of fh = matrix * X

  ConstructPerspectiveTransformMatrices(fundamental_verts, ordered_faces,
                                        face2related_verts, vert2direction,
                                        v2matrix, f2matrix);

  // use LM algorithm
  auto energy_wrapper_fun = [&f2matrix, &v2matrix,
                             &energy_fun](const VectorXd &X) -> VectorXd {
    return energy_fun(
        [&v2matrix, &X](const VertT &vert) -> double {
          return 1.0 / (v2matrix.at(vert) * X.normalized().cwiseAbs());
        },
        [&f2matrix, &X](const FaceT &face) -> Vector3d {
          return f2matrix.at(face) * X.normalized().cwiseAbs();
        });
  };

  auto initial_energy_terms = energy_wrapper_fun(inversed_depths);

  double initial_energy = 0.0;
  for (int k = 0; k < initial_energy_terms.size(); k++) {
    initial_energy += Square(initial_energy_terms[k]);
  }
  Println("initial energy = ", initial_energy);

  auto functor = misc::MakeGenericNumericDiffFunctor<double>(
      [&energy_wrapper_fun](const VectorXd &curX, VectorXd &e) {
        e = energy_wrapper_fun(curX);
      },
      inversed_depths.size(), initial_energy_terms.size());

  LevenbergMarquardt<decltype(functor)> lm(functor);
  lm.parameters.maxfev = 10000;
  lm.minimize(inversed_depths);
  inversed_depths = inversed_depths.normalized().cwiseAbs();

  auto final_energy_terms = energy_wrapper_fun(inversed_depths);
  double final_energy = 0.0;
  for (int k = 0; k < final_energy_terms.size(); k++) {
    final_energy += Square(final_energy_terms[k]);
  }

  core::Println("final energy = ", final_energy);

  if (v2matrix_ptr) {
    *v2matrix_ptr = std::move(v2matrix);
  }
  if (f2matrix_ptr) {
    *f2matrix_ptr = std::move(f2matrix);
  }

  return inversed_depths;
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
