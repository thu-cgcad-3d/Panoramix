#pragma once

#include "mesh.hpp"
#include "mesh_util.hpp"

namespace pano {
namespace experimental {

using namespace pano::core;

// SubMesh
struct SubMesh {
  std::unordered_set<VertHandle> vhs;
  std::unordered_set<HalfHandle> hhs;
  std::unordered_set<FaceHandle> fhs;
  template <class ArchiverT> void serialize(ArchiverT &ar) {
    ar(vhs, hhs, fhs);
  }
};

inline bool Contains(const SubMesh &subMesh, VertHandle vh) {
  return Contains(subMesh.vhs, vh);
}
inline bool Contains(const SubMesh &subMesh, HalfHandle hh) {
  return Contains(subMesh.hhs, hh);
}
inline bool Contains(const SubMesh &subMesh, FaceHandle fh) {
  return Contains(subMesh.fhs, fh);
}

// ExtractSubMeshes
// - HalfEdgeColinearFunT: (HalfEdgeIterT hhsBegin, HalfEdgeIterT hhsEnd) ->
// bool
template <class VertDataT, class HalfDataT, class FaceDataT>
std::vector<SubMesh>
ExtractSubMeshes(const Mesh<VertDataT, HalfDataT, FaceDataT> &mesh);

// GuessUpperBoundOfDRF
template <class VertDataT, class HalfDataT, class FaceDataT,
          class HalfEdgeColinearFunT>
int GuessUpperBoundOfDRF(const Mesh<VertDataT, HalfDataT, FaceDataT> &mesh,
                         SubMesh &subMesh, HalfEdgeColinearFunT colinearFun,
                         int guessNum = 1000);

// Reconstruct (TPAMI 2008, Plane Based Optimization ... )
// - PlaneObjectiveFunT: ( ((FaceHandle)[int]->double) fh2planeeq) -> double
// - VertDataToDirectionFunT: (VertDataT) -> Vec3
template <class VertDataT, class HalfDataT, class FaceDataT,
          class PlaneObjectiveFunT,
          class VertDataToDirectionFunT = std::identity<VertDataT>>
std::unordered_map<FaceHandle, Plane3>
Reconstruct(const Mesh<VertDataT, HalfDataT, FaceDataT> &mesh,
            const SubMesh &sub, int drf, PlaneObjectiveFunT planes2energy,
            VertDataToDirectionFunT vert2dir = VertDataToDirectionFunT(),
            int maxIters = 10000);

// ConstructSingleViewTransformMatrices
// v2matrix: inversed_depth of vh = matrix * X
// f2matrix: plane equation of fh = matrix * X
template <class VertT, class FaceT, class Face2RelatedVertsFunT,
          class Vert2Vec3FunT>
void ConstructSingleViewTransformMatrices(
    const std::vector<VertT> &fundamental_verts,
    const std::vector<FaceT> &ordered_faces,
    Face2RelatedVertsFunT face2related_verts, Vert2Vec3FunT vert2direction,
    std::map<VertHandle, DenseMatd> &v2matrix,
    std::map<FaceHandle, DenseMatd> &f2matrix);

// EstimateSingleViewInversedDepths
// Vert2InitialDepthFunT: (VertT vert) -> double
// Face2RelatedVertsFunT: (FaceT face) -> std::vector<VertT>
// Vert2Vec3FunT: (VertT vert) -> Vec3
// EnergyFunT: (((FaceT face) -> Vec3) face2equation, ((VertT vert) -> double)
// vert2depth) -> DenseMatd
// Returns: fundamental_inversed_depths_vector
template <class VertT, class FaceT, class Vert2InitialDepthFunT,
          class EnergyFunT>
DenseMatd
EstimateSingleViewInversedDepths(const std::map<VertT, DenseMatd> &v2matrix,
                                  const std::map<FaceT, DenseMatd> &f2matrix,
                                  Vert2InitialDepthFunT vert2initial_depth,
                                  EnergyFunT energy_fun);


}
}






////////////////////////////////////////////////
//// implementations
////////////////////////////////////////////////
namespace pano {
namespace experimental {
// ExtractSubMeshes
// - HalfEdgeColinearFunT: (HalfEdgeIterT hhsBegin, HalfEdgeIterT hhsEnd) ->
// bool
template <class VertDataT, class HalfDataT, class FaceDataT>
std::vector<SubMesh>
ExtractSubMeshes(const Mesh<VertDataT, HalfDataT, FaceDataT> &mesh) {
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
  return subMeshes;
}

// GuessUpperBoundOfDRF
template <class VertDataT, class HalfDataT, class FaceDataT,
          class HalfEdgeColinearFunT>
int GuessUpperBoundOfDRF(const Mesh<VertDataT, HalfDataT, FaceDataT> &mesh,
                         SubMesh &subMesh, HalfEdgeColinearFunT colinearFun,
                         int guessNum) {
  // compute drfub
  std::default_random_engine rng;
  std::vector<FaceHandle> fhs(subMesh.fhs.begin(), subMesh.fhs.end());
  int drfub = 0;
  for (int k = 0; k < guessNum; k++) {
    std::shuffle(fhs.begin(), fhs.end(), rng);
    int ub = FindUpperBoundOfDRF(mesh, fhs.begin(), fhs.end(), colinearFun);
    if (drfub < ub) {
      drfub = ub;
    }
  }
  return drfub;
}

// ExtractDeterminacy
// VH2Vec3DirT: (VertHandle vh) -> Vec3
// HH2VPIdT: (HalfHandle hh) -> int
// FH2VPIdT: (FaceHandle fh) -> int
template <class VertDataT, class HalfDataT, class FaceDataT, class VH2Vec3DirT,
          class HH2VPIdT, class FH2VPIdT>
void ExtractSingleViewMatrices(
    const Mesh<VertDataT, HalfDataT, FaceDataT> &mesh, const SubMesh &subMesh,
    const std::vector<Vec3> &vpdirs, VH2Vec3DirT vh2dir, HH2VPIdT hh2vp,
    FH2VPIdT fh2vp, std::map<VertHandle, DenseMatd> &vh2mats,
    std::map<FaceHandle, DenseMatd> &fh2mats) {
  //

}

// Reconstruct (TPAMI 2008, Plane Based Optimization ... )
// - PlaneObjectiveFunT: ( ((FaceHandle)[int]->double) fh2planeeq) -> double
// - VertDataToDirectionFunT: (VertDataT) -> Vec3
template <class VertDataT, class HalfDataT, class FaceDataT,
          class PlaneObjectiveFunT, class VertDataToDirectionFunT>
std::unordered_map<FaceHandle, Plane3>
Reconstruct(const Mesh<VertDataT, HalfDataT, FaceDataT> &mesh,
            const SubMesh &sub, int drf, PlaneObjectiveFunT planes2energy,
            VertDataToDirectionFunT vert2dir, int maxIters) {

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
  assert(drf <= nvars);
  // see Algorithm 2 of the paper
  MatrixXd H = V.rightCols(drf);
  assert(H.cols() == drf && H.rows() == nvars);

  // initialize vars
  VectorXd X = VectorXd::Ones(H.cols());

  static const bool use_simulated_annealing = false;
  // if (use_simulated_annealing) {
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

// ConstructSingleViewTransformMatrices
// inversed_depth of vh = matrix * X
// plane equation of fh = matrix * X
template <class VertT, class FaceT, class Face2RelatedVertsFunT,
          class Vert2Vec3FunT>
void ConstructSingleViewTransformMatrices(
    const std::vector<VertT> &fundamental_verts,
    const std::vector<FaceT> &ordered_faces,
    Face2RelatedVertsFunT face2related_verts, Vert2Vec3FunT vert2direction,
    std::map<VertHandle, DenseMatd> &v2matrix,
    std::map<FaceHandle, DenseMatd> &f2matrix) {
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
    DenseMatd mat = DenseMatd::zeros(1, nvars); // RowVectorXd::Zero(1, nvars);
    mat(0, v2position[v]) = 1;
    v2matrix[v] = mat;
  }

  for (const FaceT &f : ordered_faces) {
    // pick three determined adjacent verts
    auto related_verts = face2related_verts(f);

    Mat3 best_matrix;
    double best_matrix_abs_det = 0.0;
    VertT best_dependencies[3];

    bool found = false;
    for (int i1 = 0; !found && i1 < related_verts.size(); i1++) {
      const VertT &v1 = related_verts[i1];
      if (!Contains(v2matrix, v1)) {
        continue;
      }
      Vec3 dir1 = normalize(vert2direction(v1));
      for (int i2 = i1 + 1; !found && i2 < related_verts.size(); i2++) {
        const VertT &v2 = related_verts[i2];
        if (!Contains(v2matrix, v2)) {
          continue;
        }
        Vec3 dir2 = normalize(vert2direction(v2));
        for (int i3 = i2 + 1; !found && i3 < related_verts.size(); i3++) {
          const VertT &v3 = related_verts[i3];
          if (!Contains(v2matrix, v3)) {
            continue;
          }
          Vec3 dir3 = normalize(vert2direction(v3));

          Mat3 matrix;
          for (int col = 0; col < 3; col++) {
            matrix(0, col) = dir1[col];
            matrix(1, col) = dir2[col];
            matrix(2, col) = dir3[col];
          }

          double abs_det = abs(cv::determinant(matrix));
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
    DenseMatd vmatrices = DenseMatd::zeros(3, nvars);
    for (int row = 0; row < 3; row++) {
      for (int col = 0; col < nvars; col++) {
        vmatrices(row, col) = v2matrix.at(best_dependencies[row])(0, col);
      }
    }

    // 3 x nvars
    f2matrix[f] = DenseMatd(best_matrix.inv()) * vmatrices;
    // best_matrix.fullPivLu().solve(vh_matrices);

    // update all the related vh matrices
    for (const VertT &v : related_verts) {
      if (Contains(v2matrix, v)) {
        continue;
      }
      Vec3 dir = normalize(vert2direction(v));

      // 1 x nvars
      v2matrix[v] = DenseMatd(dir).t() * f2matrix[f];
      assert(v2matrix[v].cols == nvars && v2matrix[v].rows == 1);
    }
  } //  for (FaceHandle fh_proxy : sub.ordered_fhs)
}


// EstimateSingleViewInversedDepths
// Vert2InitialDepthFunT: (VertT vert) -> double
// Face2RelatedVertsFunT: (FaceT face) -> std::vector<VertT>
// Vert2Vec3FunT: (VertT vert) -> Vec3
// EnergyFunT: (((FaceT face) -> Vec3) face2equation, ((VertT vert) -> double)
// vert2depth) -> DenseMatd
// Returns: fundamental_inversed_depths_vector
template <class VertT, class FaceT, class Vert2InitialDepthFunT,
          class EnergyFunT>
DenseMatd
EstimateSingleViewInversedDepths(const std::map<VertT, DenseMatd> &v2matrix_cv,
                                 const std::map<FaceT, DenseMatd> &f2matrix_cv,
                                 Vert2InitialDepthFunT vert2initial_depth,
                                 EnergyFunT energy_fun) {

  using namespace Eigen;

  size_t nverts = v2matrix_cv.size();
  size_t nfaces = f2matrix_cv.size();

  assert(!v2matrix_cv.empty());
  size_t nvars = v2matrix_cv.begin()->second.cols;
  assert(nvars <= nverts);
  for (auto &vm : v2matrix_cv) {
    assert(vm.second.cols == nvars);
  }
  for (auto &fm : f2matrix_cv) {
    assert(fm.second.cols == nvars);
  }

  MatrixXd all_v_matrix = MatrixXd::Zero(nverts, nvars);
  VectorXd all_v_inversed_depths = VectorXd::Zero(nverts);
  std::map<VertT, RowVectorXd> v2matrix; // inversed_depth of vh = matrix * X
  for (auto &vm : v2matrix_cv) {
    const DenseMatd &m_cv = vm.second;
    assert(m_cv.cols == nvars);
    RowVectorXd m = RowVectorXd::Zero(nvars);
    for (int i = 0; i < nvars; i++) {
      m[i] = m_cv(0, i);
    }
    all_v_matrix.row(v2matrix.size()) = m;
    all_v_inversed_depths[v2matrix.size()] = 1.0 / vert2initial_depth(vm.first);
    v2matrix[vm.first] = std::move(m);
  }
  // all_v_inversed_depths = all_v_matrix * inversed_depths
  // -> inversed_depths = 
  VectorXd inversed_depths =
      all_v_matrix.fullPivLu().solve(all_v_inversed_depths).normalized();
  assert(inversed_depths.size() == nvars);

  std::map<FaceT, Matrix3Xd> f2matrix; // plane equation of fh = matrix * X
  for (auto &fm : f2matrix_cv) {
    const DenseMatd &m_cv = fm.second;
    assert(m_cv.cols == nvars);
    Matrix3Xd m = Matrix3Xd::Zero(3, m_cv.cols);
    for (int row = 0; row < 3; row++) {
      for (int col = 0; col < nvars; col++) {
        m(row, col) = m_cv(row, col);
      }
    }
    f2matrix[fm.first] = std::move(m);
  }

  // use LM algorithm
  auto energy_wrapper_fun = [&f2matrix, &v2matrix,
                             &energy_fun](const VectorXd &X) -> VectorXd {
    std::vector<double> energy_terms = energy_fun(
        [&v2matrix, &X](const VertT &vert) -> double {
          return 1.0 / (v2matrix.at(vert) * X.normalized().cwiseAbs());
        },
        [&f2matrix, &X](const FaceT &face) -> Vec3 {
          RowVector3d eq = f2matrix.at(face) * X.normalized().cwiseAbs();
          return Vec3(eq[0], eq[1], eq[2]);
        });
    /*Println("cur energy = ",
            std::accumulate(energy_terms.begin(), energy_terms.end(), 0.0,
                            [](double e, double t) { return e + t * t; }));*/
    return VectorXd::Map(energy_terms.data(), energy_terms.size());
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

  DenseMatd inversed_depths_cv = DenseMatd::zeros(nvars, 1);
  for (int i = 0; i < nvars; i++) {
    inversed_depths_cv(i, 0) = inversed_depths[i];
  }
  return inversed_depths_cv;
}
}
}