#pragma once

#include "basic_types.hpp"
#include "cameras.hpp"
#include "iterators.hpp"
#include "line_drawing.hpp"
#include "line_drawing_widget.hpp"
#include "utility.hpp"

namespace pano {
namespace experimental {

using namespace ::pano::core;

// DecomposeFaces
// assume all internal faces are already collected in face2verts
std::vector<std::set<int>>
DecomposeFaces(const std::vector<std::vector<int>> &face2verts,
               const std::vector<Point2> &vert2pos);

// CameraParam
struct CameraParam {
  Point2 pp;
  double focal;
  template <class ArchiverT> void serialize(ArchiverT &ar) { ar(pp, focal); }
};

// CalibrateCamera
std::vector<CameraParam>
CalibrateCamera(const Box2 &box, const std::vector<std::set<int>> &face_groups,
                std::function<std::vector<Chain2>(int face)> face2chain_fun,
                int k = std::numeric_limits<int>::max());

// BindPointsToLines
std::vector<std::set<int>> BindPointsToLines(const std::vector<Point2> &points,
                                             const std::vector<Line2> &lines,
                                             double angle_thres);
std::vector<std::set<int>> BindPointsToLines(const std::vector<Vec3> &points,
                                             const std::vector<Line3> &lines,
                                             double angle_thres);
std::vector<std::set<int>> BindPointsToLines(const std::vector<Point2> &points,
                                             const std::vector<Line2> &lines,
                                             const CameraParam &cam_param,
                                             double angle_thres);

// CollectLineIntersections
std::vector<Vec3>
CollectLineIntersections(const std::vector<Line3> &lines,
                         std::vector<std::pair<int, int>> *line_ids = nullptr);

// CollectVanishingPoints
struct CollectVanishingPointsParam { // best params so far
  double angle_thres_phase1 = DegreesToRadians(2);
  double angle_thres_phase2 = DegreesToRadians(0.1);
  double angle_thres_phase3 = DegreesToRadians(8);
  int max_iters = std::numeric_limits<int>::max();
  bool use_mean_shift_merge_phase1 = false;
};
// std::vector<Point2> CollectVanishingPoints(
//    const std::vector<Line2> &lines, const CameraParam &cam_param,
//    const CollectVanishingPointsParam &param = CollectVanishingPointsParam());
std::vector<Vec3> CollectVanishingPoints(
    const std::vector<Line3> &lines,
    const CollectVanishingPointsParam &param = CollectVanishingPointsParam());

// MergeColinearLines
std::vector<Line2>
MergeColinearLines(const std::vector<Line2> &lines,
                   const CameraParam &cam_param, double angle_thres,
                   std::vector<int> *oldline2newline = nullptr);

//// EstimateParallelism
// std::map<std::pair<int, int>, double>
// EstimateParallelism(const std::vector<Line2> &lines,
//                    const CameraParam &cam_param, double angle_thres);

// EstimateEdgeOrientations
struct EstimateEdgeOrientationsParam { // best params so far
  double angle_thres_allowed_vp_line_deviation = DegreesToRadians(10);
  double angle_thres_judging_colinearility = DegreesToRadians(1);
  double angle_thres_distinguishing_vps = DegreesToRadians(2);
  double angle_thres_juding_orthogonality = DegreesToRadians(10);
  double angle_thres_juding_coplanarity = DegreesToRadians(10);
  double coeff_vp_line_fitness = 50.0;
  double coeff_noncolinear_adj_line_exlusiveness = 10.0;
  double coeff_line_pair_orthogonality = 20.0;
  double coeff_line_triplet_coplanar = 30.0;
  int vp_min_degree = 3;
  int solve_max_iter = 5;
};
// std::vector<int> EstimateEdgeOrientations(
//    const std::vector<Line2> &lines, const std::vector<Point2> &vps,
//    const std::vector<std::vector<int>> &face2ordered_lines, double focal,
//    const Point2 &pp, const EstimateEdgeOrientationsParam &param =
//                          EstimateEdgeOrientationsParam());
std::vector<int> EstimateEdgeOrientations(
    const std::vector<Line3> &lines, const std::vector<Vec3> &vps,
    const std::vector<std::pair<int, int>> &adjacent_line_pairs,
    const std::vector<std::vector<int>> &coplanar_ordered_lines,
    const EstimateEdgeOrientationsParam &param =
        EstimateEdgeOrientationsParam());

// PlaneConstraint
struct PlaneConstraint {
  std::vector<int> verts;
  DenseMatd P; // the matrix P_i in my cvpr16 paper
  template <class ArchiverT> void serialize(ArchiverT &ar) { ar(verts, P); }
};
DenseMatd MakePlaneMatrix();
DenseMatd MakePlaneMatrixAlongDirection(const Vec3 &dir);
DenseMatd MakePlaneMatrixTowardDirection(const Vec3 &dir);

// InferenceFunctors
struct Inferencer {
  virtual size_t nvars() const = 0;
  virtual Vec3 getPlaneEquation(int cons, const DenseMatd &variables) const = 0;
  virtual double getInversedDepth(int vert,
                                  const DenseMatd &variables) const = 0;
  virtual DenseMatd
  recoverVariables(const std::vector<double> &vert2inversed_depths) const = 0;
};
// GenerateInferenceFunctors
std::unique_ptr<Inferencer>
GenerateInferenceFunctors(const std::vector<PlaneConstraint> &constraints,
                          const std::vector<Vec3> &vert2dir, int root_vert = 0,
                          std::vector<int> *fundamental_verts = nullptr);

// The Energy Terms
std::vector<double> AnglesBetweenAdjacentEdges(
    const std::vector<Vec3> &vert2dir,
    const std::vector<std::vector<int>> &face2verts, const DenseMatd &variables,
    const Inferencer &infer,
    std::function<bool(int v1, int v2)> edge_selected = nullptr,
    std::function<bool(int face)> face_selected = nullptr);

std::vector<double> AnglesBetweenAdjacentFaces(
    size_t nfaces, const std::vector<std::set<int>> &edge2faces,
    const DenseMatd &variables, const Inferencer &infer,
    std::function<bool(int face)> face_selected = nullptr);

std::vector<double> AnglesBetweenAdjacentFaces(
    size_t nfaces, const std::vector<std::set<int>> &edge2faces,
    const DenseMatd &variables, const Inferencer &infer,
    const std::map<std::pair<int, int>, bool> &faces_overlap,
    std::function<bool(int face)> face_selected = nullptr);

// PerformReconstruction
struct PerformReconstructionParam {
  int max_iters = 100;
};
double PerformReconstruction(
    const std::vector<PlaneConstraint> &constraints,
    const std::vector<Vec3> &vert2dir, int root_vert,
    std::function<double(const Inferencer &infer, const DenseMatd &variables,
                         const std::vector<Vec3> &vert2dir)>
        energy_fun,
    std::default_random_engine &rng, std::vector<Point3> &vert2pos,
    std::vector<int> *fundamental_verts_ptr = nullptr,
    const PerformReconstructionParam &param = PerformReconstructionParam());
}
}