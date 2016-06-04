#pragma once

#include "cache.hpp"
#include "clock.hpp"
#include "factor_graph.hpp"
#include "parallel.hpp"

#include "canvas.hpp"
#include "gui_util.hpp"
#include "qttools.hpp"
#include "scene.hpp"
#include "singleton.hpp"

#include "line_drawing.hpp"
#include "mesh_advanced_util.hpp"
#include "pi_graph_annotation.hpp"
#include "pi_graph_cg.hpp"
#include "pi_graph_control.hpp"
#include "pi_graph_occlusion.hpp"
#include "pi_graph_optimize.hpp"
#include "pi_graph_vis.hpp"

namespace pano {
namespace experimental {

// DecomposeFaces
// assume all internal faces are already collected in face2verts
std::vector<std::set<int>>
DecomposeFaces(const std::vector<std::vector<int>> &face2verts,
               const std::vector<Point2> &vert2pos);

// CameraParam
struct CameraParam {
  Point2 pp;
  double focal;
};

// CalibrateCamera
std::vector<CameraParam>
CalibrateCamera(const Box2 &box, const std::vector<std::set<int>> &face_groups,
                std::function<Chain2(int face)> face2chain_fun,
                int k = std::numeric_limits<int>::max());

// BindPointsToLines
std::vector<std::set<int>> BindPointsToLines(const std::vector<Point2> &points,
                                             const std::vector<Line2> &lines,
                                             double angle_thres);

// CollectVanishingPoints
struct CollectVanishingPointsParam { // best params so far
  double angle_thres_phase1 = DegreesToRadians(2);
  double angle_thres_phase2 = DegreesToRadians(0.1);
  double angle_thres_phase3 = DegreesToRadians(8);
  int max_iters = std::numeric_limits<int>::max();
};
std::vector<Point2> CollectVanishingPoints(
    const std::vector<Line2> &lines, double focal, const Point2 &pp,
    const CollectVanishingPointsParam &param = CollectVanishingPointsParam());

// MergeColinearLines
std::vector<Line2>
MergeColinearLines(const std::vector<Line2> &lines,
                   const CameraParam &cam_param, double angle_thres,
                   std::vector<int> *oldline2newline = nullptr);

// EstimateEdgeOrientations
struct EstimateEdgeOrientationsParam { // best params so far
  double angle_thres_allowed_vp_line_deviation = DegreesToRadians(10);
  double angle_thres_judging_colinearility = DegreesToRadians(1);
  double angle_thres_distinguishing_vps = DegreesToRadians(2);
  double angle_thres_juding_coplanarity = DegreesToRadians(10);
  double coeff_vp_line_fitness = 50.0;
  double coeff_noncolinear_adj_line_exlusiveness = 10.0;
  double coeff_line_triplet_coplanar = 30.0;
  int vp_min_degree = 3;
  int solve_max_iter = 5;
};
std::vector<int> EstimateEdgeOrientations(
    const std::vector<Line2> &lines, const std::vector<Point2> &vps,
    const std::vector<std::vector<int>> &face2ordered_lines, double focal,
    const Point2 &pp, const EstimateEdgeOrientationsParam &param =
                          EstimateEdgeOrientationsParam());

// PlaneConstraint
struct PlaneConstraint {
  std::vector<int> verts;
  DenseMatd P; // the matrix P_i in my cvpr16 paper
};
DenseMatd MakePlaneMatrix();
DenseMatd MakePlaneMatrixAlongDirection(const Vec3 & dir);
DenseMatd MakePlaneMatrixTowardDirection(const Vec3 & dir);

// InferenceFunctors
struct InferenceFunctors {
  DenseMatd variables;
  std::function<Vec3(int plane)> getPlaneEquation;
  std::function<double(int vert)> getInversedDepth;
};
// GenerateInferenceFunctors
InferenceFunctors
GenerateInferenceFunctors(int nverts,
                          const std::vector<PlaneConstraint> &constraints);
}
}