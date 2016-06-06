#include <set>

extern "C" {
#include <wclique.h>
}

#include "algorithms.hpp"
#include "containers.hpp"
#include "eigen.hpp"
#include "factor_graph.hpp"
#include "iterators.hpp"
#include "manhattan.hpp"
#include "math.hpp"
#include "tools.hpp"
#include "utility.hpp"

namespace pano {
namespace experimental {

using namespace ::pano::core;

// ICCV 2013 Complex 3D General Object Reconstruction from Line Drawings
// Algorithm 1
std::set<int> ExtendFaces(const std::vector<std::vector<int>> &face2verts,
                          const std::pair<int, int> &initial_faces,
                          const std::vector<Point2> &vert2pos) {
  std::set<int> f_fixed = {initial_faces.first, initial_faces.second};
  std::set<int> f_unfixed;
  for (int f = 0; f < face2verts.size(); f++) {
    if (f != initial_faces.first && f != initial_faces.second) {
      f_unfixed.insert(f);
    }
  }

  auto initial_verts = MakeConcatedRange(face2verts[initial_faces.first],
                                         face2verts[initial_faces.second]);
  std::set<int> v_fixed(initial_verts.begin(), initial_verts.end());

  while (true) {
    int insertible_face = -1;
    for (int f : f_unfixed) {
      std::vector<Point2> anchor_vert_pos;
      for (int v : face2verts[f]) {
        if (core::Contains(v_fixed, v)) {
          anchor_vert_pos.push_back(vert2pos[v]);
        }
      }
      if (!IsFuzzyColinear(anchor_vert_pos.begin(), anchor_vert_pos.end(),
                           0.01)) {
        insertible_face = f;
        break;
      }
    }
    if (insertible_face != -1) {
      f_unfixed.erase(insertible_face);
      f_fixed.insert(insertible_face);
      for (int v : face2verts[insertible_face]) {
        v_fixed.insert(v);
      }
    } else {
      break;
    }
  }
  return f_fixed;
}

// ICCV 2013 Complex 3D General Object Reconstruction from Line Drawings
// Algorithm 2
std::vector<std::set<int>>
DecomposeFaces(const std::vector<std::vector<int>> &face2verts,
               const std::vector<Point2> &vert2pos) {
  std::map<std::pair<int, int>, std::set<int>> edge2faces;
  for (int f = 0; f < face2verts.size(); f++) {
    auto &verts = face2verts[f];
    for (int i = 0; i < verts.size(); i++) {
      int v1 = verts[i];
      int v2 = verts[(i + 1) % verts.size()];
      edge2faces[MakeOrderedPair(v1, v2)].insert(f);
    }
  }
  // for each adj face pair collect mefs
  std::set<std::set<int>> mefs_set;
  for (auto &p : edge2faces) {
    auto &faces = p.second;
    if (faces.size() < 2) {
      continue;
    }
    for (auto it1 = faces.begin(); it1 != faces.end(); ++it1) {
      for (auto it2 = std::next(it1); it2 != faces.end(); ++it2) {
        auto f_fixed =
            ExtendFaces(face2verts, std::make_pair(*it1, *it2), vert2pos);
        if (f_fixed.size() > 2) {
          mefs_set.insert(std::move(f_fixed));
        }
      }
    }
  }

  //// find max weight cliques from mefs
  std::vector<std::set<int>> mefs;
  mefs.reserve(mefs_set.size());
  for (auto &s : mefs_set) {
    mefs.push_back(std::move(s));
  }
  // w[1...n] stores the weight, start from 1 because wclique.c requires that
  std::vector<int> w(mefs.size() + 1);
  for (int i = 0; i < mefs.size(); i++) {
    w[i + 1] = mefs[i].size();
  }
  std::vector<unsigned char> adjmat(
      mefs.size() * (mefs.size() + 1) / 2 / CHAR_BIT + 1, 0);
  // std::set<std::pair<int, int>> adj;
  { // construct the ugly adjmat according to wclique.c
    for (int i = 0; i < mefs.size(); i++) {
      for (int j = i + 1; j < mefs.size(); j++) {
        std::vector<int> common;
        std::set_intersection(mefs[i].begin(), mefs[i].end(), mefs[j].begin(),
                              mefs[j].end(), std::back_inserter(common));
        if (common.size() <= 1) {
          // adj.emplace(i, j);
          // see the is_edge macro defined in wclique.c, line: 124
          int offset = j * (j - 1) / 2 + i;
          adjmat[offset / CHAR_BIT] |=
              (unsigned char)(1 << ((CHAR_BIT - 1) - (offset) % CHAR_BIT));
        }
      }
    }
  }
  for (auto &mask : adjmat) {
    Println(std::bitset<8>(mask).to_string());
  }
  std::vector<int> mwc_inds_plus_one(mefs.size() + 1);
  int mwc_size =
      wclique(mefs.size(), w.data(), adjmat.data(), mwc_inds_plus_one.data());
  assert(mwc_size > 0);
  // get mwcs
  std::vector<std::set<int>> mwcs(mwc_size);
  for (int i = 1; i <= mwc_size; i++) {
    mwcs[i - 1] = std::move(mefs[mwc_inds_plus_one[i] - 1]);
  }

  // add the rest faces to current mwcs
  std::vector<std::set<int>> face2mwcs(face2verts.size());
  for (int m = 0; m < mwcs.size(); m++) {
    for (int f : mwcs[m]) {
      face2mwcs[f].insert(m);
    }
  }
  for (int f = 0; f < face2verts.size(); f++) {
    if (face2mwcs[f].empty()) {
      auto &verts = face2verts[f];
      std::vector<int> mwc2adjcount(mwcs.size(), 0);
      for (int i = 0; i < verts.size(); i++) {
        int v1 = verts[i];
        int v2 = verts[(i + 1) % verts.size()];
        auto &adjfs = edge2faces.at(MakeOrderedPair(v1, v2));
        for (int adjf : adjfs) {
          for (int mwc : face2mwcs[adjf]) {
            mwc2adjcount[mwc]++;
          }
        }
      }
      // get the max connected mwc
      mwcs[*std::max_element(mwc2adjcount.begin(), mwc2adjcount.end())].insert(
          f);
    }
  }

  return mwcs;
}

// PossibleKeyVanishingPoints
std::vector<Point2> ProposeVanishingPoints(const Chain2 &chain) {
  assert(chain.size() > 2);
  if (chain.size() == 3) {
    return {};
  }
  if (chain.size() == 4) {
    return {Intersection(chain.edge(0).ray(), chain.edge(2).ray()),
            Intersection(chain.edge(1).ray(), chain.edge(3).ray())};
  }
  if (chain.size() % 2 == 0) {
    std::vector<Point2> vp_positions;
    vp_positions.reserve(chain.size() * chain.size() / 4);
    for (int i = 0; i < chain.size(); i++) {
      int j = (i + chain.size() / 2) % chain.size();
      vp_positions.push_back(
          Intersection(chain.edge(i).ray(), chain.edge(j).ray()));
    }
    return vp_positions;
  } else {
    std::vector<Point2> vp_positions;
    for (int i = 0; i < chain.size(); i++) {
      int j = (i + chain.size() / 2) % chain.size();
      vp_positions.push_back(
          Intersection(chain.edge(i).ray(), chain.edge(j).ray()));
      vp_positions.push_back(
          Intersection(chain.edge(i).ray(), chain.edge(j + 1).ray()));
    }
    return vp_positions;
  }
}

// CalibrateCamera
std::vector<CameraParam>
CalibrateCamera(const Box2 &box, const std::vector<std::set<int>> &face_groups,
                std::function<Chain2(int face)> face2chain_fun, int k) {

  double scale = box.outerSphere().radius;

  // collect pp focal candidates
  std::vector<CameraParam> pp_focal_candidates;
  pp_focal_candidates.reserve(face_groups.size() * 3);

  for (auto &group : face_groups) {
    // collect edge intersections in each face
    std::vector<Point2> interps;
    for (int face : group) {
      Chain2 corners = face2chain_fun(face);
      auto key_vps = ProposeVanishingPoints(corners);
      // remove all vps that lie on certain lines
      for (auto &key_vp : key_vps) {
        bool overlaped = false;
        for (int i = 0; i < corners.size(); i++) {
          Line2 edge = corners.edge(i);
          if (Distance(key_vp, edge) < 1e-3) {
            overlaped = true;
            break;
          }
        }
        if (!overlaped) {
          interps.push_back(key_vp);
        }
      }
    }

    for (int i = 0; i < interps.size(); i++) {
      const Point2 &p1 = interps[i];
      for (int j = i + 1; j < interps.size(); j++) {
        const Point2 &p2 = interps[j];
        for (int k = j + 1; k < interps.size(); k++) {
          const Point2 &p3 = interps[k];
          // compute pp and focal
          Point2 pp;
          double focal = 0.0;
          std::tie(pp, focal) = ComputePrinciplePointAndFocalLength(p1, p2, p3);
          if (HasValue(pp, IsInfOrNaN<double>) || IsInfOrNaN(focal)) {
            continue;
          }
          if (!IsBetween(focal, scale / 5.0, scale * 5.0) ||
              Distance(pp, box.center()) > scale) {
            continue;
          }
          pp_focal_candidates.push_back(CameraParam{pp, focal});
        }
      }
    }
  }

  std::sort(pp_focal_candidates.begin(), pp_focal_candidates.end(),
            [](auto &a, auto &b) { return a.focal < b.focal; });
  // naive clustering
  std::vector<std::pair<std::set<int>, CameraParam>> pp_focal_groups;
  {
    std::vector<int> ppFocalId2group(pp_focal_candidates.size(), -1);
    int ngroups = 0;
    RTreeMap<Vec3, int> ppFocalIdTree;
    for (int i = 0; i < pp_focal_candidates.size(); i++) {
      Vec3 coordinate =
          cat(pp_focal_candidates[i].pp, pp_focal_candidates[i].focal);
      const double thres = scale / 50.0;
      // find the nearest ppFocal sample point
      int nearestPPFocalCandId = -1;
      double min_dist = thres;
      ppFocalIdTree.search(BoundingBox(coordinate).expand(thres * 2),
                           [&nearestPPFocalCandId, &min_dist,
                            &coordinate](const std::pair<Vec3, int> &cand) {
                             double dist = Distance(cand.first, coordinate);
                             if (dist < min_dist) {
                               min_dist = dist;
                               nearestPPFocalCandId = cand.second;
                             }
                             return true;
                           });
      if (nearestPPFocalCandId != -1) { // if found, assign to the same group
        ppFocalId2group[i] = ppFocalId2group[nearestPPFocalCandId];
      } else { // otherwise, create a new group
        ppFocalId2group[i] = ngroups++;
      }
      ppFocalIdTree.emplace(coordinate, i);
    }

    pp_focal_groups.resize(ngroups);
    for (auto &g : pp_focal_groups) {
      g.second.focal = 0.0;
      g.second.pp = Point2();
    }
    for (int i = 0; i < ppFocalId2group.size(); i++) {
      auto &g = pp_focal_groups[ppFocalId2group[i]];
      g.first.insert(i);
      g.second.focal += pp_focal_candidates[i].focal;
      g.second.pp += pp_focal_candidates[i].pp;
    }
    for (auto &g : pp_focal_groups) {
      g.second.focal /= g.first.size();
      g.second.pp /= double(g.first.size());
    }

    std::sort(
        pp_focal_groups.begin(), pp_focal_groups.end(),
        [](auto &g1, auto &g2) { return g1.first.size() > g2.first.size(); });
  }

  pp_focal_candidates.clear();
  for (int i = 0; i < k && i < pp_focal_groups.size(); i++) {
    pp_focal_candidates.push_back(std::move(pp_focal_groups[i].second));
  }
  return pp_focal_candidates;
}

template <class T> struct BinaryRelationTable {
  std::vector<T> relations;
  size_t nelements;
  explicit BinaryRelationTable(size_t n, const T &v)
      : nelements(n), relations(n * (n - 1) / 2, v) {}
  decltype(auto) operator()(int i, int j) const {
    if (i == j) {
      return T();
    }
    int offset = i < j ? (j * (j - 1) / 2 + i) : (i * (i - 1) / 2 + j);
    return relations[offset];
  }
  decltype(auto) operator()(int i, int j) {
    assert(i != j);
    int offset = i < j ? (j * (j - 1) / 2 + i) : (i * (i - 1) / 2 + j);
    return relations[offset];
  }
  constexpr auto neighbors(int i) const {
    return MakeConditionalRange(MakeIotaRange<int>(nelements),
                                [this, i](int ind) { return (*this)(i, ind); });
  }
};

inline Vec3 NormalizedSpatialDirection(const Point2 &p2, double focal,
                                       const Point2 &pp) {
  return normalize(Vec3(p2[0] - pp[0], p2[1] - pp[1], focal));
}
inline Point2 PlanarPoint(const Vec3 &dir, double focal, const Point2 &pp) {
  double dir2 = dir[2];
  if (dir2 == 0) {
    dir2 = 1e-10;
  }
  return focal * Point2(dir[0] / dir2, dir[1] / dir2) + pp;
}

inline void AddToDirection(Vec3 &dir, const Vec3 &add) {
  if (dir.dot(add) < 0) {
    dir -= normalize(add);
  } else {
    dir += normalize(add);
  }
}

std::vector<Point2> MergePoints(const std::vector<Point2> &intersections,
                                double focal, const Point2 &pp,
                                double angle_thres) {
  // unproject to space
  std::vector<Vec3> directions(intersections.size());
  for (int i = 0; i < intersections.size(); i++) {
    directions[i] = NormalizedSpatialDirection(intersections[i], focal, pp);
  }
  // merge
  std::vector<Vec3> merged_directions;
  merged_directions.reserve(intersections.size());
  for (auto &dir : directions) {
    int where_to_merge = -1;
    double min_angle = angle_thres;
    for (int i = 0; i < merged_directions.size(); i++) {
      double angle = AngleBetweenUndirected(dir, merged_directions[i]);
      if (angle <= min_angle) {
        min_angle = angle;
        where_to_merge = i;
      }
    }
    if (where_to_merge != -1) {
      auto &mdir = merged_directions[where_to_merge];
      AddToDirection(mdir, dir);
    } else {
      merged_directions.push_back(dir);
    }
  }
  std::vector<Point2> merged_points(merged_directions.size());
  for (int i = 0; i < merged_points.size(); i++) {
    merged_points[i] = PlanarPoint(merged_directions[i], focal, pp);
  }
  return merged_points;
}

std::vector<Point2> MergePoints(const std::vector<Point2> &intersections,
                                double focal, const Point2 &pp,
                                const BinaryRelationTable<bool> &adj_relation) {
  // unproject to space
  std::vector<Vec3> directions(intersections.size());
  for (int i = 0; i < intersections.size(); i++) {
    directions[i] = NormalizedSpatialDirection(intersections[i], focal, pp);
  }
  // merge
  std::vector<Vec3> merged_directions(directions.size(), Vec3());
  int ccs = ConnectedComponents(
      MakeIotaIterator<int>(0), MakeIotaIterator<int>(directions.size()),
      [&adj_relation](int id) { return adj_relation.neighbors(id); },
      [&merged_directions, &directions](int id, int group) {
        AddToDirection(merged_directions[group], directions[id]);
      });
  std::vector<Point2> merged_points(ccs);
  for (int i = 0; i < ccs; i++) {
    merged_points[i] = PlanarPoint(merged_directions[i], focal, pp);
  }
  return merged_points;
}

std::vector<std::set<int>> BindPointsToLines(const std::vector<Point2> &points,
                                             const std::vector<Line2> &lines,
                                             double angle_thres) {
  std::vector<std::set<int>> point2lines(points.size());
  for (int inter = 0; inter < points.size(); inter++) {
    auto &interp = points[inter];
    for (int l = 0; l < lines.size(); l++) {
      auto &line = lines[l];
      double lambda = ProjectionOfPointOnLine(interp, line).ratio;
      static const double thres = 0.1;
      // if the intersection projects on the body of the line, pass
      if (lambda >= -thres && lambda <= 1.0 + thres) {
        continue;
      }
      double angle =
          AngleBetweenUndirected(line.direction(), interp - line.center());
      assert(!IsInfOrNaN(angle));
      if (angle >= angle_thres) {
        continue;
      }
      point2lines[inter].insert(l);
    }
  }
  return point2lines;
}

// CollectVanishingPoints
std::vector<Point2>
CollectVanishingPoints(const std::vector<Line2> &lines, double focal,
                       const Point2 &pp,
                       const CollectVanishingPointsParam &param) {

  Box2 box = BoundingBoxOfContainer(lines);
  double scale = box.outerSphere().radius;

  // collect intersections
  std::vector<Point2> intersections;
  intersections.reserve(lines.size() * (lines.size() - 1) / 2);
  for (int i = 0; i < lines.size(); i++) {
    const Line2 &linei = lines[i];
    for (int j = i + 1; j < lines.size(); j++) {
      const Line2 &linej = lines[j];
      Point2 interp = Intersection(linei.ray(), linej.ray());
      if (IsFuzzyParallel(linei.direction(), linej.direction(), 1e-8)) {
        interp = normalize(linei.direction()) * 1e8;
      }
      if (std::min(Distance(interp, linei), Distance(interp, linej)) <=
          scale / 10.0) {
        continue;
      }
      intersections.push_back(interp);
    }
  }

  for (int iteration = 0; iteration < param.max_iters; iteration++) {
    Println("intersecitons num = ", intersections.size());
    size_t prev_intersections_num = intersections.size();

    // phase 1: directly merge intersections
    intersections =
        MergePoints(intersections, focal, pp, param.angle_thres_phase1);
    Println("1: after directly merge, intersecitons num = ",
            intersections.size());

    // phase 2: merge intersections who share two or more lines
    // first get adjacency relations
    auto intersection2lines =
        BindPointsToLines(intersections, lines, param.angle_thres_phase2);
    BinaryRelationTable<bool> intersection_sharing_many_lines(
        intersections.size(), false);
    for (int i = 0; i < intersections.size(); i++) {
      auto &linesi = intersection2lines[i];
      for (int j = i + 1; j < intersections.size(); j++) {
        auto &linesj = intersection2lines[j];
        // get common lines
        std::vector<int> common_lines;
        common_lines.reserve(std::min(linesi.size(), linesj.size()));
        std::set_intersection(linesi.begin(), linesi.end(), linesj.begin(),
                              linesj.end(), std::back_inserter(common_lines));
        if (common_lines.size() >= 2) { // should merge i and j
          intersection_sharing_many_lines(i, j) = true;
        }
      }
    }
    // merge using connected components
    auto group_centers =
        MergePoints(intersections, focal, pp, intersection_sharing_many_lines);
    for (auto &pos : group_centers) {
      assert(!(IsInfOrNaN(pos[0]) || IsInfOrNaN(pos[1])));
    }
    Println("2: after CC, group nums = ", group_centers.size());

    // phase 3: remove group centers that bind to too few lines
    auto group_center2lines =
        BindPointsToLines(group_centers, lines, param.angle_thres_phase3);
    intersections.clear();
    for (int group = 0; group < group_centers.size(); group++) {
      if (group_center2lines[group].size() >= 3) {
        intersections.push_back(group_centers[group]);
      }
    }
    Println("3: after removing lonely intersections, intersection nums = ",
            intersections.size());

    {
      double min_dist_in_vps = std::numeric_limits<double>::infinity();
      for (int i = 0; i < intersections.size(); i++) {
        for (int j = i + 1; j < intersections.size(); j++) {
          double dist = Distance(intersections[i], intersections[j]);
          if (dist < min_dist_in_vps) {
            min_dist_in_vps = dist;
          }
        }
      }
      Println("the min distance between two different vps is ",
              min_dist_in_vps);
    }

    assert(intersections.size() <= prev_intersections_num);
    if (intersections.size() == prev_intersections_num) {
      break;
    }
  }

  return intersections;
}

inline bool LinesAreColinear(const Line2 &line1, const Line2 &line2,
                             double focal, const Point2 &pp,
                             double angle_thres) {
  Vec3 dir1 = NormalizedSpatialDirection(line1.first, focal, pp)
                  .cross(NormalizedSpatialDirection(line1.second, focal, pp));
  Vec3 dir2 = NormalizedSpatialDirection(line2.first, focal, pp)
                  .cross(NormalizedSpatialDirection(line2.second, focal, pp));
  return AngleBetweenUndirected(dir1, dir2) < angle_thres;
}

std::vector<Line2> MergeColinearLines(const std::vector<Line2> &lines,
                                      const CameraParam &cam_param,
                                      double angle_thres,
                                      std::vector<int> *oldline2newline) {

  BinaryRelationTable<bool> lines_colinear(lines.size(), false);
  for (int i = 0; i < lines.size(); i++) {
    for (int j = i + 1; j < lines.size(); j++) {
      lines_colinear(i, j) = LinesAreColinear(
          lines[i], lines[j], cam_param.focal, cam_param.pp, angle_thres);
    }
  }
  std::vector<std::vector<int>> groups(lines.size());
  if (oldline2newline) {
    oldline2newline->resize(lines.size(), -1);
  }
  int ccs = ConnectedComponents(
      MakeIotaIterator<int>(0), MakeIotaIterator<int>(lines.size()),
      [&lines_colinear](int lineid) {
        return lines_colinear.neighbors(lineid);
      },
      [oldline2newline, &groups](int oldline, int newline) {
        groups[newline].push_back(oldline);
        if (oldline2newline) {
          (*oldline2newline)[oldline] = newline;
        }
      });
  groups.resize(ccs);

  std::vector<Line2> merged_lines(groups.size());
  for (int i = 0; i < groups.size(); i++) {
    assert(groups[i].size() > 0);
    const Line2 &first_line = lines[groups[i][0]];
    Ray2 ray = first_line.ray();
    assert(ray.direction != Vec2());
    merged_lines[i] = first_line;
    double minproj = 0.0;
    double maxproj = first_line.length();
    for (int lineid : groups[i]) {
      for (auto p : {lines[lineid].first, lines[lineid].second}) {
        double proj =
            (p - first_line.first).dot(normalize(first_line.direction()));
        if (proj < minproj) {
          minproj = proj;
          merged_lines[i].first = p;
        }
        if (proj > maxproj) {
          maxproj = proj;
          merged_lines[i].second = p;
        }
      }
    }

    assert(first_line.length() <= merged_lines[i].length());
    assert(IsFuzzyParallel(merged_lines[i].direction(), first_line.direction(),
                           0.1));
  }

  return merged_lines;
}

// EstimateEdgeOrientations
std::vector<int> EstimateEdgeOrientations(
    const std::vector<Line2> &lines, const std::vector<Point2> &vps,
    const std::vector<std::vector<int>> &face2ordered_lines, double focal,
    const Point2 &pp, const EstimateEdgeOrientationsParam &param) {

  size_t nlines = lines.size();
  size_t nvps = vps.size();
  std::vector<std::vector<Scored<int>>> line2ordered_vp_angles(nlines);

  for (int l = 0; l < nlines; l++) {
    auto &line = lines[l];
    for (int vp = 0; vp < nvps; vp++) {
      auto &vp_pos = vps[vp];
      double lambda = ProjectionOfPointOnLine(vp_pos, line).ratio;
      static const double thres = 0.1;
      if (lambda >= -thres &&
          lambda <= 1.0 + thres) { // if vp's projection lies on the line, pass
        continue;
      }
      double angle =
          AngleBetweenUndirected(line.direction(), vp_pos - line.center());
      assert(!IsInfOrNaN(angle));
      if (angle >= param.angle_thres_allowed_vp_line_deviation) { // deviation
                                                                  // too big,
                                                                  // pass //
        // DegreesToRadians(10)
        continue;
      }
      line2ordered_vp_angles[l].push_back(ScoreAs(vp, angle));
    }
    std::sort(line2ordered_vp_angles[l].begin(),
              line2ordered_vp_angles[l].end());
  }

  // construct a factor graph to optimize line-vp bindings
  FactorGraph fg;
  std::vector<int> result;
  {
    // add lines as vars
    for (int l = 0; l < nlines; l++) {
      fg.addVar(fg.addVarCategory(line2ordered_vp_angles[l].size() + 1, 1.0));
    }

    // potential 1: the edge should bind to some vp with better
    // scored
    for (int vh = 0; vh < nlines; vh++) {
      auto &related_vp_angles = line2ordered_vp_angles[vh];
      auto fc = fg.addFactorCategory(
          [&related_vp_angles, nlines, &param](
              const std::vector<int> &varlabels, void *givenData) -> double {
            assert(varlabels.size() == 1);
            int label = varlabels[0];
            assert(label <= related_vp_angles.size());
            const double K = param.coeff_vp_line_fitness / nlines; // 50.0
            if (label == related_vp_angles.size()) { // not bind to any vp
              return K;
            }
            double angle = related_vp_angles[label].score;
            assert(!IsInfOrNaN(angle));
            return (1.0 - Gaussian(angle, DegreesToRadians(3))) * K;
          },
          1.0);
      fg.addFactor(fc, {vh});
    }

    // potential 2: two adjacent but not colinear lines should not bind to a
    // near vp
    // first collect all adjacent vh(line) pairs and record whether they are
    // colinear in 2d
    std::set<std::pair<int, int>> noncolineear_adj_vh_pairs;
    for (int face = 0; face < face2ordered_lines.size(); face++) {
      auto &face_lines = face2ordered_lines[face];
      for (int i = 0; i < face_lines.size(); i++) {
        auto vh1 = face_lines[i];
        auto vh2 = face_lines[(i + 1) % face_lines.size()];
        auto &line1 = lines[vh1];
        auto &line2 = lines[vh2];
        if (!LinesAreColinear(
                line1, line2, focal, pp,
                param.angle_thres_judging_colinearility /* / 5.0*/)) {
          noncolineear_adj_vh_pairs.insert(MakeOrderedPair(vh1, vh2));
        }
      }
    }
    double scale = BoundingBoxOfContainer(lines).outerSphere().radius;
    // second add factors
    for (auto &p : noncolineear_adj_vh_pairs) {
      auto vh1 = p.first;
      auto vh2 = p.second;
      auto fc = fg.addFactorCategory(
          [vh1, vh2, &line2ordered_vp_angles, &vps, scale,
           &noncolineear_adj_vh_pairs, focal, pp, &param](
              const std::vector<int> &varlabels, void *givenData) -> double {
            auto &related_vp_angles1 = line2ordered_vp_angles[vh1];
            auto &related_vp_angles2 = line2ordered_vp_angles[vh2];
            assert(varlabels.size() == 2);
            int bindedVP1 = varlabels[0] == related_vp_angles1.size()
                                ? -1
                                : (related_vp_angles1[varlabels[0]].component);
            int bindedVP2 = varlabels[1] == related_vp_angles2.size()
                                ? -1
                                : (related_vp_angles2[varlabels[1]].component);
            if (bindedVP1 == -1 || bindedVP2 == -1) {
              return 0;
            }
            auto &vp_pos1 = vps[bindedVP1];
            auto &vp_pos2 = vps[bindedVP2];
            // const double thres = scale / 10.0;
            const double K = param.coeff_noncolinear_adj_line_exlusiveness /
                             noncolineear_adj_vh_pairs.size(); // 10.0
            // if (Distance(vp_pos1, vp_pos2) < thres) { // todo
            //    return K;
            //  }
            if (AngleBetweenUndirected(
                    NormalizedSpatialDirection(vp_pos1, focal, pp),
                    NormalizedSpatialDirection(vp_pos2, focal, pp)) <
                param.angle_thres_distinguishing_vps) { // todo //
                                                        // DegreesToRadians(2)
              return K;
            }
            return 0.0;
          },
          1.0);
      fg.addFactor(fc, {vh1, vh2});
    }

    // potential 3: the vps of edges sharing a same face should lie on
    // the same line (the vanishing line of the face)
    // first collect all line triplets that should be coplanar in space
    std::vector<std::array<int, 3>> coplanar_line_triplets;
    for (int face = 0; face < face2ordered_lines.size(); face++) {
      auto &face_lines = face2ordered_lines[face];
      assert(face_lines.size() >= 3);
      for (int gap = 1; gap <= (face_lines.size() > 4 ? 2 : 1); gap++) {
        for (int i = 0; i < face_lines.size(); i++) {
          int prev_l =
              face_lines[(i + face_lines.size() - gap) % face_lines.size()];
          int l = face_lines[i];
          int next_l = face_lines[(i + gap) % face_lines.size()];
          coplanar_line_triplets.push_back({{prev_l, l, next_l}});
        }
      }
    }
    // second add factors
    for (auto &triplet : coplanar_line_triplets) {
      auto fc = fg.addFactorCategory(
          [&line2ordered_vp_angles, &vps, &triplet, &coplanar_line_triplets,
           focal, pp, &param](const std::vector<int> &varlabels,
                              void *givenData) -> double {
            assert(varlabels.size() == 3);
            int vp_ids[3] = {-1, -1, -1};
            for (int i = 0; i < 3; i++) {
              if (varlabels[i] == line2ordered_vp_angles[triplet[i]].size()) {
                return 0.0;
              }
              vp_ids[i] =
                  line2ordered_vp_angles[triplet[i]][varlabels[i]].component;
            }
            if (vp_ids[0] == vp_ids[1] || vp_ids[1] == vp_ids[2] ||
                vp_ids[2] == vp_ids[0]) {
              return 0.0;
            }
            Vec3 vpdirs[3];
            for (int i = 0; i < 3; i++) {
              vpdirs[i] = NormalizedSpatialDirection(vps[vp_ids[i]], focal, pp);
            }
            double angle = AngleBetweenUndirected(vpdirs[0].cross(vpdirs[1]),
                                                  vpdirs[1].cross(vpdirs[2]));
            assert(!IsInfOrNaN(angle));
            const double K = param.coeff_line_triplet_coplanar /
                             coplanar_line_triplets.size(); // 30.0
            return (1.0 -
                    Gaussian(angle, param.angle_thres_juding_coplanarity)) *
                   K; // DegreesToRadians(10)
          },
          1.0);
      fg.addFactor(fc, triplet.begin(), triplet.end());
    }

    // solve the factor graph
    result = fg.solve(param.solve_max_iter, 1,
                      [](int epoch, double energy, double denergy,
                         const std::vector<int> &results) -> bool {
                        Println("epoch: ", epoch, "  energy: ", energy);
                        return true;
                      });
  }

  // install the vp labels to line2vp and vp2lines
  std::vector<int> line2vp(nlines, -1);
  std::vector<std::set<int>> vp2lines(nvps);
  for (int l = 0; l < nlines; l++) {
    int id = result[l];
    if (id == line2ordered_vp_angles[l].size()) {
      continue;
    }
    line2vp[l] = line2ordered_vp_angles[l][id].component;
    vp2lines[line2vp[l]].insert(l);
  }

  for (int vp = 0; vp < nvps; vp++) {
    if (vp2lines[vp].size() < param.vp_min_degree) {
      for (int l : vp2lines[vp]) {
        line2vp[l] = -1;
      }
      vp2lines[vp].clear();
    }
  }
  return line2vp;
}

DenseMatd MakePlaneMatrix() { return DenseMatd::eye(3, 3); }
DenseMatd MakePlaneMatrixAlongDirection(const Vec3 &dir) {
  Vec3 u = normalize(dir);
  DenseMatd mat = DenseMatd::zeros(3, 2);
  Vec3 abs_u(abs(u[0]), abs(u[1]), abs(u[2]));
  if (abs_u[2] >= abs_u[1] && abs_u[2] >= abs_u[0]) {
    mat(0, 0) = mat(1, 1) = 1;
    mat(2, 0) = -u[0] / u[2];
    mat(2, 1) = -u[1] / u[2];
  } else if (abs_u[1] >= abs_u[2] && abs_u[1] >= abs_u[0]) {
    mat(0, 0) = mat(2, 1) = 1;
    mat(1, 0) = -u[0] / u[1];
    mat(1, 1) = -u[2] / u[1];
  } else {
    mat(1, 0) = mat(2, 1) = 1;
    mat(0, 0) = -u[1] / u[0];
    mat(0, 1) = -u[2] / u[0];
  }
  return mat;
}
DenseMatd MakePlaneMatrixTowardDirection(const Vec3 &dir) {
  return DenseMatd(normalize(dir));
}

// GenerateInferenceFunctors
std::unique_ptr<Inferencer>
GenerateInferenceFunctors(const std::vector<PlaneConstraint> &constraints,
                          const std::vector<Vec3> &vert2dir,
                          std::vector<int> *fundamental_verts) {
  size_t nverts = vert2dir.size();
  assert(nverts > 0);
  size_t ncons = constraints.size();

  using Eigen::Matrix3Xd;
  using Eigen::MatrixX3d;
  using Eigen::RowVectorXd;
  using Eigen::RowVector3d;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  std::vector<Matrix3Xd> Ps(ncons);
  for (int c = 0; c < ncons; c++) {
    misc::ToEigenMat(constraints[c].P, Ps[c]);
  }

  std::vector<int> cons_order;
  cons_order.reserve(ncons);
  std::unordered_map<int, MatrixXd> cons_fixed2DP;
  std::unordered_map<int, std::vector<int>> cons_fixed2parent_verts;
  std::unordered_map<int, int> vert_fixed2parent_cons;
  vert_fixed2parent_cons[0] = -1; // -1 indicates this is a fundamental vert

  while (vert_fixed2parent_cons.size() < nverts || cons_order.size() < ncons) {
    bool more_cons_fixed = false;
    // find fixable constraints and compute its DP matrix
    for (int c = 0; c < constraints.size(); c++) {
      if (Contains(cons_fixed2DP, c)) {
        continue;
      }
      // collect all related verts fixed
      std::vector<int> related_verts_fixed;
      related_verts_fixed.reserve(constraints[c].verts.size());
      for (int v : constraints[c].verts) {
        if (Contains(vert_fixed2parent_cons, v)) {
          related_verts_fixed.push_back(v);
        }
      }
      if (related_verts_fixed.empty()) {
        continue;
      }
      // build the D matrix from vert dirs
      MatrixX3d D = MatrixX3d::Zero(related_verts_fixed.size(), 3);
      for (int r = 0; r < related_verts_fixed.size(); r++) {
        Vec3 dir = normalize(vert2dir[related_verts_fixed[r]]);
        for (int c = 0; c < 3; c++) {
          D(r, c) = dir[c];
        }
      }

      // D_i * P_i * \\pi_i = \\inv_depths_i
      MatrixXd DP = D * Ps[c];
      assert(DP.cols() == Ps[c].cols());
      Eigen::FullPivLU<MatrixXd> lu_decomp(DP);
      lu_decomp.setThreshold(0.001);
      if (lu_decomp.rank() == Ps[c].cols()) { // full rank, fixable
        cons_fixed2DP[c] = DP;
        cons_fixed2parent_verts[c] = std::move(related_verts_fixed);
        cons_order.push_back(c);
        more_cons_fixed = true;
        // update all its related unfixed verts
        for (int v : constraints[c].verts) {
          if (!Contains(vert_fixed2parent_cons, v)) {
            vert_fixed2parent_cons[v] = c;
          }
        }
      }
    }
    // no more cons are fixed
    if (!more_cons_fixed) {
      // but still has unfixed verts
      if (vert_fixed2parent_cons.size() < nverts) {
        // then make a new fundamental vert
        for (int v = 0; v < nverts; v++) {
          if (!Contains(vert_fixed2parent_cons, v)) {
            vert_fixed2parent_cons[v] = -1;
            break;
          }
        }
      } else {
        break;
      }
    }
  }

  assert(cons_order.size() == ncons);

  if (fundamental_verts) {
    fundamental_verts->clear();
    for (auto &p : vert_fixed2parent_cons) {
      if (p.second == -1) {
        fundamental_verts->push_back(p.first);
      }
    }
  }

  // now let's build the matrices
  std::vector<RowVectorXd> v2matrix(nverts);
  std::vector<MatrixXd> c2matrix(ncons);

  // collect vars
  size_t nvars = 0;
  std::map<int, int> root_vert2pos_in_vars;
  for (auto &p : vert_fixed2parent_cons) {
    if (p.second == -1) {
      root_vert2pos_in_vars[p.first] = nvars;
      nvars++;
    }
  }

  // initialize the matrices of fundamental verts
  for (auto &p : vert_fixed2parent_cons) {
    if (p.second == -1) {
      v2matrix[p.first] = RowVectorXd::Zero(nvars);
      v2matrix[p.first][root_vert2pos_in_vars.at(p.first)] = 1.0;
    }
  }

  // follow the cons order
  for (int cons : cons_order) {
    auto &parent_verts = cons_fixed2parent_verts.at(cons);
    auto &DP = cons_fixed2DP.at(cons);
    // construct the inversed depth matrix of parent verts
    MatrixXd verts_mat = MatrixXd::Zero(parent_verts.size(), nvars);
    for (int r = 0; r < parent_verts.size(); r++) {
      verts_mat.row(r) = v2matrix.at(parent_verts.at(r));
      assert(verts_mat.row(r).norm() != 0.0);
    }
    // D * P * \pi = verts_mat * X
    // \pi = (DP)^-1 * verts_mat * X
    // equation = P * (DP)^-1 * verts_mat * X
    c2matrix[cons] = Ps[cons] * DP.fullPivLu().solve(verts_mat);

    // update the matrices of nonparent verts of this cons
    for (int v : constraints[cons].verts) {
      if (Contains(parent_verts, v)) {
        continue;
      }
      assert(vert_fixed2parent_cons.at(v) == cons);
      // inv_depth = dir^T * equation
      RowVector3d dir(vert2dir.at(v)[0], vert2dir.at(v)[1], vert2dir.at(v)[2]);
      v2matrix[v] = dir.normalized() * c2matrix[cons];
    }
  }

  struct InferImpl : public Inferencer {
    size_t nvariables;
    std::vector<DenseMatd> v2matrix;
    std::vector<DenseMatd> c2matrix;
    virtual size_t nvars() const override { return nvariables; }
    virtual Vec3 getPlaneEquation(int cons,
                                  const DenseMatd &variables) const override {
      DenseMatd eq = c2matrix.at(cons) * variables;
      return Vec3(eq(0, 0), eq(1, 0), eq(2, 0));
    }
    virtual double getInversedDepth(int vert,
                                    const DenseMatd &variables) const override {
      DenseMatd inv_depth = v2matrix.at(vert) * variables;
      return inv_depth(0, 0);
    }
    virtual DenseMatd recoverVariables(
        const std::vector<double> &vert2inversed_depths) const override {
      MatrixXd vmat =
          MatrixXd::Zero(vert2inversed_depths.size(), this->nvars());
      for (int v = 0; v < vert2inversed_depths.size(); v++) {
        RowVectorXd r;
        misc::ToEigenMat(v2matrix.at(v), r);
        vmat.row(v) = r;
      }
      VectorXd vars = vmat.fullPivLu().solve(VectorXd::Map(
          vert2inversed_depths.data(), vert2inversed_depths.size()));
      return misc::ToCVMat(vars);
    }
  };

  auto infer_ptr = std::make_unique<InferImpl>();
  InferImpl &infer = *infer_ptr;
  infer.v2matrix.resize(nverts);
  for (int v = 0; v < nverts; v++) {
    infer.v2matrix[v] = misc::ToCVMat(v2matrix.at(v));
    assert(!infer.v2matrix.at(v).empty());
  }
  infer.c2matrix.resize(ncons);
  for (int c = 0; c < ncons; c++) {
    infer.c2matrix[c] = misc::ToCVMat(c2matrix.at(c));
    assert(!infer.c2matrix.at(c).empty());
  }
  infer.nvariables = nvars;

  return infer_ptr;
}

std::vector<double>
AnglesBetweenAdjacentEdges(const std::vector<Vec3> &vert2dir,
                           const std::vector<std::vector<int>> &face2verts,
                           const DenseMatd &variables, const Inferencer &infer,
                           std::function<bool(int face)> face_selected) {
  std::vector<double> edge_angles;
  std::vector<Point3> vert2pos(vert2dir.size());
  for (int v = 0; v < vert2dir.size(); v++) {
    vert2pos[v] = normalize(vert2dir[v]) / infer.getInversedDepth(v, variables);
  }
  for (int face = 0; face < face2verts.size(); face++) {
    if (face_selected && !face_selected(face)) {
      continue;
    }
    auto & vs = face2verts[face];
    assert(vs.size() >= 3);
    for (int i = 0; i < vs.size(); i++) {
      int v1 = vs[i];
      int v2 = vs[(i + 1) % vs.size()];
      int v3 = vs[(i + 2) % vs.size()];
      double angle = AngleBetweenDirected(vert2pos[v1] - vert2pos[v2],
                                          vert2pos[v3] - vert2pos[v2]);
      edge_angles.push_back(angle);
    }
  }
  return edge_angles;
}

std::vector<double>
AnglesBetweenAdjacentFaces(size_t nfaces,
                           const std::vector<std::vector<int>> &edge2faces,
                           const DenseMatd &variables, const Inferencer &infer,
                           std::function<bool(int face)> face_selected) {
  std::vector<double> face_angles;
  std::vector<Vec3> face_equations(nfaces);
  for (int face = 0; face < nfaces; face++) {
    face_equations[face] = infer.getPlaneEquation(face, variables);
  }
  for (auto &fs : edge2faces) {
    if (fs.size() < 2) {
      continue;
    }
    for (int i = 0; i < fs.size(); i++) {
      if (face_selected && !face_selected(fs[i])) {
        continue;
      }
      for (int j = i + 1; j < fs.size(); j++) {
        if (face_selected && !face_selected(fs[j])) {
          continue;
        }
        double angle = AngleBetweenUndirected(face_equations[fs[i]],
                                              face_equations[fs[j]]);
        face_angles.push_back(angle);
      }
    }
  }
  return face_angles;
}
}
}