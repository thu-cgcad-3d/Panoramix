#include <QtWidgets>

#include "../../src/core/factor_graph.hpp"
#include "../../src/core/parallel.hpp"
#include "../../src/misc/cache.hpp"
#include "../../src/misc/clock.hpp"

#include "../../src/gui/canvas.hpp"
#include "../../src/gui/qttools.hpp"
#include "../../src/gui/scene.hpp"
#include "../../src/gui/singleton.hpp"
#include "../../src/gui/utility.hpp"

#include "../../src/experimental/line_drawing.hpp"
#include "../../src/experimental/pi_graph_annotation.hpp"
#include "../../src/experimental/pi_graph_cg.hpp"
#include "../../src/experimental/pi_graph_control.hpp"
#include "../../src/experimental/pi_graph_occlusion.hpp"
#include "../../src/experimental/pi_graph_optimize.hpp"
#include "../../src/experimental/pi_graph_vis.hpp"

using namespace pano;
using namespace pano::core;
using namespace pano::experimental;

#define DISABLED_main MACRO_CONCAT(main_, __COUNTER__)

// VHFSelectedFunT: (VertHandle/HalfHandle/FaceHandle) -> bool
// VertPositionFunT: (VertHandle)->Point3
// HalfColorerFunT: (HalfHandle)->gui::Color
// FaceColorerFunT: (FaceHandle)->gui::Color
template <class VertDataT, class HalfDataT, class FaceDataT,
          class VHFSelectedFunT, class VertPositionFunT, class HalfColorerFunT,
          class FaceColorerFunT>
void AddToScene(gui::SceneBuilder &sb,
                const Mesh<VertDataT, HalfDataT, FaceDataT> &m,
                VHFSelectedFunT selected, VertPositionFunT vertPosFun,
                HalfColorerFunT colorHalf, FaceColorerFunT colorFace) {
  sb.installingOptions().lineWidth = 1.0;
  sb.installingOptions().defaultShaderSource =
      gui::OpenGLShaderSourceDescriptor::XLines;
  HandledTable<HalfHandle, int> added(m.internalHalfEdges().size(), false);
  for (auto &h : m.halfedges()) {
    if (!selected(h.topo.hd) || !selected(h.topo.from()) ||
        !selected(h.topo.to())) {
      continue;
    }
    Line3 line(vertPosFun(m.data(h.topo.from())),
               vertPosFun(m.data(h.topo.to())));
    auto hh = h.topo.hd;
    auto fh = h.topo.face;
    auto oppohh = h.topo.opposite;
    auto oppofh = oppohh.valid() ? m.topo(oppohh).face : FaceHandle();
    bool hasFace = fh.valid();
    bool hasOppo = oppohh.valid();
    gui::Color color = colorHalf(hh);
    if (color.isTransparent()) {
      continue;
    }
    if (!hasFace && hasOppo) {
      color = gui::Red;
    } else if (hasFace && !hasOppo) {
      color = gui::Blue;
    } else if (!hasFace && !hasOppo) {
      color = gui::Yellow;
    }
    if (added[oppohh]) {
      continue;
    }
    sb.add(gui::ColorAs(line, color), [hh, fh, oppohh, oppofh](auto &...) {
      std::cout << "halfedge id: " << hh.id
                << ", opposite halfedge id: " << oppohh.id
                << ", face id: " << fh.id << ", opposite face id: " << oppofh.id
                << '\n';
    });
    added[hh] = true;
  }
  sb.installingOptions().defaultShaderSource =
      gui::OpenGLShaderSourceDescriptor::XTriangles;
  for (auto &f : m.faces()) {
    if (!selected(f.topo.hd)) {
      continue;
    }
    Polygon3 poly;
    for (auto h : f.topo.halfedges) {
      auto v = m.topo(h).to();
      poly.corners.push_back(vertPosFun(m.data(v)));
    }
    assert(poly.corners.size() > 2);
    poly.normal = (poly.corners[0] - poly.corners[1])
                      .cross(poly.corners[0] - poly.corners[2]);
    auto fh = f.topo.hd;
    gui::Color color = colorFace(f.topo.hd);
    if (color.isTransparent()) {
      continue;
    }
    sb.add(gui::ColorAs(poly, color),
           [fh](auto &...) { std::cout << "face id: " << fh.id << '\n'; });
  }
}

inline std::vector<Point2> PossibleKeyVanishingPoints(const Chain2 &chain) {
  assert(chain.size() > 2);
  if (chain.size() == 3) {
    return {};
  }
  if (chain.size() == 4) {
    return {Intersection(chain.edge(0).ray(), chain.edge(2).ray()),
            Intersection(chain.edge(1).ray(), chain.edge(3).ray())};
  }
  if (chain.size() % 2 == 0) {
    std::vector<Point2> vpPositions;
    vpPositions.reserve(chain.size() * chain.size() / 4);
    for (int i = 0; i < chain.size() / 2; i++) {
      for (int j = i + 1; j < i + chain.size() / 2; j++) {
        vpPositions.push_back(
            Intersection(Line2(chain.at(i), chain.at(j)).ray(),
                         Line2(chain.at(i + chain.size() / 2),
                               chain.at(j + chain.size() / 2))
                             .ray()));
      }
    }
    return vpPositions;
  } else {
    std::vector<Point2> vpPositions;
    for (int i = 0; i < chain.size(); i++) {
      for (int j = i + 1; j < chain.size(); j++) {
        vpPositions.push_back(
            Intersection(chain.edge(i).ray(), chain.edge(j).ray()));
      }
    }
    return vpPositions;
  }
}


struct EnergyWeights {
  double vertMSDAWeight;
  double faceMSDAWeight;
  double faceAngleWeight;
  double allMSDAWeight;
};

template <class VT, class HT, class FT, class FaceHandleIterT,
          class VertHandleToDirFunT, class FaceHandleToPlaneEqFunT>
std::vector<double>
ComputeEnergy(const Mesh<VT, HT, FT> &mesh, FaceHandleIterT fhs_begin,
              FaceHandleIterT fhs_end, const EnergyWeights &weights,
              VertHandleToDirFunT vh2dir, FaceHandleToPlaneEqFunT fh2planeEq,
              const Point3 &eye = Point3()) {
  using namespace Eigen;
  std::vector<double> energyTerms;

  // fh2planeEq: (FaceHandle)[i] -> double
  std::map<VertHandle, Point3> vpositions;
  std::map<VertHandle, int> vfacedegrees;
  //for (VertHandle vh : sub.vhs) {
  //  vpositions[vh] = Vec3();
  //  vfacedegrees[vh] = 0;
  //}
  for (FaceHandle fh : MakeRange(fhs_begin, fhs_end)) {
    decltype(auto) planeEq = fh2planeEq(fh);
    Plane3 plane = Plane3FromEquation(planeEq[0], planeEq[1], planeEq[2]);
    if(HasValue(plane, IsInfOrNaN<double>)) {
		return {};
    }
    for (HalfHandle hh : mesh.topo(fh).halfedges) {
      VertHandle vh = mesh.topo(hh).to();
	  Vec3 dir = normalize(vh2dir(vh));
      Point3 p3d = Intersection(Ray3(eye, dir), plane);
      vpositions[vh] += p3d;
      vfacedegrees[vh]++;
    }
  }
  for (auto vhpos : vpositions) {
    vhpos.second /= vfacedegrees.at(vhpos.first);
  }

  std::vector<double> allAngles;
  std::map<VertHandle, std::vector<double>> vh2angles;
  for (FaceHandle fh : MakeRange(fhs_begin, fhs_end)) {
    auto &loophhs = mesh.topo(fh).halfedges;
    std::vector<double> faceAngles(loophhs.size(), 0.0);
    for (int k = 0; k < loophhs.size(); k++) {
      int knext = (k + 1) % loophhs.size();
      HalfHandle hh1 = loophhs[k];
      HalfHandle hh2 = loophhs[knext];
      assert(mesh.topo(hh1).to() == mesh.topo(hh2).from());
      VertHandle v1 = mesh.topo(hh1).from();
      VertHandle v2 = mesh.topo(hh1).to();
      VertHandle v3 = mesh.topo(hh2).to();
      auto &p1 = vpositions.at(v1);
      auto &p2 = vpositions.at(v2);
      auto &p3 = vpositions.at(v3);
      double angle = AngleBetweenDirected(p1 - p2, p3 - p2);
      faceAngles[k] = angle;
      vh2angles[v2].push_back(angle);
	  allAngles.push_back(angle);
    }
    double faceMeanAngle =
        std::accumulate(faceAngles.begin(), faceAngles.end(), 0.0) /
        faceAngles.size();
    for (double a : faceAngles) {
      energyTerms.push_back((a - faceMeanAngle) * weights.faceMSDAWeight / faceAngles.size());
    }
  }

  for (auto &vhangles : vh2angles) {
    VertHandle vh = vhangles.first;
    auto &angles = vhangles.second;
    double vertMeanAngle =
        std::accumulate(angles.begin(), angles.end(), 0.0) / angles.size();
    for (double a : angles) {
      energyTerms.push_back((a - vertMeanAngle) * weights.vertMSDAWeight /
                            angles.size());
    }
  }

  double allMeanAngle =
      std::accumulate(allAngles.begin(), allAngles.end(), 0.0) /
      allAngles.size();
  for (double a : allAngles) {
    energyTerms.push_back((a - allMeanAngle) * weights.allMSDAWeight /
                          allAngles.size());
  }

  std::vector<HalfHandle> validHhs;
  for (auto &h : mesh.halfedges()) {
    auto hh = h.topo.hd;
    auto fh1 = mesh.topo(hh).face;
    auto fh2 = mesh.topo(mesh.topo(hh).opposite).face;
    if (std::find(fhs_begin, fhs_end, fh1) != fhs_end &&
        std::find(fhs_begin, fhs_end, fh2) != fhs_end) {
      validHhs.push_back(hh);
    }
  }

  for (auto &hh : validHhs) {
    auto fh1 = mesh.topo(hh).face;
    auto planeEq1 = fh2planeEq(fh1);
    auto fh2 = mesh.topo(mesh.topo(hh).opposite).face;
    auto planeEq2 = fh2planeEq(fh2);
    double angle =
        AngleBetweenUndirected(Vec3(planeEq1[0], planeEq1[1], planeEq1[2]),
                               Vec3(planeEq2[0], planeEq2[1], planeEq2[2]));
    assert(!IsInfOrNaN(angle));
    energyTerms.push_back((angle - M_PI_2) * weights.faceAngleWeight /
                          validHhs.size());
  }

  return energyTerms;
}

// ReconstructWithOrientations
Mesh3 ReconstructWithOrientations(
    const std::vector<Line2> &edge2line, const std::vector<int> &edge2vp,
    const std::vector<Vec3> &vp2dir,
    const HandledTable<HalfHandle, int> &hh2d2edge,
    const PerspectiveCamera &curCam,
    const Mesh<VertHandle, HalfHandle, FaceHandle> &meshProxy,
    const Mesh2 &mesh2d, misc::Matlab &matlab) {

  int nedges = edge2line.size();
  assert(edge2vp.size() == nedges);

  // get the vh2d2reconstructedPosition data
  // entity
  using SupportingPlane = PIConstraintGraph::Entity::SupportingPlane;
  struct EntityBase {
    SupportingPlane supportingPlane;
    EntityBase(const SupportingPlane &sp) : supportingPlane(sp) {}
    virtual ~EntityBase() {}
    virtual FaceHandle fhProxy() const { return FaceHandle(); }
    virtual int edge() const { return -1; }
    bool isFace() const { return fhProxy().valid(); }
    bool isEdge() const { return edge() != -1; }
  };
  struct FaceEntity : EntityBase {
    FaceHandle fh;
    FaceEntity(FaceHandle fh, const Vec3 &center)
        : EntityBase(SupportingPlane(SegControl{-1, -1}, center, {})), fh(fh) {}
    virtual ~FaceEntity() {}
    virtual FaceHandle fhProxy() const override { return fh; }
  };
  struct EdgeEntity : EntityBase {
    int e;
    EdgeEntity(int e, const Classified<Line3> &line,
               const std::vector<Vec3> &vp2dir)
        : EntityBase(SupportingPlane(line, vp2dir)), e(e) {}
    virtual ~EdgeEntity() {}
    virtual int edge() const override { return e; }
  };

  std::vector<std::unique_ptr<EntityBase>> entities;

  std::vector<int> edge2ent(nedges, -1);
  HandledTable<FaceHandle, int> fhProxy2ent(meshProxy.internalFaces().size(),
                                            -1);

  // install supporting planes from edges and faces
  // from edges
  for (int edge = 0; edge < nedges; edge++) {
    auto &line2 = edge2line[edge];
    Line3 line3(normalize(curCam.direction(line2.first)),
                normalize(curCam.direction(line2.second)));
    int vp = edge2vp[edge];
    entities.push_back(
        std::make_unique<EdgeEntity>(edge, ClassifyAs(line3, vp), vp2dir));
    int ent = entities.size() - 1;
    edge2ent[edge] = ent;
  }
  // from faces
  HandledTable<FaceHandle, Point2> fhProxy2center2d(
      meshProxy.internalFaces().size());
  for (auto &f : meshProxy.faces()) {
    Point2 center2 = Origin<2>();
    for (auto hh : f.topo.halfedges) {
      center2 += mesh2d.data(meshProxy.data(meshProxy.topo(hh).to()));
    }
    center2 /= double(f.topo.halfedges.size());
    fhProxy2center2d[f.topo.hd] = center2;
    entities.push_back(std::make_unique<FaceEntity>(
        f.topo.hd, normalize(curCam.direction(center2))));
    int ent = entities.size() - 1;
    fhProxy2ent[f.topo.hd] = ent;
  }

  //// install connections
  // std::map<std::pair<int, int>, std::vector<Vec3>> ents2anchors;
  // for (auto &f : meshProxy.faces()) {
  //  for (auto hhProxy : f.topo.halfedges) {
  //    HalfHandle hh2d = meshProxy.data(hhProxy);
  //    if (hh2d.invalid()) {
  //      hh2d = meshProxy.data(meshProxy.topo(hhProxy).opposite);
  //    }
  //    assert(hh2d.valid());
  //    int edge = hh2d2edge[hh2d];

  //    auto &line2 = edge2line[edge];
  //    ents2anchors[std::make_pair(edge2ent[edge], fhProxy2ent[f.topo.hd])]
  //    =
  //        {normalize(curCam.direction(line2.first)),
  //         normalize(curCam.direction(line2.second))};
  //  }
  //}

  //// install face angles constraints
  // std::vector<std::pair<int, int>> adjFaceEnts;
  // std::vector<int> adjFace2SubMeshId;
  // std::vector<bool> adjFaceOverlapInView;
  // for (auto &half : meshProxy.halfedges()) {
  //  FaceHandle fh1 = half.topo.face;
  //  FaceHandle fh2 = meshProxy.topo(half.topo.opposite).face;

  //  bool onSameSide = IsOnLeftSide(faceCorners[0], line.first,
  //  line.second) ==
  //                    IsOnLeftSide(faceCorners[1], line.first,
  //                    line.second);

  //  adjFaceOverlapInView.push_back(onSameSide);
  //}
  // int nadjFaces = adjFaceEnts.size();

  // Start Building Matrices
  // vert start position in variable vector
  std::vector<int> ent2varPosition(entities.size(), -1);
  std::vector<int> ent2nvar(entities.size(), -1);
  std::vector<DenseMatd> ent2matFromVarToPlaneCoeffs(entities.size());
  std::vector<int> var2ent;
  int nvars = 0;
  for (int ent = 0; ent < entities.size(); ent++) {
    auto &e = *entities[ent];
    int nvar = e.supportingPlane.dof;
    ent2nvar[ent] = nvar;
    ent2varPosition[ent] = nvars;
    ent2matFromVarToPlaneCoeffs[ent] =
        e.supportingPlane.matFromVarsToPlaneCoeffs();
    var2ent.insert(var2ent.end(), (size_t)nvar, ent);
    nvars += nvar;
  }

  // P : [(3 * nents) x nvars]
  // P * X -> plane coefficients
  int nents = entities.size();
  std::vector<SparseMatElementd> Ptriplets;
  for (int ent = 0; ent < nents; ent++) {
    int varpos = ent2varPosition[ent];
    // [3 x k]
    auto &matFromVarsToPlaneCoeffs = ent2matFromVarToPlaneCoeffs[ent];
    assert(matFromVarsToPlaneCoeffs.rows == 3 &&
           matFromVarsToPlaneCoeffs.cols == ent2nvar[ent]);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < ent2nvar[ent]; j++) {
        assert(varpos + j < nvars);
        Ptriplets.emplace_back(i + ent * 3, varpos + j,
                               matFromVarsToPlaneCoeffs(i, j));
      }
    }
  }

  // anchorPositions: [nanchors x 3]
  int nanchors = 0;
  std::vector<SparseMatElementd> anchorPositions;
  HandledTable<VertHandle, int> vhProxy2anchor(
      meshProxy.internalVertices().size(), -1);
  for (auto &vert : meshProxy.vertices()) {
    auto vh2d = vert.data;
    Vec3 anchor = normalize(curCam.direction(mesh2d.data(vh2d)));
    for (int i = 0; i < 3; i++) {
      anchorPositions.emplace_back(nanchors, i, anchor[i]);
    }
    vhProxy2anchor[vert.topo.hd] = nanchors;
    nanchors++;
  }

  // connection2anchor: [nconnection -> nanchors],
  // connection2leftEnt, connection2rightEnt:
  // [nconnection -> nents]
  std::vector<double> connection2anchor, connection2leftEnt,
      connection2rightEnt;
  // angleHead2anchor, angleTailA2anchor, angleTailB2anchor
  std::vector<double> angleHead2anchor, angleTailA2anchor, angleTailB2anchor;
  // angle2faceEnt
  std::vector<double> angle2faceEnt;

  int nconnections = 0;
  int nangles = 0;
  for (auto &f : meshProxy.faces()) {
    auto &hhProxies = f.topo.halfedges;
    for (int i = 0; i < hhProxies.size(); i++) {
      auto hhProxy = hhProxies[i];
      HalfHandle hh2d = meshProxy.data(hhProxy);
      if (hh2d.invalid()) {
        hh2d = meshProxy.data(meshProxy.topo(hhProxy).opposite);
      }
      assert(hh2d.valid());
      int edge = hh2d2edge[hh2d];
      int edgeEnt = edge2ent[edge];
      int faceEnt = fhProxy2ent[f.topo.hd];

      VertHandle curVhProxy = meshProxy.topo(hhProxy).to();
      int curAnchorId = vhProxy2anchor[curVhProxy];
      int prevAnchorId = vhProxy2anchor[meshProxy.topo(hhProxy).from()];
      assert(meshProxy
                 .topo(hhProxies[(i + hhProxies.size() - 1) % hhProxies.size()])
                 .to() == meshProxy.topo(hhProxy).from());
      int nextAnchorId =
          vhProxy2anchor[meshProxy.topo(hhProxies[(i + 1) % hhProxies.size()])
                             .to()];

      connection2leftEnt.push_back(edgeEnt);
      connection2rightEnt.push_back(faceEnt);
      connection2anchor.push_back(curAnchorId);
      connection2leftEnt.push_back(edgeEnt);
      connection2rightEnt.push_back(faceEnt);
      connection2anchor.push_back(prevAnchorId);

      angleHead2anchor.push_back(curAnchorId);
      angleTailA2anchor.push_back(prevAnchorId);
      angleTailB2anchor.push_back(nextAnchorId);
      angle2faceEnt.push_back(faceEnt);

      nconnections += 2;
      nangles++;
    }
  }
  assert(angle2faceEnt.size() == nangles);

  // solve!
  matlab << "clear()";

  // matrices
  matlab.setVar("P", MakeSparseMatFromElements(
                         3 * nents, nvars, Ptriplets.begin(), Ptriplets.end()));
  matlab.setVar("AnchorPositions",
                MakeSparseMatFromElements(nanchors, 3, anchorPositions.begin(),
                                          anchorPositions.end()));
  // index mappings
  matlab.setVar("connection2anchor", DenseMatd(connection2anchor));
  matlab.setVar("connection2leftEnt", DenseMatd(connection2leftEnt));
  matlab.setVar("connection2rightEnt", DenseMatd(connection2rightEnt));

  matlab.setVar("angleHead2anchor", DenseMatd(angleHead2anchor));
  matlab.setVar("angleTailA2anchor", DenseMatd(angleTailA2anchor));
  matlab.setVar("angleTailB2anchor", DenseMatd(angleTailB2anchor));
  matlab.setVar("angle2faceEnt", DenseMatd(angle2faceEnt));

  // sizes
  matlab.setVar("nvars", nvars);
  matlab.setVar("nents", nents);
  matlab.setVar("nanchors", nanchors);
  matlab.setVar("nconnections", nconnections);
  matlab.setVar("nangles", nangles);

  double minE = std::numeric_limits<double>::infinity();
  matlab << "FinalX = zeros(nvars, 1);";
  static const int maxIter = 10;

  matlab << "Prev_InverseDepthsOnLeftOfConnections = ones(nconnections, 1);";
  matlab << "Prev_InverseDepthsOnRightOfConnections = ones(nconnections, 1);";

  matlab << "Prev_Scalar1 = ones(nangles, 1);";
  matlab << "Prev_Scalar2 = ones(nangles, 1);";
  matlab << "Prev_Scalar3 = ones(nangles, 1);";
  // matlab << "Prev_AvgAngleCos = zeros(nangles, 1);";

  matlab.setPrintMessage(false);
  for (int t = 0; t < maxIter; t++) {
    matlab << "cvx_begin quiet";
    matlab << "variable X(nvars);";

    auto constructProblem = [&matlab](const std::string &objectiveVarName) {
      matlab << "PlaneEqs = reshape(P * X, [3 nents])';"; // [nents x 3]

      matlab << "PlaneEqsOnLeftOfConnections = "
                "PlaneEqs(connection2leftEnt + 1, :);";
      matlab << "PlaneEqsOnRightOfConnections = "
                "PlaneEqs(connection2rightEnt + 1, :);";

      matlab << "ConnectionPositions = "
                "AnchorPositions(connection2anchor + 1, :);";
      matlab << "InverseDepthsOnLeftOfConnections = "
                "dot(PlaneEqsOnLeftOfConnections, ConnectionPositions, 2);";
      matlab << "InverseDepthsOnRightOfConnections = "
                "dot(PlaneEqsOnRightOfConnections, ConnectionPositions, 2);";

      // energy of connections
      matlab << "EnergyOfConnections = "
                "sum_square((InverseDepthsOnLeftOfConnections - "
                "InverseDepthsOnRightOfConnections) ./ "
                "Prev_InverseDepthsOnLeftOfConnections ./ "
                "Prev_InverseDepthsOnRightOfConnections);";
      const std::string trueEnergyOfConnectionsExpr =
          "sum_square((InverseDepthsOnLeftOfConnections - "
          "InverseDepthsOnRightOfConnections) ./ "
          "InverseDepthsOnLeftOfConnections ./ "
          "InverseDepthsOnRightOfConnections)";

      matlab << "PlaneEqsOfAngles = PlaneEqs(angle2faceEnt + 1, :);";

      matlab << "AngleHeadPositions = "
                "AnchorPositions(angleHead2anchor + 1, :);";
      matlab << "AngleTailAPositions = "
                "AnchorPositions(angleTailA2anchor + 1, :);";
      matlab << "AngleTailBPositions = "
                "AnchorPositions(angleTailB2anchor + 1, :);";

      // m : AngleTailAPositions
      // n : AngleTailBPositions
      matlab << "Scalar1 = dot(AngleHeadPositions, PlaneEqsOfAngles, 2);";
      matlab << "Scalar2 = dot(AngleTailBPositions, PlaneEqsOfAngles, 2);";
      matlab << "Scalar3 = dot(AngleTailAPositions, PlaneEqsOfAngles, 2);";

      // px : AngleHeadPositions(:, 1)
      // py : AngleHeadPositions(:, 2)
      // pz : AngleHeadPositions(:, 3)

      // mx : AngleTailAPositions(:, 1)
      // my : AngleTailAPositions(:, 2)
      // mz : AngleTailAPositions(:, 3)

      // nx : AngleTailBPositions(:, 1)
      // ny : AngleTailBPositions(:, 2)
      // nz : AngleTailBPositions(:, 3)

      // / mx   px \ / nx   px \   / my   py \ / ny   py \   / mz   pz \ / nz   pz \
        // | -- - -- | | -- - -- | + | -- - -- | | -- - -- | + | -- - -- | | -- - -- |
      // \ #3   #1 / \ #2   #1 /   \ #3   #1 / \ #2   #1 /   \ #3   #1 / \
          // #2
      // #1 /
      matlab << "AngleCosineNumerator1 = "
                "(Scalar1 .* AngleTailAPositions(:, 1) - "
                " Scalar3 .* AngleHeadPositions(:, 1)) .*"
                "(Prev_Scalar1 .* AngleTailBPositions(:, 1) - "
                "Prev_Scalar2 .* AngleHeadPositions(:, 1));";
      matlab << "AngleCosineNumerator2 = "
                "(Scalar1 .* AngleTailAPositions(:, 2) - "
                " Scalar3 .* AngleHeadPositions(:, 2)) .*"
                "(Prev_Scalar1 .* AngleTailBPositions(:, 2) - "
                "Prev_Scalar2 .* AngleHeadPositions(:, 2));";
      matlab << "AngleCosineNumerator3 = "
                "(Scalar1 .* AngleTailAPositions(:, 3) - "
                " Scalar3 .* AngleHeadPositions(:, 3)) .*"
                "(Prev_Scalar1 .* AngleTailBPositions(:, 3) - "
                "Prev_Scalar2 .* AngleHeadPositions(:, 3));";

      matlab << "AngleCosines = "
                "(AngleCosineNumerator1 + AngleCosineNumerator2 + "
                "AngleCosineNumerator3) "
                "./ Prev_Scalar1 ./ Prev_Scalar1 "
                "./ Prev_Scalar2 ./ Prev_Scalar3;";

      // restrict face angles to be orhtogonal, aka, cosine values to be
      // zero
      matlab << "EnergyOfOrthogonalFaceAngles = sum_square(AngleCosines);";

      // energy
      matlab << (objectiveVarName + " = 1e3 * EnergyOfConnections + "
                                    "EnergyOfOrthogonalFaceAngles;");
    };

    constructProblem("ObjectiveEnergy");
    matlab << "minimize 1e3 * ObjectiveEnergy";
    matlab << "subject to";
    matlab << "   ones(nconnections, 1) <= "
              "InverseDepthsOnLeftOfConnections;";
    matlab << "   ones(nconnections, 1) <= "
              "InverseDepthsOnRightOfConnections;";
    matlab << "cvx_end";

    matlab << "Prev_InverseDepthsOnLeftOfConnections = "
              "InverseDepthsOnLeftOfConnections ./ "
              "norm(InverseDepthsOnLeftOfConnections);";
    matlab << "Prev_InverseDepthsOnRightOfConnections = "
              "InverseDepthsOnRightOfConnections ./ "
              "norm(InverseDepthsOnRightOfConnections);";
    matlab << "Prev_Scalar1 = Scalar1;";
    matlab << "Prev_Scalar2 = Scalar2;";
    matlab << "Prev_Scalar3 = Scalar3;";

    constructProblem("TrueObjectiveEnergy");

    double curEnergy = matlab.var("TrueObjectiveEnergy");
    Println("cur energy - ", curEnergy);
    if (IsInfOrNaN(curEnergy)) {
      break;
    }
    if (curEnergy < minE) {
      minE = curEnergy;
      matlab << "FinalX = 2 * X ./ median(InverseDepthsOnLeftOfConnections + "
                "InverseDepthsOnRightOfConnections);";
      matlab << "norm(FinalX)";
    }
  }
  matlab.setPrintMessage(true);

  // use the solved X to recover supporting planes in each ent
  matlab << "PlaneEqsVector = P * FinalX;"; // [nents*3 x 1]
  DenseMatd planeEqsVector = matlab.var("PlaneEqsVector");
  for (int ent = 0; ent < nents; ent++) {
    auto &entity = *entities[ent];
    entity.supportingPlane.reconstructed = Plane3FromEquation(
        planeEqsVector(ent * 3 + 0), planeEqsVector(ent * 3 + 1),
        planeEqsVector(ent * 3 + 2));
  }

  if (true) {
    // show reconstructed edge lines directly
    gui::SceneBuilder sb;
    sb.installingOptions().defaultShaderSource =
        gui::OpenGLShaderSourceDescriptor::XLines;
    sb.installingOptions().discretizeOptions.color(gui::Black);
    sb.installingOptions().lineWidth = 10;
    for (int edge = 0; edge < nedges; edge++) {
      auto &line2 = edge2line[edge];
      auto &plane = entities[edge2ent[edge]]->supportingPlane.reconstructed;
      Line3 line3(
          Intersection(Ray3(curCam.eye(), curCam.direction(line2.first)),
                       plane),
          Intersection(Ray3(curCam.eye(), curCam.direction(line2.second)),
                       plane));
      sb.add(line3);
    }
    sb.show(true, true);
  }

  // reconstruct vertex positions
  Mesh3 curReconstruction =
      Transform(mesh2d, [](const Point2 &) { return Origin(); });
  HandledTable<VertHandle, int> vh2d2faceProxyCount(
      mesh2d.internalVertices().size(), 0);
  for (auto &faceProxy : meshProxy.faces()) {
    auto fhProxy = faceProxy.topo.hd;
    int ent = fhProxy2ent[fhProxy];
    const Plane3 &plane = entities[ent]->supportingPlane.reconstructed;
    for (auto hhProxy : faceProxy.topo.halfedges) {
      VertHandle vhProxy = meshProxy.topo(hhProxy).to();
      VertHandle vh2d = meshProxy.data(vhProxy);
      Vec3 direction = normalize(curCam.direction(mesh2d.data(vh2d)));
      Point3 position = Intersection(Ray3(curCam.eye(), direction), plane);
      curReconstruction.data(vh2d) += position;
      vh2d2faceProxyCount[vh2d]++;
    }
  }
  for (auto &vert : curReconstruction.vertices()) {
    vert.data /= double(vh2d2faceProxyCount[vert.topo.hd]);
  }

  if (true) { // show current reconstruction
    gui::SceneBuilder sb;
    sb.installingOptions().defaultShaderSource =
        gui::OpenGLShaderSourceDescriptor::XTriangles;
    sb.installingOptions().discretizeOptions.color(gui::Black);
    sb.installingOptions().lineWidth = 10;
    AddToScene(sb, curReconstruction, [](auto) { return true; },
               [](const Point3 &pos) { return pos; },
               [](HalfHandle hh) { return gui::Black; },
               [](FaceHandle fh) {
                 return gui::Color(
                     Vec3i(rand() % 256, rand() % 256, rand() % 256));
               });
    sb.show(true, true, gui::RenderOptions()
                            .backgroundColor(gui::White)
                            .renderMode(gui::All)
                            .bwTexColor(0.0)
                            .bwColor(1.0)
                            .fixUpDirectionInCameraMove(false)
                            .cullBackFace(false)
                            .cullFrontFace(false));
  }
  return curReconstruction;
}

// OptimizeWithoutOrientations
void OptimizeWithoutOrientations(
    Mesh3 &cur_reconstruction, const PerspectiveCamera &cur_cam,
    const Mesh<VertHandle, HalfHandle, FaceHandle> &mesh_proxy,
	const std::vector<SubMesh> &sub_meshs) {
  using namespace Eigen;

   { // show together
    auto color_table = gui::CreateRandomColorTableWithSize(sub_meshs.size());
    gui::SceneBuilder sb;
    sb.installingOptions().defaultShaderSource =
        gui::OpenGLShaderSourceDescriptor::XLines;
    sb.installingOptions().lineWidth = 10;
    for (int i = 0; i < sub_meshs.size(); i++) {
      sb.installingOptions().discretizeOptions.color(color_table[i]);
      for (HalfHandle hh : sub_meshs[i].hhs) {
        VertHandle vh1 = mesh_proxy.topo(hh).from();
        VertHandle vh2 = mesh_proxy.topo(hh).to();
        Point3 p1 = cur_reconstruction.data(mesh_proxy.data(vh1));
        Point3 p2 = cur_reconstruction.data(mesh_proxy.data(vh2));
        sb.add(Line3(p1, p2));
      }
    }
    sb.show(true, true, gui::RenderOptions()
                            .winName("Together Before Optimization")
                            .fixUpDirectionInCameraMove(false));
  }


  std::vector<std::map<VertHandle, double>> vh_proxy2inversed_depth(
      sub_meshs.size());

  // optimize in each submesh first
  for (int i = 0; i < sub_meshs.size(); i++) {
    auto &sub = sub_meshs[i];
    std::vector<VertHandle> fundamental_vhs(sub.fundamental_vhs.begin(),
                                            sub.fundamental_vhs.end());
    std::map<VertHandle, RowVectorXd> v2matrix;
    std::map<FaceHandle, Matrix3Xd> f2matrix;

    auto inversed_depths = EstimatePerspectiveDepths(
        fundamental_vhs, sub.ordered_fhs,
        // vert2initial_depth
        [&mesh_proxy, &cur_reconstruction,
         &cur_cam](VertHandle vh_proxy) -> double {
          return Distance(cur_reconstruction.data(mesh_proxy.data(vh_proxy)),
                          cur_cam.eye());
        },
        // face2related_verts
        [&mesh_proxy](FaceHandle fh_proxy) -> std::vector<VertHandle> {
          std::vector<VertHandle> related_vh_proxies;
          for (HalfHandle hh : mesh_proxy.topo(fh_proxy).halfedges) {
            related_vh_proxies.push_back(mesh_proxy.topo(hh).from());
          }
          return related_vh_proxies;
        },
        // vert2direction
        [&cur_reconstruction, &mesh_proxy,
         &cur_cam](VertHandle vh_proxy) -> RowVector3d {
          Vec3 dir = cur_reconstruction.data(mesh_proxy.data(vh_proxy));
          return RowVector3d(dir[0], dir[1], dir[2]).normalized();
        },
        // energy_fun
        [&mesh_proxy, &sub, &cur_reconstruction, &cur_cam](
            auto &&face2equation, auto &&vert2depth) -> VectorXd {
          EnergyWeights weights;
          weights.faceAngleWeight = 1e3;
          weights.faceMSDAWeight = 0;
          weights.vertMSDAWeight = 0;
          weights.allMSDAWeight = 1;
          auto energy_terms = ComputeEnergy(
              mesh_proxy, sub.fhs.begin(), sub.fhs.end(), weights,
              [&mesh_proxy, &cur_reconstruction](VertHandle vh_proxy) -> Vec3 {
                return normalize(
                    cur_reconstruction.data(mesh_proxy.data(vh_proxy)));
              },
              face2equation, cur_cam.eye());
          return VectorXd::Map(energy_terms.data(), energy_terms.size());
        },
        &v2matrix, &f2matrix);

    if (inversed_depths.mean() < 0) {
      inversed_depths = -inversed_depths;
    }

    if (true) {
      gui::SceneBuilder sb;
      sb.installingOptions().defaultShaderSource =
          gui::OpenGLShaderSourceDescriptor::XTriangles;
      sb.installingOptions().discretizeOptions.color(gui::Black);
      sb.installingOptions().lineWidth = 10;

      for (HalfHandle hh : sub.hhs) {
        VertHandle vh1 = mesh_proxy.topo(hh).from();
        VertHandle vh2 = mesh_proxy.topo(hh).to();
        Point3 p1 = normalize(cur_reconstruction.data(mesh_proxy.data(vh1))) /
                    (v2matrix.at(vh1) * inversed_depths);
        Point3 p2 = normalize(cur_reconstruction.data(mesh_proxy.data(vh2))) /
                    (v2matrix.at(vh2) * inversed_depths);
        sb.add(Line3(p1, p2));
      }
      sb.show(true, false, gui::RenderOptions().camera(cur_cam).winName(
                               "After Optimization"));
    }

    // install
    for (VertHandle vh_proxy : sub.vhs) {
      vh_proxy2inversed_depth[i][vh_proxy] =
          v2matrix.at(vh_proxy) * inversed_depths;
    }
  }

  { // show together
    auto color_table = gui::CreateRandomColorTableWithSize(sub_meshs.size());
    gui::SceneBuilder sb;
    sb.installingOptions().defaultShaderSource =
        gui::OpenGLShaderSourceDescriptor::XLines;
    sb.installingOptions().lineWidth = 10;
    for (int i = 0; i < vh_proxy2inversed_depth.size(); i++) {
      sb.installingOptions().discretizeOptions.color(color_table[i]);
      for (HalfHandle hh : sub_meshs[i].hhs) {
        VertHandle vh1 = mesh_proxy.topo(hh).from();
        VertHandle vh2 = mesh_proxy.topo(hh).to();
        Point3 p1 = normalize(cur_reconstruction.data(mesh_proxy.data(vh1))) /
                    vh_proxy2inversed_depth[i].at(vh1);
        Point3 p2 = normalize(cur_reconstruction.data(mesh_proxy.data(vh2))) /
                    vh_proxy2inversed_depth[i].at(vh2);
        sb.add(Line3(p1, p2));
      }
    }
    sb.show(true, false, gui::RenderOptions().camera(cur_cam).winName(
                             "Together After Optimization"));
  }

  // merge them together
  std::sort(vh_proxy2inversed_depth.begin(), vh_proxy2inversed_depth.end(),
            [](auto &&a, auto &&b) { return a.size() > b.size(); });
  while (vh_proxy2inversed_depth.size() > 1) {
    auto &group1 = vh_proxy2inversed_depth[0];
    auto &group2 = vh_proxy2inversed_depth[1];
    // find common vhs
    std::vector<std::pair<VertHandle, VertHandle>> common_vh_proxies;
    for (auto &v_id1 : group1) {
      VertHandle vh_proxy1 = v_id1.first;
      for (auto &v_id2 : group2) {
        VertHandle vh_proxy2 = v_id2.first;
        if (mesh_proxy.data(vh_proxy1) == mesh_proxy.data(vh_proxy2)) {
          common_vh_proxies.emplace_back(vh_proxy1, vh_proxy2);
        }
      }
    }
    double inversed_depths_mean1 = 0.0;
    double inversed_depths_mean2 = 0.0;
    for (auto &p : common_vh_proxies) {
      inversed_depths_mean1 += group1.at(p.first);
      inversed_depths_mean2 += group2.at(p.second);
    }
    inversed_depths_mean1 /= common_vh_proxies.size();
    inversed_depths_mean2 /= common_vh_proxies.size();

    double inversed_depths_sigma1 = 0.0;
    double inversed_depths_sigma2 = 0.0;
    for (auto &p : common_vh_proxies) {
      inversed_depths_sigma1 +=
          Square(group1.at(p.first) - inversed_depths_mean1);
      inversed_depths_sigma2 +=
          Square(group2.at(p.second) - inversed_depths_mean2);
    }
    inversed_depths_sigma1 =
        sqrt(inversed_depths_sigma1 / common_vh_proxies.size());
    inversed_depths_sigma2 =
        sqrt(inversed_depths_sigma2 / common_vh_proxies.size());

    // merge group2 into group1
    for (auto &v_id2 : group2) {
      assert(!Contains(group1, v_id2.first));
      group1[v_id2.first] = inversed_depths_mean1 +
                            (v_id2.second - inversed_depths_mean2) *
                                inversed_depths_sigma1 / inversed_depths_sigma2;
    }
	// erase group2
	std::swap(group2, vh_proxy2inversed_depth.back());
	vh_proxy2inversed_depth.pop_back();
  }

  HandledTable<VertHandle, double> vh2inversed_depth =
      cur_reconstruction.createVertexTable(-1.0);
  HandledTable<VertHandle, int> vh2vh_proxy_count =
      cur_reconstruction.createVertexTable(0);
  for (auto & v_id : vh_proxy2inversed_depth.front()) {
	  VertHandle vh = mesh_proxy.data(v_id.first);
	  double inv_depth = v_id.second;
	  vh2inversed_depth[vh] += inv_depth;
	  vh2vh_proxy_count[vh] ++;
  }
  for (auto & v : cur_reconstruction.vertices()) {
    vh2inversed_depth[v.topo.hd] /= vh2vh_proxy_count[v.topo.hd];
  }
  {
      gui::SceneBuilder sb;
      sb.installingOptions().defaultShaderSource =
          gui::OpenGLShaderSourceDescriptor::XTriangles;
      sb.installingOptions().discretizeOptions.color(gui::Black);
      sb.installingOptions().lineWidth = 10;

      for (auto &h : cur_reconstruction.halfedges()) {
        HalfHandle hh = h.topo.hd;
        VertHandle vh1 = cur_reconstruction.topo(hh).from();
        VertHandle vh2 = cur_reconstruction.topo(hh).to();
        Point3 p1 = normalize(cur_reconstruction.data(vh1)) / vh2inversed_depth[vh1];
        Point3 p2 = normalize(cur_reconstruction.data(vh2)) / vh2inversed_depth[vh2];
        sb.add(Line3(p1, p2));
      }
      sb.show(true, true,
              gui::RenderOptions().winName("Before Holistic Optimization"));
    }


  // holistic optimization
  {
    std::vector<FaceHandle> fhs;
    for (auto &f : cur_reconstruction.faces()) {
      fhs.push_back(f.topo.hd);
    }
    auto fundamental_vhs_set = SortFacesWithPlanarityDependency(
        cur_reconstruction, fhs.begin(), fhs.end(),
        [&cur_reconstruction](auto vhs_begin, auto vhs_end) -> bool {
          // whether these vhs are colinear?
          if (std::distance(vhs_begin, vhs_end) < 3) {
            return true;
          }
          const Point3 &p0 = cur_reconstruction.data(*vhs_begin++);
          const Point3 &p1 = cur_reconstruction.data(*vhs_begin++);
          Vec3 dir = normalize(p1 - p0);
          while (vhs_begin != vhs_end) {
            const Point3 &pi = cur_reconstruction.data(*vhs_begin++);
            Vec3 diri = normalize(pi - p0);
            if (!IsFuzzyParallel(dir, diri, 1e-6)) {
              return false;
            }
          }
          return true;
        });

    std::vector<VertHandle> fundamental_vhs(fundamental_vhs_set.begin(),
                                            fundamental_vhs_set.end());

    std::map<VertHandle, RowVectorXd> v2matrix;
    std::map<FaceHandle, Matrix3Xd> f2matrix;

    auto inversed_depths = EstimatePerspectiveDepths(
        fundamental_vhs, fhs,
        // vert2initial_depth
        [&vh2inversed_depth](VertHandle vh) -> double {
          return 1.0 / vh2inversed_depth[vh];
        },
        // face2related_verts
        [&cur_reconstruction](FaceHandle fh) -> std::vector<VertHandle> {
          std::vector<VertHandle> related_vhs;
          for (HalfHandle hh : cur_reconstruction.topo(fh).halfedges) {
            related_vhs.push_back(cur_reconstruction.topo(hh).from());
          }
          return related_vhs;
        },
        // vert2direction
        [&cur_reconstruction, &cur_cam](VertHandle vh) -> RowVector3d {
          Vec3 dir = cur_reconstruction.data(vh);
          return RowVector3d(dir[0], dir[1], dir[2]).normalized();
        },
        // energy_fun
        [&cur_reconstruction, &cur_cam, &fhs](auto &&face2equation,
                                              auto &&vert2depth) -> VectorXd {
          EnergyWeights weights;
          weights.faceAngleWeight = 0.1;
          weights.faceMSDAWeight = 1e4;
          weights.vertMSDAWeight = 1e4;
          weights.allMSDAWeight = 1e3;

          auto energy_terms =
              ComputeEnergy(cur_reconstruction, fhs.begin(), fhs.end(), weights,
                            [&cur_reconstruction](VertHandle vh) -> Vec3 {
                              return normalize(cur_reconstruction.data(vh));
                            },
                            face2equation, cur_cam.eye());
          return VectorXd::Map(energy_terms.data(), energy_terms.size());
        },
        &v2matrix, &f2matrix);

    {
      gui::SceneBuilder sb;
      sb.installingOptions().defaultShaderSource =
          gui::OpenGLShaderSourceDescriptor::XTriangles;
      sb.installingOptions().discretizeOptions.color(gui::Black);
      sb.installingOptions().lineWidth = 10;

      for (auto &h : cur_reconstruction.halfedges()) {
        HalfHandle hh = h.topo.hd;
        VertHandle vh1 = cur_reconstruction.topo(hh).from();
        VertHandle vh2 = cur_reconstruction.topo(hh).to();
        Point3 p1 = normalize(cur_reconstruction.data(vh1)) /
                    (v2matrix.at(vh1) * inversed_depths);
        Point3 p2 = normalize(cur_reconstruction.data(vh2)) /
                    (v2matrix.at(vh2) * inversed_depths);
        sb.add(Line3(p1, p2));
      }
      sb.show(true, true,
              gui::RenderOptions().winName("After Holistic Optimization"));
    }
  }
}

auto main(int argc, char **argv, char **env) -> int {
  gui::Singleton::SetCmdArgs(argc, argv, env);
  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath("D:\\Panoramix\\LineDrawing\\");
  misc::Matlab matlab;

  std::string name = "tower";
  //std::string name = "hex";
  std::string camName = "cam1";
  bool resetCam = false;

  std::string objFile = "F:\\LineDrawings\\manifold\\" + name +
                        "\\" + name + ".obj";
  std::string camFile = "F:\\LineDrawings\\manifold\\" + name +
                        "\\" + name + ".obj." + camName + ".cereal";

  //// [Load Mesh]
  auto mesh = LoadFromObjFile(objFile);
  auto meshProxy = MakeMeshProxy(mesh);

  //// [Decompose]
  auto cutFacePairs = DecomposeAll(
      meshProxy, [](HalfHandle hh1, HalfHandle hh2) -> bool { return false; });
  // check validity
  for (auto &hProxy : meshProxy.halfedges()) {
    HalfHandle hh = hProxy.data;
    if (hh.invalid()) {
      hh = meshProxy.data(hProxy.topo.opposite);
    }
    assert(hh.valid());
  }
  std::unordered_map<FaceHandle, FaceHandle> cutFace2Another;
  for (auto &cutFacePair : cutFacePairs) {
    cutFace2Another[cutFacePair.first] = cutFacePair.second;
    cutFace2Another[cutFacePair.second] = cutFacePair.first;
  }
  auto subMeshes = ExtractSubMeshes(meshProxy,
                                    [](auto hhbegin, auto hhend) -> bool {
                                      return std::distance(hhbegin, hhend) <= 1;
                                    },
                                    10);
  Println("found ", subMeshes.size(), " subMeshes");
  HandledTable<FaceHandle, int> fhProxy2subMeshId(
      meshProxy.internalFaces().size(), -1);
  for (int i = 0; i < subMeshes.size(); i++) {
    for (auto fh : subMeshes[i].fhs) {
      fhProxy2subMeshId[fh] = i;
    }
  }

  //// [Load Camera]
  PerspectiveCamera cam;
  if (!LoadFromDisk(camFile, cam) || resetCam) {
    auto sphere = BoundingBoxOfContainer(mesh.vertices()).outerSphere();
    PerspectiveCamera projCam(500, 500, Point2(250, 250), 200,
                              sphere.center + Vec3(1, 2, 3) * sphere.radius * 2,
                              sphere.center);

    const gui::ColorTable ctable =
        gui::CreateRandomColorTableWithSize(subMeshes.size());

    HandledTable<HalfHandle, int> hhProxy2subMeshId(
        meshProxy.internalHalfEdges().size(), -1);
    for (auto &h : meshProxy.halfedges()) {
      int &id = hhProxy2subMeshId[h.topo.hd];
      for (int i = 0; i < subMeshes.size(); i++) {
        if (subMeshes[i].contains(h.topo.hd)) {
          id = i;
          break;
        }
      }
    }
    HandledTable<FaceHandle, int> fhProxy2subMeshId(
        meshProxy.internalFaces().size(), -1);
    for (auto &f : meshProxy.faces()) {
      int &id = fhProxy2subMeshId[f.topo.hd];
      for (int i = 0; i < subMeshes.size(); i++) {
        if (subMeshes[i].contains(f.topo.hd)) {
          id = i;
          break;
        }
      }
    }

    // show each subMesh
    if (true) {
      for (int i = 0; i < subMeshes.size(); i++) {
        Println("subMesh - ", i);
        gui::SceneBuilder sb;
        sb.installingOptions().defaultShaderSource =
            gui::OpenGLShaderSourceDescriptor::XTriangles;
        sb.installingOptions().discretizeOptions.color(gui::Black);
        sb.installingOptions().lineWidth = 0.03;
        AddToScene(
            sb, meshProxy,
            [&subMeshes, i](auto h) { return subMeshes[i].contains(h); },
            [&mesh](VertHandle vh) { return mesh.data(vh); },
            [&hhProxy2subMeshId, &ctable](HalfHandle hh) { return gui::Black; },
            [&fhProxy2subMeshId, &ctable](FaceHandle fh) {
              return ctable[fhProxy2subMeshId[fh]];
            });
        sb.show(true, false, gui::RenderOptions()
                                 .camera(projCam)
                                 .backgroundColor(gui::White)
                                 .renderMode(gui::All)
                                 .bwTexColor(0.0)
                                 .bwColor(1.0)
                                 .fixUpDirectionInCameraMove(false)
                                 .cullBackFace(false)
                                 .cullFrontFace(false));
      }
    }

    { // show together
      gui::SceneBuilder sb;
      sb.installingOptions().defaultShaderSource =
          gui::OpenGLShaderSourceDescriptor::XTriangles;
      sb.installingOptions().discretizeOptions.color(gui::Black);
      sb.installingOptions().lineWidth = 0.03;
      AddToScene(
          sb, meshProxy, [](auto) { return true; },
          [&mesh](VertHandle vh) { return mesh.data(vh); },
          [&hhProxy2subMeshId, &ctable](HalfHandle hh) { return gui::Black; },
          [&fhProxy2subMeshId, &ctable](FaceHandle fh) {
            return ctable[fhProxy2subMeshId[fh]];
          });
      cam = sb.show(true, false, gui::RenderOptions()
                                     .camera(projCam)
                                     .backgroundColor(gui::White)
                                     .renderMode(gui::All)
                                     .bwTexColor(0.0)
                                     .bwColor(1.0)
                                     .fixUpDirectionInCameraMove(false)
                                     .cullBackFace(false)
                                     .cullFrontFace(false))
                .camera();
    }
    SaveToDisk(camFile, cam);
  }

  //// [Make 2D Mesh]
  // convert to 2d
  auto mesh2d = Transform(
      mesh, [&cam](const Point3 &p) -> Point2 { return cam.toScreen(p); });

  // add offset noise
  Vec2 offsetNoise = Vec2(20, -20);
  for (auto &v : mesh2d.vertices()) {
    v.data += offsetNoise;
  }

  if (true) {
    Image3ub im(cam.screenSize(), Vec3ub(255, 255, 255));
    auto canvas = gui::MakeCanvas(im);
    canvas.color(gui::Black);
    canvas.thickness(2);
    for (auto &h : mesh2d.halfedges()) {
      auto &p1 = mesh2d.data(h.topo.from());
      auto &p2 = mesh2d.data(h.topo.to());
      canvas.add(Line2(p1, p2));
    }
    canvas.show(0, "mesh2d");
  }

  // analyze the fh/vh dependencies in each subMesh
  for (int subMeshId = 0; subMeshId < subMeshes.size(); subMeshId++) {
    // find the vhs that should serve as anchoring variables
    // and sort the order of the fhs to be determined
    subMeshes[subMeshId].ComputeVertexFaceDependencies(
        meshProxy, [&meshProxy, &mesh2d](auto vhs_begin, auto vhs_end) -> bool {
          // whether these vhs are colinear?
          if (std::distance(vhs_begin, vhs_end) < 3) {
            return true;
          }
          const Point2 &p0 = mesh2d.data(meshProxy.data(*vhs_begin++));
          const Point2 &p1 = mesh2d.data(meshProxy.data(*vhs_begin++));
          Vec2 dir = normalize(p1 - p0);
          while (vhs_begin != vhs_end) {
            const Point2 &pi = mesh2d.data(meshProxy.data(*vhs_begin++));
            Vec2 diri = normalize(pi - p0);
            if (!IsFuzzyParallel(dir, diri, 1e-6)) {
              return false;
            }
          }
          return true;
        });
  }

  //// [Estimate PP & Focal Candidates from 2D Mesh]
  auto point2dAtProxy = [&mesh2d, &meshProxy](VertHandle vhInProxy) -> Point2 {
    return mesh2d.data(meshProxy.data(vhInProxy));
  };
  auto line2dAtProxy = [&mesh2d, &meshProxy,
                        point2dAtProxy](HalfHandle hhInProxy) -> Line2 {
    return Line2(point2dAtProxy(meshProxy.topo(hhInProxy).from()),
                 point2dAtProxy(meshProxy.topo(hhInProxy).to()));
  };
  Box2 box = BoundingBoxOfContainer(mesh2d.vertices());
  double scale = box.outerSphere().radius;

  struct PPFocalCandidate {
    Point2 pp;
    double focal;
  };

  // collect pp focal candidates
  std::vector<PPFocalCandidate> ppFocalCandidates;
  ppFocalCandidates.reserve(subMeshes.size() * 3);

  for (int subMeshId = 0; subMeshId < subMeshes.size(); subMeshId++) {
    // collect edge intersections in each face
    std::vector<Point2> interps;
    for (auto fhProxy : subMeshes[subMeshId].fhs) {
      auto &hhProxys = meshProxy.topo(fhProxy).halfedges;
      Chain2 corners;
      for (auto hhProxy : hhProxys) {
        corners.append(point2dAtProxy(meshProxy.topo(hhProxy).to()));
      }
      auto keyVPs = PossibleKeyVanishingPoints(corners);
      interps.insert(interps.end(), keyVPs.begin(), keyVPs.end());
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
              Distance(pp, box.center()) > scale * 2.0) {
            continue;
          }
          ppFocalCandidates.push_back(PPFocalCandidate{pp, focal});
        }
      }
    }
  }

  std::sort(ppFocalCandidates.begin(), ppFocalCandidates.end(),
            [](auto &a, auto &b) { return a.focal < b.focal; });
  // naive clustering
  std::vector<std::pair<std::set<int>, PPFocalCandidate>> ppFocalGroups;
  {
    std::vector<int> ppFocalId2group(ppFocalCandidates.size(), -1);
    int ngroups = 0;
    RTreeMap<Vec3, int> ppFocalIdTree;
    for (int i = 0; i < ppFocalCandidates.size(); i++) {
      Vec3 coordinate =
          cat(ppFocalCandidates[i].pp, ppFocalCandidates[i].focal);
      const double thres = scale / 50.0;
      // find the nearest ppFocal sample point
      int nearestPPFocalCandId = -1;
      double minDist = thres;
      ppFocalIdTree.search(BoundingBox(coordinate).expand(thres * 2),
                           [&nearestPPFocalCandId, &minDist,
                            &coordinate](const std::pair<Vec3, int> &cand) {
                             double dist = Distance(cand.first, coordinate);
                             if (dist < minDist) {
                               minDist = dist;
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

    ppFocalGroups.resize(ngroups);
    for (auto &g : ppFocalGroups) {
      g.second.focal = 0.0;
      g.second.pp = Point2();
    }
    for (int i = 0; i < ppFocalId2group.size(); i++) {
      auto &g = ppFocalGroups[ppFocalId2group[i]];
      g.first.insert(i);
      g.second.focal += ppFocalCandidates[i].focal;
      g.second.pp += ppFocalCandidates[i].pp;
    }
    for (auto &g : ppFocalGroups) {
      g.second.focal /= g.first.size();
      g.second.pp /= double(g.first.size());
    }

    std::sort(
        ppFocalGroups.begin(), ppFocalGroups.end(),
        [](auto &g1, auto &g2) { return g1.first.size() > g2.first.size(); });
  }

  //// [Orient Edges]
  // record edges
  std::vector<std::pair<HalfHandle, HalfHandle>> edge2hh2ds;
  std::vector<Line2> edge2line;
  HandledTable<HalfHandle, int> hh2d2edge(mesh2d.internalHalfEdges().size(),
                                          -1);
  int nedges = 0;
  {
    for (auto &h : mesh2d.halfedges()) {
      auto hh = h.topo.hd;
      auto oppohh = h.topo.opposite;
      if (hh2d2edge[hh] == -1 && hh2d2edge[oppohh] == -1) {
        hh2d2edge[hh] = hh2d2edge[oppohh] = nedges;
        nedges++;
        edge2hh2ds.push_back(MakeOrderedPair(hh, oppohh));
        edge2line.push_back(Line2(mesh2d.data(mesh2d.topo(hh).from()),
                                  mesh2d.data(mesh2d.topo(hh).to())));
      }
    }
    assert(edge2hh2ds.size() == nedges && edge2line.size() == nedges);
  }

  // collect edge intersections and
  // get vpPositions from the intersections
  std::vector<Point2> vpPositions;
  int nvps = 0;
  std::vector<std::vector<Scored<int>>> edge2OrderedVPAndAngles(nedges);
  {
    std::vector<Point2> intersections;
    std::vector<std::pair<int, int>> intersection2edges;
    intersections.reserve(nedges * (nedges - 1) / 2);
    intersection2edges.reserve(nedges * (nedges - 1) / 2);
    for (int i = 0; i < nedges; i++) {
      const Line2 &linei = edge2line[i];
      for (int j = i + 1; j < nedges; j++) {
        const Line2 &linej = edge2line[j];
        Point2 interp = Intersection(linei.ray(), linej.ray());
        if (std::min(Distance(interp, linei), Distance(interp, linej)) <=
            scale / 10.0) {
          continue;
        }
        intersections.push_back(interp);
        intersection2edges.emplace_back(i, j);
      }
    }
    assert(intersections.size() == intersection2edges.size());

    std::vector<int> intersection2Rawvp(intersections.size(), -1);
    RTreeMap<Point2, int> intersectionTree;
    for (int i = 0; i < intersections.size(); i++) {
      const double thres = scale / 30.0;
      auto &p = intersections[i];
      int nearestIntersectionId = -1;
      double minDist = thres;
      intersectionTree.search(
          BoundingBox(p).expand(thres * 2),
          [&nearestIntersectionId, &minDist,
           &p](const std::pair<Point2, int> &locationAndIntersectionId) {
            double dist = Distance(locationAndIntersectionId.first, p);
            if (dist < minDist) {
              minDist = dist;
              nearestIntersectionId = locationAndIntersectionId.second;
            }
            return true;
          });
      if (nearestIntersectionId != -1) {
        intersection2Rawvp[i] = intersection2Rawvp[nearestIntersectionId];
      } else {
        intersection2Rawvp[i] = nvps++;
      }
      intersectionTree.emplace(p, i);
    }

    // further merge vps
    // if any two vps share two or more edges, merge them
    std::vector<std::set<int>> rawVP2edges(nvps);
    for (int i = 0; i < intersection2Rawvp.size(); i++) {
      int vpid = intersection2Rawvp[i];
      rawVP2edges[vpid].insert(i);
    }
    std::vector<std::set<int>> rawVp2rawVpShouldMerge(nvps);
    for (int vp1 = 0; vp1 < nvps; vp1++) {
      auto &edges1 = rawVP2edges[vp1];
      for (int vp2 = vp1 + 1; vp2 < nvps; vp2++) {
        auto &edges2 = rawVP2edges[vp2];
        std::set<int> commonEdges;
        std::set_intersection(edges1.begin(), edges1.end(), edges2.begin(),
                              edges2.end(),
                              std::inserter(commonEdges, commonEdges.begin()));
        if (commonEdges.size() >= 2) {
          rawVp2rawVpShouldMerge[vp1].insert(vp2);
          rawVp2rawVpShouldMerge[vp2].insert(vp1);
        }
      }
    }
    std::map<int, std::set<int>> newVP2rawVPs;
    std::vector<int> rawVP2newVP(nvps, -1);
    std::vector<int> rawVPIds(nvps);
    std::iota(rawVPIds.begin(), rawVPIds.end(), 0);
    nvps = ConnectedComponents(
        rawVPIds.begin(), rawVPIds.end(),
        [&rawVp2rawVpShouldMerge](int vp) -> const std::set<int> & {
          return rawVp2rawVpShouldMerge.at(vp);
        },
        [&newVP2rawVPs, &rawVP2newVP](int rawVP, int newVP) {
          newVP2rawVPs[newVP].insert(rawVP);
          rawVP2newVP[rawVP] = newVP;
        });

    vpPositions.resize(nvps, Origin<2>());
    std::vector<std::set<int>> vp2intersections(nvps);
    for (int i = 0; i < intersection2Rawvp.size(); i++) {
      int rawVP = intersection2Rawvp[i];
      int newVP = rawVP2newVP[rawVP];
      vpPositions[newVP] += intersections[i]; // TODO: what if some
                                              // intersections are oppsite far
                                              // points?
      vp2intersections[newVP].insert(i);
    }
    for (int i = 0; i < nvps; i++) {
      vpPositions[i] /= double(vp2intersections[i].size());
    }

    // initial edge vp bindings
    std::vector<std::map<int, double>> vp2edgeWithAngles(nvps);
    std::vector<bool> vpIsGood(nvps, true);
    for (int vp = 0; vp < nvps; vp++) {
      auto &vpPos = vpPositions[vp];
      for (int edge = 0; edge < nedges; edge++) {
        auto &line = edge2line[edge];
        double lambda = ProjectionOfPointOnLine(vpPos, line).ratio;
        static const double thres = 0.1;
        if (lambda >= -thres && lambda <= 1.0 + thres) {
          continue;
        }
        double angle =
            AngleBetweenUndirected(line.direction(), vpPos - line.center());
        static const double theta = DegreesToRadians(5); ////// TODO
        if (angle >= theta) {
          continue;
        }
        vp2edgeWithAngles[vp][edge] = angle;
      }
      vpIsGood[vp] = vp2edgeWithAngles[vp].size() >= 3;
    }

    if (false) {
      for (int i = 0; i < std::min(10, nvps); i++) {
        Image3ub im(cam.screenSize(), Vec3ub(255, 255, 255));
        auto canvas = gui::MakeCanvas(im);
        canvas.color(gui::LightGray);
        canvas.thickness(2);
        for (auto &line : edge2line) {
          canvas.add(line);
        }
        canvas.color(gui::Gray);
        canvas.thickness(1);
        for (auto &edgeWithAngle : vp2edgeWithAngles[i]) {
          canvas.add(edge2line[edgeWithAngle.first].ray());
        }
        canvas.color(gui::Black);
        for (auto &edgeWithAngle : vp2edgeWithAngles[i]) {
          canvas.add(edge2line[edgeWithAngle.first]);
        }
        canvas.show(0, "before removing bad vps: raw vp_" + std::to_string(i));
      }
    }

    // remove bad vps
    int newvp = 0;
    for (int oldvp = 0; oldvp < nvps; oldvp++) {
      if (vpIsGood[oldvp]) {
        // update vpPositions
        // update vp2edgeWithAngles
        vpPositions[newvp] = vpPositions[oldvp];
        vp2edgeWithAngles[newvp] = std::move(vp2edgeWithAngles[oldvp]);
        newvp++;
      }
    }
    nvps = newvp;
    vpPositions.resize(nvps);
    vp2edgeWithAngles.resize(nvps);

    // construct edge2OrderedVPs
    for (int vp = 0; vp < nvps; vp++) {
      for (auto &edgeAndAngle : vp2edgeWithAngles[vp]) {
        int edge = edgeAndAngle.first;
        double angle = edgeAndAngle.second;
        edge2OrderedVPAndAngles[edge].push_back(ScoreAs(vp, angle));
      }
    }
    for (auto &vpAndAngles : edge2OrderedVPAndAngles) {
      std::sort(vpAndAngles.begin(), vpAndAngles.end());
    }

    if (false) {
      for (int i = 0; i < std::min(10, nvps); i++) {
        Image3ub im(cam.screenSize(), Vec3ub(255, 255, 255));
        auto canvas = gui::MakeCanvas(im);
        canvas.color(gui::LightGray);
        canvas.thickness(2);
        for (auto &line : edge2line) {
          canvas.add(line);
        }
        canvas.color(gui::Gray);
        canvas.thickness(2);
        for (auto &edgeWithAngle : vp2edgeWithAngles[i]) {
          canvas.add(edge2line[edgeWithAngle.first].ray());
        }
        canvas.color(gui::Black);
        for (auto &edgeWithAngle : vp2edgeWithAngles[i]) {
          canvas.add(edge2line[edgeWithAngle.first]);
        }
        canvas.show(0, "raw vp_" + std::to_string(i));
      }
    }
  }

  std::vector<int> edge2vp;
  std::vector<std::vector<int>> vp2edges;

  // factor graph
  if (false ||
      !misc::LoadCache(objFile + "-" + camName, "edge2vp_vp2edges", edge2vp,
                       vp2edges)) {
    edge2vp.resize(nedges, -1);
    vp2edges.resize(nvps);

    // construct a factor graph to optimize edge-vp bindings
    FactorGraph fg;
    std::vector<FactorGraph::VarHandle> edge2vh(nedges);
    {
      for (int edge = 0; edge < nedges; edge++) {
        auto vc =
            fg.addVarCategory(edge2OrderedVPAndAngles[edge].size() + 1, 1.0);
        edge2vh[edge] = fg.addVar(vc);
      }

      // potential 1: the edge should bind to some vp, should prefer better
      // scored
      for (int edge = 0; edge < nedges; edge++) {
        auto vh = edge2vh[edge];
        auto &relatedVPAndAngles = edge2OrderedVPAndAngles[edge];
        auto fc = fg.addFactorCategory(
            [&relatedVPAndAngles, nedges](const int *varlabels, size_t nvar,
                                          FactorGraph::FactorCategoryId fcid,
                                          void *givenData) -> double {
              assert(nvar == 1);
              int label = varlabels[0];
              assert(label <= relatedVPAndAngles.size());
              const double K = 50.0 / nedges;
              if (label == relatedVPAndAngles.size()) { // not bind to any vp
                return K;
              }
              double angle = relatedVPAndAngles[label].score;
              assert(!IsInfOrNaN(angle));
              return (1.0 - Gaussian(angle, DegreesToRadians(3))) * K;
            },
            1.0);
        fg.addFactor({vh}, fc);
      }

      // potential 2: two adjacent edges should not bind to a near vp
      int ncorners = 0;
      for (auto &f : mesh2d.faces()) {
        ncorners += f.topo.halfedges.size();
      }
      for (auto &f : mesh2d.faces()) {
        auto &hhs = f.topo.halfedges;
        for (int i = 0; i < hhs.size(); i++) {
          auto hh1 = hhs[i];
          auto hh2 = hhs[(i + 1) % hhs.size()];
          int edge1 = hh2d2edge[hh1];
          int edge2 = hh2d2edge[hh2];
          auto &relatedVPAndAngles1 = edge2OrderedVPAndAngles[edge1];
          auto &relatedVPAndAngles2 = edge2OrderedVPAndAngles[edge2];
          auto vh1 = edge2vh[edge1];
          auto vh2 = edge2vh[edge2];
          auto fc = fg.addFactorCategory(
              [edge1, edge2, &relatedVPAndAngles1, &relatedVPAndAngles2,
               &vpPositions, ncorners,
               scale](const int *varlabels, size_t nvar,
                      FactorGraph::FactorCategoryId fcid,
                      void *givenData) -> double {
                assert(nvar == 2);
                int bindedVP1 =
                    varlabels[0] == relatedVPAndAngles1.size()
                        ? -1
                        : (relatedVPAndAngles1[varlabels[0]].component);
                int bindedVP2 =
                    varlabels[1] == relatedVPAndAngles2.size()
                        ? -1
                        : (relatedVPAndAngles2[varlabels[1]].component);
                if (bindedVP1 == -1 || bindedVP2 == -1) {
                  return 0;
                }
                auto &vpPos1 = vpPositions[bindedVP1];
                auto &vpPos2 = vpPositions[bindedVP2];
                const double thres = scale / 10.0;
                const double K = 10.0 / ncorners;
                if (Distance(vpPos1, vpPos2) < thres) { // todo
                  return K;
                }
                return 0.0;
              },
              1.0);
          fg.addFactor({vh1, vh2}, fc);
        }
      }

      // potential 3: the vpPositions of edges sharing a same face should lie on
      // the same
      // line (the vanishing line of the face)
      int ntris = 0;
      for (auto &f : mesh2d.faces()) {
        auto &hhs = f.topo.halfedges;
        if (hhs.size() <= 3) {
          continue;
        }
        ntris = hhs.size() > 4 ? (2 * hhs.size()) : hhs.size();
      }
      for (auto &f : mesh2d.faces()) {
        auto &hhs = f.topo.halfedges;
        if (hhs.size() <= 3) {
          continue;
        }
        for (int i = 0; i < hhs.size(); i++) {
          int maxGap = hhs.size() > 4 ? 2 : 1;
          for (int gap = 1; gap <= maxGap; gap++) {
            int prevEdge = hh2d2edge[hhs[(i + hhs.size() - gap) % hhs.size()]];
            int edge = hh2d2edge[hhs[i]];
            int nextEdge = hh2d2edge[hhs[(i + gap) % hhs.size()]];
            auto fc = fg.addFactorCategory(
                [&edge2OrderedVPAndAngles, &vpPositions, prevEdge, edge,
                 nextEdge, ntris](const int *varlabels, size_t nvar,
                                  FactorGraph::FactorCategoryId fcid,
                                  void *givenData) -> double {
                  assert(nvar == 3);
                  int vp1 =
                      varlabels[0] == edge2OrderedVPAndAngles[prevEdge].size()
                          ? -1
                          : (edge2OrderedVPAndAngles[prevEdge][varlabels[0]]
                                 .component);
                  if (vp1 == -1) {
                    return 0.0;
                  }
                  int vp2 = varlabels[1] == edge2OrderedVPAndAngles[edge].size()
                                ? -1
                                : (edge2OrderedVPAndAngles[edge][varlabels[1]]
                                       .component);
                  if (vp2 == -1) {
                    return 0.0;
                  }
                  int vp3 =
                      varlabels[2] == edge2OrderedVPAndAngles[nextEdge].size()
                          ? -1
                          : (edge2OrderedVPAndAngles[nextEdge][varlabels[2]]
                                 .component);
                  if (vp3 == -1) {
                    return 0.0;
                  }
                  if (vp1 == vp2 || vp2 == vp3 || vp1 == vp3) {
                    return 0.0;
                  }
                  double angle = AngleBetweenUndirected(
                      vpPositions[vp1] - vpPositions[vp2],
                      vpPositions[vp3] - vpPositions[vp2]);
                  assert(!IsInfOrNaN(angle));
                  const double K = 30.0 / ntris;
                  return (1.0 - Gaussian(angle, DegreesToRadians(10))) * K;
                },
                1.0);
            fg.addFactor({edge2vh[prevEdge], edge2vh[edge], edge2vh[nextEdge]},
                         fc);
          }
        }
      }

      //// potential 4: a vp is good only when edges are bound to it!
      //// a vp should be with either no edges or >= 3 edges
      // for (int vpid = 0; vpid < nvps; vpid++) {
      //  std::vector<FactorGraph::VarHandle> relatedVhs;
      //  relatedVhs.reserve(vp2edges[vpid].size());
      //  for (int edge : vp2edges[vpid]) {
      //    relatedVhs.push_back(edge2vh[edge]);
      //  }
      //  auto fc = fg.addFactorCategory(
      //      [&edge2OrderedVPAndAngles, vpid, &vp2edges](
      //          const int *varlabels, size_t nvar,
      //          FactorGraph::FactorCategoryId fcid, void *givenData) -> double
      //          {
      //        int bindedEdges = 0;
      //        auto &relatedEdges = vp2edges[vpid];
      //        assert(nvar == relatedEdges.size());
      //        for (int i = 0; i < relatedEdges.size(); i++) {
      //          int edge = relatedEdges[i];
      //          int edgeLabel = varlabels[i];
      //          if (edgeLabel == edge2OrderedVPAndAngles[edge].size()) {
      //            continue; // not bound to any vp
      //          }
      //          int edgeBindedVPId =
      //          edge2OrderedVPAndAngles[edge][edgeLabel].component;
      //          if (edgeBindedVPId == vpid) {
      //            bindedEdges++;
      //          }
      //        }
      //        // todo
      //        if (bindedEdges == 1 || bindedEdges == 2) {
      //          return 10.0;
      //        }
      //        return 0.0;
      //      },
      //      1.0);
      //  fg.addFactor(relatedVhs.begin(), relatedVhs.end(), fc);
      //}
    }

    // solve the factor graph
    auto result =
        fg.solve(5, 1, [](int epoch, double energy, double denergy,
                          const FactorGraph::ResultTable &results) -> bool {
          Println("epoch: ", epoch, "  energy: ", energy);
          return true;
        });

    for (int edge = 0; edge < nedges; edge++) {
      int id = result[edge2vh[edge]];
      if (id == edge2OrderedVPAndAngles[edge].size()) {
        continue;
      }
      edge2vp[edge] = edge2OrderedVPAndAngles[edge][id].component;
      vp2edges[edge2vp[edge]].push_back(edge);
    }

    // invalidate the edge bindings for vps who have <= 3 edges
    for (int vp = 0; vp < nvps; vp++) {
      if (vp2edges[vp].size() <= 3) { //
        for (int edge : vp2edges[vp]) {
          edge2vp[edge] = -1;
        }
        vp2edges[vp].clear();
      }
    }

    misc::SaveCache(objFile + "-" + camName, "edge2vp_vp2edges", edge2vp,
                    vp2edges);
  }

  if (true) { // show line classification results
    for (int i = 0; i < nvps; i++) {
      if (vp2edges[i].empty()) {
        continue;
      }
      Image3ub im(cam.screenSize(), Vec3ub(255, 255, 255));
      auto canvas = gui::MakeCanvas(im);
      canvas.color(gui::LightGray);
      canvas.thickness(2);
      for (auto &line : edge2line) {
        canvas.add(line);
      }
      canvas.color(gui::Gray);
      canvas.thickness(2);
      for (int edge : vp2edges[i]) {
        canvas.add(edge2line[edge].ray());
      }
      canvas.color(gui::Black);
      for (int edge : vp2edges[i]) {
        canvas.add(edge2line[edge]);
      }
      canvas.show(0, "optimized vp_" + std::to_string(i));
    }
  }

  assert(edge2vp.size() == nedges);
  assert(vp2edges.size() == nvps);

  // for each pp focal candidate
  std::vector<Mesh3> reconstructions;
  for (int configId = 0;
       configId < std::min(5ull, ppFocalGroups.size()) &&
       ppFocalGroups[configId].first.size() * 30 >= ppFocalCandidates.size();
       configId++) {

    double focal = ppFocalGroups[configId].second.focal;
    auto &pp = ppFocalGroups[configId].second.pp;

    PerspectiveCamera curCam(cam.screenWidth(), cam.screenHeight(), pp, focal);
    std::vector<Vec3> vp2dir(nvps);
    for (int i = 0; i < nvps; i++) {
      vp2dir[i] = curCam.direction(vpPositions[i]);
    }

    auto meshReconstructed =
        ReconstructWithOrientations(edge2line, edge2vp, vp2dir, hh2d2edge,
                                    curCam, meshProxy, mesh2d, matlab);

    // optimize without orientations
    OptimizeWithoutOrientations(meshReconstructed, curCam, meshProxy,
                                subMeshes);

    reconstructions.push_back(std::move(meshReconstructed));
  }

  return 0;
}
