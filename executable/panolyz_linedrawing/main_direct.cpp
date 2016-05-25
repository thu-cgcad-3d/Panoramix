#include <QtWidgets>

#include "../../src/core/parallel.hpp"
#include "../../src/misc/cache.hpp"
#include "../../src/misc/clock.hpp"

#include "../../src/gui/canvas.hpp"
#include "../../src/gui/qttools.hpp"
#include "../../src/gui/scene.hpp"
#include "../../src/gui/singleton.hpp"
#include "../../src/gui/gui_util.hpp"

#include "../../src/experimental/line_drawing.hpp"
#include "../../src/experimental/mesh_advanced_util.hpp"
#include "../../src/experimental/pi_graph_annotation.hpp"
#include "../../src/experimental/pi_graph_cg.hpp"
#include "../../src/experimental/pi_graph_control.hpp"
#include "../../src/experimental/pi_graph_occlusion.hpp"
#include "../../src/experimental/pi_graph_optimize.hpp"
#include "../../src/experimental/pi_graph_vis.hpp"

using namespace pano;
using namespace pano::core;
using namespace pano::experimental;

template <class HalfColorerFunT, class FaceColorerFunT>
void AddToScene(gui::SceneBuilder &sb, const Mesh<Point3> &m,
                HalfColorerFunT colorHalf, FaceColorerFunT colorFace) {
  sb.installingOptions().lineWidth = 10.0;
  sb.installingOptions().defaultShaderSource =
      gui::OpenGLShaderSourceDescriptor::XLines;
  for (auto &h : m.halfedges()) {
    Line3 line(m.data(h.topo.from()), m.data(h.topo.to()));
    auto hh = h.topo.hd;
    auto fh = h.topo.face;
    auto oppohh = h.topo.opposite;
    auto oppofh = oppohh.valid() ? m.topo(oppohh).face : FaceHandle();
    bool hasFace = fh.valid();
    bool hasOppo = oppohh.valid();
    gui::Color color = colorHalf(hh);
    if (!hasFace && hasOppo) {
      color = gui::Red;
    } else if (hasFace && !hasOppo) {
      color = gui::Blue;
    } else if (!hasFace && !hasOppo) {
      color = gui::Yellow;
    }
    sb.add(gui::ColorAs(line, color), [hh, fh, oppohh, oppofh](auto &...) {
      std::cout << "halfedge id: " << hh.id
                << ", opposite halfedge id: " << oppohh.id
                << ", face id: " << fh.id << ", opposite face id: " << oppofh.id
                << '\n';
    });
  }
  sb.installingOptions().defaultShaderSource =
      gui::OpenGLShaderSourceDescriptor::XTriangles;
  for (auto &f : m.faces()) {
    Polygon3 poly;
    for (auto h : f.topo.halfedges) {
      auto v = m.topo(h).to();
      poly.corners.push_back(m.data(v));
    }
    assert(poly.corners.size() > 2);
    poly.normal = (poly.corners[0] - poly.corners[1])
                      .cross(poly.corners[0] - poly.corners[2]);
    auto fh = f.topo.hd;
    sb.add(gui::ColorAs(poly, colorFace(f.topo.hd)),
           [fh](auto &...) { std::cout << "face id: " << fh.id << '\n'; });
  }
}

template <class HalfColorerFunT, class FaceColorerFunT>
void AddToScene(gui::SceneBuilder &sb, const Mesh<Point3> &m,
                const SubMesh &sub, HalfColorerFunT colorHalf,
                FaceColorerFunT colorFace) {
  sb.installingOptions().lineWidth = 10.0;
  sb.installingOptions().defaultShaderSource =
      gui::OpenGLShaderSourceDescriptor::XLines;
  for (auto &h : m.halfedges()) {
    if (!Contains(sub.hhs, h.topo.hd))
      continue;
    Line3 line(m.data(h.topo.from()), m.data(h.topo.to()));
    auto hh = h.topo.hd;
    auto fh = h.topo.face;
    auto oppohh = h.topo.opposite;
    auto oppofh = oppohh.valid() ? m.topo(oppohh).face : FaceHandle();
    bool hasFace = fh.valid();
    bool hasOppo = oppohh.valid();
    gui::Color color = colorHalf(hh);
    if (!hasFace && hasOppo) {
      color = gui::Red;
    } else if (hasFace && !hasOppo) {
      color = gui::Blue;
    } else if (!hasFace && !hasOppo) {
      color = gui::Yellow;
    }
    sb.add(gui::ColorAs(line, color), [hh, fh, oppohh, oppofh](auto &...) {
      std::cout << "halfedge id: " << hh.id
                << ", opposite halfedge id: " << oppohh.id
                << ", face id: " << fh.id << ", opposite face id: " << oppofh.id
                << '\n';
    });
  }
  sb.installingOptions().defaultShaderSource =
      gui::OpenGLShaderSourceDescriptor::XTriangles;
  for (auto &f : m.faces()) {
    if (!Contains(sub.fhs, f.topo.hd))
      continue;
    Polygon3 poly;
    for (auto h : f.topo.halfedges) {
      auto v = m.topo(h).to();
      poly.corners.push_back(m.data(v));
    }
    assert(poly.corners.size() > 2);
    poly.normal = (poly.corners[0] - poly.corners[1])
                      .cross(poly.corners[0] - poly.corners[2]);
    auto fh = f.topo.hd;
    sb.add(gui::ColorAs(poly, colorFace(f.topo.hd)),
           [fh](auto &...) { std::cout << "face id: " << fh.id << '\n'; });
  }
}

struct EnergyWeights {
  double vertMSDAWeight;
  double faceMSDAWeight;
  double faceAngleWeight;
};

template <class FaceHandleToPlaneEqFunT>
double ComputeEnergy(const Mesh2 &mesh, const SubMesh &sub,
                     const PerspectiveCamera &cam, const EnergyWeights &weights,
                     FaceHandleToPlaneEqFunT fh2planeEq) {
  int angleNum = sub.hhs.size();
  // fh2planeEq: (FaceHandle)[i] -> double
  std::map<VertHandle, Point3> vpositions;
  std::map<VertHandle, int> vfacedegrees;
  for (VertHandle vh : sub.vhs) {
    vpositions[vh] = Vec3();
    vfacedegrees[vh] = 0;
  }
  for (FaceHandle fh : sub.fhs) {
    decltype(auto) planeEq = fh2planeEq(fh);
    Plane3 plane = Plane3FromEquation(planeEq[0], planeEq[1], planeEq[2]);
    if (HasValue(plane, IsInfOrNaN<double>)) {
      return std::numeric_limits<double>::infinity();
    }
    for (HalfHandle hh : mesh.topo(fh).halfedges) {
      VertHandle vh = mesh.topo(hh).to();
      const Point2 &p2d = mesh.data(vh);
      auto dir = normalize(cam.direction(p2d));
      Point3 p3d = Intersection(Ray3(cam.eye(), dir), plane);
      vpositions[vh] += p3d;
      vfacedegrees[vh]++;
    }
  }
  for (auto vhpos : vpositions) {
    vhpos.second /= vfacedegrees.at(vhpos.first);
  }

  double faceMSDA = 0.0;
  int faceMSDANum = 0;

  double vertexMSDA = 0.0;
  int vertexMSDANum = 0;

  std::map<VertHandle, std::vector<double>> vh2angles;
  for (FaceHandle fh : sub.fhs) {
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
    }
    double faceMeanAngle =
        std::accumulate(faceAngles.begin(), faceAngles.end(), 0.0) /
        faceAngles.size();
    for (double a : faceAngles) {
      faceMSDA += Square(a - faceMeanAngle);
      faceMSDANum++;
    }
  }
  faceMSDA /= faceMSDANum;

  for (auto &vhangles : vh2angles) {
    VertHandle vh = vhangles.first;
    auto &angles = vhangles.second;
    double vertMeanAngle =
        std::accumulate(angles.begin(), angles.end(), 0.0) / angles.size();
    for (double a : angles) {
      vertexMSDA += Square(a - vertMeanAngle);
      vertexMSDANum++;
    }
  }
  vertexMSDA /= vertexMSDANum;

  double faceAngleEnergy = 0.0;
  for (auto &hh : sub.hhs) {
    auto fh1 = mesh.topo(hh).face;
    auto planeEq1 = fh2planeEq(fh1);
    auto fh2 = mesh.topo(mesh.topo(hh).opposite).face;
    auto planeEq2 = fh2planeEq(fh2);
    double angle = AngleBetweenUndirected(
        Vec3(planeEq1[0], planeEq1[1], planeEq1[2]),
        Vec3(planeEq2[0], planeEq2[1], planeEq2[2]));
    assert(!IsInfOrNaN(angle));
    faceAngleEnergy += Square(angle - M_PI_2);
  }
  faceAngleEnergy /= sub.hhs.size();

  double energy = faceMSDA * weights.faceMSDAWeight +
                  vertexMSDA * weights.vertMSDAWeight +
                  faceAngleEnergy * weights.faceAngleWeight;
  // std::cout << "energy: " << energy << "\n";
  return energy;
}

void TestOnMesh(const Mesh3 &meshGT, int maxIters,
                const EnergyWeights &weights) {
  auto sphere = BoundingBoxOfContainer(meshGT.vertices()).outerSphere();
  PerspectiveCamera projCam(500, 500, Point2(250, 250), 200,
                            sphere.center + Vec3(1, 2, 3) * sphere.radius * 2,
                            sphere.center);
  {
    gui::SceneBuilder sb;
    sb.installingOptions().defaultShaderSource =
        gui::OpenGLShaderSourceDescriptor::XTriangles;
    sb.installingOptions().discretizeOptions.color(gui::Black);
    sb.installingOptions().lineWidth = 0.1;
    const gui::ColorTable ctable =
        gui::CreateRandomColorTableWithSize(meshGT.internalFaces().size());
    AddToScene(sb, meshGT, ConstantFunctor<gui::Color>(gui::Black),
               [&ctable](FaceHandle fh) { return ctable[fh.id]; });
    projCam = sb.show(true, false, gui::RenderOptions()
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

  Mesh2 meshProjected =
      Transform(meshGT, [&projCam](const Point3 &p) -> Point2 {
        return projCam.toScreen(p);
      });

  PerspectiveCamera cam = projCam;
  cam.setEye(Origin(), false);
  cam.setCenter(X(), true);

  auto subs = ExtractSubMeshes(meshProjected,
                               [](auto hhsBegin, auto hhsEnd) -> bool {
                                 return std::distance(hhsBegin, hhsEnd) <= 1;
                               },
                               10);
  auto &sub = subs.front();

  int count = 0;
  auto energyFun = [&meshProjected, &cam, &sub, &weights,
                    &count](auto &&fh2planeEq) -> double {
    double e = ComputeEnergy(meshProjected, sub, cam, weights, fh2planeEq);
    if ((count++) % 30000 == 1) {
      std::cout << "energy: " << e << '\n';
    }
    return e;
  };

  double gtEnergy = energyFun([&meshGT, &meshProjected, &projCam,
                               &cam](FaceHandle fh) -> Vec3 {
    Point3 ps[3];
    for (int i = 0; i < 3; i++) {
      auto p3 = meshGT.data(meshGT.topo(meshGT.topo(fh).halfedges[i]).to());
      double depthGT = Distance(p3, projCam.eye());

      auto dir = normalize(cam.direction(meshProjected.data(
          meshProjected.topo(meshProjected.topo(fh).halfedges[i]).to())));
      ps[i] = dir * depthGT;
    }
    auto plane = Plane3From3Points(ps[0], ps[1], ps[2]);
    auto eq = Plane3ToEquation(plane);
    assert(!HasValue(eq, IsInfOrNaN<double>));
    return eq;
  });
  std::cout << "GT Energy: " << gtEnergy << '\n';

  auto planes = Reconstruct(
      meshProjected, sub, energyFun,
      [&cam](const Vec2 &p2d) -> Vec3 { return normalize(cam.direction(p2d)); },
      0, maxIters);
  HandledTable<FaceHandle, Plane3> facePlanes(
      meshProjected.internalFaces().size());
  auto meshReconstructed = Transform(
      meshProjected, [](const Point2 &) -> Point3 { return Point3(); });
  for (auto &fhplane : planes) {
    facePlanes[fhplane.first] = fhplane.second;
  }

  // compute corners using planes
  for (auto &vh : sub.vhs) {
    for (auto hh : meshProjected.topo(vh).halfedges) {
      FaceHandle fh = meshProjected.topo(hh).face;
      auto &plane = facePlanes[fh];
      auto p = Intersection(
          Ray3(cam.eye(), cam.direction(meshProjected.data(vh))), plane);
      meshReconstructed.data(vh) += p;
    }
    meshReconstructed.data(vh) /=
        (double)meshProjected.topo(vh).halfedges.size();
  }

  {
    gui::SceneBuilder sb;
    sb.installingOptions().defaultShaderSource =
        gui::OpenGLShaderSourceDescriptor::XTriangles;
    sb.installingOptions().discretizeOptions.color(gui::Black);
    sb.installingOptions().lineWidth = 0.1;
    const gui::ColorTable ctable = gui::CreateRandomColorTableWithSize(
        meshReconstructed.internalFaces().size());
    AddToScene(sb, meshReconstructed, ConstantFunctor<gui::Color>(gui::Black),
               [&ctable](FaceHandle fh) { return ctable[fh.id]; });
    sb.show(true, true, gui::RenderOptions()
                            .backgroundColor(gui::White)
                            .renderMode(gui::All)
                            .bwTexColor(0.0)
                            .bwColor(1.0)
                            .fixUpDirectionInCameraMove(false)
                            .cullBackFace(false)
                            .cullFrontFace(false));
  }
}

int DISABLED_main(int argc, char **argv) {
  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath("D:\\Panoramix\\LineDrawing\\");

  Mesh3 meshGT;
  //MakeTetrahedron(meshGT);
  //MakeQuadFacedCube(meshGT);
  //MakePrism(meshGT, 5, 2.5);
  MakeCone(meshGT, 5, 1);
  //MakeStarPrism(meshGT, 6, 0.5, 1.0, 2.0);

  EnergyWeights weights;
  weights.vertMSDAWeight = 1.0;
  weights.faceMSDAWeight = 1.0;
  weights.faceAngleWeight = 1.0;

  TestOnMesh(meshGT, 30000, weights);

  return 0;
}
