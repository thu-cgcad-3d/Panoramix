#include "containers.hpp"

#include "panoramix.hpp"
#include <ctime>
#include <tuple>

template <class CameraT>
std::vector<Imagei> GTFaceLabels(const PILayoutAnnotation &anno,
                                 const std::vector<CameraT> &testCams) {
  auto impath = anno.impath;
  std::vector<Polygon3> polygons(anno.nfaces());
  for (int i = 0; i < anno.nfaces(); i++) {
    auto &plane = anno.face2plane[i];
    auto &poly = polygons[i];
    poly.normal = plane.normal;
    for (int c : anno.face2corners[i]) {
      Ray3 ray(Origin(), anno.corners[c]);
      poly.corners.push_back(Intersection(ray, plane));
    }
  }

  std::vector<Imagei> faceLabels(testCams.size());
  ParallelRun(
      testCams.size(), std::thread::hardware_concurrency() - 1, [&](int i) {
        std::cout << "computing gt face label map - " << i << std::endl;
        auto &cam = testCams[i];
        Imagei faceLabelMap(cam.screenSize(), -1);
        for (auto it = faceLabelMap.begin(); it != faceLabelMap.end(); ++it) {
          Vec3 direction = normalize(cam.toSpace(it.pos()));
          Ray3 ray(Origin(), direction);
          double depth = std::numeric_limits<double>::infinity();
          for (int i = 0; i < anno.nfaces(); i++) {
            auto &poly = polygons[i];
            auto inter = IntersectionOfLineAndPolygon(ray, poly);
            if (!inter.failed()) {
              *it = i;
              break;
            }
          }
        }
        for (auto it = faceLabelMap.begin(); it != faceLabelMap.end(); ++it) {
          if (*it == -1) {
            for (int dx = -1; dx <= 1; dx++) {
              if (*it != -1) {
                break;
              }
              for (int dy = -1; dy <= 1; dy++) {
                auto p = it.pos() + Pixel(dx, dy);
                if (Contains(faceLabelMap.size(), p) && faceLabelMap(p) != -1) {
                  *it = faceLabelMap(p);
                  break;
                }
              }
            }
          }
        }
        faceLabels[i] = faceLabelMap;
      });

  return faceLabels;
}

std::vector<Imagei>
OrientationMaps(const PILayoutAnnotation &anno,
                const std::vector<PerspectiveCamera> &testCams) {
  std::vector<Imagei> oms(testCams.size());
  ParallelRun(
      testCams.size(), std::thread::hardware_concurrency() - 1, [&](int i) {
        std::vector<HPoint2> hvps(anno.vps.size());
        for (int k = 0; k < anno.vps.size(); k++) {
          hvps[k] = testCams[i].toScreenInHPoint(anno.vps[k]);
        }
        // collect line 2ds in this view
        auto pim = anno.view.sampled(testCams[i]).image;
        LineSegmentExtractor lineExtractor;
        lineExtractor.params().algorithm = LineSegmentExtractor::LSD;
        auto ls = lineExtractor(pim); // use pyramid

        auto rawLines2 = ClassifyEachAs(ls, -1);
        ClassifyLines(rawLines2, hvps);

        Imagei omap =
            ComputeOrientationMaps(rawLines2, hvps, testCams[i].screenSize());
        oms[i] = omap;
      });
  return oms;
}

std::vector<Image7d>
GeometricContext(const PILayoutAnnotation &anno,
                 const std::vector<PerspectiveCamera> &testCams,
                 misc::Matlab &matlab) {
  auto up = normalize(anno.vps[anno.vertVPId]);
  if (up.dot(-anno.view.camera.up()) < 0) {
    up = -up;
  }
  std::vector<Image7d> rawgcs(testCams.size());
  for (int i = 0; i < testCams.size(); i++) {
    auto pim = anno.view.sampled(testCams[i]).image;
    rawgcs[i] = ComputeRawIndoorGeometricContextHedau(matlab, pim);
  }
  return rawgcs;
}

// evaluation helpers
std::vector<int>
FilterHorizontalCameras(const PILayoutAnnotation &anno,
                        const std::vector<PerspectiveCamera> &cams) {
  std::vector<int> selected;
  for (int i = 0; i < cams.size(); i++) {
    if (AngleBetweenUndirected(
            anno.view.camera.up(), cams[i].forward()) < DegreesToRadians(60)) {
      continue;
    }
    selected.push_back(i);
  }
  return selected;
}

double ErrorOfSurfaceNormal(const Vec3 &gt, const Vec3 &cand) {
  if (norm(cand) == 0) {
    return M_PI_2;
  }
  return AngleBetweenUndirected(gt, cand);
}

int GeometricContextToLabel(const Vec7 &labelValues, double clutterThres) {
  // 0: front, 1: left, 2: right, 3: floor, 4: ceiling, 5: clutter, 6: unknown
  if (labelValues[5] >= clutterThres) {
    return -1;
  }
  int labels[5];
  std::iota(labels, labels + 5, 0);
  std::sort(labels, labels + 5, [&labelValues](int a, int b) {
    return labelValues[a] > labelValues[b];
  });
  int maxid = labels[0];
  // if (labelValues[maxid] < failIfLabelScoreLessThan) {
  //    return -1;
  //}
  return maxid;
}

int SurfaceNormalToLabel(const Vec3 &normal, const Vec3 &up,
                         const Vec3 &forward, const Vec3 &left) {
  if (norm(normal) == 0) {
    return -1;
  }
  // 0: front, 1: left, 2: right, 3: floor, 4: ceiling, 5: clutter, 6: unknown
  Vec3 dirs[] = {forward, left, -left, -up, up};
  double angleToDirs[5];
  for (int i = 0; i < 5; i++) {
    angleToDirs[i] = AngleBetweenDirected(dirs[i], normal);
  }
  return std::min_element(std::begin(angleToDirs), std::end(angleToDirs)) -
         std::begin(angleToDirs);
}

int SurfaceNormalToLabel(const Vec3 &normal, const Vec3 &up,
                         const Vec3 &forward, const Vec3 &left,
                         bool pixelOnLeft, bool pixelOnTop) {
  if (norm(normal) == 0) {
    return -1;
  }
  // 0: front, 1: left, 2: right, 3: floor, 4: ceiling, 5: clutter, 6: unknown
  Vec3 dirs[] = {forward, left, -left, -up, up};
  double angleToDirs[5];
  for (int i = 0; i < 5; i++) {
    angleToDirs[i] = AngleBetweenUndirected(dirs[i], normal);
  }
  int result =
      std::min_element(std::begin(angleToDirs), std::end(angleToDirs)) -
      std::begin(angleToDirs);
  if (result == 1 || result == 2) {
    result = pixelOnLeft ? 1 : 2;
  }
  if (result == 3 || result == 4) {
    result = pixelOnTop ? 4 : 3;
  }
  return result;
}

int main(int argc, char **argv) {

  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath("E:\\STORAGE\\CACHE\\Panoramix\\");
  misc::Matlab matlab;

  using TaskQueue =
      std::vector<std::function<misc::MXA(const std::string &impath)>>;
  TaskQueue activeQ;

  TaskQueue Qtest;
  Qtest.push_back([&matlab](const std::string &impath) -> misc::MXA {
    auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    anno.impath = impath;

    PanoramixOptions options;
    options.useWallPrior = true;
    options.usePrincipleDirectionPrior = true;
    options.useGeometricContextPrior = true;

    options.useGTOcclusions = false;
    options.looseLinesSecondTime = false;
    options.looseSegsSecondTime = false;
    options.restrictSegsSecondTime = false;

    options.notUseOcclusions = false;
    options.notUseCoplanarity = false;

    options.refresh_preparation = true;
    options.refresh_mg_init = options.refresh_preparation || false;
    options.refresh_mg_oriented = options.refresh_mg_init || false;
    options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
    options.refresh_lsw = options.refresh_mg_oriented || false;
    options.refresh_mg_occdetected =
        options.refresh_lsw || options.refresh_line2leftRightSegs || false;
    options.refresh_mg_reconstructed = options.refresh_mg_occdetected || false;

    RunPanoramix(anno, options, matlab, false).print();
    return misc::MXA();
  });

  TaskQueue Qstore;

  // set tasks
  // default
  Qstore.push_back([&matlab](const std::string &impath) -> misc::MXA {
    auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    anno.impath = impath;

    PanoramixOptions options;
    options.useWallPrior = true;
    options.usePrincipleDirectionPrior = true;
    options.useGeometricContextPrior = true;

    options.useGTOcclusions = false;
    options.looseLinesSecondTime = false;
    options.looseSegsSecondTime = false;
    options.restrictSegsSecondTime = false;

    options.notUseOcclusions = false;
    options.notUseCoplanarity = false;

    options.refresh_preparation = false;
    options.refresh_mg_init = options.refresh_preparation || true;
    options.refresh_mg_oriented = options.refresh_mg_init || false;
    options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
    options.refresh_lsw = options.refresh_mg_oriented || false;
    options.refresh_mg_occdetected =
        options.refresh_lsw || options.refresh_line2leftRightSegs || false;
    options.refresh_mg_reconstructed = options.refresh_mg_occdetected || false;

    RunPanoramix(anno, options, matlab, false);
    return misc::MXA();
  });
  // use gt occ
  Qstore.push_back([&matlab](const std::string &impath) -> misc::MXA {
    auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    anno.impath = impath;

    PanoramixOptions options;
    options.useWallPrior = true;
    options.usePrincipleDirectionPrior = true;
    options.useGeometricContextPrior = true;

    options.useGTOcclusions = true;
    options.looseLinesSecondTime = false;
    options.looseSegsSecondTime = false;
    options.restrictSegsSecondTime = false;

    options.notUseOcclusions = false;
    options.notUseCoplanarity = false;

    options.refresh_preparation = false;
    options.refresh_mg_init = options.refresh_preparation || true;
    options.refresh_mg_oriented = options.refresh_mg_init || false;
    options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
    options.refresh_lsw = options.refresh_mg_oriented || false;
    options.refresh_mg_occdetected =
        options.refresh_lsw || options.refresh_line2leftRightSegs || false;
    options.refresh_mg_reconstructed = options.refresh_mg_occdetected || false;

    RunPanoramix(anno, options, matlab, false);
    return misc::MXA();
  });
  // without wall
  Qstore.push_back([&matlab](const std::string &impath) -> misc::MXA {
    auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    anno.impath = impath;

    PanoramixOptions options;
    options.useWallPrior = false;
    options.usePrincipleDirectionPrior = true;
    options.useGeometricContextPrior = true;

    options.useGTOcclusions = false;
    options.looseLinesSecondTime = false;
    options.looseSegsSecondTime = false;
    options.restrictSegsSecondTime = false;

    options.notUseOcclusions = false;
    options.notUseCoplanarity = false;

    options.refresh_preparation = false;
    options.refresh_mg_init = options.refresh_preparation || true;
    options.refresh_mg_oriented = options.refresh_mg_init || false;
    options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
    options.refresh_lsw = options.refresh_mg_oriented || false;
    options.refresh_mg_occdetected =
        options.refresh_lsw || options.refresh_line2leftRightSegs || false;
    options.refresh_mg_reconstructed = options.refresh_mg_occdetected || false;

    RunPanoramix(anno, options, matlab, false);
    return misc::MXA();
  });
  // without principle direction
  Qstore.push_back([&matlab](const std::string &impath) -> misc::MXA {
    auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    anno.impath = impath;

    PanoramixOptions options;
    options.useWallPrior = true;
    options.usePrincipleDirectionPrior = false;
    options.useGeometricContextPrior = true;

    options.useGTOcclusions = false;
    options.looseLinesSecondTime = false;
    options.looseSegsSecondTime = false;
    options.restrictSegsSecondTime = false;

    options.notUseOcclusions = false;
    options.notUseCoplanarity = false;

    options.refresh_preparation = false;
    options.refresh_mg_init = options.refresh_preparation || true;
    options.refresh_mg_oriented = options.refresh_mg_init || false;
    options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
    options.refresh_lsw = options.refresh_mg_oriented || false;
    options.refresh_mg_occdetected =
        options.refresh_lsw || options.refresh_line2leftRightSegs || false;
    options.refresh_mg_reconstructed = options.refresh_mg_occdetected || false;

    RunPanoramix(anno, options, matlab, false);
    return misc::MXA();
  });
  // without gc!
  Qstore.push_back([&matlab](const std::string &impath) -> misc::MXA {
    auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    anno.impath = impath;

    PanoramixOptions options;
    options.useWallPrior = true;
    options.usePrincipleDirectionPrior = true;
    options.useGeometricContextPrior = false;

    options.useGTOcclusions = false;
    options.looseLinesSecondTime = false;
    options.looseSegsSecondTime = false;
    options.restrictSegsSecondTime = false;

    options.notUseOcclusions = false;
    options.notUseCoplanarity = false;

    options.refresh_preparation = false;
    options.refresh_mg_init = options.refresh_preparation || true;
    options.refresh_mg_oriented = options.refresh_mg_init || false;
    options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
    options.refresh_lsw = options.refresh_mg_oriented || false;
    options.refresh_mg_occdetected =
        options.refresh_lsw || options.refresh_line2leftRightSegs || false;
    options.refresh_mg_reconstructed = options.refresh_mg_occdetected || false;

    RunPanoramix(anno, options, matlab, false);
    return misc::MXA();
  });
  // without gc + gt occ !
  Qstore.push_back([&matlab](const std::string &impath) -> misc::MXA {
    auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    anno.impath = impath;

    PanoramixOptions options;
    options.useWallPrior = true;
    options.usePrincipleDirectionPrior = true;
    options.useGeometricContextPrior = false;

    options.useGTOcclusions = true;
    options.looseLinesSecondTime = false;
    options.looseSegsSecondTime = false;
    options.restrictSegsSecondTime = false;

    options.notUseOcclusions = false;
    options.notUseCoplanarity = false;

    options.refresh_preparation = false;
    options.refresh_mg_init = options.refresh_preparation || true;
    options.refresh_mg_oriented = options.refresh_mg_init || false;
    options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
    options.refresh_lsw = options.refresh_mg_oriented || false;
    options.refresh_mg_occdetected =
        options.refresh_lsw || options.refresh_line2leftRightSegs || false;
    options.refresh_mg_reconstructed = options.refresh_mg_occdetected || false;

    RunPanoramix(anno, options, matlab, false);
    return misc::MXA();
  });
  // only gc!
  Qstore.push_back([&matlab](const std::string &impath) -> misc::MXA {
    auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    anno.impath = impath;

    PanoramixOptions options;
    options.useWallPrior = false;
    options.usePrincipleDirectionPrior = false;
    options.useGeometricContextPrior = true;

    options.useGTOcclusions = false;
    options.looseLinesSecondTime = false;
    options.looseSegsSecondTime = false;
    options.restrictSegsSecondTime = false;

    options.notUseOcclusions = false;
    options.notUseCoplanarity = false;

    options.refresh_preparation = false;
    options.refresh_mg_init = options.refresh_preparation || true;
    options.refresh_mg_oriented = options.refresh_mg_init || false;
    options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
    options.refresh_lsw = options.refresh_mg_oriented || false;
    options.refresh_mg_occdetected =
        options.refresh_lsw || options.refresh_line2leftRightSegs || false;
    options.refresh_mg_reconstructed = options.refresh_mg_occdetected || false;

    RunPanoramix(anno, options, matlab, false);
    return misc::MXA();
  });
  // only gc + gt occ!
  Qstore.push_back([&matlab](const std::string &impath) -> misc::MXA {
    auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    anno.impath = impath;

    PanoramixOptions options;
    options.useWallPrior = false;
    options.usePrincipleDirectionPrior = false;
    options.useGeometricContextPrior = true;

    options.useGTOcclusions = true;
    options.looseLinesSecondTime = false;
    options.looseSegsSecondTime = false;
    options.restrictSegsSecondTime = false;

    options.notUseOcclusions = false;
    options.notUseCoplanarity = false;

    options.refresh_preparation = false;
    options.refresh_mg_init = options.refresh_preparation || true;
    options.refresh_mg_oriented = options.refresh_mg_init || false;
    options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
    options.refresh_lsw = options.refresh_mg_oriented || false;
    options.refresh_mg_occdetected =
        options.refresh_lsw || options.refresh_line2leftRightSegs || false;
    options.refresh_mg_reconstructed = options.refresh_mg_occdetected || false;

    RunPanoramix(anno, options, matlab, false);
    return misc::MXA();
  });

  // no coplanarity
  TaskQueue QstoreNoCop;
  QstoreNoCop.push_back([&matlab](const std::string &impath) -> misc::MXA {
    auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    anno.impath = impath;

    PanoramixOptions options;
    options.useWallPrior = true;
    options.usePrincipleDirectionPrior = true;
    options.useGeometricContextPrior = true;

    options.useGTOcclusions = false;
    options.looseLinesSecondTime = false;
    options.looseSegsSecondTime = false;
    options.restrictSegsSecondTime = false;

    options.notUseOcclusions = false;
    options.notUseCoplanarity = true;

    options.refresh_preparation = false;
    options.refresh_mg_init = options.refresh_preparation || true;
    options.refresh_mg_oriented = options.refresh_mg_init || false;
    options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
    options.refresh_lsw = options.refresh_mg_oriented || false;
    options.refresh_mg_occdetected =
        options.refresh_lsw || options.refresh_line2leftRightSegs || false;
    options.refresh_mg_reconstructed = options.refresh_mg_occdetected || false;

    RunPanoramix(anno, options, matlab, false);
    return misc::MXA();
  });


  TaskQueue QshowModel;
  QshowModel.push_back([&](const std::string &impath) -> misc::MXA {
    auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    anno.impath = impath;

    PanoramixOptions options;
    options.useWallPrior = true;
    options.usePrincipleDirectionPrior = true;
    options.useGeometricContextPrior = true;

    options.useGTOcclusions = false;
    options.looseLinesSecondTime = false;
    options.looseSegsSecondTime = false;
    options.restrictSegsSecondTime = false;

    options.notUseOcclusions = false;
    options.notUseCoplanarity = false;

    options.refresh_preparation = false;
    options.refresh_mg_init = options.refresh_preparation || false;
    options.refresh_mg_oriented = options.refresh_mg_init || false;
    options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
    options.refresh_lsw = options.refresh_mg_oriented || false;
    options.refresh_mg_occdetected =
        options.refresh_lsw || options.refresh_line2leftRightSegs || false;
    options.refresh_mg_reconstructed = options.refresh_mg_occdetected || true;

    RunPanoramix(anno, options, matlab, true, true);
    return misc::MXA();
  });

  TaskQueue Qdata2;
  Qdata2.push_back([&](const std::string &impath) -> misc::MXA {
    auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    // gt1 -> gt2
    {
      std::vector<PerspectiveCamera> testCams1;
      std::vector<Imagei> gt1;
      misc::LoadCache(impath, "testCams1_gt1", testCams1, gt1);

      std::vector<PerspectiveCamera> testCams2;
      std::vector<Imagei> gt2;
      auto selected = FilterHorizontalCameras(anno, testCams1);
      for (int i : selected) {
        testCams2.push_back(testCams1[i]);
        gt2.push_back(gt1[i]);
      }
      std::cout << "cam num: " << selected.size() << std::endl;
      misc::SaveCache(impath, "testCams2_gt2", testCams2, gt2);
    }

    // gc1 -> gc2
    {
      std::vector<PerspectiveCamera> testCams1;
      std::vector<Image7d> gc1;
      misc::LoadCache(impath, "testCams1_gc1", testCams1, gc1);

      std::vector<PerspectiveCamera> testCams2;
      std::vector<Image7d> gc2;
      auto selected = FilterHorizontalCameras(anno, testCams1);
      for (int i : selected) {
        testCams2.push_back(testCams1[i]);
        gc2.push_back(gc1[i]);
      }
      std::cout << "cam num: " << selected.size() << std::endl;
      misc::SaveCache(impath, "testCams2_gc2", testCams2, gc2);
    }

    // om1 -> om2
    {
      std::vector<PerspectiveCamera> testCams1;
      std::vector<Imagei> om1;
      misc::LoadCache(impath, "testCams1_om1", testCams1, om1);

      std::vector<PerspectiveCamera> testCams2;
      std::vector<Imagei> om2;
      auto selected = FilterHorizontalCameras(anno, testCams1);
      for (int i : selected) {
        testCams2.push_back(testCams1[i]);
        om2.push_back(om1[i]);
      }
      std::cout << "cam num: " << selected.size() << std::endl;
      misc::SaveCache(impath, "testCams2_om2", testCams2, om2);
    }

    return misc::MXA();
  });

  // cached gt
  auto getCachedGT = [](const std::string &impath) {
    std::vector<PerspectiveCamera> testCams2;
    std::vector<Imagei> gt2;
    misc::LoadCache(impath, "testCams2_gt2", testCams2, gt2);
    return gt2;
  };

  // cache gc
  auto getCachedGC = [&matlab](const std::string &impath) {
    std::vector<PerspectiveCamera> testCams2;
    std::vector<Image7d> gc2;
    misc::LoadCache(impath, "testCams2_gc2", testCams2, gc2);
    return gc2;
  };

  // cache om
  auto getCachedOM = [](const std::string &impath) {
    std::vector<PerspectiveCamera> testCams2;
    std::vector<Imagei> om2;
    misc::LoadCache(impath, "testCams2_om2", testCams2, om2);
    return om2;
  };

  TaskQueue QcompareLabels;
  QcompareLabels.push_back([&](const std::string &impath) -> misc::MXA {
    auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    anno.impath = impath;

    bool consider_horizontal_cams_only = true;

    PanoramixOptions options;
    options.useWallPrior = true;
    options.usePrincipleDirectionPrior = true;
    options.useGeometricContextPrior = false;

    options.useGTOcclusions = false;
    options.looseLinesSecondTime = false;
    options.looseSegsSecondTime = false;
    options.restrictSegsSecondTime = false;

    options.notUseOcclusions = false;
    options.notUseCoplanarity = false;

    options.refresh_preparation = false;
    options.refresh_mg_init = options.refresh_preparation || false;
    options.refresh_mg_oriented = options.refresh_mg_init || false;
    options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
    options.refresh_lsw = options.refresh_mg_oriented || false;
    options.refresh_mg_occdetected =
        options.refresh_lsw || options.refresh_line2leftRightSegs || false;
    options.refresh_mg_reconstructed = options.refresh_mg_occdetected || true;

    std::vector<PerspectiveCamera> testCams;
    std::vector<Imagei> gtData;
    misc::LoadCache(impath, "testCams2_gt2", testCams, gtData);
    std::vector<Image7d> gcData = getCachedGC(impath);
    if (gcData.size() == testCams.size() && testCams.size() == gtData.size()) {
      std::cout << "size matched!!!" << std::endl;
    }

    auto panoramixResults =
        GetSurfaceNormalMapsOfPanoramix(testCams, anno, options, matlab);

    double error_pn = 0.0;
    double error_gc = 0.0;

    int npixels = 0;
    std::vector<double> completenessTable(testCams.size(), 0);
    for (int i = 0; i < testCams.size(); i++) {
      auto &pn = panoramixResults[i];

      auto &gc = gcData[i];
      auto &gt = gtData[i];

      auto &cam = testCams[i];

      Imagei gcLabels(cam.screenSize(), -1);
      Imagei pnLabels(cam.screenSize(), -1);
      Imagei gtLabels(cam.screenSize(), -1);

      double &completeness = completenessTable[i];
      for (auto it = gt.begin(); it != gt.end(); ++it) {
        auto p = it.pos();
        int gtFaceId = *it;
        if (gtFaceId == -1) {
          std::cout << "!!!gt face id is -1!!!!" << std::endl;
          continue;
        }

        Vec3 gtNormal = anno.face2plane[gtFaceId].normal;
        if (norm(gtNormal) == 0) {
          continue;
        }

        int gtLabel = SurfaceNormalToLabel(gtNormal, cam.up(), cam.forward(),
                                           cam.leftward());
        gtLabels(p) = gtLabel;

        Vec3 pnNormal = pn(p);
        int pnLabel = SurfaceNormalToLabel(pnNormal, cam.up(), cam.forward(),
                                           cam.leftward());

        Vec7 gcScores = gc(p);
        int gcLabel = GeometricContextToLabel(gcScores, 1.0);
        if (gcLabel == -1) {
          continue;
        }

        pnLabels(p) = pnLabel;
        gcLabels(p) = gcLabel;

        error_pn += pnLabel != gtLabel;
        error_gc += gcLabel != gtLabel;

        npixels++;
        completeness += 1.0 / (gt.cols * gt.rows);
      }

      std::cout << "completeness:" << completeness << std::endl;

      if (false) {
        // 0: front, 1: left, 2: right, 3: floor, 4: ceiling, 5: clutter, 6:
        // unknown
        std::vector<gui::Color> colors = {gui::Green, gui::Blue, gui::Blue,
                                          gui::Red,   gui::Red,  gui::White};
        gui::ColorTable ctable(colors);
        ctable.exceptionalColor() = gui::Black;
        const std::string path = "F:\\GitHub\\write-"
                                 "papers\\papers\\a\\figure\\experiments\\compa"
                                 "regc\\";
        // gui::MakeCanvas((Image3ub)anno.view.sampled(cam).image).saveAs(path +
        // "im" + std::to_string(i) + ".png");
        gui::MakeCanvas(ctable(gcLabels))
            .show(1, "gcLabels"); //.saveAs(path + "gc" + std::to_string(i) +
                                  //".png");
        gui::MakeCanvas(ctable(pnLabels))
            .show(1, "pnLabels"); // .saveAs(path + "pn" + std::to_string(i) +
                                  // ".png");
        gui::MakeCanvas(ctable(gtLabels))
            .show(0, "gtLabels"); // .saveAs(path + "gt" + std::to_string(i) +
                                  // ".png");
      }
    }

    error_pn /= npixels;
    error_gc /= npixels;

    std::cout << std::endl << std::endl;
    std::cout << "{{{{{{{{{  PN[" << error_pn << "] vs GC[" << error_gc
              << "]  }}}}}}}}}}}} on image\"" << impath << "\"" << std::endl;
    std::cout << std::endl << std::endl;

    DenseMatd pn_gc_error(3, 1, 0.0);
    pn_gc_error(0, 0) = error_pn;
    pn_gc_error(1, 0) = error_gc;
    pn_gc_error(2, 0) = std::accumulate(completenessTable.begin(),
                                        completenessTable.end(), 0.0) /
                        testCams.size();
    return misc::MXA(pn_gc_error, false);
  });

  TaskQueue QcompareOMLabels;
  QcompareOMLabels.push_back([&](const std::string &impath) -> misc::MXA {
    auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    anno.impath = impath;

    PanoramixOptions options;
    options.useWallPrior = true;
    options.usePrincipleDirectionPrior = true;
    options.useGeometricContextPrior = true;

    options.useGTOcclusions = false;
    options.looseLinesSecondTime = false;
    options.looseSegsSecondTime = false;
    options.restrictSegsSecondTime = false;

    options.notUseOcclusions = false;
    options.notUseCoplanarity = false;

    options.refresh_preparation = false;
    options.refresh_mg_init = options.refresh_preparation || false;
    options.refresh_mg_oriented = options.refresh_mg_init || false;
    options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
    options.refresh_lsw = options.refresh_mg_oriented || false;
    options.refresh_mg_occdetected =
        options.refresh_lsw || options.refresh_line2leftRightSegs || false;
    options.refresh_mg_reconstructed = options.refresh_mg_occdetected || false;

    std::vector<PerspectiveCamera> testCams;
    std::vector<Imagei> gtData;
    misc::LoadCache(impath, "testCams2_gt2", testCams, gtData);
    auto omData = getCachedOM(impath);
    if (omData.size() == testCams.size() && testCams.size() == gtData.size()) {
      std::cout << "size matched!!!" << std::endl;
    }

    auto panoramixResults =
        GetSurfaceNormalMapsOfPanoramix(testCams, anno, options, matlab);

    double error_pn = 0.0;
    double error_om = 0.0;

    int npixels = 0;
    std::vector<double> completenessTable(testCams.size(), 0);
    for (int i = 0; i < testCams.size(); i++) {
      auto &pn = panoramixResults[i];

      auto &om = omData[i];
      auto &gt = gtData[i];

      auto &cam = testCams[i];

      Imagei omLabels(cam.screenSize(), -1);
      Imagei pnLabels(cam.screenSize(), -1);
      Imagei gtLabels(cam.screenSize(), -1);

      double &completeness = completenessTable[i];
      for (auto it = gt.begin(); it != gt.end(); ++it) {
        auto p = it.pos();
        bool onLeft = p.x < gt.cols / 2;
        bool onTop = p.y < gt.rows / 2;

        int gtFaceId = *it;
        if (gtFaceId == -1) {
          std::cout << "!!!gt face id is -1!!!!" << std::endl;
          continue;
        }

        Vec3 gtNormal = anno.face2plane[gtFaceId].normal;
        if (norm(gtNormal) == 0) {
          continue;
        }

        int gtLabel = SurfaceNormalToLabel(gtNormal, cam.up(), cam.forward(),
                                           cam.leftward());
        gtLabels(p) = gtLabel;

        Vec3 pnNormal = pn(p);
        int pnLabel = SurfaceNormalToLabel(pnNormal, cam.up(), cam.forward(),
                                           cam.leftward());
        pnLabels(p) = pnLabel;

        int omvpid = om(p);
        if (omvpid == -1) {
          continue;
        }
        // om -> surface label
        int omLabel =
            SurfaceNormalToLabel(anno.vps[omvpid], cam.up(), cam.forward(),
                                 cam.leftward(), onLeft, onTop);
        omLabels(p) = omLabel;

        error_pn += pnLabel != gtLabel;
        error_om += omLabel != gtLabel;

        npixels++;
        completeness += 1.0 / (gt.cols * gt.rows);
      }

      std::cout << "completeness:" << completeness << std::endl;

      if (false) {
        // 0: front, 1: left, 2: right, 3: floor, 4: ceiling, 5: clutter, 6:
        // unknown
        std::vector<gui::Color> colors = {gui::Green,     gui::DarkGray,
                                          gui::LightGray, gui::Red,
                                          gui::Red,       gui::White};
        gui::ColorTable ctable(colors);
        ctable.exceptionalColor() = gui::Black;
        // const std::string path =
        // "F:\\GitHub\\write-papers\\papers\\a\\figure\\experiments\\comparegc\\";
        // gui::MakeCanvas((Image3ub)anno.view.sampled(cam).image).saveAs(path +
        // "im" + std::to_string(i) + ".png");
        gui::MakeCanvas(ctable(omLabels))
            .show(1, "omLabels"); //.saveAs(path + "gc" + std::to_string(i) +
                                  //".png");
        gui::MakeCanvas(ctable(pnLabels))
            .show(1, "pnLabels"); // .saveAs(path + "pn" + std::to_string(i) +
                                  // ".png");
        gui::MakeCanvas(ctable(gtLabels))
            .show(0, "gtLabels"); // .saveAs(path + "gt" + std::to_string(i) +
                                  // ".png");
      }
    }

    error_pn /= npixels;
    error_om /= npixels;

    std::cout << std::endl << std::endl;
    std::cout << "{{{{{{{{{  PN[" << error_pn << "] vs OM[" << error_om
              << "]  }}}}}}}}}}}} on image\"" << impath << "\"" << std::endl;
    std::cout << std::endl << std::endl;

    auto results = misc::MXA::createStructMatrix(
        1, 1, {"error_pn", "error_om", "completeness"});
    results.setField("error_pn", 0, error_pn);
    results.setField("error_om", 0, error_om);
    results.setField("completeness", 0,
                     std::accumulate(completenessTable.begin(),
                                     completenessTable.end(), 0.0) /
                         testCams.size());
    return results;
  });

  TaskQueue Qstat;
  Qstat.push_back([&](const std::string &impath) -> misc::MXA {
    auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    anno.impath = impath;

    PanoramixOptions options;
    options.useWallPrior = true;
    options.usePrincipleDirectionPrior = true;
    options.useGeometricContextPrior = true;

    options.useGTOcclusions = false;
    options.looseLinesSecondTime = false;
    options.looseSegsSecondTime = false;
    options.restrictSegsSecondTime = false;

    options.notUseOcclusions = false;
    options.notUseCoplanarity = true;

    options.refresh_preparation = false;
    options.refresh_mg_init = options.refresh_preparation || false;
    options.refresh_mg_oriented = options.refresh_mg_init || false;
    options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
    options.refresh_lsw = options.refresh_mg_oriented || false;
    options.refresh_mg_occdetected =
        options.refresh_lsw || options.refresh_line2leftRightSegs || false;
    options.refresh_mg_reconstructed = options.refresh_mg_occdetected || false;

    PIGraph<PanoramicCamera> mg;
    PIConstraintGraph cg;
    PICGDeterminablePart dp;

   
    bool succ = GetPanoramixResult(anno, options, mg, cg, dp);
    if (!succ) {
        std::cout << "failed in " << impath << std::endl;
        return misc::MXA();
    }

    std::map<std::string, double> resultData;

    resultData["linesNum"] = mg.nlines();
    for (int i = 0; i < mg.nlines(); i++) {
      if (mg.lines[i].claz == -1) {
        resultData["line2sNum"]++;
      } else {
        resultData["line1sNum"]++;
      }
    }
    resultData["segsNum"] = mg.nsegs;
    for (int i = 0; i < mg.nsegs; i++) {
      int dof = mg.seg2control[i].dof();
      if (dof == 1) {
        resultData["seg1sNum"]++;
      } else if (dof == 2) {
        resultData["seg2sNum"]++;
      } else {
        resultData["seg3sNum"]++;
      }
    }

    resultData["entsNum"] = cg.entities.size();
    resultData["consNum"] = cg.constraints.size();

    resultData["determinableEntNum"] = dp.determinableEnts.size();
    resultData["consNumBetweenDeterminableEnts"] =
        dp.consBetweenDeterminableEnts.size();

    for (int ent : dp.determinableEnts) {
      if (cg.entities[ent].isLine()) {
        if (cg.entities[ent].supportingPlane.dof == 1) {
          resultData["determinableLine1sNum"]++;
        } else if (cg.entities[ent].supportingPlane.dof == 2) {
          resultData["determinableLine2sNum"]++;
        }
      } else {
        if (cg.entities[ent].supportingPlane.dof == 1) {
          resultData["determinableSeg1sNum"]++;
        } else if (cg.entities[ent].supportingPlane.dof == 2) {
          resultData["determinableSeg2sNum"]++;
        } else if (cg.entities[ent].supportingPlane.dof == 3) {
          resultData["determinableSeg3sNum"]++;
        }
        resultData["determinableSegAreaRatio"] +=
            mg.seg2areaRatio[cg.entities[ent].id];
      }
    }

    for (int cons : dp.consBetweenDeterminableEnts) {
      auto &c = cg.constraints[cons];
      if (c.isConnection()) {
        resultData["consConNumBetweenDeterminableEnts"]++;
      } else {
        resultData["consCopNumBetweenDeterminableEnts"]++;
      }
    }

    std::vector<std::string> fieldNames;
    for (auto &p : resultData) {
      fieldNames.push_back(p.first);
    }
    auto result = misc::MXA::createStructMatrix(1, 1, fieldNames, false);
    for (auto &p : resultData) {
      result.setField(p.first, 0, misc::MXA(p.second));
    }
    return result;
  });

  std::vector<std::pair<Point3, Point3>> cuboidCands;
  std::vector<Imagei> cuboidCandsFaceMap;

  if (true) {
    core::LoadFromDisk(
        "F:\\GitHub\\write-papers\\papers\\a\\data\\cuboidcands_5.cereal",
        cuboidCands);
    core::LoadFromDisk("F:\\GitHub\\write-papers\\papers\\a\\data\\cuboidcands_"
                       "facemap_5.cereal",
                       cuboidCandsFaceMap);
  } else {
    std::cout << "CUBOID CANDS NOT LOADED" << std::endl;
  }

  std::cout << "size of cands: " << cuboidCands.size() << std::endl;

  TaskQueue QfindBestCuboids;
  QfindBestCuboids.push_back([&](const std::string &impath) -> misc::MXA {
    misc::Clock clock = "processing" + impath;

    auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    anno.impath = impath;

    // load gt and test cams
    std::vector<PerspectiveCamera> testCams;
    std::vector<Imagei> gtData;
    misc::LoadCache(impath, "testCams2_gt2", testCams, gtData);

    // create gtSurfaceLabelMaps
    std::vector<Imagei> gtSurfaceLabelMaps(gtData.size());
    for (int c = 0; c < testCams.size(); c++) {
      auto &cam = testCams[c];
      auto &gt = gtData[c];
      Imagei gtSurfaceLabelMap(gt.size(), -1);
      for (auto it = gt.begin(); it != gt.end(); ++it) {
        auto p = it.pos();
        int gtFaceId = *it;
        if (gtFaceId == -1) {
          std::cout << "!!!gt face id is -1!!!!" << std::endl;
          continue;
        }
        Vec3 gtNormal = anno.face2plane[gtFaceId].normal;
        if (norm(gtNormal) == 0) {
          continue;
        }

        int gtLabel = SurfaceNormalToLabel(gtNormal, cam.up(), cam.forward(),
                                           cam.leftward());
        gtSurfaceLabelMap(p) = gtLabel;
      }
      gtSurfaceLabelMaps[c] = gtSurfaceLabelMap;
    }

    auto &vps = anno.vps;
    std::vector<Vec3> faceId2VP = {-vps[1], -vps[0], -vps[2],
                                   vps[0],  vps[2],  vps[1]};

    DenseMati faceId_camId_to_label(6, testCams.size(), -1);
    for (int faceid = 0; faceid < 6; faceid++) {
      for (int camid = 0; camid < testCams.size(); camid++) {
        auto &cam = testCams[camid];
        faceId_camId_to_label(faceid, camid) = SurfaceNormalToLabel(
            faceId2VP[faceid], cam.up(), cam.forward(), cam.leftward());
      }
    }

    // select best cuboid
    std::vector<double> cuboidErrorRatios(cuboidCands.size(), 0.0);
    PanoramicCamera cuboidCandCam =
        CreatePanoramicCamera(cuboidCandsFaceMap[0]);
    ParallelRun(
        cuboidCands.size(), std::thread::hardware_concurrency() - 1, 100,
        [&cuboidCands, &cuboidCandsFaceMap, &faceId_camId_to_label,
         &cuboidErrorRatios, &anno, &testCams, &gtSurfaceLabelMaps,
         &cuboidCandCam](int i) {
          auto &cuboidCandFaceMap = cuboidCandsFaceMap[i];

          auto vps = anno.vps;
          for (auto &vp : vps) {
            vp /= norm(vp);
          }

          int errornum = 0;
          int allnum = 0;
          for (int c = 0; c < testCams.size(); c++) {
            auto &cam = testCams[c];
            auto &gtSurfaceLabelMap = gtSurfaceLabelMaps[c];

            static const int sampleStep = 10;
            for (int y = 0; y < gtSurfaceLabelMap.rows; y += sampleStep) {
              for (int x = 0; x < gtSurfaceLabelMap.cols; x += sampleStep) {

                auto p = Pixel(x, y);
                auto dir = cam.toSpace(p);
                Vec3 dirInVPs;
                for (int k = 0; k < 3; k++) {
                  dirInVPs[k] = dir.dot(vps[k]);
                }

                // get the face of cuboid cand
                auto pixel = ToPixel(cuboidCandCam.toScreen(dirInVPs));
                pixel.x = WrapBetween(pixel.x, 0, cuboidCandFaceMap.cols);
                pixel.y = BoundBetween(pixel.y, 0, cuboidCandFaceMap.rows - 1);
                int faceid = cuboidCandFaceMap(pixel);
                int cuboidSurfaceLabel = -1;
                if (faceid != -1) {
                  cuboidSurfaceLabel = faceId_camId_to_label(faceid, c);
                }
                int gtSurfaceLabel = gtSurfaceLabelMap(p);
                if (gtSurfaceLabel != cuboidSurfaceLabel) {
                  errornum++;
                }
                allnum++;
              }
            }
          }

          cuboidErrorRatios[i] = errornum / double(allnum);
        });

    int bestCuboidId =
        std::min_element(cuboidErrorRatios.begin(), cuboidErrorRatios.end()) -
        cuboidErrorRatios.begin();
    std::cout << "error ratio: " << cuboidErrorRatios[bestCuboidId]
              << std::endl;
    misc::SaveCache(impath, "bestCuboidId_errorRatio", bestCuboidId,
                    cuboidErrorRatios[bestCuboidId]);
    auto result = misc::MXA::createStructMatrix(
        2, 1, {"bestCuboidId", "errorRatio"}, false);
    result.setField("bestCuboidId", 0, bestCuboidId);
    result.setField("errorRatio", 0, cuboidErrorRatios[bestCuboidId]);
    return result;
  });

  TaskQueue QsaveCuboids;
  QsaveCuboids.push_back([&](const std::string &impath) -> misc::MXA {
    int bestCuboidId = -1;
    double bestErrorRatio = 0.0;
    misc::LoadCache(impath, "bestCuboidId_errorRatio", bestCuboidId,
                    bestErrorRatio);
    auto result = misc::MXA::createStructMatrix(
        1, 1, {"bestCuboidId", "errorRatio"}, false);
    result.setField("bestCuboidId", 0, bestCuboidId);
    result.setField("errorRatio", 0, bestErrorRatio);
    return result;
  });

  TaskQueue QshowBestCuboid;
  QshowBestCuboid.push_back([&](const std::string &impath) -> misc::MXA {
    auto anno = LoadOrInitializeNewLayoutAnnotation(impath);

    // load gt and test cams
    std::vector<PerspectiveCamera> testCams;
    std::vector<Imagei> gtData;
    misc::LoadCache(impath, "testCams2_gt2", testCams, gtData);

    // create gtSurfaceLabelMaps
    std::vector<Imagei> gtSurfaceLabelMaps(gtData.size());
    for (int c = 0; c < testCams.size(); c++) {
      auto &cam = testCams[c];
      auto &gt = gtData[c];
      Imagei gtSurfaceLabelMap(gt.size(), -1);
      for (auto it = gt.begin(); it != gt.end(); ++it) {
        auto p = it.pos();
        int gtFaceId = *it;
        if (gtFaceId == -1) {
          std::cout << "!!!gt face id is -1!!!!" << std::endl;
          continue;
        }
        Vec3 gtNormal = anno.face2plane[gtFaceId].normal;
        if (norm(gtNormal) == 0) {
          continue;
        }

        int gtLabel = SurfaceNormalToLabel(gtNormal, cam.up(), cam.forward(),
                                           cam.leftward());
        gtSurfaceLabelMap(p) = gtLabel;
      }
      gtSurfaceLabelMaps[c] = gtSurfaceLabelMap;
    }

    auto &vps = anno.vps;
    std::vector<Vec3> faceId2VP = {-vps[1], -vps[0], -vps[2],
                                   vps[0],  vps[2],  vps[1]};

    int bestCuboidId = -1;
    double errorRatio = 0.0;
    misc::LoadCache(impath, "bestCuboidId_errorRatio", bestCuboidId,
                    errorRatio);

    std::cout << "bestCuboidId = " << bestCuboidId << std::endl;

    auto cuboidCandFaceView =
        CreatePanoramicView(cuboidCandsFaceMap[bestCuboidId]);

    std::vector<Imagei> surfaceLabelMaps(testCams.size());
    for (int c = 0; c < testCams.size(); c++) {
      auto &cam = testCams[c];
      auto &surfaceLabelMap = surfaceLabelMaps[c];
      auto &gtSurfaceLabelMap = gtSurfaceLabelMaps[c];
      surfaceLabelMap = Imagei(gtSurfaceLabelMap.size(), -1);
      for (auto it = surfaceLabelMap.begin(); it != surfaceLabelMap.end();
           ++it) {
        auto p = it.pos();
        auto dir = cam.toSpace(p);
        Vec3 dirInVPs;
        for (int k = 0; k < 3; k++) {
          dirInVPs[k] = dir.dot(normalize(vps[k]));
        }
        // get the face of cuboid cand
        int faceid = cuboidCandFaceView.image(
            ToPixel(cuboidCandFaceView.camera.toScreen(dirInVPs)));
        int cuboidSurfaceLabel = -1;
        if (faceid != -1) {
          cuboidSurfaceLabel = SurfaceNormalToLabel(
              faceId2VP[faceid], cam.up(), cam.forward(), cam.leftward());
        }
        *it = cuboidSurfaceLabel;
      }

      if (true) {
        // 0: front, 1: left, 2: right, 3: floor, 4: ceiling, 5: clutter, 6:
        // unknown
        std::vector<gui::Color> colors = {gui::Green, gui::Blue, gui::Blue,
                                          gui::Red,   gui::Red,  gui::White};
        gui::ColorTable ctable(colors);
        ctable.exceptionalColor() = gui::Black;
        gui::MakeCanvas(ctable(surfaceLabelMap))
            .show(1, "cuboidsurfaceLabelMap");
        gui::MakeCanvas(ctable(gtSurfaceLabelMap)).show(0, "gtSurfaceLabelMap");
      }
    }

    return misc::MXA();

  });

  TaskQueue QstoreGTPanoFaceIds;
  QstoreGTPanoFaceIds.push_back([&](const std::string &impath) -> misc::MXA {
    auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    PanoramicCamera cam = anno.view.camera;

    Imagei gtPanoFaceIds(cam.screenSize(), -1);
    // misc::LoadCache(impath, "gtpanofaceids", gtPanoFaceIds);
    // gt pano face ids
    {
      std::vector<Polygon3> polygons(anno.nfaces());
      for (int i = 0; i < anno.nfaces(); i++) {
        auto &plane = anno.face2plane[i];
        auto &poly = polygons[i];
        poly.normal = plane.normal;
        for (int c : anno.face2corners[i]) {
          Ray3 ray(Origin(), anno.corners[c]);
          poly.corners.push_back(Intersection(ray, plane));
        }
      }
      for (auto it = gtPanoFaceIds.begin(); it != gtPanoFaceIds.end(); ++it) {
        Vec3 direction = normalize(cam.toSpace(it.pos()));
        Ray3 ray(Origin(), direction);
        double depth = std::numeric_limits<double>::infinity();
        for (int i = 0; i < anno.nfaces(); i++) {
          auto &poly = polygons[i];
          auto inter = IntersectionOfLineAndPolygon(ray, poly);
          if (!inter.failed()) {
            *it = i;
            break;
          }
        }
      }
      for (auto it = gtPanoFaceIds.begin(); it != gtPanoFaceIds.end(); ++it) {
        if (*it == -1) {
          for (int dx = -1; dx <= 1; dx++) {
            if (*it != -1) {
              break;
            }
            for (int dy = -1; dy <= 1; dy++) {
              auto p = it.pos() + Pixel(dx, dy);
              if (Contains(gtPanoFaceIds.size(), p) && gtPanoFaceIds(p) != -1) {
                *it = gtPanoFaceIds(p);
                break;
              }
            }
          }
        }
      }
      misc::SaveCache(impath, "gtpanofaceids", gtPanoFaceIds);
    }

    return misc::MXA();
  });

  TaskQueue QcomputePanoDepth;
  QcomputePanoDepth.push_back([&](const std::string &impath) -> misc::MXA {
    auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    PanoramicCamera cam = anno.view.camera;

    Imaged gtDepths, bestCuboidDepths, pnDepths, pnDepthNoOcc, pnDepthGTOcc;

    // gt depths
    {
      Imagei gtPanoFaceIds(cam.screenSize(), -1);
      misc::LoadCache(impath, "gtpanofaceids", gtPanoFaceIds);
      gtDepths = Imaged(cam.screenSize(), 0.0);
      for (auto it = gtDepths.begin(); it != gtDepths.end(); ++it) {
        int faceId = gtPanoFaceIds(it.pos());
        auto &plane = anno.face2plane[faceId];
        auto dir = normalize(cam.toSpace(it.pos()));
        Ray3 ray(Origin(), dir);
        double depth = norm(Intersection(ray, plane));
        *it = depth;
      }
    }

    // best cuboid depths
    {
      int bestCuboidId = -1;
      double bestErrorRatio = 0.0;
      misc::LoadCache(impath, "bestCuboidId_errorRatio", bestCuboidId,
                      bestErrorRatio);
      auto &p1 = cuboidCands[bestCuboidId].second;
      auto &p2 = cuboidCands[bestCuboidId].first;

      // corners
      Vec3 corners[8] = {{p1[0], p1[1], p1[2]}, {p1[0], p2[1], p1[2]},
                         {p2[0], p2[1], p1[2]}, {p2[0], p1[1], p1[2]},
                         {p1[0], p1[1], p2[2]}, {p1[0], p2[1], p2[2]},
                         {p2[0], p2[1], p2[2]}, {p2[0], p1[1], p2[2]}};
      // faces
      int faces[6][4] = {{1, 4, 8, 5}, {1, 5, 6, 2}, {1, 2, 3, 4},
                         {3, 7, 8, 4}, {5, 8, 7, 6}, {2, 6, 7, 3}};
      Polygon3 polys[6];
      for (int k = 0; k < 6; k++) {
        polys[k].corners.resize(4);
        for (int m = 0; m < 4; m++) {
          polys[k].corners[m] = corners[faces[k][m] - 1];
        }
        polys[k].normal =
            -(polys[k].corners[0] - polys[k].corners[1])
                 .cross(polys[k].corners[1] - polys[k].corners[2]);
        polys[k].normal /= norm(polys[k].normal);
      }

      bestCuboidDepths = Imaged(cam.screenSize(), 0.0);
      for (auto it = bestCuboidDepths.begin(); it != bestCuboidDepths.end();
           ++it) {
        auto p = it.pos();
        auto dir = normalize(cam.toSpace(it.pos()));
        Vec3 dirInVPs;
        for (int k = 0; k < 3; k++) {
          dirInVPs[k] = dir.dot(normalize(anno.vps[k]));
        }

        Ray3 ray(Origin(), dirInVPs);
        *it = std::numeric_limits<double>::infinity();
        for (int faceid = 0; faceid < 6; faceid++) {
          auto test = IntersectionOfLineAndPlane(ray, polys[faceid].plane());
          if (test.ratio > 0 && test.ratio < *it) {
            *it = test.ratio;
          }
        }
      }
    }

    // pn depths
    {
      PanoramixOptions options;
      options.useWallPrior = true;
      options.usePrincipleDirectionPrior = true;
      options.useGeometricContextPrior = true;

      options.useGTOcclusions = false;
      options.looseLinesSecondTime = false;
      options.looseSegsSecondTime = false;
      options.restrictSegsSecondTime = false;

      options.notUseOcclusions = false;
      options.notUseCoplanarity = false;

      options.refresh_preparation = false;
      options.refresh_mg_init = options.refresh_preparation || false;
      options.refresh_mg_oriented = options.refresh_mg_init || false;
      options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
      options.refresh_lsw = options.refresh_mg_oriented || false;
      options.refresh_mg_occdetected =
          options.refresh_lsw || options.refresh_line2leftRightSegs || false;
      options.refresh_mg_reconstructed =
          options.refresh_mg_occdetected || false;

      pnDepths = GetSurfaceDepthMapsOfPanoramix(
                     std::vector<PanoramicCamera>{cam}, anno, options, matlab)
                     .front();
    }

    // pn depth no occ
    {
      PanoramixOptions options;
      options.useWallPrior = true;
      options.usePrincipleDirectionPrior = true;
      options.useGeometricContextPrior = true;

      options.useGTOcclusions = false;
      options.looseLinesSecondTime = false;
      options.looseSegsSecondTime = false;
      options.restrictSegsSecondTime = false;

      options.notUseOcclusions = true;
      options.notUseCoplanarity = false;

      options.refresh_preparation = false;
      options.refresh_mg_init = options.refresh_preparation || false;
      options.refresh_mg_oriented = options.refresh_mg_init || false;
      options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
      options.refresh_lsw = options.refresh_mg_oriented || false;
      options.refresh_mg_occdetected =
          options.refresh_lsw || options.refresh_line2leftRightSegs || false;
      options.refresh_mg_reconstructed =
          options.refresh_mg_occdetected || false;

      pnDepthNoOcc =
          GetSurfaceDepthMapsOfPanoramix(std::vector<PanoramicCamera>{cam},
                                         anno, options, matlab)
              .front();
    }

    // pn depth gt occ
    {
      PanoramixOptions options;
      options.useWallPrior = true;
      options.usePrincipleDirectionPrior = true;
      options.useGeometricContextPrior = true;

      options.useGTOcclusions = true;
      options.looseLinesSecondTime = false;
      options.looseSegsSecondTime = false;
      options.restrictSegsSecondTime = false;

      options.notUseOcclusions = false;
      options.notUseCoplanarity = false;

      options.refresh_preparation = false;
      options.refresh_mg_init = options.refresh_preparation || false;
      options.refresh_mg_oriented = options.refresh_mg_init || false;
      options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
      options.refresh_lsw = options.refresh_mg_oriented || false;
      options.refresh_mg_occdetected =
          options.refresh_lsw || options.refresh_line2leftRightSegs || false;
      options.refresh_mg_reconstructed =
          options.refresh_mg_occdetected || false;

      pnDepthGTOcc =
          GetSurfaceDepthMapsOfPanoramix(std::vector<PanoramicCamera>{cam},
                                         anno, options, matlab)
              .front();
    }

    auto results = misc::MXA::createStructMatrix(
        1, 1, {"gtDepths", "bestCuboidDepths", "pnDepths", "pnDepthNoOcc"},
        false);
    results.setField("gtDepths", 0, gtDepths);
    results.setField("bestCuboidDepths", 0, bestCuboidDepths);
    results.setField("pnDepths", 0, pnDepths);
    results.setField("pnDepthNoOcc", 0, pnDepthNoOcc);

    misc::SaveCache(impath, "pano_gtDepths", gtDepths);
    misc::SaveCache(impath, "pano_bestCuboidDepths", bestCuboidDepths);
    misc::SaveCache(impath, "pano_pnDepths", pnDepths);
    misc::SaveCache(impath, "pano_pnDepthNoOcc", pnDepthNoOcc);

    return results;
  });

  TaskQueue QcomputePanoGTOccDepth;
  QcomputePanoGTOccDepth.push_back([&](const std::string &impath) -> misc::MXA {
    auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    PanoramicCamera cam = anno.view.camera;

    Imaged pnDepthGTOcc;
    // pn depth gt occ
    {
      PanoramixOptions options;
      options.useWallPrior = true;
      options.usePrincipleDirectionPrior = true;
      options.useGeometricContextPrior = true;

      options.useGTOcclusions = true;
      options.looseLinesSecondTime = false;
      options.looseSegsSecondTime = false;
      options.restrictSegsSecondTime = false;

      options.notUseOcclusions = false;
      options.notUseCoplanarity = false;

      options.refresh_preparation = false;
      options.refresh_mg_init = options.refresh_preparation || false;
      options.refresh_mg_oriented = options.refresh_mg_init || false;
      options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
      options.refresh_lsw = options.refresh_mg_oriented || false;
      options.refresh_mg_occdetected =
          options.refresh_lsw || options.refresh_line2leftRightSegs || false;
      options.refresh_mg_reconstructed =
          options.refresh_mg_occdetected || false;

      pnDepthGTOcc =
          GetSurfaceDepthMapsOfPanoramix(std::vector<PanoramicCamera>{cam},
                                         anno, options, matlab)
              .front();
    }

    auto results = misc::MXA::createStructMatrix(1, 1, {"pnDepthGTOcc"}, false);
    results.setField("pnDepthGTOcc", 0, pnDepthGTOcc);
    misc::SaveCache(impath, "pano_pnDepthGTOcc", pnDepthGTOcc);

    return results;
  });

  TaskQueue QcomputePanoNoCopDepth;
  QcomputePanoNoCopDepth.push_back([&](const std::string &impath) -> misc::MXA {
    auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    PanoramicCamera cam = anno.view.camera;

    Imaged pnDepthNoCop;
    // pn depth gt occ
    {
      PanoramixOptions options;
      options.useWallPrior = true;
      options.usePrincipleDirectionPrior = true;
      options.useGeometricContextPrior = true;

      options.useGTOcclusions = false;
      options.looseLinesSecondTime = false;
      options.looseSegsSecondTime = false;
      options.restrictSegsSecondTime = false;

      options.notUseOcclusions = false;
      options.notUseCoplanarity = true;

      options.refresh_preparation = false;
      options.refresh_mg_init = options.refresh_preparation || false;
      options.refresh_mg_oriented = options.refresh_mg_init || false;
      options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
      options.refresh_lsw = options.refresh_mg_oriented || false;
      options.refresh_mg_occdetected =
          options.refresh_lsw || options.refresh_line2leftRightSegs || false;
      options.refresh_mg_reconstructed =
          options.refresh_mg_occdetected || false;

      pnDepthNoCop =
          GetSurfaceDepthMapsOfPanoramix(std::vector<PanoramicCamera>{cam},
                                         anno, options, matlab)
              .front();
    }

    auto results = misc::MXA::createStructMatrix(1, 1, {"pnDepthNoCop"}, false);
    results.setField("pnDepthNoCop", 0, pnDepthNoCop);
    misc::SaveCache(impath, "pano_pnDepthNoCop", pnDepthNoCop);

    return results;
  });

  TaskQueue QwritePanoMasks;
  QwritePanoMasks.push_back([&](const std::string &impath) -> misc::MXA {
    auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    // estimate a mask representing the original panorama embeded in the
    // rectified panorama
    int w = anno.rectifiedImage.cols;
    int h = anno.rectifiedImage.rows;
    Image3ub im = anno.rectifiedImage;
    int topH = 0;
    if (anno.extendedOnTop) {
      for (; topH < h / 2; topH++) {
        bool hasContent = false;
        for (int i = 0; i < w; i++) {
          if (im(topH, i) != Vec3ub(0, 0, 0)) {
            hasContent = true;
            break;
          }
        }
        if (hasContent) {
          break;
        }
      }
    }
    int bottomH = h - 1;
    if (anno.extendedOnBottom) {
      for (; bottomH > h / 2; bottomH--) {
        bool hasContent = false;
        for (int i = 0; i < w; i++) {
          if (im(bottomH, i) != Vec3ub(0, 0, 0)) {
            hasContent = true;
            break;
          }
        }
        if (hasContent) {
          break;
        }
      }
    }

    Imageb mask(im.size(), false);
    for (int y = topH; y <= bottomH; y++) {
      for (int x = 0; x < w; x++) {
        mask(y, x) = true;
      }
    }

    auto results = misc::MXA::createStructMatrix(1, 1, {"rectified", "mask"});
    results.setField("rectified", 0, im);
    results.setField("mask", 0, mask);
    return results;
  });

  TaskQueue QwriteGT;
  QwriteGT.push_back([&](const std::string &impath) -> misc::MXA {
    auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    auto annoFaceIds =
        GTFaceLabels(anno, std::vector<PanoramicCamera>{anno.view.camera})
            .front();
    Image3ub originalIm = anno.view.image;
    Image3ub pim(annoFaceIds.size());
    for (auto it = pim.begin(); it != pim.end(); ++it) {
      Vec3ub impixel = originalIm(it.pos());
      int faceid = annoFaceIds(it.pos());
      Vec3 normal = normalize(anno.face2plane[faceid].normal);
      Vec3 dirInVPs;
      for (int k = 0; k < 3; k++) {
        dirInVPs[k] = abs(normal.dot(normalize(anno.vps[k])));
      }
      dirInVPs /= norm(dirInVPs);
      static const double imweight = 0.1;
      *it = (Vec3)(impixel)*imweight + (dirInVPs * 255 * (1 - imweight));
    }
    // draw disconnected borders
    for (int i = 0; i < anno.nborders(); i++) {
      if (anno.border2connected[i]) {
        continue;
      }
      Line3 line;
      line.first = anno.corners[anno.border2corners[i].first];
      line.second = anno.corners[anno.border2corners[i].second];
      double angle = AngleBetweenDirected(line.first, line.second);
      std::vector<Pixel> ps;
      double stepAngle = 0.001;
      for (double a = 0.0; a <= angle; a += stepAngle) {
        ps.push_back(ToPixel(anno.view.camera.toScreen(
            RotateDirection(line.first, line.second, a))));
      }
      for (int i = 1; i < ps.size(); i++) {
        auto p1 = ps[i - 1];
        auto p2 = ps[i];
        if (Distance(p1, p2) >= pim.cols / 2) {
          continue;
        }
        cv::clipLine(cv::Rect(0, 0, pim.cols, pim.rows), p1, p2);
        cv::line(pim, p1, p2, cv::Scalar(0, 0, 0), 5);
      }
    }

    gui::AsCanvas(pim).show();
    return misc::MXA(pim);
  });

  activeQ = QshowModel;

  if (true) {
    std::vector<std::string> impaths;
    gui::PickImages("F:\\PanoContext\\", &impaths);
    for (int i = 0; i < activeQ.size(); i++) {
      auto &task = activeQ[i];
      std::cout << "[[[[[[[[[ TASK " << i << "]]]]]]]]" << std::endl;

      auto timeTag = misc::CurrentTimeString(true);
      /*misc::MAT dataFile("F:\\GitHub\\write-papers\\papers\\a\\data\\task_" +
                             std::to_string(i) + timeTag + ".mat",
                         misc::MAT::Write_7_3);*/

      misc::MXA impathsForTask =
          misc::MXA::createCellMatrix(impaths.size(), 1, true);
      misc::MXA resultsForTask =
          misc::MXA::createCellMatrix(impaths.size(), 1, true);
      for (int j = 0; j < impaths.size(); j++) {
        auto &impath = impaths[j];
        misc::Clock clock =
            "Task " + std::to_string(i) + " on \"" + impath + "\"";
        try {
          misc::MXA result = task(impath);
          impathsForTask.setCell(j, misc::MXA::createString(impath));
          resultsForTask.setCell(j, std::move(result));
        } catch (...) {
          std::cout << "############### ERROR #############" << std::endl;
        }
      }

      /*dataFile.setVar("results", resultsForTask, false);
      dataFile.setVar("impaths", impathsForTask, false);*/
    }
  }

  if (false) {
    RTreeMap<Vec<double, 6>, int> cuboidTree;
    std::vector<std::pair<Point3, Point3>> cuboidCands;
    cuboidCands.reserve(1e5);
    static const int N = 10;
    static const double distThres = 0.01;
    for (int x1 = 1; x1 <= N; x1++) {
      for (int y1 = 1; y1 <= N; y1++) {
        for (int z1 = 1; z1 <= N; z1++) {

          for (int x2 = 1; x2 <= N; x2++) {
            for (int y2 = 1; y2 <= N; y2++) {
              for (int z2 = 1; z2 <= N; z2++) {

                Point3 p1(x1, y1, z1);
                Point3 p2(-x2, -y2, -z2);

                double scale = Distance(p1, p2);
                p1 /= scale;
                p2 /= scale;

                bool duplicated = false;
                cuboidTree.search(
                    BoundingBox(cat(p1, p2)).expand(distThres * 2),
                    [&duplicated, &cuboidCands, &p1,
                     &p2](const std::pair<Vec<double, 6>, int> &cbid) {
                      auto &cuboid = cuboidCands[cbid.second];
                      if (Distance(cuboid.first, p1) <= distThres &&
                          Distance(cuboid.second, p2) <= distThres) {
                        duplicated = true;
                        return false;
                      }
                      return true;
                    });
                if (!duplicated) {
                  cuboidCands.emplace_back(p1, p2);
                  cuboidTree.emplace(cat(p1, p2), cuboidCands.size() - 1);
                }
              }
            }
          }
        }
      }
    }

    std::cout << "size of cands: " << cuboidCands.size() << std::endl;
    core::SaveToDisk("F:\\GitHub\\write-papers\\papers\\a\\data\\cuboidcands_" +
                         std::to_string(N) + ".cereal",
                     cuboidCands);

    gui::SceneBuilder sb;
    std::vector<Box3> boxes;
    for (int i = 0; i < 1000; i++) {
      boxes.push_back(Box3(cuboidCands[i].first, cuboidCands[i].second));
    }

    sb.installingOptions().lineWidth = 1.0;
    sb.begin(boxes)
        .shaderSource(gui::OpenGLShaderSourceDescriptor::XLines)
        .end();
    sb.show();
  }

  if (false) {
    std::vector<std::pair<Point3, Point3>> cuboidCands;
    core::LoadFromDisk(
        "F:\\GitHub\\write-papers\\papers\\a\\data\\cuboidcands_10.cereal",
        cuboidCands);
    std::vector<Imagei> cuboidCandsFaceMap(cuboidCands.size());
    ParallelRun(
        cuboidCands.size() / 1000 + 1, std::thread::hardware_concurrency() - 1,
        [&cuboidCands, &cuboidCandsFaceMap](const int ss) {
          std::cout << "start processing " << ss * 1000 << "-"
                    << ((ss + 1) * 1000) << " /" << cuboidCandsFaceMap.size()
                    << std::endl;
          for (int i = ss * 1000;
               i < std::min((int)cuboidCands.size(), (ss + 1) * 1000); i++) {
            int W = 100;
            auto &p2 = cuboidCands[i].first;
            auto &p1 = cuboidCands[i].second;
            // corners
            Vec3 corners[8] = {{p1[0], p1[1], p1[2]}, {p1[0], p2[1], p1[2]},
                               {p2[0], p2[1], p1[2]}, {p2[0], p1[1], p1[2]},
                               {p1[0], p1[1], p2[2]}, {p1[0], p2[1], p2[2]},
                               {p2[0], p2[1], p2[2]}, {p2[0], p1[1], p2[2]}};
            // faces
            int faces[6][4] = {{1, 4, 8, 5}, {1, 5, 6, 2}, {1, 2, 3, 4},
                               {3, 7, 8, 4}, {5, 8, 7, 6}, {2, 6, 7, 3}};
            Polygon3 polys[6];
            for (int k = 0; k < 6; k++) {
              polys[k].corners.resize(4);
              for (int m = 0; m < 4; m++) {
                polys[k].corners[m] = corners[faces[k][m] - 1];
              }
              polys[k].normal =
                  -(polys[k].corners[0] - polys[k].corners[1])
                       .cross(polys[k].corners[1] - polys[k].corners[2]);
              polys[k].normal /= norm(polys[k].normal);
            }

            Imagei faceIdDist(W / 2, W, -1);
            auto view = CreatePanoramicView(faceIdDist);
            for (auto it = faceIdDist.begin(); it != faceIdDist.end(); ++it) {
              auto dir = normalize(view.camera.toSpace(it.pos()));
              Ray3 ray(Origin(), dir);
              double d = std::numeric_limits<double>::infinity();
              for (int faceid = 0; faceid < 6; faceid++) {
                auto test =
                    IntersectionOfLineAndPlane(ray, polys[faceid].plane());
                if (test.ratio > 0) {
                  if (test.ratio < d) {
                    *it = faceid;
                    d = test.ratio;
                  }
                }
              }
            }

            cuboidCandsFaceMap[i] = faceIdDist;
          }
        });

    core::SaveToDisk("F:\\GitHub\\write-papers\\papers\\a\\data\\cuboidcands_"
                     "facemap_5.cereal",
                     cuboidCandsFaceMap);
  }

  return 0;
}