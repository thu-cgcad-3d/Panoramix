#pragma once

#include "parallel.hpp"
#include "cache.hpp"
#include "clock.hpp"

#include "canvas.hpp"
#include "singleton.hpp"
#include "gui_util.hpp"

#include "pi_graph_annotation.hpp"
#include "pi_graph_cg.hpp"
#include "pi_graph_control.hpp"
#include "pi_graph_occlusion.hpp"
#include "pi_graph_optimize.hpp"
#include "pi_graph_vis.hpp"

using namespace pano;
using namespace pano::core;
using namespace pano::experimental;

// options
struct PanoramixOptions {
  // algorithm options
  static const int LayoutVersion = 0;
  bool useWallPrior;
  bool usePrincipleDirectionPrior;
  bool useGeometricContextPrior;
  bool useGTOcclusions;
  bool looseLinesSecondTime;
  bool looseSegsSecondTime;
  bool restrictSegsSecondTime;
  bool notUseOcclusions;

  bool notUseCoplanarity;

  static const std::string parseOption(bool b);
  std::string algorithmOptionsTag() const;
  std::string identityOfImage(const std::string &impath) const;

  // cache options
  bool refresh_preparation;
  bool refresh_mg_init;
  bool refresh_line2leftRightSegs;
  bool refresh_mg_oriented;
  bool refresh_lsw;
  bool refresh_mg_occdetected;
  bool refresh_mg_reconstructed;

  // print options out
  void print() const;

  template <class Archiver> void serialize(Archiver &ar) {
    ar(useWallPrior, usePrincipleDirectionPrior, useGeometricContextPrior,
       useGTOcclusions, looseLinesSecondTime, looseSegsSecondTime,
       restrictSegsSecondTime, notUseOcclusions, notUseCoplanarity);
    ar(refresh_preparation, refresh_mg_init, refresh_line2leftRightSegs,
       refresh_mg_oriented, refresh_lsw, refresh_mg_occdetected,
       refresh_mg_reconstructed);
  }
};

// report
struct PanoramixReport {
  double time_preparation;
  double time_mg_init;
  double time_line2leftRightSegs;
  double time_mg_oriented;
  double time_lsw;
  double time_mg_occdetected;
  double time_mg_reconstructed;

  double time_solve_lp;

  bool succeeded;

  PanoramixReport();

  void print() const;

  template <class Archiver> void serialize(Archiver &ar) {
    ar(time_preparation, time_mg_init, time_line2leftRightSegs,
       time_mg_oriented, time_lsw, time_mg_occdetected, time_mg_reconstructed,
       time_solve_lp, succeeded);
  }
};

// run the main algorithm
PanoramixReport RunPanoramix(const PILayoutAnnotation &anno,
                             const PanoramixOptions &options,
                             misc::Matlab &matlab, bool showGUI,
                             bool writeToFile = false);

// get result
bool GetPanoramixResult(const PILayoutAnnotation &anno,
                        const PanoramixOptions &options,
                        PIGraph<PanoramicCamera> &mg, PIConstraintGraph &cg,
                        PICGDeterminablePart &dp);
std::vector<LineSidingWeight>
GetPanoramixOcclusionResult(const PILayoutAnnotation &anno,
                            const PanoramixOptions &options);

// get surface normal maps
template <class CameraT>
std::vector<Image3d> GetSurfaceNormalMapsOfPanoramix(
    const std::vector<CameraT> &testCams, const PILayoutAnnotation &anno,
    const PanoramixOptions &options, misc::Matlab &matlab) {

  PIGraph<PanoramicCamera> mg;
  PIConstraintGraph cg;
  PICGDeterminablePart dp;
  if (!GetPanoramixResult(anno, options, mg, cg, dp)) {
    std::cout << "failed to load panoramix result, performing RunPanoramix ..."
              << std::endl;
    RunPanoramix(anno, options, matlab, false);
    GetPanoramixResult(anno, options, mg, cg, dp);
  }

  std::vector<Image3d> surfaceNormalMaps(testCams.size());
  ParallelRun(
      testCams.size(), std::thread::hardware_concurrency() - 1, [&](int i) {
        std::cout
            << "computing surface normal map for panoramix on testCamera - "
            << i << std::endl;
        auto &map = surfaceNormalMaps[i];
        auto &cam = testCams[i];
        map = SurfaceNormalMap(cam, dp, cg, mg, true);
      });
  return surfaceNormalMaps;
}

// get surface depth maps
template <class CameraT>
std::vector<Imaged> GetSurfaceDepthMapsOfPanoramix(
    const std::vector<CameraT> &testCams, const PILayoutAnnotation &anno,
    const PanoramixOptions &options, misc::Matlab &matlab) {

  PIGraph<PanoramicCamera> mg;
  PIConstraintGraph cg;
  PICGDeterminablePart dp;
  if (!GetPanoramixResult(anno, options, mg, cg, dp)) {
    std::cout << "failed to load panoramix result, performing RunPanoramix ..."
              << std::endl;
    RunPanoramix(anno, options, matlab, false);
    GetPanoramixResult(anno, options, mg, cg, dp);
  }

  std::vector<Imaged> surfaceDepthMaps(testCams.size());
  ParallelRun(
      testCams.size(), std::thread::hardware_concurrency() - 1, [&](int i) {
        std::cout
            << "computing surface normal map for panoramix on testCamera - "
            << i << std::endl;
        auto &map = surfaceDepthMaps[i];
        auto &cam = testCams[i];
        map = SurfaceDepthMap(cam, dp, cg, mg, true);
      });
  return surfaceDepthMaps;
}
