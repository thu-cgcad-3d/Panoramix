#include "panorama_reconstruction.hpp"

int main_run(int argc, char **argv,
             const std::map<std::string, std::string> &options) {
  bool showUI = true;
  if (options.count("ui")) {
    showUI = options.at("ui") == "on";
  }
  bool refreshAll = false;
  if (options.count("refresh")) {
    refreshAll = options.at("refresh") == "on";
  }
  std::string dataDir = PANORAMIX_TEST_DATA_DIR_STR;
  if (options.count("data")) {
    dataDir = options.at("data");
  }
  bool saveResults = false;
  if (options.count("save")) {
    saveResults = options.at("save") == "on";
  }

  gui::UI::InitGui(argc, argv);
  misc::SetCachePath(PANORAMIX_CACHE_DATA_DIR_STR "\\Panorama\\");
  misc::Matlab matlab;

  std::vector<std::string> impaths;
  gui::FileDialog::PickImages(dataDir, &impaths);

  for (auto &&impath : impaths) {
    auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    anno.impath = impath;

    PanoramaReconstructionOptions options;
    options.useWallPrior = true;
    options.usePrincipleDirectionPrior = true;
    options.useGeometricContextPrior = true;

    options.useGTOcclusions = false;
    options.looseLinesSecondTime = false;
    options.looseSegsSecondTime = false;
    options.restrictSegsSecondTime = false;

    options.notUseOcclusions = false;
    options.notUseCoplanarity = false;

    options.refresh_preparation = refreshAll;
    options.refresh_mg_init = options.refresh_preparation || false;
    options.refresh_mg_oriented = options.refresh_mg_init || false;
    options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
    options.refresh_lsw = options.refresh_mg_oriented || false;
    options.refresh_mg_occdetected =
        options.refresh_lsw || options.refresh_line2leftRightSegs || false;
    options.refresh_mg_reconstructed = options.refresh_mg_occdetected || false;

    RunPanoramaReconstruction(anno, options, matlab, showUI, false);

    if (saveResults) {
      SaveMatlabResultsOfPanoramaReconstruction(anno, options, matlab,
                                                impath + ".result.mat");
      SaveObjModelResultsOfPanoramaReconstruction(anno, options, matlab,
                                                  impath + ".result.obj");
    }
  }

  return 0;
}